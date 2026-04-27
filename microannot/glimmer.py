"""
Glimmer3 gene prediction runner and start-codon refinement.

Workflow
--------
1. If ≥ 50 high-confidence CDS were found by homology, train a new Glimmer ICM
   from those CDS (``build-icm``); otherwise use the supplied pre-trained ICM.
2. Run ``glimmer3`` on the preprocessed FASTA.
3. Parse ``res_glimmer.predict`` and, for each predicted ORF, scan upstream for
   the best-scoring Met codon (using the microsporidia start-signal heuristics).

Returns a dict with the structure::

    {seq_id: {start: {end: {'c': complement, 's': note_signal, 'o': orf_range}}}}

where start/end are 1-based +strand coordinates, ``c`` is 1 (forward) or -1
(reverse), and ``s`` / ``o`` are diagnostic strings.
"""

import os
import re
import subprocess
from Bio.Seq import Seq

from microannot.signal import check_signal
from microannot.find_orf import find_orf


def _trunc(seq_str: str, a: int, b: int) -> str:
    return seq_str[a - 1: b]


def _translate_fwd(seq_str: str, a: int, b: int) -> str:
    s = _trunc(seq_str, a, b)
    return str(Seq(s).translate()) if s else ""


def _translate_rev(seq_str: str, a: int, b: int) -> str:
    s = _trunc(seq_str, a, b)
    return str(Seq(s).reverse_complement().translate()) if s else ""


# ---------------------------------------------------------------------------
# Build ICM + run Glimmer
# ---------------------------------------------------------------------------

def build_icm(fasta_file: str, icm_file: str, dir_glimmer: str = "") -> None:
    """Build a Glimmer ICM from a FASTA file of CDS sequences."""
    build_icm_bin = (
        os.path.join(dir_glimmer, "build-icm") if dir_glimmer else "build-icm"
    )
    with open(fasta_file) as stdin_fh:
        subprocess.run(
            [build_icm_bin, "-r", icm_file],
            stdin=stdin_fh,
            check=True,
        )


def run_glimmer(
    input_fasta: str,
    icm_file: str,
    out_prefix: str,
    glimmer_size: int = 300,
    dir_glimmer: str = "",
) -> None:
    """Run glimmer3 on the preprocessed FASTA.

    Parameters
    ----------
    input_fasta:
        Preprocessed FASTA file.
    icm_file:
        Path to a pre-built Glimmer ICM.
    out_prefix:
        Prefix for Glimmer output files (``{out_prefix}.predict``,
        ``{out_prefix}.detail``).
    glimmer_size:
        Minimum CDS size passed to ``-g``.
    dir_glimmer:
        Directory containing Glimmer binaries (empty → use PATH).
    """
    glimmer3 = (
        os.path.join(dir_glimmer, "glimmer3") if dir_glimmer else "glimmer3"
    )
    subprocess.run(
        [
            glimmer3,
            "-X",
            "-A", "atg",
            "-l",
            "-o10",
            "-g", str(glimmer_size),
            "-t", "30",
            input_fasta,
            icm_file,
            out_prefix,
        ],
        check=True,
    )


# ---------------------------------------------------------------------------
# Parse + refine
# ---------------------------------------------------------------------------

def parse_and_refine_glimmer(
    predict_file: str,
    annotation: dict,
    debug: bool = False,
) -> dict:
    """Parse ``glimmer3`` ``.predict`` output and refine start codons.

    For each predicted ORF the function:

    * Recovers the full ORF boundaries via :func:`~microannot.find_orf.find_orf`.
    * Scans from the ORF start toward the Glimmer-predicted start for a Met
      codon preceded by a microsporidia start signal.
    * Records the refined start (or the Glimmer start if no signal is found).

    Parameters
    ----------
    predict_file:
        Path to ``res_glimmer.predict``.
    annotation:
        Master annotation dict (needs ``bioseq`` entries).
    debug:
        Enable extra diagnostic prints.

    Returns
    -------
    glim : dict
        ``{seq_id: {start: {end: {'c': int, 's': str, 'o': str}}}}``
    """
    glim: dict = {}
    SIZE_CHECK_AA = 30

    with open(predict_file) as fh:
        display_id = ""
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue

            if line.startswith(">"):
                display_id = line[1:].split()[0]
                continue

            cols = line.split()
            if len(cols) < 4:
                continue

            # Glimmer columns: name, start, end, frame_string, [score]
            raw_s = int(cols[1])
            raw_e = int(cols[2])
            frame_str = cols[3]          # e.g. "+1" or "-2"

            bioseq_str = str(annotation[display_id]["bioseq"].seq).upper()
            seq_len = len(bioseq_str)

            # --- Boundary sanity checks (matching original Perl) ----------
            if raw_s > seq_len:
                raw_s -= 3
            if raw_e > seq_len:
                raw_e -= 3
            if raw_s <= 0:
                raw_s += 3
            if raw_e <= 0:
                raw_e += 3

            strand = -1 if "-" in frame_str else 1
            frame_val = int(frame_str.replace("-", "")) - 1  # 0-indexed

            if strand == 1:
                orf_start, orf_stop = find_orf(
                    bioseq_str, strand, frame_val, raw_s, raw_e
                )
            else:
                neg_frame = (seq_len - raw_e + 1) % 3  # BioPerl recalc
                orf_start, orf_stop = find_orf(
                    bioseq_str, strand, neg_frame, raw_e, raw_s
                )

            note_signal = ""
            bool_signal = False
            complement = 1

            if raw_e < raw_s:
                # --- Reverse strand ---
                complement = -1
                a = raw_e
                b = orf_start
                prot = _translate_rev(bioseq_str, a, b)

                # Limit search to Glimmer start + SIZE_CHECK_AA AAs
                glimmer_offset = (b - raw_s) // 3
                check_len = glimmer_offset + SIZE_CHECK_AA + 1
                check_AA = prot[:check_len] if len(prot) >= check_len else prot

                i = 0
                for mobj in re.finditer(r"(.*?)M", check_AA):
                    i += len(mobj.group(1)) * 3
                    b2 = b - i
                    end2 = b2 + 21
                    if end2 > seq_len + 1:
                        end2 = seq_len + 1
                    if b2 + 1 >= end2 - 1:
                        amont_M_20nt = ""
                    else:
                        amont_M_20nt = str(
                            Seq(_trunc(bioseq_str, b2 + 1, end2 - 1))
                            .reverse_complement()
                        )
                    note_signal = ""
                    if amont_M_20nt:
                        bool_signal, note_signal, _ = check_signal(amont_M_20nt)
                        if i < orf_start - raw_s and "weak" in note_signal:
                            bool_signal = False

                    if bool_signal:
                        note_signal = "M with a signal " + note_signal
                        break
                    else:
                        note_signal = "M without a signal"
                    i += 3

                if bool_signal:
                    rel = (i - (orf_start - raw_s)) // 3
                    if rel != 0:
                        if i < orf_start - raw_s:
                            note_signal += (
                                f". The selected M is at {-rel} AA before the M "
                                f"selected by Glimmer because it had a signal"
                            )
                        else:
                            note_signal += (
                                f". The selected M is at {rel} AA after the M "
                                f"selected by Glimmer because it had a signal"
                            )
                        raw_s = orf_start - i
                # Swap so that start < end for storage
                raw_s, raw_e = raw_e, raw_s
                raw_s += 3

            else:
                # --- Forward strand ---
                a = orf_start
                b = raw_e
                prot = _translate_fwd(bioseq_str, a, b)

                glimmer_offset = (raw_s - a) // 3
                check_len = glimmer_offset + SIZE_CHECK_AA + 1
                check_AA = prot[:check_len] if len(prot) >= check_len else prot

                i = 0
                for mobj in re.finditer(r"(.*?)M", check_AA):
                    i += len(mobj.group(1)) * 3
                    a2 = orf_start + i
                    s2 = a2 - 21
                    if s2 < 1:
                        s2 = 1
                    if s2 + 1 >= a2 - 1:
                        amont_M_20nt = ""
                    else:
                        amont_M_20nt = _trunc(bioseq_str, s2 + 1, a2 - 1)
                    note_signal = ""
                    if amont_M_20nt:
                        bool_signal, note_signal, _ = check_signal(amont_M_20nt)
                        if i < raw_s - orf_start and "weak" in note_signal:
                            bool_signal = False

                    if bool_signal:
                        note_signal = "M with a signal " + note_signal
                        break
                    else:
                        note_signal = "M without a signal"
                    i += 3

                if bool_signal:
                    rel = (i - (raw_s - orf_start)) // 3
                    if rel != 0:
                        if i < raw_s - orf_start:
                            note_signal += (
                                f". The selected M is at {-rel} AA before the M "
                                f"selected by Glimmer because it had a signal"
                            )
                        else:
                            note_signal += (
                                f". The selected M is at {rel} AA after the M "
                                f"selected by Glimmer because it had a signal"
                            )
                    raw_s = orf_start + i
                raw_e -= 3

            if not note_signal:
                note_signal = "Prediction without M"

            glim.setdefault(display_id, {}).setdefault(
                raw_s, {}
            )[raw_e] = {
                "c": complement,
                "s": note_signal,
                "o": f"{orf_start}-{orf_stop}",
            }

    return glim
