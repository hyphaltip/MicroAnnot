"""
Ribosomal RNA detection via blastn against a 16S rRNA database.

Returns a nested dict of rRNA locus coordinates::

    {seq_id: {start: {end: {'c': strand}}}}

where start/end are 1-based +strand coordinates (extended by ±2500 nt to
capture the full rRNA unit) and ``c`` is the subject strand (1 or -1).
"""

import os
import subprocess
from Bio.Blast import NCBIXML


def run_rrna_blast(
    input_fasta: str,
    db: str,
    out_file: str,
    dir_blast: str = "",
    evalue: float = 1e-50,
) -> None:
    """Run blastn to detect 16S rRNA loci.

    Parameters
    ----------
    input_fasta:
        Preprocessed FASTA file.
    db:
        Path to the 16S rRNA BLAST database (without extension).
    out_file:
        Destination for BLAST XML output (``-outfmt 5``).
    dir_blast:
        Directory containing BLAST+ binaries (may be empty if they are on PATH).
    evalue:
        E-value cutoff.
    """
    blastn = os.path.join(dir_blast, "blastn") if dir_blast else "blastn"
    db_nin = db + ".nin"
    if not os.path.exists(db_nin):
        makeblastdb = (
            os.path.join(dir_blast, "makeblastdb") if dir_blast else "makeblastdb"
        )
        subprocess.run(
            [makeblastdb, "-in", db, "-out", db, "-dbtype", "nucl"],
            check=True,
        )
    subprocess.run(
        [
            blastn,
            "-num_alignments", "500",
            "-word_size", "7",
            "-gapopen", "4",
            "-gapextend", "2",
            "-penalty", "-1",
            "-reward", "1",
            "-evalue", str(evalue),
            "-db", db,
            "-query", input_fasta,
            "-out", out_file,
            "-outfmt", "5",
        ],
        check=True,
    )


def parse_rrna_blast(out_file: str, annotation: dict) -> dict:
    """Parse the blastn XML output to build the rRNA dict.

    Only the *best* hit for each query is used (first alignment), matching
    the original Perl ``last;`` after the inner hit loop.

    Parameters
    ----------
    out_file:
        BLAST XML file produced by :func:`run_rrna_blast`.
    annotation:
        The master annotation dict (needed for sequence lengths).

    Returns
    -------
    rrna : dict
        ``{seq_id: {start: {end: {'c': strand}}}}``
    """
    rrna: dict = {}
    if not os.path.exists(out_file):
        return rrna

    with open(out_file) as fh:
        for record in NCBIXML.parse(fh):
            if not record.alignments:
                continue
            display_id = record.query.split()[0]
            bioseq_len = len(annotation[display_id]["bioseq"].seq)

            # Only process the best hit (first alignment)
            alignment = record.alignments[0]
            for hsp in alignment.hsps:
                strand = hsp.strand[1]   # subject strand: ('Plus','Plus')…
                strand_int = -1 if strand == "Minus" else 1
                start = hsp.query_start
                end = hsp.query_end

                if strand_int == -1:
                    start = start - 2500
                else:
                    end = end + 2500

                start = max(1, start)
                end = min(bioseq_len, end)
                rrna.setdefault(display_id, {}).setdefault(
                    start, {}
                )[end] = {"c": strand_int}

    return rrna
