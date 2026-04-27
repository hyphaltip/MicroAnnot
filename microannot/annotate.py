"""
Annotation merge, overlap resolution and output writers.

This module is the Python equivalent of the large output-generation block in the
original Perl script.  It:

1. Merges all evidence sources (blastp CDS, small CDS blastx, tRNA, rRNA,
   Glimmer predictions) into a per-sequence ordered feature list.
2. Resolves overlaps (small CDS are suppressed when they overlap a homology CDS;
   Glimmer predictions are filtered against all other features).
3. Computes start/stop completeness for each feature and sets the appropriate
   partial-feature flags (``<1``, ``>end``).
4. Writes per-sequence GenBank, EMBL, and GFF3 files to three result directories.
5. Archives those directories as ``.tar.gz`` files.

All coordinate arithmetic uses 1-based inclusive coordinates (BioPerl convention)
except when constructing Biopython ``FeatureLocation`` objects, which use 0-based
half-open coordinates.
"""

import copy
import os
import re
import tarfile
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    AfterPosition,
    BeforePosition,
    ExactPosition,
)
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _trunc(seq_str: str, a: int, b: int) -> str:
    return seq_str[a - 1: b]


def _translate_fwd(seq_str: str, a: int, b: int) -> str:
    s = _trunc(seq_str, a, b)
    return str(Seq(s).translate()) if s else ""


def _translate_rev(seq_str: str, a: int, b: int) -> str:
    s = _trunc(seq_str, a, b)
    return str(Seq(s).reverse_complement().translate()) if s else ""


def _feat_location(start: int, end: int, strand: int,
                   bool_start: bool, bool_stop: bool,
                   seq_len: int) -> FeatureLocation:
    """Build a Biopython ``FeatureLocation`` with partial-feature handling.

    All input coordinates are 1-based inclusive.
    Biopython uses 0-based half-open, so start_0 = start-1, end_0 = end.
    """
    start_0 = start - 1

    if strand == -1:
        # For revcomp features the "start" (5′-most in +strand) may be partial
        if not bool_stop:
            bp_start: object = BeforePosition(0)
        else:
            bp_start = ExactPosition(start_0)
        if end >= seq_len - 3 and not bool_start:
            bp_end: object = AfterPosition(end)
        else:
            bp_end = ExactPosition(end)
    else:
        if start <= 3 and not bool_start:
            bp_start = BeforePosition(0)
        else:
            bp_start = ExactPosition(start_0)
        if not bool_stop:
            bp_end = AfterPosition(seq_len)
        else:
            bp_end = ExactPosition(end)

    return FeatureLocation(bp_start, bp_end, strand=strand)


# ---------------------------------------------------------------------------
# Small-CDS overlap filter
# ---------------------------------------------------------------------------

def _filter_small_cds(small_cds: dict, cds_start_end: dict) -> dict:
    """Return small CDS entries that do not overlap any homology-based CDS."""
    filtered = {}
    # Sort existing CDS by start for efficient early-exit
    cds_list = sorted(
        ((s, e) for s, ends in cds_start_end.items() for e in ends),
        key=lambda x: x[0],
    )
    for anno_start, entry in sorted(small_cds.items()):
        anno_end = entry["end"]
        overlaps = False
        for (cs, ce) in cds_list:
            if anno_end < cs:
                break
            if anno_start <= ce and anno_end >= cs:
                overlaps = True
                break
        if not overlaps:
            filtered[anno_start] = entry
    return filtered


# ---------------------------------------------------------------------------
# Glimmer overlap filter
# ---------------------------------------------------------------------------

def _select_glimmer(glim_seq: dict, cds_start_end: dict,
                    others_start_end: dict) -> dict:
    """Select non-overlapping Glimmer predictions.

    Mirrors the complex Glimmer-selection block in the original Perl script.
    Returns a dict ``{start: {end: glim_entry}}`` of kept predictions.
    """
    # Pool all Glimmer predictions that don't overlap tRNA/rRNA
    candidates = {}   # (start, end) -> glim entry
    for start, ends in sorted(glim_seq.items()):
        for end, entry in sorted(ends.items()):
            # Check against tRNA / rRNA
            overlap = False
            for (os_, oe) in sorted(
                (s, e)
                for s, es in others_start_end.items()
                for e in es
            ):
                if end < os_:
                    break
                if start <= oe and end >= os_:
                    overlap = True
                    break
            if not overlap:
                candidates[(start, end)] = entry

    # Group overlapping candidates
    groups = []
    group_bounds = []
    for (start, end), entry in sorted(candidates.items()):
        placed = False
        for gi, (gs, ge) in enumerate(group_bounds):
            if start <= ge and end >= gs:
                groups[gi].append(((start, end), entry))
                group_bounds[gi] = (
                    min(gs, start), max(ge, end)
                )
                placed = True
                break
        if not placed:
            groups.append([((start, end), entry)])
            group_bounds.append((start, end))

    selected = {}
    for group in groups:
        # Sort by size descending (largest first) to pick the representative
        group_sorted = sorted(group, key=lambda x: x[0][1] - x[0][0],
                              reverse=True)

        low_start = None
        max_end = 0
        first_big = None

        for ((start, end), entry) in group_sorted:
            size = end - start + 1
            big = size > 500
            align_val = entry.get("alignment")

            if low_start is None:
                low_start = start
                max_end = end
                first_big = big

            if start <= max_end and end >= low_start:
                # Check if this overlaps any homology CDS
                overlaps_cds = False
                for (cs, ce) in sorted(
                    (s, e)
                    for s, es in cds_start_end.items()
                    for e in es
                ):
                    if end < cs:
                        break
                    if start <= ce and end >= cs:
                        overlaps_cds = True
                        break

                if overlaps_cds and (
                    (big and first_big)
                    and (
                        (end >= low_start and end <= max_end
                         and end - low_start + 1 <= 50)
                        or (start >= low_start and start <= max_end
                            and max_end - start + 1 <= 50)
                    )
                ):
                    selected.setdefault(start, {})[end] = entry
                else:
                    selected.setdefault(start, {})[end] = entry

                if start < low_start:
                    low_start = start
                if end > max_end:
                    max_end = end

    return selected


# ---------------------------------------------------------------------------
# start_end_extremity  (port of the Perl sub of the same name)
# ---------------------------------------------------------------------------

def _start_end_extremity(
    translation: str,
    strand: int,
    start: int,
    end: int,
    seq_len: int,
    mode_translate: str,
    start_orf: int,
    end_orf: int,
    note: str = "",
):
    """Determine start/stop completeness and produce feature info.

    Parameters
    ----------
    translation:
        Amino-acid translation of the CDS region (may end with ``*``).
    strand:
        1 or -1.
    start, end:
        1-based +strand coordinates of the translated region
        (includes the stop codon).
    seq_len:
        Length of the parent sequence.
    mode_translate:
        ``"all"`` (normal CDS), ``"extremity"`` (partial/no-M CDS),
        ``"no"`` (tRNA/rRNA – no translation check),
        ``"frameshift_stop"`` / ``"frameshift_nostop"``.
    start_orf, end_orf:
        1-based +strand boundaries of the enclosing ORF.
    note:
        Annotation note (for warning generation).

    Returns
    -------
    feat_start : int or str  – the 1-based feature start (or "<1")
    feat_end   : int or str  – the 1-based feature end   (or ">len")
    bool_start : bool
    bool_stop  : bool
    warning_note : str  – non-empty if the start is far from the ORF start
    """
    bool_start = False
    bool_stop = False
    warning_note = ""

    if mode_translate == "no":
        return start, end, False, False, ""

    # Determine bool_start / bool_stop from translation
    tmp = translation
    if mode_translate == "frameshift_stop":
        bool_stop = True
    elif mode_translate == "frameshift_nostop":
        bool_stop = False
    else:
        if tmp.endswith("*"):
            tmp = tmp[:-1]
            bool_stop = True
        if tmp.startswith("M"):
            bool_start = True

    feat_start = start
    feat_end = end

    if strand == -1:
        if not bool_stop:
            feat_start = "<1"
        if end >= seq_len - 3 and not bool_start:
            feat_end = f">{end}"
        # Extremity check
        if start_orf > seq_len - 3:
            warning_note = ""   # extremity – no warning
        if start_orf - end > 300:
            warning_note = (
                f"Warning: the start is at {start_orf - end} nucleotides "
                f"from the start of the ORF"
            )
    else:
        if start <= 3 and not bool_start:
            feat_start = "<1"
        if not bool_stop:
            feat_end = f">{seq_len}"
        # Extremity check
        if start_orf <= 3:
            warning_note = ""   # extremity – no warning
        if start - start_orf > 300:
            warning_note = (
                f"Warning: the start is at {start - start_orf} nucleotides "
                f"from the start of the ORF"
            )

    return feat_start, feat_end, bool_start, bool_stop, warning_note


# ---------------------------------------------------------------------------
# GFF3 writer helper
# ---------------------------------------------------------------------------

def _gff3_line(seq_id: str, source: str, feat_type: str,
               start, end, strand: int, phase: int,
               attributes: dict) -> str:
    """Format one GFF3 data line."""
    strand_ch = "+" if strand >= 0 else "-"
    start_s = str(start).lstrip("<>")
    end_s = str(end).lstrip("<>")
    attr_str = ";".join(
        f"{k}={v}" for k, v in attributes.items() if v
    )
    return (
        f"{seq_id}\tMicroAnnot\t{feat_type}\t{start_s}\t{end_s}\t.\t"
        f"{strand_ch}\t{phase}\t{attr_str}"
    )


# ---------------------------------------------------------------------------
# Core output writer
# ---------------------------------------------------------------------------

def write_annotations(
    annotation: dict,
    small_cds: dict,
    trna: dict,
    rrna: dict,
    glim: dict,
    id_2_name: dict,
    out_dir: str,
    cds_gene_nt_path: str,
    debug: bool = False,
) -> tuple:
    """Write per-sequence GenBank, EMBL, and GFF3 files (pre-TE/InterproScan).

    Returns
    -------
    list_file_cds_aa : dict
        Maps ``"start-end"`` strings to amino-acid sequences (used later for
        TE / InterproScan annotation).
    """
    dir_gb = os.path.join(out_dir, "Res_gb")
    dir_embl = os.path.join(out_dir, "Res_embl")
    dir_warnings = os.path.join(out_dir, "Res_warnings")
    for d in (dir_gb, dir_embl, dir_warnings):
        os.makedirs(d, exist_ok=True)

    list_file_cds_aa: dict = {}
    list_name: dict = {}

    with open(cds_gene_nt_path, "w") as cds_nt_out:
        for display_id in sorted(annotation, key=lambda x: int(x)):
            seq_rec = annotation[display_id]["bioseq"]
            seq_str = str(seq_rec.seq).upper()
            seq_len = len(seq_str)

            # Build a clean SeqRecord for output
            out_rec = SeqRecord(
                seq_rec.seq,
                id=id_2_name.get(display_id, display_id).lstrip(">"),
                description="",
            )
            out_rec.annotations["molecule_type"] = "DNA"

            # Sanitise the file-system name
            raw_name = id_2_name.get(display_id, display_id)
            raw_name = raw_name.lstrip(">").replace(" ", "_").replace(">", "")
            tmp_name = raw_name[:20]
            uniq_name = tmp_name
            suffix = 2
            while uniq_name in list_name:
                uniq_name = f"{tmp_name}_{suffix}"
                suffix += 1
            list_name[uniq_name] = 1

            gb_path = os.path.join(dir_gb, f"{uniq_name}.gb")
            embl_path = os.path.join(dir_embl, f"{uniq_name}.embl")
            warnings_path = os.path.join(dir_warnings, f"{uniq_name}.txt")

            # ---- Build ordered feature list --------------------------------
            order = {}   # start_pos (int) → label
            cds_start_end: dict = {}  # {start: {end: 1}}

            # 1. Homology-based CDS
            for orf_id, orf_data in annotation[display_id].items():
                if orf_id == "bioseq":
                    continue
                cds = orf_data.get("CDS")
                if cds:
                    order[cds["start"]] = orf_id
                    cds_start_end.setdefault(
                        cds["start"], {}
                    )[cds["end"]] = 1

            # 2. Small CDS (only if not overlapping homology CDS)
            seq_small = small_cds.get(display_id, {})
            kept_small = _filter_small_cds(seq_small, cds_start_end)
            for anno_start in sorted(kept_small):
                order[anno_start] = "small_cds"
                cds_start_end.setdefault(
                    anno_start, {}
                )[kept_small[anno_start]["end"]] = 1

            # 3. tRNA
            others_start_end: dict = {}
            for t_start in sorted(trna.get(display_id, {})):
                for t_end in sorted(trna[display_id][t_start]):
                    order_s = min(t_start, t_end)
                    order_e = max(t_start, t_end)
                    key = f"trna:{t_start}-{t_end}"
                    order[order_s] = key
                    others_start_end.setdefault(order_s, {})[order_e] = 1

            # 4. rRNA
            for r_start in sorted(rrna.get(display_id, {})):
                for r_end in sorted(rrna[display_id][r_start]):
                    order_s = min(r_start, r_end)
                    order_e = max(r_start, r_end)
                    key = f"rrna:{r_start}-{r_end}"
                    order[order_s] = key
                    others_start_end.setdefault(order_s, {})[order_e] = 1

            # 5. Glimmer (filtered)
            glim_sel = _select_glimmer(
                glim.get(display_id, {}), cds_start_end, others_start_end
            )
            for g_start in sorted(glim_sel):
                for g_end in sorted(glim_sel[g_start]):
                    entry = glim_sel[g_start][g_end]
                    key = (
                        f"Glimmer:{g_start}-{g_end}:"
                        f"{entry['s']}"
                    )
                    order[g_start] = key

            # ---- Write features -------------------------------------------
            gene_counter = 0
            with open(warnings_path, "w") as warn_out:
                warned_header = False

                for a_order in sorted(order):
                    label = order[a_order]
                    qualifiers: dict = {}
                    feat_type = "CDS"
                    mode_translate = "all"
                    strand_f = 1
                    a = 0
                    b = 0
                    start_orf = 0
                    end_orf = 0

                    if label == "small_cds":
                        entry = kept_small[a_order]
                        end_val = entry["end"]
                        orf_parts = entry["orf"].split("-")
                        start_orf = int(orf_parts[0])
                        end_orf = int(orf_parts[1])
                        strand_f = entry["strand"]

                        sc_note = entry["note"]
                        if ("Alignment without M and there is no M at all "
                                "upstream of the alignment") in sc_note:
                            feat_type = "gene"
                            sc_note += (" Possible frameshift or intron or "
                                        "N-terminal truncated protein.")
                            mode_translate = "extremity"
                        else:
                            feat_type = "CDS"
                            mode_translate = "all"

                        if strand_f == -1:
                            a = a_order - 3
                            if a <= 0:
                                a += 3
                            b = end_val
                            nt = _trunc(seq_str, a, b)
                            tmp_prot = _translate_rev(seq_str, a, b)
                        else:
                            a = a_order
                            b = end_val + 3
                            if b > seq_len:
                                b -= 3
                            nt = _trunc(seq_str, a, b)
                            tmp_prot = _translate_fwd(seq_str, a, b)

                        f_start, f_end, bs, bst, w_note = (
                            _start_end_extremity(
                                tmp_prot if mode_translate != "no" else "",
                                strand_f, a, b, seq_len, mode_translate,
                                start_orf, end_orf, sc_note,
                            )
                        )
                        qualifiers["hit"] = entry["hit"]
                        qualifiers["evalue"] = str(entry["evalue"])
                        qualifiers["note"] = sc_note
                        if entry["warning"] == 1:
                            if not warned_header:
                                warn_out.write(f">{uniq_name}\n")
                                warned_header = True
                            warn_out.write(
                                f"{feat_type}-{a}-{b}\n{sc_note}\n"
                            )
                        if mode_translate == "all" and tmp_prot:
                            tp = tmp_prot.rstrip("*")
                            qualifiers["translation"] = tp

                    elif label.startswith("Glimmer:"):
                        feat_type = "CDS"
                        parts = label.split(":")
                        pos_parts = parts[1].split("-")
                        g_start_pos = int(pos_parts[0])
                        g_end_pos = int(pos_parts[1])
                        g_entry = glim_sel.get(g_start_pos, {}).get(g_end_pos)
                        if g_entry is None:
                            continue
                        strand_f = g_entry["c"]
                        orf_range = g_entry["o"].split("-")
                        start_orf = int(orf_range[0])
                        end_orf = int(orf_range[1])
                        g_note = ":".join(parts[2:])
                        qualifiers["note"] = f"CDS found by Glimmer. {g_note}"

                        if strand_f == -1:
                            a = g_start_pos - 3
                            if a <= 0:
                                a += 3
                            b = g_end_pos
                            nt = _trunc(seq_str, a, b)
                            tmp_prot = _translate_rev(seq_str, a, b)
                        else:
                            a = g_start_pos
                            b = g_end_pos + 3
                            if b > seq_len:
                                b -= 3
                            nt = _trunc(seq_str, a, b)
                            tmp_prot = _translate_fwd(seq_str, a, b)

                        f_start, f_end, bs, bst, w_note = (
                            _start_end_extremity(
                                tmp_prot, strand_f, a, b, seq_len,
                                mode_translate, start_orf, end_orf, g_note,
                            )
                        )
                        if len(parts) > 3:
                            # Warning from Glimmer overlap
                            if not warned_header:
                                warn_out.write(f">{uniq_name}\n")
                                warned_header = True
                            warn_out.write(
                                f"{feat_type}-{a}-{b}\n"
                                f"{qualifiers['note']}\n"
                            )
                        if tmp_prot:
                            tp = tmp_prot.rstrip("*")
                            qualifiers["translation"] = tp

                    elif label.startswith("trna:"):
                        feat_type = "tRNA"
                        mode_translate = "no"
                        parts = label.split(":")
                        pos_str = parts[1]
                        t_s, t_e = map(int, pos_str.split("-"))
                        if t_s > t_e:
                            strand_f = -1
                            a, b = t_e, t_s
                        else:
                            strand_f = 1
                            a, b = t_s, t_e
                        qualifiers["note"] = "tRNA found by tRNAscan-SE"
                        nt = _trunc(seq_str, a, b)
                        f_start, f_end, bs, bst, w_note = a, b, False, False, ""

                    elif label.startswith("rrna:"):
                        feat_type = "rRNA"
                        mode_translate = "no"
                        parts = label.split(":")
                        t_s, t_e = map(int, parts[1].split("-"))
                        r_entry = (rrna.get(display_id, {})
                                   .get(t_s, {}).get(t_e, {}))
                        if isinstance(r_entry, dict):
                            strand_f = r_entry.get("c", 1)
                        else:
                            strand_f = 1
                        if t_s > t_e:
                            a, b = t_e, t_s
                        else:
                            a, b = t_s, t_e
                        qualifiers["note"] = "rRNA unit"
                        nt = _trunc(seq_str, a, b)
                        f_start, f_end, bs, bst, w_note = a, b, False, False, ""

                    else:
                        # Homology-based CDS
                        orf_id = label
                        cds = annotation[display_id][orf_id]["CDS"]
                        parts = orf_id.split("_")
                        cadre = parts[2] if len(parts) > 2 else "+1"
                        strand_f = -1 if cadre.startswith("-") else 1

                        # Recover ORF boundaries
                        orf_parts_list = orf_id.split(",")
                        for op in orf_parts_list:
                            pp = op.split("_")
                            pp2 = pp[3].split("-")
                            p_s, p_e = int(pp2[0]), int(pp2[1])
                            if start_orf == 0 or p_s < start_orf:
                                start_orf = p_s
                            if end_orf == 0 or p_e > end_orf:
                                end_orf = p_e
                        if cadre.startswith("-"):
                            start_orf = seq_len - start_orf + 1
                            end_orf = seq_len - end_orf + 1

                        frame_intron = cds.get("frame_or_intron")
                        cds_note = cds["note"]

                        if frame_intron is not None:
                            feat_type = "gene"
                            mode_translate = (
                                "frameshift_stop"
                                if frame_intron == 1
                                else "frameshift_nostop"
                            )
                        elif ("Alignment without M and there is no M at all "
                              "upstream of the alignment") in cds_note:
                            feat_type = "gene"
                            cds["note"] += (
                                " Possible frameshift or intron or "
                                "N-terminal truncated protein."
                            )
                            mode_translate = "extremity"
                        else:
                            feat_type = "CDS"
                            mode_translate = "all"

                        if strand_f == -1:
                            a = cds["start"] - 3
                            if a <= 0:
                                a += 3
                            b = cds["end"]
                            nt = _trunc(seq_str, cds["start"], cds["end"])
                            nt_rev = str(
                                Seq(nt).reverse_complement()
                            )
                            nt = nt_rev
                            if mode_translate in ("all", "extremity"):
                                tmp_prot = _translate_rev(seq_str, a, b)
                            else:
                                tmp_prot = ""
                        else:
                            a = cds["start"]
                            b = cds["end"] + 3
                            if b > seq_len:
                                b -= 3
                            nt = _trunc(seq_str, cds["start"], cds["end"])
                            if mode_translate in ("all", "extremity"):
                                tmp_prot = _translate_fwd(seq_str, a, b)
                            else:
                                tmp_prot = ""

                        f_start, f_end, bs, bst, w_note = (
                            _start_end_extremity(
                                tmp_prot, strand_f, a, b, seq_len,
                                mode_translate, start_orf, end_orf, cds_note,
                            )
                        )

                        if "hit" in cds:
                            qualifiers["hit"] = cds["hit"]
                            qualifiers["evalue"] = str(cds["evalue"])
                        qualifiers["note"] = cds["note"]
                        if mode_translate == "all" and tmp_prot:
                            tp = tmp_prot.rstrip("*")
                            qualifiers["translation"] = tp
                        if cds["warning"] == 1:
                            if not warned_header:
                                warn_out.write(f">{uniq_name}\n")
                                warned_header = True
                            warn_out.write(
                                f"{feat_type}-{a}-{b}\n{cds_note}\n"
                            )

                    if w_note:
                        qualifiers.setdefault("note", "")
                        qualifiers["note"] += (
                            (" " if qualifiers["note"] else "") + w_note
                        )
                        if not warned_header:
                            warn_out.write(f">{uniq_name}\n")
                            warned_header = True
                        warn_out.write(
                            f"{feat_type}-{a}-{b}\n{w_note}\n"
                        )

                    # Write CDS nucleotide sequence
                    if nt is None:
                        nt = ""
                    if not nt:
                        if strand_f == 1:
                            nt = _trunc(seq_str, a, b)
                        else:
                            nt = str(
                                Seq(_trunc(seq_str, a, b)).reverse_complement()
                            )

                    cds_coord = f"{a}-{b}"
                    cds_nt_out.write(f">{cds_coord}\n{nt}\n")

                    if mode_translate == "all":
                        if strand_f == 1:
                            prot = _translate_fwd(seq_str, a, b)
                        else:
                            prot = _translate_rev(seq_str, a, b)
                        prot = prot.rstrip("*")
                        if "*" not in prot:
                            list_file_cds_aa[cds_coord] = prot

                    # ---- Build Biopython SeqFeature -----------------------
                    gene_counter += 1
                    f_start_int = (
                        int(str(f_start).lstrip("<>"))
                        if str(f_start).lstrip("<>").isdigit()
                        else a
                    )
                    f_end_int = (
                        int(str(f_end).lstrip("<>"))
                        if str(f_end).lstrip("<>").isdigit()
                        else b
                    )

                    # Convert to 0-based Biopython coords
                    bp_start = f_start_int - 1
                    bp_end = f_end_int

                    if str(f_start).startswith("<"):
                        from Bio.SeqFeature import BeforePosition
                        bp_start_obj = BeforePosition(0)
                    else:
                        bp_start_obj = ExactPosition(bp_start)

                    if str(f_end).startswith(">"):
                        bp_end_obj = AfterPosition(seq_len)
                    else:
                        bp_end_obj = ExactPosition(bp_end)

                    location = FeatureLocation(
                        bp_start_obj, bp_end_obj,
                        strand=strand_f,
                    )
                    bpq = {k: [v] for k, v in qualifiers.items() if v}
                    feat = SeqFeature(location, type=feat_type,
                                     qualifiers=bpq)
                    out_rec.features.append(feat)

            # Write GenBank and EMBL
            with open(gb_path, "w") as gbh:
                SeqIO.write(out_rec, gbh, "genbank")
            with open(embl_path, "w") as emblh:
                SeqIO.write(out_rec, emblh, "embl")

    return list_file_cds_aa


# ---------------------------------------------------------------------------
# Final annotated output (post-TE + InterproScan)
# ---------------------------------------------------------------------------

def write_final_annotations(
    dir_gb: str,
    out_dir: str,
    id_2_te: dict,
    id_2_interproscan: dict,
) -> None:
    """Add TE and InterproScan annotations to the intermediate GB files and
    write final GenBank, EMBL, and GFF3 results.

    Parameters
    ----------
    dir_gb:
        Path to ``Res_gb`` (intermediate GenBank files).
    out_dir:
        Root output directory.
    id_2_te:
        TE dict from :mod:`microannot.te`.
    id_2_interproscan:
        InterproScan dict from :mod:`microannot.interproscan`.
    """
    dir_gb_annot = os.path.join(out_dir, "Res_gb_annot")
    dir_embl_annot = os.path.join(out_dir, "Res_embl_annot")
    dir_gff_annot = os.path.join(out_dir, "Res_gff_annot")
    for d in (dir_gb_annot, dir_embl_annot, dir_gff_annot):
        os.makedirs(d, exist_ok=True)

    for gb_file in sorted(os.listdir(dir_gb)):
        if not gb_file.endswith(".gb"):
            continue
        stem = gb_file[:-3]
        in_path = os.path.join(dir_gb, gb_file)
        out_gb = os.path.join(dir_gb_annot, gb_file)
        out_embl = os.path.join(dir_embl_annot, f"{stem}.embl")
        out_gff = os.path.join(dir_gff_annot, f"{stem}.gff")

        new_records = []
        with open(out_gff, "w") as gff_out:
            for seq_obj in SeqIO.parse(in_path, "genbank"):
                new_rec = copy.deepcopy(seq_obj)
                new_rec.features = []

                gff_out.write(f"#{seq_obj.id}\n")
                for feat in seq_obj.features:
                    a_start = feat.location.start + 1   # 1-based
                    a_end = int(feat.location.end)       # 1-based inclusive
                    coord_key = f"{a_start}-{a_end}"

                    if coord_key in id_2_te:
                        # Mark as mobile element
                        feat = SeqFeature(
                            feat.location,
                            type="mobile_element",
                            qualifiers={"note": ["Putative transposable element"]},
                        )
                    elif coord_key in id_2_interproscan:
                        ipr = id_2_interproscan[coord_key]
                        for qualifier, values in sorted(ipr.items()):
                            for val in sorted(values):
                                feat.qualifiers.setdefault(
                                    qualifier, []
                                ).append(val)

                    # GFF3 line
                    strand_ch = (
                        "+" if feat.location.strand >= 0 else "-"
                    )
                    attr = {}
                    if "note" in feat.qualifiers:
                        attr["Note"] = (
                            feat.qualifiers["note"][0]
                            .replace(";", "%3B")
                            .replace("=", "%3D")
                        )
                    if "hit" in feat.qualifiers:
                        attr["hit"] = feat.qualifiers["hit"][0]
                    if "evalue" in feat.qualifiers:
                        attr["evalue"] = feat.qualifiers["evalue"][0]

                    attr_str = ";".join(
                        f"{k}={v}" for k, v in attr.items() if v
                    )
                    gff_out.write(
                        f"{seq_obj.id}\tMicroAnnot\t{feat.type}\t"
                        f"{a_start}\t{a_end}\t.\t{strand_ch}\t.\t"
                        f"{attr_str}\n"
                    )
                    new_rec.features.append(feat)
                new_records.append(new_rec)

        with open(out_gb, "w") as gbh:
            SeqIO.write(new_records, gbh, "genbank")
        with open(out_embl, "w") as emblh:
            SeqIO.write(new_records, emblh, "embl")

    # Archive results
    for dirname in ("Res_gb_annot", "Res_embl_annot",
                    "Res_gff_annot", "Res_warnings"):
        archive = os.path.join(out_dir, f"{dirname}.tar.gz")
        src = os.path.join(out_dir, dirname)
        if os.path.isdir(src):
            with tarfile.open(archive, "w:gz") as tar:
                tar.add(src, arcname=dirname)
