"""
BLAST-based CDS annotation.

Contains three main operations (matching the original Perl functions):

annot_alignment(qresult, hit, hsp, annotation, hit_2_frameshift, debug)
    Parse a blastp HSP against the ORF amino-acid sequences and annotate the
    CDS start/end coordinates.

annot_alignment_small(qresult, hit, hsp, annotation, small_cds,
                      orf_hit_small, debug)
    Parse a blastx HSP (small CDS database, proteins < 80 aa) and annotate
    the CDS coordinates directly on the nucleotide sequence.

group_frameshifts(annotation, hit_2_frameshift)
    Merge adjacent ORFs that share BLAST hits (indicating frameshifts or
    introns) into a single compound feature.

Coordinate convention: 1-based inclusive, matching BioPerl.
"""

import re
from Bio.Seq import Seq

from pymicrosporidiaannot.signal import check_signal
from pymicrosporidiaannot.find_orf import find_orf


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _trunc_str(seq_str: str, a: int, b: int) -> str:
    """BioPerl-style trunc: return seq_str[a-1:b] (1-based inclusive)."""
    if b < a:
        return ""
    return seq_str[a - 1: b]


def _translate_fwd(seq_str: str, a: int, b: int) -> str:
    s = _trunc_str(seq_str, a, b)
    if not s:
        return ""
    return str(Seq(s).translate())


def _translate_rev(seq_str: str, a: int, b: int) -> str:
    s = _trunc_str(seq_str, a, b)
    if not s:
        return ""
    return str(Seq(s).reverse_complement().translate())


# ---------------------------------------------------------------------------
# annot_alignment  (blastp – ORF amino-acid query)
# ---------------------------------------------------------------------------

def annot_alignment(qresult, hit, hsp, annotation, hit_2_frameshift,
                    bool_hit: bool, debug: bool = False):
    """Annotate a CDS from one blastp HSP.

    Modifies *annotation* and *hit_2_frameshift* in place.

    Parameters
    ----------
    qresult:
        ``Bio.Blast.Record.Blast`` record (NCBIXML) for the query.
    hit:
        ``Bio.Blast.Record.Alignment`` alignment object.
    hsp:
        ``Bio.Blast.Record.HSP`` object.
    annotation:
        Master annotation dict (keyed by seq_id → orf_id).
    hit_2_frameshift:
        Dict used to detect potential frameshifts between ORFs.
    bool_hit:
        ``True`` when this is not the best (first) hit for this query,
        so only the frameshift register is updated.
    debug:
        Emit extra diagnostic notes when ``True``.
    """
    # --- Parse the ORF ID to recover the contig id, frame, start, end ------
    query_name = qresult.query.split()[0]   # e.g. "1_orf_+1_100-300_extr_5prim"
    parts = query_name.split("_")
    # parts = [display_id, "orf", cadre, "start-end", ...]
    display_id = parts[0]
    cadre = parts[2]                        # e.g. "+1" or "-2"
    pos_parts = parts[3].split("-")
    orf_start_str = pos_parts[0]
    orf_end_str = pos_parts[1]

    # Register in hit_2_frameshift
    if cadre.startswith("-"):
        hit_2_frameshift.setdefault(display_id, {}).setdefault("-", {})
        hit_2_frameshift[display_id]["-"].setdefault(orf_end_str, {})
        hit_2_frameshift[display_id]["-"][orf_end_str]["orf"] = query_name
        hit_2_frameshift[display_id]["-"][orf_end_str].setdefault(
            "hit", {})[hit.hit_id] = 1
    else:
        hit_2_frameshift.setdefault(display_id, {}).setdefault("+", {})
        hit_2_frameshift[display_id]["+"].setdefault(orf_start_str, {})
        hit_2_frameshift[display_id]["+"][orf_start_str]["orf"] = query_name
        hit_2_frameshift[display_id]["+"][orf_start_str].setdefault(
            "hit", {})[hit.hit_id] = 1

    if bool_hit:
        return

    # --- Determine anno_start in the ORF AA frame -------------------------
    # All BioPerl NCBI BLAST XML uses 1-based coordinates:
    # hsp.query_start / hsp.sbjct_start are 1-based inclusive
    hsp_query_start = hsp.query_start   # 1-based AA
    hsp_sbjct_start = hsp.sbjct_start   # 1-based AA
    hit_seq = hsp.sbjct                 # aligned hit (subject) string
    query_seq = hsp.query               # aligned query string

    seq_aa = annotation[display_id][query_name]["seqAA"]
    seq_nt_add = annotation[display_id][query_name]["seqNT-add"]
    seq_nt_add_n = annotation[display_id][query_name]["seqNT-add-n"]

    anno_start = ""
    note = ""
    warning = 0
    amont_M_20nt = ""
    amont10 = ""

    if debug:
        print(query_name)

    # Case 1: alignment starts at M in both hit and query, at N-terminus
    if (hit_seq and hit_seq[0] == "M"
            and query_seq and query_seq[0] == "M"
            and hsp_sbjct_start == 1):
        anno_start = hsp_query_start - 1   # 0-based within the AA seq
        note = "Alignment starting with M"
        warning = 0

    # Case 2: 5′-extremity ORF – start is undetermined
    elif (query_name.endswith("_extr_5prim")
          and ((cadre.startswith("-") and "_5prim" in query_name)
               or (not cadre.startswith("-") and "_5prim" in query_name))):
        note = "Start undetermined"
        anno_start = 0

    # Case 3: search for an M in the ±5 AA window around the alignment start
    else:
        substr_a = (hsp_query_start - hsp_sbjct_start) - 5
        substr_b = 10
        if substr_a < 0:
            substr_b += substr_a
            substr_a = 0
        if (substr_a + substr_b) >= (hsp_query_start + 5):
            substr_b = hsp_query_start - 1

        if substr_b > 0:
            amont10 = seq_aa[substr_a: substr_a + substr_b]
        else:
            amont10 = ""

        m10 = re.search(r"(.*)M", amont10)
        if m10:
            anno_start = substr_a + len(m10.group(1))
            note = ("Alignment with a M into the 5AA before or after "
                    "the start of the hit")
            warning = 0
        else:
            # Search for an upstream M with a start signal
            offset = -3
            amont = seq_aa[: hsp_query_start - 1]
            found_M = False
            first_M_pos = None

            for mobj in re.finditer(r"(.*?)M", amont):
                amont_M = mobj.group(1)
                offset += 3
                offset += len(amont_M) * 3

                # 20 nt upstream in the NT sequence
                nt_pos = offset + seq_nt_add_n - 1   # 0-based in seq_nt_add
                upstream_start = max(0, nt_pos - 20)
                amont_M_20nt = seq_nt_add[upstream_start: nt_pos]

                bool_signal, note_signal, _ = check_signal(amont_M_20nt)

                if bool_signal:
                    anno_start = offset // 3
                    note = ("Alignment without M but in upstream there is a M "
                            "with a signal " + note_signal)
                    warning = 0
                    found_M = True
                    break

                if first_M_pos is None:
                    first_M_pos = offset // 3

                offset += 3   # skip past the M itself

            if not found_M:
                if first_M_pos is not None:
                    anno_start = first_M_pos
                    note = ("Alignment without M but in upstream there is a M "
                            "without a signal.")
                    warning = 1
                else:
                    anno_start = hsp_query_start - 1
                    note = ("Alignment without M and there is no M at all "
                            "upstream of the alignment.")
                    warning = 1

    # Debug extra notes
    if debug:
        if amont10:
            note += (f"    The 5AA +/- start(query - hit):{amont10}.")
        if amont_M_20nt:
            note += f"    The 20nt upstream:{amont_M_20nt}."

    # --- Convert AA offset → nucleotide coordinates on the contig ----------
    anno_start_nt = anno_start * 3
    split2 = [int(orf_start_str), int(orf_end_str)]

    if cadre.startswith("-"):
        bioseq_len = len(annotation[display_id]["bioseq"].seq)
        anno_end_out = bioseq_len - split2[0] - anno_start_nt + 1
        anno_start_out = bioseq_len - split2[1] + 1
    else:
        anno_start_out = anno_start_nt + split2[0]
        anno_end_out = split2[1]

    annotation[display_id][query_name]["CDS"] = {
        "start":   anno_start_out,
        "end":     anno_end_out,
        "note":    note,
        "warning": warning,
        "hit":     hit.hit_id,
        "evalue":  hsp.expect,
    }


# ---------------------------------------------------------------------------
# annot_alignment_small  (blastx – nucleotide query vs protein DB)
# ---------------------------------------------------------------------------

def annot_alignment_small(qresult, hit, hsp, annotation, small_cds,
                          orf_hit_small, debug: bool = False):
    """Annotate a small CDS (<80 aa) from one blastx HSP.

    Modifies *small_cds* and *orf_hit_small* in place.
    """
    display_id = qresult.query.split()[0]
    bioseq_str = str(annotation[display_id]["bioseq"].seq).upper()
    bioseq_len = len(bioseq_str)

    # Biopython NCBIXML: coordinates are 1-based inclusive
    hsp_strand = 1 if hsp.frame[0] > 0 else -1
    frame = abs(hsp.frame[0]) - 1          # BioPerl-style frame (0, 1, or 2)
    hsp_q_start = hsp.query_start          # 1-based nt on + strand
    hsp_q_end = hsp.query_end              # 1-based nt on + strand
    hsp_h_start = hsp.sbjct_start          # 1-based aa (subject)
    hit_seq = hsp.sbjct                    # aligned hit string
    query_seq = hsp.query                  # translated query string
    hit_length = hit.length                # total subject protein length

    # --- Find the enclosing ORF boundaries --------------------------------
    orf_start, orf_stop = find_orf(
        bioseq_str, hsp_strand, frame, hsp_q_start, hsp_q_end
    )

    # Deduplication: skip if we already processed an HSP from this ORF
    dedup_key = f"{orf_stop}-{hsp_strand}"
    if dedup_key in orf_hit_small:
        return

    anno_start = ""
    anno_end = ""
    note = ""
    warning = 0
    amont10 = ""
    amont_M_20nt = ""
    bool_small = False

    # Case 1: alignment starts with M at the N-terminus of the subject ------
    if (hit_seq and hit_seq[0] == "M"
            and query_seq and query_seq[0] == "M"
            and hsp_h_start == 1):
        if hsp_strand == 1:
            anno_start = hsp_q_start
        else:
            anno_end = hsp_q_end
        bool_small = True
        note = "Alignment starting with M."
        warning = 0

    else:
        # ---- Look for M in the ±5 AA upstream window ---------------------
        if hsp_strand == 1:
            substr_a = (hsp_q_start - hsp_h_start * 3) + 3 - (5 * 3)
            substr_b = 10 * 3
            if substr_a <= 0:
                substr_b = substr_b + substr_a - 3 - frame - 1
                substr_a = frame + 1
            if substr_b > 0:
                a = substr_a
                b = (substr_a + substr_b) - 1
                amont10 = _translate_fwd(bioseq_str, a, b)
        else:
            substr_a = (hsp_q_end + hsp_h_start * 3) - 3 + (5 * 3)
            substr_b = 10 * 3
            if substr_a > bioseq_len:
                substr_b = substr_b - (substr_a + 3 - bioseq_len + frame)
                substr_a = bioseq_len - frame
            if substr_b > 0:
                a = substr_a - substr_b + 1
                b = substr_a
                amont10 = _translate_rev(bioseq_str, a, b)

        # Strip everything up to and including the last stop codon
        if amont10:
            m_stop = re.search(r"(.*\*)(.*)", amont10)
            if m_stop:
                prefix_len = len(m_stop.group(1))
                amont10 = m_stop.group(2)
                if hsp_strand == 1:
                    substr_a += prefix_len * 3
                else:
                    substr_a -= prefix_len * 3

        if (amont10
                and ((hsp_strand == 1 and orf_start <= substr_a)
                     or (hsp_strand == -1 and orf_start >= substr_a))):
            m_M = re.search(r"(.*)M", amont10)
            if m_M:
                if hsp_strand == 1:
                    anno_start = substr_a + len(m_M.group(1)) * 3
                else:
                    anno_end = substr_a - len(m_M.group(1)) * 3
                note = ("Alignment with a M into the 5AA before or after "
                        "the start of the hit.")
                warning = 0
            else:
                amont10 = ""  # not found, fall through

        if not amont10 or (not anno_start and not anno_end):
            # ---- Search the full ORF upstream region for M ----------------
            pos_small = 0
            if hsp_strand == 1:
                a = frame + 1
                b = hsp_q_start - 1
                amont_full = ""
                if a < b:
                    amont_full = _translate_fwd(bioseq_str, a, b)
                m_stop2 = re.search(r"(.*\*)(.*)", amont_full)
                if m_stop2:
                    pos_small = (len(m_stop2.group(1)) * 3) + frame + 1
                    amont_full = m_stop2.group(2)
                else:
                    pos_small = frame + 1
            else:
                a = hsp_q_end + 1
                b = bioseq_len - frame
                amont_full = ""
                if a < b:
                    amont_full = _translate_rev(bioseq_str, a, b)
                m_stop2 = re.search(r"(.*\*)(.*)", amont_full)
                if m_stop2:
                    pos_small = bioseq_len - (len(m_stop2.group(1)) * 3) - frame
                    amont_full = m_stop2.group(2)
                else:
                    pos_small = bioseq_len - frame

            # Scan for M with upstream signal
            start_off = 0
            bool_signal_found = False
            first_M = None

            for mobj in re.finditer(r"(.*?)M", amont_full):
                amont_M = mobj.group(1)
                if amont_M:   # skip zero-length prefix (i.e. M at very start)
                    start_off += len(amont_M) * 3
                    add = 20
                    if hsp_strand == 1:
                        pos = pos_small + start_off
                        if pos - add < 1:
                            add = pos - 1
                        if add > 0:
                            amont_M_20nt = _trunc_str(
                                bioseq_str, pos - add, pos - 1
                            )
                        else:
                            amont_M_20nt = ""
                    else:
                        pos = pos_small - start_off
                        if pos + add > bioseq_len:
                            add = bioseq_len - pos
                        if add > 0:
                            amont_M_20nt = str(
                                Seq(_trunc_str(
                                    bioseq_str, pos + 1, pos + add
                                )).reverse_complement()
                            )
                        else:
                            amont_M_20nt = ""

                    bool_sig, note_sig, _ = check_signal(amont_M_20nt)

                    if bool_sig:
                        if hsp_strand == 1:
                            anno_start = pos_small + start_off
                        else:
                            anno_end = pos_small - start_off
                        note = ("Alignment without M but in upstream there is "
                                "a M with a signal. " + note_sig)
                        warning = 0
                        bool_signal_found = True
                        break

                    if first_M is None:
                        if hsp_strand == 1:
                            first_M = pos_small + start_off
                        else:
                            first_M = pos_small - start_off

                start_off += 3   # skip past the M codon

            if not bool_signal_found:
                if first_M is not None:
                    if hsp_strand == 1:
                        anno_start = first_M
                    else:
                        anno_end = first_M
                    note = ("Alignment without M but in upstream there is a M "
                            "without a signal.")
                    warning = 1
                else:
                    if hsp_strand == 1:
                        anno_start = hsp_q_start
                    else:
                        anno_end = hsp_q_end
                    note = ("Alignment without M and there is no M at all "
                            "upstream of the alignment.")
                    warning = 1

    # --- Adjust for stop codon inclusion ----------------------------------
    if hsp_strand == -1:
        anno_start = (anno_start if anno_start else 0) + 3
    else:
        anno_end = (anno_end if anno_end else 0) - 3

    # --- Verify size is within 80–120 % of the subject protein length ------
    end_val = anno_end if anno_end else 0
    start_val = anno_start if anno_start else 0

    seq_obj = bioseq_str
    a_chk = start_val
    b_chk = end_val
    if hsp_strand == -1:
        a_chk -= 3
        if a_chk <= 0:
            a_chk += 3
    else:
        b_chk += 3
        if b_chk > bioseq_len:
            b_chk -= 3

    if hsp_strand == -1:
        tmp_prot = _translate_rev(seq_obj, a_chk, b_chk)
    else:
        tmp_prot = _translate_fwd(seq_obj, a_chk, b_chk)

    bool_stop_prot = tmp_prot.endswith("*")
    bool_start_prot = tmp_prot.startswith("M") and warning == 0

    # Check if at least one extremity is truncated
    bool_check_extremity = True
    if hsp_strand == -1:
        if not bool_stop_prot or (orf_start >= (bioseq_len - 3)
                                   and not bool_start_prot):
            bool_check_extremity = False
    else:
        if not bool_stop_prot or (orf_start <= 3
                                   and not bool_start_prot):
            bool_check_extremity = False

    cds_len_aa = (end_val - start_val + 1) / 3

    if (
        (not bool_check_extremity and cds_len_aa <= 1.20 * hit_length)
        or (0.80 * hit_length <= cds_len_aa <= 1.20 * hit_length)
    ):
        entry = {
            "end":     end_val,
            "note":    note,
            "warning": warning,
            "hit":     hit.hit_id,
            "evalue":  hsp.expect,
            "strand":  hsp_strand,
            "orf":     f"{orf_start}-{orf_stop}",
        }
        if debug:
            entry["note"] += (
                f"      length : hit80%<CDS>hit80%)"
                f":{hit_length}  "
                f"{0.80 * hit_length:.1f}<={cds_len_aa:.1f}"
                f"<={1.20 * hit_length:.1f}."
            )
            if amont10:
                entry["note"] += (
                    f"    The 5AA +/- start(query - hit):{amont10}."
                )
            if amont_M_20nt:
                entry["note"] += f"    The 20nt upstream:{amont_M_20nt}."
            entry["note"] += f"    orf_start:{orf_start}."
            entry["note"] += f"    orf_stop:{orf_stop}."

        small_cds.setdefault(display_id, {})[start_val] = entry
        orf_hit_small[dedup_key] = 1


# ---------------------------------------------------------------------------
# group_frameshifts
# ---------------------------------------------------------------------------

def group_frameshifts(annotation: dict, hit_2_frameshift: dict):
    """Merge ORFs that share BLAST hits and are close in genomic space.

    These likely represent frameshifted genes or genes interrupted by an
    intron. The merged entry replaces the individual ORF entries in
    *annotation*.
    """
    for display_id in sorted(hit_2_frameshift):
        already_done = {}
        for strand in sorted(hit_2_frameshift[display_id]):
            for order in sorted(
                hit_2_frameshift[display_id][strand],
                key=lambda x: int(x) if x.lstrip("-").isdigit() else 0,
            ):
                orf = hit_2_frameshift[display_id][strand][order]["orf"]
                if orf in already_done:
                    continue
                if orf not in annotation.get(display_id, {}):
                    continue
                orf_cds = annotation[display_id][orf].get("CDS")
                if not orf_cds:
                    continue

                already_done[orf] = 1
                list_hit = dict(
                    hit_2_frameshift[display_id][strand][order]["hit"]
                )

                parts = orf.split("_")
                pos_parts = parts[3].split("-")
                s = int(pos_parts[0])
                e = int(pos_parts[1])
                group = {orf: 1}

                changed = True
                while changed:
                    changed = False
                    for order_b in sorted(
                        hit_2_frameshift[display_id][strand],
                        key=lambda x: int(x) if x.lstrip("-").isdigit() else 0,
                    ):
                        orf_b = hit_2_frameshift[display_id][strand][order_b][
                            "orf"
                        ]
                        if orf_b in already_done:
                            continue
                        cds_b = annotation.get(display_id, {}).get(
                            orf_b, {}
                        ).get("CDS", {})
                        if not cds_b:
                            continue
                        # Skip ORFs whose alignment clearly identified M start
                        note_b = cds_b.get("note", "")
                        if ("Alignment starting with M" in note_b
                                or "Alignment with a M into the 5AA" in note_b):
                            continue
                        if orf_b in group:
                            continue

                        parts_b = orf_b.split("_")
                        pos_b = parts_b[3].split("-")
                        for_s = int(pos_b[0])
                        for_e = int(pos_b[1])

                        for h in hit_2_frameshift[display_id][strand][order_b][
                            "hit"
                        ]:
                            if h in list_hit and for_s <= (e + 50):
                                if s > for_s:
                                    s = for_s
                                    changed = True
                                if e < for_e:
                                    e = for_e
                                    changed = True
                                group[orf_b] = 1
                                for h2 in hit_2_frameshift[display_id][strand][
                                    order_b
                                ]["hit"]:
                                    if h2 not in list_hit:
                                        list_hit[h2] = 1
                                        changed = True
                                break

                if len(group) <= 1:
                    continue

                # Merge group into a compound annotation
                merged_start = 0
                merged_end = 0
                merged_note = ""
                merged_hits = ""
                merged_evalues = ""
                merged_orf_key = ""
                bool_stop = ""
                bioseq_str = str(
                    annotation[display_id]["bioseq"].seq
                ).upper()
                bioseq_len = len(bioseq_str)

                for g_orf in sorted(group):
                    g_cds = annotation[display_id][g_orf]["CDS"]
                    merged_orf_key += g_orf + ","
                    g_parts = g_orf.split("_")
                    g_cadre = g_parts[2]

                    if merged_start == 0 or merged_start > g_cds["start"]:
                        merged_start = g_cds["start"]
                        if not g_cadre.startswith("-"):
                            merged_note = (
                                "Frameshift or intron. " + g_cds["note"]
                            )
                        else:
                            a = g_cds["start"] - 3
                            b = g_cds["end"]
                            if a <= 0:
                                a += 3
                            tmp = str(
                                Seq(bioseq_str[a - 1: b])
                                .reverse_complement()
                                .translate()
                            )
                            bool_stop = 1 if tmp.endswith("*") else 0

                    if merged_end == 0 or merged_end < g_cds["end"]:
                        merged_end = g_cds["end"]
                        if g_cadre.startswith("-"):
                            merged_note = (
                                "Frameshift or intron. " + g_cds["note"]
                            )
                        else:
                            a = g_cds["start"]
                            b = g_cds["end"] + 3
                            if b > bioseq_len:
                                b -= 3
                            tmp = str(
                                Seq(bioseq_str[a - 1: b]).translate()
                            )
                            bool_stop = 1 if tmp.endswith("*") else 0

                    merged_evalues += str(g_cds["evalue"]) + " "
                    merged_hits += g_cds.get("hit", "") + " "
                    del annotation[display_id][g_orf]
                    already_done[g_orf] = 1

                merged_orf_key = merged_orf_key.rstrip(",")
                merged_evalues = merged_evalues.rstrip()
                merged_hits = merged_hits.rstrip()

                annotation[display_id][merged_orf_key] = {
                    "CDS": {
                        "frame_or_intron": bool_stop,
                        "start":   merged_start,
                        "end":     merged_end,
                        "note":    merged_note,
                        "hit":     "multiple: " + merged_hits,
                        "evalue":  "multiple: " + merged_evalues,
                        "warning": 1,
                    }
                }
