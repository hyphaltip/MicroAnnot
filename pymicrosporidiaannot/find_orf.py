"""
ORF boundary extension: given the nucleotide coordinates of a BLAST HSP, extend
to the nearest upstream and downstream stop codons to recover the full containing
reading frame.

All coordinates follow the 1-based *inclusive* convention used by BioPerl (and
the original Perl code), so::

    trunc(a, b)  ≡  seq_str[a-1 : b]   (Python 0-based slicing)

For the **forward strand** the function returns:
  ``orf_start`` – first nucleotide of the ORF (after the upstream stop codon, or
                  the start of the reading frame if there is none).
  ``orf_stop``  – last nucleotide of the downstream stop codon (or end of the
                  reading frame if no stop codon is found downstream).

For the **reverse strand** (+strand coordinates are used throughout):
  ``orf_start`` – the higher +strand position (far boundary, start-codon side on
                  the − strand).
  ``orf_stop``  – the lower +strand position (first nucleotide of the stop codon
                  on the + strand; corresponds to the 3′-end of the − strand ORF).

The calling code (``annot_alignment_small``, Glimmer refinement) applies an
additional ``+3`` or ``−3`` adjustment to ``orf_stop`` when storing coordinates,
matching the original Perl behaviour where ``$anno_start`` was set to
``$orf_stop`` then incremented by 3.
"""

import re
from Bio.Seq import Seq


def _trunc(seq_str: str, a: int, b: int) -> str:
    """Return the 1-based inclusive subsequence seq_str[a-1:b]."""
    return seq_str[a - 1: b]


def find_orf(seq_str: str, strand: int, frame: int,
             start_query: int, end_query: int) -> tuple:
    """Find the ORF boundaries enclosing a BLAST HSP region.

    Parameters
    ----------
    seq_str:
        Full nucleotide sequence of the contig (uppercase, + strand).
    strand:
        ``1`` for forward, ``-1`` for reverse.
    frame:
        Reading-frame offset (0, 1, or 2) in BioPerl style:

        * forward strand – nt to skip from the 5′ end of the sequence
          before the reading frame begins.
        * reverse strand – nt to skip from the 3′ end of the sequence
          before the reading frame begins (counting from the end).
    start_query:
        1-based start position of the BLAST HSP on the + strand.
    end_query:
        1-based end position of the BLAST HSP on the + strand.

    Returns
    -------
    orf_start : int
        First nt of the ORF (+ strand coordinate, 1-based).
    orf_stop : int
        Last nt of the stop codon (or reading-frame endpoint if none found),
        1-based + strand coordinate.
    """
    seq_len = len(seq_str)

    def translate_fwd(a: int, b: int) -> str:
        if b < a:
            return ""
        return str(Seq(_trunc(seq_str, a, b)).translate())

    def translate_rev(a: int, b: int) -> str:
        if b < a:
            return ""
        return str(Seq(_trunc(seq_str, a, b)).reverse_complement().translate())

    if strand == -1:
        # ---- Reverse strand ----
        # Reading frame on +strand: runs from
        #   ((seq_len - frame) % 3) + 1  ..  seq_len - frame
        # (the revcom is read right-to-left on the +strand)

        # Find the upstream (−strand direction) stop codon:
        # look left of start_query in the +strand, translate as revcom.
        a = ((seq_len - frame) % 3) + 1   # 1-based frame start on + strand
        b = start_query - 1               # 1-based end of upstream region

        if start_query > 3 and b >= a:
            prot = translate_rev(a, b)
            m = re.search(r"(.*?)\*", prot)
            if m:
                neg_anno_start = start_query - (len(m.group(1)) * 3 + 3)
            else:
                neg_anno_start = ((seq_len - frame) % 3) + 1
        else:
            neg_anno_start = ((seq_len - frame) % 3) + 1

        orf_stop = neg_anno_start  # lower +strand coord (stop-codon side)

        # Find the downstream (−strand direction) stop codon:
        # look right of end_query in the +strand, translate as revcom.
        a = end_query + 1       # 1-based first nt after the HSP
        b = seq_len - frame     # 1-based last nt of the reading frame

        if a < b:
            prot = translate_rev(a, b)
            m = re.search(r".*\*(.*)", prot)
            if m:
                orf_start = len(m.group(1)) * 3 + end_query
            else:
                orf_start = seq_len - frame
        else:
            orf_start = seq_len - frame

    else:
        # ---- Forward strand ----

        # Find the downstream stop codon:
        a = end_query + 1                         # 1-based first nt after HSP
        b = seq_len - ((seq_len - frame) % 3)     # last nt of reading frame

        if end_query + 1 < seq_len - 3 and b >= a:
            prot = translate_fwd(a, b)
            m = re.search(r"(.*?)\*", prot)
            if m:
                orf_stop = end_query + (len(m.group(1)) * 3 + 3)
            else:
                orf_stop = seq_len - ((seq_len - frame) % 3)
        else:
            orf_stop = seq_len - ((seq_len - frame) % 3)

        # Find the upstream stop codon:
        a = frame + 1        # 1-based first nt of reading frame
        b = start_query - 1  # 1-based end of upstream region

        if frame + 1 < start_query - 1 and b >= a:
            prot = translate_fwd(a, b)
            m = re.search(r"(.*)\*", prot)
            if m:
                orf_start = len(m.group(1)) * 3 + 3 + frame + 1
            else:
                orf_start = frame + 1
        else:
            orf_start = frame + 1

    return orf_start, orf_stop
