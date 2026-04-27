"""
Detection of microsporidia translation-start signals in the 20 nt upstream of a
candidate Met codon.

The original MicroAnnot criteria (in decreasing priority):
  Strong signals  – AT-content >= 80 %, or CCC / GGG motif in the upstream window
  Weak signals    – CC or GG motif (only evaluated when no strong signal was found)
"""


def check_signal(upstream_20nt: str) -> tuple:
    """Analyse up to 20 nt immediately upstream of a candidate Met.

    Parameters
    ----------
    upstream_20nt:
        Nucleotide string in 5'→3' orientation of the ORF (may be shorter than
        20 nt when the candidate Met is near the sequence edge).

    Returns
    -------
    bool_signal : bool
        ``True`` if any signal was found.
    note_signal : str
        Human-readable description of the detected signal(s).
    bool_strong : bool
        ``True`` if at least one *strong* signal was found.
    """
    seq = upstream_20nt.upper()
    note_signal = ""
    bool_signal = False
    bool_strong = False

    if not seq:
        return bool_signal, note_signal, bool_strong

    # AT-content
    at_count = seq.count("A") + seq.count("T")
    perc_AT = at_count / len(seq)

    if perc_AT >= 0.8:
        note_signal += "(percent of AT >= 80% into the 20nt before the M)"
        bool_signal = True
        bool_strong = True

    if "CCC" in seq:
        note_signal += "(CCC into the 20nt before the M)"
        bool_signal = True
        bool_strong = True

    if "GGG" in seq:
        note_signal += "(GGG into the 20nt before the M)"
        bool_signal = True
        bool_strong = True

    if bool_strong:
        note_signal = "-> strong signal " + note_signal
    else:
        if "CC" in seq:
            note_signal += "(CC into the 20nt before the M)"
            bool_signal = True
        if "GG" in seq:
            note_signal += "(GG into the 20nt before the M)"
            bool_signal = True
        if bool_signal:
            note_signal = "-> weak signal " + note_signal

    return bool_signal, note_signal, bool_strong
