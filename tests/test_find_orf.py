"""
Unit tests for ``pymicrosporidiaannot``.find_orf.find_orf.
"""

import pytest
from pymicrosporidiaannot.find_orf import find_orf


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_seq(length: int, frame: int = 0, strand: int = 1,
              stop_before: bool = True, stop_after: bool = True) -> str:
    """Build a synthetic nucleotide sequence for testing.

    The central region (positions 31..60, 1-based) acts as the 'HSP'.
    We embed stop codons before and after the HSP (in the given frame).
    """
    seq = ["A"] * length
    # Convert to mutable
    seq = list("A" * length)

    # Place an in-frame stop codon at position 16 (1-based, frame 0 → 0-based 15)
    if stop_before and strand == 1:
        stop_pos = 15 + frame    # 0-based
        if stop_pos + 2 < length:
            for i, c in enumerate("TAA"):
                seq[stop_pos + i] = c

    # Place an in-frame stop codon at position 64 (1-based, frame 0 → 0-based 63)
    if stop_after and strand == 1:
        stop_pos = 63 + frame    # 0-based
        if stop_pos + 2 < length:
            for i, c in enumerate("TAA"):
                seq[stop_pos + i] = c

    return "".join(seq)


# ---------------------------------------------------------------------------
# Forward strand
# ---------------------------------------------------------------------------

class TestFindOrfForward:
    SEQ_LEN = 100

    def test_returns_tuple_of_two(self):
        seq = "A" * self.SEQ_LEN
        result = find_orf(seq, 1, 0, 30, 60)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_orf_start_leq_start_query(self):
        seq = "A" * self.SEQ_LEN
        orf_start, orf_stop = find_orf(seq, 1, 0, 30, 60)
        assert orf_start <= 30

    def test_orf_stop_geq_end_query(self):
        seq = "A" * self.SEQ_LEN
        orf_start, orf_stop = find_orf(seq, 1, 0, 30, 60)
        assert orf_stop >= 60

    def test_stops_at_upstream_stop_codon(self):
        # Insert stop codon at position 16 in frame 0
        seq = list("A" * self.SEQ_LEN)
        seq[15] = "T"; seq[16] = "A"; seq[17] = "A"  # TAA at pos 16-18 (1-based)
        seq = "".join(seq)
        orf_start, _ = find_orf(seq, 1, 0, 40, 60)
        # ORF should start after the stop codon (at position 19)
        assert orf_start == 19

    def test_stops_at_downstream_stop_codon(self):
        # Insert stop codon at position 64 in frame 0
        seq = list("A" * self.SEQ_LEN)
        seq[63] = "T"; seq[64] = "A"; seq[65] = "A"  # TAA at pos 64-66 (1-based)
        seq = "".join(seq)
        _, orf_stop = find_orf(seq, 1, 0, 30, 60)
        # orf_stop should be at the end of that stop codon (position 66)
        assert orf_stop == 66

    def test_fallback_when_no_upstream_stop(self):
        seq = "A" * self.SEQ_LEN          # no stop codons at all
        orf_start, _ = find_orf(seq, 1, 0, 40, 60)
        # Fallback: first nt of the reading frame
        assert orf_start == 1             # frame 0 → start at 1

    def test_fallback_frame1(self):
        seq = "A" * self.SEQ_LEN
        orf_start, _ = find_orf(seq, 1, 1, 40, 60)
        assert orf_start == 2             # frame 1 → start at 2

    def test_coordinates_within_bounds(self):
        seq = "A" * self.SEQ_LEN
        orf_start, orf_stop = find_orf(seq, 1, 0, 30, 60)
        assert 1 <= orf_start <= self.SEQ_LEN
        assert 1 <= orf_stop <= self.SEQ_LEN


# ---------------------------------------------------------------------------
# Reverse strand
# ---------------------------------------------------------------------------

class TestFindOrfReverse:
    SEQ_LEN = 100

    def test_returns_tuple_of_two(self):
        seq = "A" * self.SEQ_LEN
        result = find_orf(seq, -1, 0, 40, 70)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_orf_start_geq_end_query(self):
        """For -strand, orf_start (higher coord) should be ≥ end_query."""
        seq = "A" * self.SEQ_LEN
        orf_start, orf_stop = find_orf(seq, -1, 0, 40, 70)
        assert orf_start >= 70

    def test_orf_stop_leq_start_query(self):
        """For -strand, orf_stop (lower coord) should be ≤ start_query."""
        seq = "A" * self.SEQ_LEN
        orf_start, orf_stop = find_orf(seq, -1, 0, 40, 70)
        assert orf_stop <= 40

    def test_coordinates_within_bounds(self):
        seq = "A" * self.SEQ_LEN
        orf_start, orf_stop = find_orf(seq, -1, 0, 40, 70)
        assert 0 <= orf_stop <= self.SEQ_LEN
        assert 0 <= orf_start <= self.SEQ_LEN
