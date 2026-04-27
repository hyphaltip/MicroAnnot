"""
Unit tests for microannot.signal.check_signal.
"""

import pytest
from microannot.signal import check_signal


class TestCheckSignalEmpty:
    def test_empty_string(self):
        sig, note, strong = check_signal("")
        assert sig is False
        assert note == ""
        assert strong is False


class TestStrongSignal:
    def test_at_rich_sequence(self):
        # 17/20 = 85 % AT
        sig, note, strong = check_signal("AAATTTAAATTTAAATTTAA")
        assert sig is True
        assert strong is True
        assert "80%" in note
        assert "strong" in note

    def test_ccc_motif(self):
        sig, note, strong = check_signal("GGGCCCGGGCCCGGGCCCGG")
        assert sig is True
        assert strong is True
        assert "CCC" in note

    def test_ggg_motif(self):
        sig, note, strong = check_signal("TTGGGTTGGGTTGGGTTGGG")
        assert sig is True
        assert strong is True
        assert "GGG" in note

    def test_strong_prefix(self):
        sig, note, strong = check_signal("AAAAAAAAAAAAAAAAAACCC")
        assert "strong signal" in note

    def test_multiple_strong_signals(self):
        # High AT AND CCC
        sig, note, strong = check_signal("AAATTTAAATTTCCCAAATT")
        assert sig is True
        assert strong is True
        assert "80%" in note
        assert "CCC" in note


class TestWeakSignal:
    def test_cc_motif(self):
        sig, note, strong = check_signal("GCGCCCGCGCGCGCCGCGCG"
                                         [:-1] + "A")  # ~50% AT, no CCC/GGG
        # This sequence has CCC so it's strong; use a CC-only sequence:
        sig, note, strong = check_signal("GCGCGCGCGCGCGCGCCCGC")
        # contains CCC → strong
        # Use a clean CC-only sequence
        sig, note, strong = check_signal("GCGCGCGCGCGCGCGCGCCG")
        assert sig is True
        assert strong is False
        assert "weak signal" in note
        assert "CC" in note

    def test_gg_motif(self):
        sig, note, strong = check_signal("GCGCGCGCGCGCGCGGCGCG")
        assert sig is True
        assert strong is False
        assert "weak signal" in note
        assert "GG" in note

    def test_no_signal(self):
        # Low AT, no special motif
        sig, note, strong = check_signal("GCGCGCGCGCGCGCGCGCGC")
        assert sig is False
        assert note == ""
        assert strong is False


class TestCaseInsensitive:
    def test_lowercase_input(self):
        sig, note, strong = check_signal("aaatttaaatttaaatttaa")
        assert sig is True
        assert strong is True
