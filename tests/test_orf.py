"""
Unit tests for microannot.orf.extract_orfs.
"""

import os
import tempfile
import pytest
from microannot.orf import extract_orfs


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

def _write_fasta(path: str, sequences: list) -> None:
    """Write a FASTA file.  *sequences* is a list of (header, seq) tuples."""
    with open(path, "w") as fh:
        for header, seq in sequences:
            fh.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# Basic extraction
# ---------------------------------------------------------------------------

class TestExtractOrfs:
    def test_single_forward_orf(self, tmp_path):
        # ATG...stop in frame +1, ≥240 nt
        orf_nt = "ATG" + "AAA" * 79 + "TAA"   # 3 + 237 + 3 = 243 nt
        seq = "NNN" + orf_nt + "NNN"           # pad with rubbish
        fasta = str(tmp_path / "in.fa")
        out_fa = str(tmp_path / "orfs.fa")
        _write_fasta(fasta, [("1", seq)])

        annotation = extract_orfs(fasta, out_fa, min_orf_size=240)
        assert "1" in annotation
        assert "bioseq" in annotation["1"]
        # At least one ORF should have been found
        orf_ids = [k for k in annotation["1"] if k != "bioseq"]
        assert len(orf_ids) >= 1

    def test_short_orf_excluded(self, tmp_path):
        # ORF of only 60 nt – below the 200 nt minimum
        orf_nt = "ATG" + "AAA" * 17 + "TAA"   # 60 nt
        seq = orf_nt
        fasta = str(tmp_path / "in.fa")
        out_fa = str(tmp_path / "orfs.fa")
        _write_fasta(fasta, [("1", seq)])

        annotation = extract_orfs(fasta, out_fa, min_orf_size=240)
        orf_ids = [k for k in annotation["1"] if k != "bioseq"]
        assert len(orf_ids) == 0

    def test_reverse_orf_found(self, tmp_path):
        # Build a reverse-strand ORF that is ≥240 nt
        orf_aa_len = 80                           # 80 aa × 3 = 240 nt
        orf_fwd = "ATG" + "GCC" * orf_aa_len + "TAA"
        from Bio.Seq import Seq
        revcom_orf = str(Seq(orf_fwd).reverse_complement())
        fasta = str(tmp_path / "in.fa")
        out_fa = str(tmp_path / "orfs.fa")
        _write_fasta(fasta, [("1", revcom_orf)])

        annotation = extract_orfs(fasta, out_fa, min_orf_size=240)
        orf_ids = [k for k in annotation["1"] if k != "bioseq"]
        assert len(orf_ids) >= 1
        # At least one must be on the minus strand
        assert any("-" in oid for oid in orf_ids)

    def test_orf_id_format(self, tmp_path):
        orf_nt = "ATG" + "GCC" * 79 + "TAA"   # 240 nt coding
        fasta = str(tmp_path / "in.fa")
        out_fa = str(tmp_path / "orfs.fa")
        _write_fasta(fasta, [("1", orf_nt)])

        annotation = extract_orfs(fasta, out_fa, min_orf_size=240)
        orf_ids = [k for k in annotation["1"] if k != "bioseq"]
        for oid in orf_ids:
            parts = oid.split("_")
            # parts: [seq_id, "orf", cadre, "start-end", ...]
            assert parts[1] == "orf"
            assert parts[2][0] in ("+", "-")
            assert "-" in parts[3]

    def test_seqAA_stored(self, tmp_path):
        orf_nt = "ATG" + "GCC" * 79 + "TAA"
        fasta = str(tmp_path / "in.fa")
        out_fa = str(tmp_path / "orfs.fa")
        _write_fasta(fasta, [("1", orf_nt)])

        annotation = extract_orfs(fasta, out_fa, min_orf_size=240)
        orf_ids = [k for k in annotation["1"] if k != "bioseq"]
        for oid in orf_ids:
            assert "seqAA" in annotation["1"][oid]
            # seqAA should NOT contain stop codon
            assert "*" not in annotation["1"][oid]["seqAA"]

    def test_multiple_sequences(self, tmp_path):
        orf_nt = "ATG" + "GCC" * 79 + "TAA"
        fasta = str(tmp_path / "in.fa")
        out_fa = str(tmp_path / "orfs.fa")
        _write_fasta(fasta, [("1", orf_nt), ("2", orf_nt)])

        annotation = extract_orfs(fasta, out_fa, min_orf_size=240)
        assert "1" in annotation
        assert "2" in annotation

    def test_min_orf_clamp(self, tmp_path):
        """min_orf_size below 200 should be clamped to 200."""
        orf_nt = "ATG" + "AAA" * 66 + "TAA"   # 201 nt
        fasta = str(tmp_path / "in.fa")
        out_fa = str(tmp_path / "orfs.fa")
        _write_fasta(fasta, [("1", orf_nt)])

        annotation = extract_orfs(fasta, out_fa, min_orf_size=10)
        orf_ids = [k for k in annotation["1"] if k != "bioseq"]
        # ORF ≥ 200 nt so it should be found even with clamped threshold
        assert len(orf_ids) >= 1
