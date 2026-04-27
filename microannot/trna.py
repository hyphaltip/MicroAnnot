"""
tRNAscan-SE runner and result parser.

Runs ``tRNAscan-SE`` on the preprocessed FASTA file and returns a nested dict::

    {seq_id: {start: {end: 1}}}

where start/end are the 1-based coordinates reported by tRNAscan-SE.
"""

import os
import subprocess


def run_trnascan(input_fasta: str, out_file: str,
                 trnascan_bin: str = "tRNAscan-SE") -> None:
    """Execute tRNAscan-SE.

    Parameters
    ----------
    input_fasta:
        Path to the preprocessed FASTA file.
    out_file:
        Destination for tRNAscan-SE plain-text output.
    trnascan_bin:
        Name or path of the ``tRNAscan-SE`` executable.
    """
    if os.path.exists(out_file):
        os.remove(out_file)
    subprocess.run(
        [trnascan_bin, "-q", "--score", "50", "-o", out_file, input_fasta],
        check=True,
    )


def parse_trnascan(out_file: str) -> dict:
    """Parse tRNAscan-SE tabular output.

    The first three lines are a header; subsequent lines are whitespace-
    separated with sequence ID in column 0, tRNA start in column 2 and
    tRNA end in column 3 (1-based).

    Returns
    -------
    trna : dict
        ``{seq_id: {start: {end: 1}}}``
    """
    trna: dict = {}
    if not os.path.exists(out_file):
        return trna
    with open(out_file) as fh:
        for _ in range(3):          # skip header lines
            fh.readline()
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split()
            if len(cols) < 4:
                continue
            seq_id = cols[0]
            start = int(cols[2])
            end = int(cols[3])
            trna.setdefault(seq_id, {}).setdefault(start, {})[end] = 1
    return trna
