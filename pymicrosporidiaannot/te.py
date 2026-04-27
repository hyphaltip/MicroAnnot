"""
Transposable element detection via tblastx.

Runs ``tblastx`` against the consensus TE nucleotide database and returns a
dict mapping ``start-end`` coordinate strings to TE descriptions::

    id_2_te = {"100-500": "TE found by alignment: Tc1-Mariner", ...}
"""

import os
import subprocess
from Bio.Blast import NCBIXML


def run_te_blast(
    cds_nt_fasta: str,
    db: str,
    out_file: str,
    dir_blast: str = "",
    evalue: float = 1e-10,
    matrix: str = "BLOSUM45",
) -> None:
    """Run tblastx to screen CDS sequences for transposable elements.

    Parameters
    ----------
    cds_nt_fasta:
        FASTA of CDS nucleotide sequences (``CDS_gene_nt.fa``).
    db:
        Path to the TE consensus nucleotide BLAST database.
    out_file:
        Destination for BLAST XML output.
    dir_blast:
        Directory of BLAST+ binaries (empty → use PATH).
    evalue:
        E-value threshold.
    matrix:
        Substitution matrix name.
    """
    tblastx = os.path.join(dir_blast, "tblastx") if dir_blast else "tblastx"
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
            tblastx,
            "-matrix", matrix,
            "-word_size", "3",
            "-num_alignments", "10",
            "-evalue", str(evalue),
            "-db", db,
            "-query", cds_nt_fasta,
            "-out", out_file,
            "-outfmt", "5",
        ],
        check=True,
    )


def parse_te_blast(out_file: str) -> dict:
    """Parse tblastx XML results.

    Returns
    -------
    id_2_te : dict
        Maps ``"start-end"`` key strings (from query FASTA headers like
        ``>100-500``) to TE name descriptions.
    """
    id_2_te: dict = {}
    if not os.path.exists(out_file):
        return id_2_te

    with open(out_file) as fh:
        for record in NCBIXML.parse(fh):
            if not record.alignments:
                continue
            query_name = record.query.split()[0]   # e.g. "100-500"
            hit_name = record.alignments[0].hit_def.split()[0]
            if query_name not in id_2_te:
                id_2_te[query_name] = f"TE found by alignment: {hit_name}"
            else:
                id_2_te[query_name] += f", {hit_name}"

    return id_2_te
