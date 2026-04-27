"""
InterproScan runner and XML result parser.

Runs ``interproscan.sh`` (if enabled) and parses the resulting XML to build a
dict of per-CDS functional annotations::

    id_2_interproscan = {
        "100-500": {
            "function": {"kinase activity": 1, ...},
            "db_xref":  {"IPR001234": 1, "GO:0006468": 1, ...},
        },
        ...
    }
"""

import os
import re
import subprocess


def run_interproscan(
    cds_aa_fasta: str,
    out_dir: str,
    interproscan_dir: str = "",
) -> None:
    """Run InterproScan on the CDS amino-acid FASTA.

    Parameters
    ----------
    cds_aa_fasta:
        Path to ``CDS_aa.fa``.
    out_dir:
        Output directory for InterproScan results.
    interproscan_dir:
        Root directory of the InterproScan installation (the script
        ``interproscan.sh`` must be at its top level).
    """
    ipr_sh = (
        os.path.join(interproscan_dir, "interproscan.sh")
        if interproscan_dir
        else "interproscan.sh"
    )
    subprocess.run(
        [
            ipr_sh,
            "--goterms",
            "--iprlookup",
            "--output-dir", out_dir,
            "-i", cds_aa_fasta,
            "-f", "xml",
        ],
        check=True,
    )


def parse_interproscan_xml(xml_file: str) -> dict:
    """Parse an InterproScan XML result file.

    The parser is a lightweight regex-based implementation that mirrors the
    original Perl line-by-line approach, avoiding a full XML library.

    Parameters
    ----------
    xml_file:
        Path to ``CDS_aa.fa.xml`` produced by InterproScan.

    Returns
    -------
    id_2_interproscan : dict
        Nested dict keyed by CDS coordinate string then by qualifier then by
        value (see module docstring).
    """
    id_2_interproscan: dict = {}
    if not os.path.exists(xml_file):
        return id_2_interproscan

    id_list: dict = {}

    with open(xml_file) as fh:
        for line in fh:
            line = line.rstrip("\n")

            # New sequence block
            if "<sequence md5=" in line:
                id_list = {}
                continue

            # CDS coordinate identifier, e.g. id="100-500"
            m = re.search(r'^\s+<xref\sid="(\d+-\d+)', line)
            if m:
                id_list[m.group(1)] = 1
                continue

            # InterPro entry (function + db_xref)
            m = re.search(
                r'^\s+<entry\s+ac="(\w+)"\s+desc="(.+?)"\s+', line
            )
            if m:
                acc = m.group(1)
                desc = m.group(2)
                for cds_id in id_list:
                    id_2_interproscan.setdefault(cds_id, {})
                    id_2_interproscan[cds_id].setdefault(
                        "function", {}
                    )[desc] = 1
                    id_2_interproscan[cds_id].setdefault(
                        "db_xref", {}
                    )[acc] = 1
                continue

            # GO cross-reference (db_xref)
            m = re.search(
                r'^\s+<go-xref\s+category="\w+"\s+db="GO"\s+id="(\w+:\w+)"',
                line,
            )
            if m:
                go_id = m.group(1)
                for cds_id in id_list:
                    id_2_interproscan.setdefault(cds_id, {})
                    id_2_interproscan[cds_id].setdefault(
                        "db_xref", {}
                    )[go_id] = 1
                continue

            # GO cross-reference with name (function)
            m = re.search(
                r'^\s+<go-xref\s+category="\w+_\w+"\s+db="GO"\s+'
                r'id="\w+:\w+"\s+name="(.+?)"',
                line,
            )
            if m:
                go_name = m.group(1)
                for cds_id in id_list:
                    id_2_interproscan.setdefault(cds_id, {})
                    id_2_interproscan[cds_id].setdefault(
                        "function", {}
                    )[go_name] = 1

    return id_2_interproscan
