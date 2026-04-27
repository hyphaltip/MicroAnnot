"""
Command-line interface for the MicroAnnot pipeline.

Usage
-----
::

    microannot \\
        --input     genome.fasta \\
        --db-dir    /opt/microannot_db \\
        --icm       microsporidia.icm \\
        --proteomes "Nosema_ceranae,Encephalitozoon_cuniculi" \\
        [options]

Options mirror the positional arguments of the original Perl script but use
named, self-documenting flags instead.
"""

import argparse
import sys

from microannot.pipeline import MicroAnnotPipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="microannot",
        description=(
            "MicroAnnot – microsporidia genome annotation pipeline.\n"
            "Python port of microannot.pl."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input", "-i", required=True, metavar="FASTA",
        help="Input genomic FASTA file.",
    )
    p.add_argument(
        "--db-dir", "-d", required=True, metavar="DIR",
        help="Root of the MicroAnnot database tree.",
    )
    p.add_argument(
        "--icm", default="", metavar="FILE",
        help=(
            "Glimmer ICM filename (relative to {db_dir}/db_ref_Glimmer/). "
            "A new ICM is trained from the blastp CDS when ≥50 high-quality "
            "CDS are found; this file is used as the fallback."
        ),
    )
    p.add_argument(
        "--proteomes", default="", metavar="NAMES",
        help=(
            "Comma-separated list of proteome names used to build the BLAST "
            "database (e.g. 'Nosema_ceranae,Encephalitozoon_cuniculi')."
        ),
    )
    p.add_argument(
        "--proteome-id", type=int, default=0, metavar="ID",
        help="Numeric proteome-DB selector (0 = use --proteomes).",
    )
    p.add_argument(
        "--min-orf", type=int, default=240, metavar="NT",
        help="Minimum ORF size in nucleotides [default: 240, min: 200].",
    )
    p.add_argument(
        "--glimmer-size", type=int, default=300, metavar="NT",
        help="Minimum gene size for Glimmer [default: 300, min: 300].",
    )
    p.add_argument(
        "--evalue-align", type=float, default=1e-15, metavar="E",
        help="E-value for blastp CDS alignment [default: 1e-15].",
    )
    p.add_argument(
        "--evalue-small", type=float, default=1e-5, metavar="E",
        help="E-value for blastx small-CDS scan [default: 1e-5].",
    )
    p.add_argument(
        "--evalue-te", type=float, default=1e-10, metavar="E",
        help="E-value for tblastx TE detection [default: 1e-10].",
    )
    p.add_argument(
        "--interproscan", action="store_true",
        help="Enable InterproScan functional annotation.",
    )
    p.add_argument(
        "--blast-dir", default="", metavar="DIR",
        help="Directory containing BLAST+ binaries (default: use PATH).",
    )
    p.add_argument(
        "--glimmer-dir", default="", metavar="DIR",
        help="Directory containing Glimmer3 binaries (default: use PATH).",
    )
    p.add_argument(
        "--interproscan-dir", default="", metavar="DIR",
        help="Root directory of the InterproScan installation.",
    )
    p.add_argument(
        "--debug", action="store_true",
        help="Emit verbose diagnostic notes in annotation qualifiers.",
    )
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    pipeline = MicroAnnotPipeline(
        input_file=args.input,
        db_dir=args.db_dir,
        icm_file=args.icm,
        list_dat_name=args.proteomes,
        list_dat_id=args.proteome_id,
        min_orf_size=args.min_orf,
        glimmer_size=args.glimmer_size,
        evalue_alignment=args.evalue_align,
        evalue_small=args.evalue_small,
        evalue_te=args.evalue_te,
        bool_interpro=args.interproscan,
        dir_blast=args.blast_dir,
        dir_glimmer=args.glimmer_dir,
        dir_interproscan=args.interproscan_dir,
        debug=args.debug,
    )
    pipeline.run()


if __name__ == "__main__":
    main()
