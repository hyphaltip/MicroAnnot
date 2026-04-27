"""
Funannotate microsporidia shunt.

``funannotate predict`` runs a standard eukaryotic gene-prediction stack
(Augustus, GeneMark, SNAP …) that is poorly suited to microsporidia because of
their highly reduced genomes, compact overlapping genes and non-standard codon
usage.

This module provides :class:`MicrosporidiaAnnotShunt`, which:

1. Detects the ``--organism microsporidia`` flag in a funannotate predict run
   (or is invoked directly via the ``pymicrosporidiaannot-funannotate`` CLI).
2. Runs the full :class:`~``pymicrosporidiaannot``.pipeline.MicrosporidiaAnnotPipeline` instead of
   the standard Augustus/GeneMark/SNAP pipeline.
3. Uses :mod:```pymicrosporidiaannot``.funannotate.adapter` to convert MicroAnnot's output
   into the file layout that funannotate's downstream ``annotate`` and
   ``compare`` commands expect.

Typical use via the CLI
-----------------------
::

    pymicrosporidiaannot-funannotate \\
        --input         genome.fasta \\
        --out           predict_results \\
        --db-dir        /opt/``pymicrosporidiaannot``_db \\
        --icm           microsporidia.icm \\
        --species       "Nosema ceranae" \\
        --proteomes     "Nosema_ceranae,Encephalitozoon_cuniculi"

Integration with existing funannotate installations
----------------------------------------------------
If funannotate is installed you can import this module from a funannotate
plugin / wrapper script and call :meth:`MicrosporidiaAnnotShunt.run`.  The resulting
``predict_results/`` directory is then passed unchanged to
``funannotate annotate``.
"""

import argparse
import os

from pymicrosporidiaannot.pipeline import MicrosporidiaAnnotPipeline
from pymicrosporidiaannot.funannotate.adapter import convert_to_funannotate


class MicrosporidiaAnnotShunt:
    """Run MicroAnnot in place of the standard funannotate predict pipeline.

    Parameters
    ----------
    input_file:
        Genomic FASTA file.
    out_dir:
        Destination directory for funannotate-compatible output
        (``predict_results/``).
    db_dir:
        pyMicrosporidiaAnnot database root.
    icm_file:
        Glimmer ICM filename within ``{db_dir}/db_ref_Glimmer/``.
    species:
        Species name (used in GFF3 ``source`` field and sequence metadata).
    proteomes:
        Comma-separated proteome names for the BLAST database.
    min_orf_size:
        Minimum ORF length (nt).
    glimmer_size:
        Minimum Glimmer gene size (nt).
    evalue_alignment:
        blastp E-value.
    evalue_small:
        blastx small-CDS E-value.
    evalue_te:
        tblastx TE E-value.
    bool_interpro:
        Enable InterproScan.
    dir_blast:
        BLAST+ binary directory.
    dir_glimmer:
        Glimmer3 binary directory.
    dir_interproscan:
        InterproScan root directory.
    debug:
        Verbose diagnostic mode.
    """

    def __init__(
        self,
        input_file: str,
        out_dir: str = "predict_results",
        db_dir: str = "",
        icm_file: str = "",
        species: str = "Microsporidia sp.",
        proteomes: str = "",
        min_orf_size: int = 240,
        glimmer_size: int = 300,
        evalue_alignment: float = 1e-15,
        evalue_small: float = 1e-5,
        evalue_te: float = 1e-10,
        bool_interpro: bool = False,
        dir_blast: str = "",
        dir_glimmer: str = "",
        dir_interproscan: str = "",
        debug: bool = False,
    ):
        self.input_file = input_file
        self.out_dir = out_dir
        self.db_dir = db_dir
        self.species = species

        self._pipeline = MicrosporidiaAnnotPipeline(
            input_file=input_file,
            db_dir=db_dir,
            icm_file=icm_file,
            list_dat_name=proteomes,
            min_orf_size=min_orf_size,
            glimmer_size=glimmer_size,
            evalue_alignment=evalue_alignment,
            evalue_small=evalue_small,
            evalue_te=evalue_te,
            bool_interpro=bool_interpro,
            dir_blast=dir_blast,
            dir_glimmer=dir_glimmer,
            dir_interproscan=dir_interproscan,
            debug=debug,
        )

    def run(self) -> str:
        """Execute the microsporidia shunt and return the output directory.

        Returns
        -------
        str
            Path to the funannotate-compatible ``predict_results/`` directory.
        """
        # 1. Run the full pyMicrosporidiaAnnot pipeline
        self._pipeline.run()

        # 2. Locate the GFF3 results
        ma_gff_dir = os.path.join(self._pipeline.out_dir, "Res_gff_annot")

        # 3. Convert to funannotate output layout
        os.makedirs(self.out_dir, exist_ok=True)
        convert_to_funannotate(
            gff_dir=ma_gff_dir,
            genome_fasta=self._pipeline.new_input_file,
            id_2_name=self._pipeline.id_2_name,
            out_dir=self.out_dir,
            species=self.species,
        )

        return self.out_dir


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pymicrosporidiaannot-funannotate",
        description=(
            "Run the pyMicrosporidiaAnnot microsporidia prediction shunt and produce\n"
            "a funannotate-compatible predict_results/ directory."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", "-i", required=True, metavar="FASTA",
                   help="Input genomic FASTA.")
    p.add_argument("--out", "-o", default="predict_results", metavar="DIR",
                   help="Output directory [default: predict_results].")
    p.add_argument("--db-dir", "-d", required=True, metavar="DIR",
                   help="pyMicrosporidiaAnnot database root directory.")
    p.add_argument("--icm", default="", metavar="FILE",
                   help="Glimmer ICM filename within {db_dir}/db_ref_Glimmer/.")
    p.add_argument("--species", default="Microsporidia sp.", metavar="NAME",
                   help="Species name [default: 'Microsporidia sp.'].")
    p.add_argument("--proteomes", default="", metavar="NAMES",
                   help="Comma-separated proteome names for BLAST DB.")
    p.add_argument("--min-orf", type=int, default=240, metavar="NT",
                   help="Minimum ORF size (nt) [default: 240].")
    p.add_argument("--glimmer-size", type=int, default=300, metavar="NT",
                   help="Minimum Glimmer gene size (nt) [default: 300].")
    p.add_argument("--evalue-align", type=float, default=1e-15, metavar="E",
                   help="E-value for blastp [default: 1e-15].")
    p.add_argument("--evalue-small", type=float, default=1e-5, metavar="E",
                   help="E-value for blastx small-CDS [default: 1e-5].")
    p.add_argument("--evalue-te", type=float, default=1e-10, metavar="E",
                   help="E-value for tblastx TE scan [default: 1e-10].")
    p.add_argument("--interproscan", action="store_true",
                   help="Enable InterproScan.")
    p.add_argument("--blast-dir", default="", metavar="DIR",
                   help="BLAST+ binary directory.")
    p.add_argument("--glimmer-dir", default="", metavar="DIR",
                   help="Glimmer3 binary directory.")
    p.add_argument("--interproscan-dir", default="", metavar="DIR",
                   help="InterproScan installation directory.")
    p.add_argument("--debug", action="store_true",
                   help="Verbose diagnostic output.")
    return p


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    shunt = MicrosporidiaAnnotShunt(
        input_file=args.input,
        out_dir=args.out,
        db_dir=args.db_dir,
        icm_file=args.icm,
        species=args.species,
        proteomes=args.proteomes,
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
    out = shunt.run()
    print(f"Done. Results in: {out}")


if __name__ == "__main__":
    main()
