"""
MicrosporidiaAnnotPipeline – top-level orchestration class.

Each public ``run_*`` method corresponds to one logical stage of the pipeline.
The stages can be called individually for testing / partial runs, or all at once
via :meth:`run`.

Example
-------
::

    from pymicrosporidiaannot import MicrosporidiaAnnotPipeline

    pipeline = MicrosporidiaAnnotPipeline(
        input_file   = "genome.fasta",
        db_dir       = "/opt/``pymicrosporidiaannot``_db",
        icm_file     = "microsporidia.icm",
        list_dat_name= "Nosema_ceranae,Encephalitozoon_cuniculi",
    )
    pipeline.run()
"""

import os
import subprocess
from Bio.Blast import NCBIXML

from pymicrosporidiaannot.io import preprocess_fasta
from pymicrosporidiaannot.orf import extract_orfs
from pymicrosporidiaannot.alignment import (
    annot_alignment,
    annot_alignment_small,
    group_frameshifts,
)
from pymicrosporidiaannot.trna import run_trnascan, parse_trnascan
from pymicrosporidiaannot.rrna import run_rrna_blast, parse_rrna_blast
from pymicrosporidiaannot.glimmer import build_icm, run_glimmer, parse_and_refine_glimmer
from pymicrosporidiaannot.te import run_te_blast, parse_te_blast
from pymicrosporidiaannot.interproscan import run_interproscan, parse_interproscan_xml
from pymicrosporidiaannot.annotate import write_annotations, write_final_annotations


class MicrosporidiaAnnotPipeline:
    """Orchestrate the full MicroAnnot gene-prediction and annotation workflow.

    Parameters
    ----------
    input_file:
        Path to the input genomic FASTA file.
    db_dir:
        Root of the pyMicrosporidiaAnnot database tree (must contain ``db_blast/``,
        ``small_db/``, ``db_ref_Glimmer/``).
    icm_file:
        Filename *within* ``{db_dir}/db_ref_Glimmer/`` for the Glimmer ICM
        training file (or path to a pre-built ``.icm``).
    list_dat_name:
        Comma-separated list of proteome names (used to select BLAST
        databases, e.g. ``"Nosema_ceranae,Encephalitozoon_cuniculi"``).
    list_dat_id:
        Numeric ID used to select the combined proteome BLAST DB
        (``0`` = use *list_dat_name* directly).
    min_orf_size:
        Minimum ORF size in nucleotides (default 240, minimum 200).
    glimmer_size:
        Minimum gene size for Glimmer (default 300, minimum 300).
    evalue_alignment:
        E-value for the blastp CDS-alignment step.
    evalue_small:
        E-value for the blastx small-CDS step.
    evalue_te:
        E-value for the tblastx TE-detection step.
    bool_interpro:
        Run InterproScan when ``True``.
    dir_blast:
        Directory containing BLAST+ executables (empty → use PATH).
    dir_glimmer:
        Directory containing Glimmer3 executables (empty → use PATH).
    dir_interproscan:
        Root directory of the InterproScan installation.
    debug:
        Enable verbose diagnostic output.
    """

    def __init__(
        self,
        input_file: str,
        db_dir: str,
        icm_file: str = "",
        list_dat_name: str = "",
        list_dat_id: int = 0,
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
        self.db_dir = db_dir
        self.icm_file = icm_file
        self.list_dat_name = list_dat_name.replace(" ", ".")
        self.list_dat_id = list_dat_id
        self.min_orf_size = max(200, min_orf_size)
        self.glimmer_size = max(300, glimmer_size)
        self.evalue_alignment = evalue_alignment
        self.evalue_small = evalue_small
        self.evalue_te = evalue_te
        self.bool_interpro = bool_interpro
        self.dir_blast = dir_blast
        self.dir_glimmer = dir_glimmer
        self.dir_interproscan = dir_interproscan
        self.debug = debug

        # Output directory
        file_name = os.path.basename(input_file)
        self.file_name = file_name
        self.out_dir = os.path.join("output", file_name)
        os.makedirs(self.out_dir, exist_ok=True)

        # Internal state (populated by run_* methods)
        self.new_input_file: str = ""
        self.id_2_name: dict = {}
        self.annotation: dict = {}
        self.small_cds: dict = {}
        self.orf_hit_small: dict = {}
        self.hit_2_frameshift: dict = {}
        self.trna: dict = {}
        self.rrna: dict = {}
        self.glim: dict = {}
        self.id_2_te: dict = {}
        self.id_2_interproscan: dict = {}
        self.list_file_cds_aa: dict = {}

        # Resolve the Glimmer ICM path
        if icm_file:
            icm_path = os.path.join(db_dir, "db_ref_Glimmer", icm_file)
            if not os.path.exists(icm_path):
                # Try to build it from the corresponding .fa file
                fa_path = icm_path.replace(".icm", "")
                if os.path.exists(fa_path):
                    build_icm(fa_path, icm_path, dir_glimmer=dir_glimmer)
            self.icm_path = icm_path
        else:
            self.icm_path = ""

        self._blast_bin = (
            lambda name: os.path.join(dir_blast, name) if dir_blast else name
        )

    # ------------------------------------------------------------------
    # Stage 1 – FASTA preprocessing
    # ------------------------------------------------------------------

    def run_preprocessing(self) -> tuple:
        """Sanitise the input FASTA and renumber sequence headers."""
        self.new_input_file, self.id_2_name = preprocess_fasta(
            self.input_file, self.out_dir
        )
        return self.new_input_file, self.id_2_name

    # ------------------------------------------------------------------
    # Stage 2 – ORF extraction
    # ------------------------------------------------------------------

    def run_orf_extraction(self) -> dict:
        """Extract 6-frame ORFs and populate *annotation*."""
        orf_fasta = os.path.join(self.out_dir, "orf.fa")
        self.annotation = extract_orfs(
            self.new_input_file, orf_fasta, self.min_orf_size
        )
        return self.annotation

    # ------------------------------------------------------------------
    # Stage 3 – blastp alignment
    # ------------------------------------------------------------------

    def run_blastp_alignment(self) -> dict:
        """Run blastp against the combined proteome DB and annotate CDS."""
        db_blast_dir = os.path.join(self.db_dir, "db_blast")
        db_name = f"Proteomes_{self.list_dat_name}.txt"
        db_path = os.path.join(db_blast_dir, db_name)

        # Build combined DB if needed
        if not os.path.exists(db_path + ".pdb"):
            with open(db_path, "w") as combined:
                combined.write(" \n")
            for name_part in self.list_dat_name.split(","):
                part_db = os.path.join(
                    db_blast_dir, f"Proteomes_{name_part}.txt"
                )
                with open(part_db) as src, open(db_path, "a") as dst:
                    dst.write(src.read())
            subprocess.run(
                [
                    self._blast_bin("makeblastdb"),
                    "-in", db_path,
                    "-out", db_path,
                    "-dbtype", "prot",
                ],
                check=True,
            )

        orf_fasta = os.path.join(self.out_dir, "orf.fa")
        orf_bls = os.path.join(self.out_dir, "orf.bls")
        subprocess.run(
            [
                self._blast_bin("blastp"),
                "-word_size", "3",
                "-num_alignments", "20",
                "-matrix", "BLOSUM45",
                "-evalue", str(self.evalue_alignment),
                "-db", db_path,
                "-query", orf_fasta,
                "-out", orf_bls,
                "-outfmt", "5",
            ],
            check=True,
        )

        with open(orf_bls) as fh:
            bool_hit = False
            for record in NCBIXML.parse(fh):
                bool_hit = False
                for alignment in record.alignments:
                    hsp = alignment.hsps[0]
                    annot_alignment(
                        record, alignment, hsp,
                        self.annotation, self.hit_2_frameshift,
                        bool_hit, self.debug,
                    )
                    bool_hit = True

        return self.annotation

    # ------------------------------------------------------------------
    # Stage 4 – blastx small CDS
    # ------------------------------------------------------------------

    def run_blastx_small_cds(self) -> dict:
        """Run blastx against the small-CDS database."""
        small_db = os.path.join(self.db_dir, "small_db", "small_db.fa")
        if not os.path.exists(small_db + ".pdb"):
            subprocess.run(
                [
                    self._blast_bin("makeblastdb"),
                    "-in", small_db,
                    "-out", small_db,
                    "-dbtype", "prot",
                ],
                check=True,
            )

        small_bls = os.path.join(self.out_dir, "small.bls")
        subprocess.run(
            [
                self._blast_bin("blastx"),
                "-word_size", "3",
                "-num_alignments", "500",
                "-seg", "no",
                "-evalue", str(self.evalue_small),
                "-matrix", "BLOSUM45",
                "-db", small_db,
                "-query", self.new_input_file,
                "-out", small_bls,
                "-outfmt", "5",
            ],
            check=True,
        )

        with open(small_bls) as fh:
            for record in NCBIXML.parse(fh):
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        # Skip if the translated query contains a stop codon
                        if "*" in hsp.query:
                            continue
                        annot_alignment_small(
                            record, alignment, hsp,
                            self.annotation, self.small_cds,
                            self.orf_hit_small, self.debug,
                        )

        return self.small_cds

    # ------------------------------------------------------------------
    # Stage 5 – frameshift grouping
    # ------------------------------------------------------------------

    def run_frameshift_grouping(self) -> None:
        """Group adjacent ORFs that are likely frameshifted or intron-split."""
        group_frameshifts(self.annotation, self.hit_2_frameshift)

    # ------------------------------------------------------------------
    # Stage 6 – tRNAscan-SE
    # ------------------------------------------------------------------

    def run_trna(self) -> dict:
        """Run tRNAscan-SE and parse results."""
        trna_file = os.path.join(self.out_dir, "trna.txt")
        run_trnascan(self.new_input_file, trna_file)
        self.trna = parse_trnascan(trna_file)
        return self.trna

    # ------------------------------------------------------------------
    # Stage 7 – rRNA blastn
    # ------------------------------------------------------------------

    def run_rrna(self) -> dict:
        """Run blastn against the 16S rRNA database."""
        rrna_db = os.path.join(self.db_dir, "db_blast", "16S_rRNA.txt")
        rrna_bls = os.path.join(self.out_dir, "rRNA.bls")
        run_rrna_blast(
            self.new_input_file, rrna_db, rrna_bls,
            dir_blast=self.dir_blast,
        )
        self.rrna = parse_rrna_blast(rrna_bls, self.annotation)
        return self.rrna

    # ------------------------------------------------------------------
    # Stage 8 – Glimmer / write pre-TE output
    # ------------------------------------------------------------------

    def run_glimmer(self) -> dict:
        """Build Glimmer ICM (if enough CDS) and run gene prediction."""
        learn_file = os.path.join(self.out_dir, "learn_Glim_orf.fa")

        # Write high-confidence CDS for Glimmer training
        nb_cds = 0
        with open(learn_file, "w") as out:
            for display_id, seq_data in self.annotation.items():
                for orf_id, orf_data in seq_data.items():
                    if orf_id == "bioseq":
                        continue
                    cds = orf_data.get("CDS")
                    if cds and cds["warning"] == 0:
                        parts = orf_id.split("_")
                        cadre = parts[2] if len(parts) > 2 else "+1"
                        from pymicrosporidiaannot.annotate import _trunc as _t
                        from Bio.Seq import Seq as _S
                        seq_str = str(
                            seq_data["bioseq"].seq
                        ).upper()
                        start = cds["start"]
                        end = cds["end"]
                        if cadre.startswith("-"):
                            nt = str(
                                _S(_t(seq_str, start, end))
                                .reverse_complement()
                            )
                        else:
                            nt = _t(seq_str, start, end)
                        out.write(
                            f">{orf_id}-{cds['hit']}\n{nt}\n"
                        )
                        nb_cds += 1

        # Retrain if sufficient CDS
        icm_to_use = self.icm_path
        if nb_cds >= 50:
            new_icm = os.path.join(self.out_dir, "learn_Glim_orf.icm")
            build_icm(learn_file, new_icm, dir_glimmer=self.dir_glimmer)
            icm_to_use = new_icm

        res_prefix = os.path.join(self.out_dir, "res_glimmer")
        run_glimmer(
            self.new_input_file, icm_to_use, res_prefix,
            glimmer_size=self.glimmer_size,
            dir_glimmer=self.dir_glimmer,
        )

        predict_file = res_prefix + ".predict"
        self.glim = parse_and_refine_glimmer(
            predict_file, self.annotation, self.debug
        )
        return self.glim

    # ------------------------------------------------------------------
    # Stage 9 – write intermediate output + TE + InterproScan
    # ------------------------------------------------------------------

    def run_write_annotations(self) -> None:
        """Write intermediate GenBank/EMBL/GFF3 and CDS-NT files."""
        cds_gene_nt = os.path.join(self.out_dir, "CDS_gene_nt.fa")
        self.list_file_cds_aa = write_annotations(
            self.annotation, self.small_cds, self.trna, self.rrna,
            self.glim, self.id_2_name, self.out_dir, cds_gene_nt,
            self.debug,
        )

    def run_te_detection(self) -> dict:
        """Screen CDS nucleotide sequences for transposable elements."""
        te_db = os.path.join(self.db_dir, "db_blast", "ConsensusTE.txt")
        cds_gene_nt = os.path.join(self.out_dir, "CDS_gene_nt.fa")
        te_bls = os.path.join(self.out_dir, "TE.bls")
        run_te_blast(
            cds_gene_nt, te_db, te_bls,
            dir_blast=self.dir_blast,
            evalue=self.evalue_te,
        )
        self.id_2_te = parse_te_blast(te_bls)
        return self.id_2_te

    def run_interproscan(self) -> dict:
        """Run InterproScan (when enabled) and parse results."""
        cds_aa_file = os.path.join(self.out_dir, "CDS_aa.fa")
        # Write CDS amino-acid FASTA
        with open(cds_aa_file, "w") as out:
            for coord, seq in sorted(self.list_file_cds_aa.items()):
                out.write(f">{coord}\n{seq}\n")

        ipr_dir = os.path.join(self.out_dir, "Res_Interproscan")
        os.makedirs(ipr_dir, exist_ok=True)

        if self.bool_interpro:
            run_interproscan(cds_aa_file, ipr_dir, self.dir_interproscan)

        xml_path = os.path.join(ipr_dir, "CDS_aa.fa.xml")
        self.id_2_interproscan = parse_interproscan_xml(xml_path)
        return self.id_2_interproscan

    def run_final_output(self) -> None:
        """Add TE/InterproScan annotations and write final outputs."""
        dir_gb = os.path.join(self.out_dir, "Res_gb")
        write_final_annotations(
            dir_gb, self.out_dir, self.id_2_te, self.id_2_interproscan
        )

    # ------------------------------------------------------------------
    # Convenience: run the full pipeline
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Execute all pipeline stages in order."""
        print("pyMicrosporidiaAnnot pipeline starting …")
        print(f"  Input file   : {self.input_file}")
        print(f"  Min ORF size : {self.min_orf_size} nt")
        print(f"  Glimmer size : {self.glimmer_size} nt")
        print(f"  ICM file     : {self.icm_path}")
        print(f"  InterproScan : {'Enabled' if self.bool_interpro else 'Disabled'}")
        print(f"  E-value (blastp): {self.evalue_alignment}")
        print(f"  E-value (blastx): {self.evalue_small}")
        print(f"  E-value (TE)    : {self.evalue_te}")
        print()

        self.run_preprocessing()
        self.run_orf_extraction()
        self.run_blastp_alignment()
        self.run_blastx_small_cds()
        self.run_frameshift_grouping()
        self.run_trna()
        self.run_rrna()
        self.run_glimmer()
        self.run_write_annotations()
        self.run_te_detection()
        self.run_interproscan()
        self.run_final_output()

        print("pyMicrosporidiaAnnot pipeline complete.")
        print(f"  Results in: {self.out_dir}/")
