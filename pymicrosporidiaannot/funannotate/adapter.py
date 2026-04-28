"""
Output adapter: convert pyMicrosporidiaAnnot GFF3 output to the file layout expected by
``funannotate annotate`` and ``funannotate compare``.

funannotate expects the ``predict_results/`` directory to contain:

==========================================  ====================================
File                                        Description
==========================================  ====================================
``genome.fasta``                            Soft-masked (or unmasked) genome
``genome.gff3``                             GFF3 with gene/mRNA/exon/CDS rows
``genome.proteins.fa``                      Protein FASTA (one entry per gene)
``genome.mrna-transcripts.fa``              mRNA-transcript FASTA
``genome.cds-transcripts.fa``               CDS nucleotide FASTA
``funannotate_train.stringtie.gtf``         (optional, left empty here)
==========================================  ====================================

The GFF3 format used by funannotate has the following mandatory hierarchy::

    gene         ID=gene_1; locus_tag=gene_1
    mRNA         ID=mRNA_1; Parent=gene_1
    exon         Parent=mRNA_1
    CDS          ID=CDS_1; Parent=mRNA_1

Genes are numbered sequentially as ``{species_code}_{XXXXXX}`` where the
species code is the first three letters of the genus + first three of the
species (e.g. ``NOSCE_`` for *Nosema ceranae*).
"""

import os
import re
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# GFF3 record container
# ---------------------------------------------------------------------------

class _GFFRecord:
    __slots__ = ("seqid", "source", "ftype", "start", "end",
                 "score", "strand", "phase", "attrs")

    def __init__(self, line: str):
        cols = line.rstrip("\n").split("\t")
        self.seqid = cols[0]
        self.source = cols[1]
        self.ftype = cols[2]
        self.start = int(cols[3])          # 1-based
        self.end = int(cols[4])            # 1-based inclusive
        self.score = cols[5]
        self.strand = cols[6]
        self.phase = cols[7]
        self.attrs: dict = {}
        for part in cols[8].split(";"):
            part = part.strip()
            if "=" in part:
                k, _, v = part.partition("=")
                self.attrs[k] = v

    def attr_str(self) -> str:
        return ";".join(f"{k}={v}" for k, v in self.attrs.items())

    def gff3_line(self) -> str:
        return (
            f"{self.seqid}\t{self.source}\t{self.ftype}\t"
            f"{self.start}\t{self.end}\t{self.score}\t"
            f"{self.strand}\t{self.phase}\t{self.attr_str()}"
        )


# ---------------------------------------------------------------------------
# Species code helper
# ---------------------------------------------------------------------------

def _species_code(species: str) -> str:
    """Generate a 6-letter species code from a binomial name."""
    parts = species.split()
    if len(parts) >= 2:
        genus = parts[0][:3].upper()
        epithet = parts[1][:3].upper()
        return genus + epithet + "_"
    return "MICRO_"


# ---------------------------------------------------------------------------
# Main converter
# ---------------------------------------------------------------------------

def convert_to_funannotate(
    gff_dir: str,
    genome_fasta: str,
    id_2_name: dict,
    out_dir: str,
    species: str = "Microsporidia sp.",
) -> None:
    """Convert pyMicrosporidiaAnnot output to the funannotate ``predict_results/`` layout.

    Parameters
    ----------
    gff_dir:
        Directory containing per-contig ``.gff`` files from MicroAnnot.
    genome_fasta:
        Path to the *original* (pre-sanitised) genomic FASTA file.  If the
        sanitised copy is used instead, supply *id_2_name* so that sequence
        IDs can be mapped back to original headers.
    id_2_name:
        Mapping of integer IDs (``"1"``, ``"2"``, …) to original FASTA
        headers.  Pass ``{}`` when *genome_fasta* already has original headers.
    out_dir:
        Output directory (funannotate's ``predict_results/``).
    species:
        Species name used to derive the locus-tag prefix.
    """
    os.makedirs(out_dir, exist_ok=True)
    sp_code = _species_code(species)

    # Read the genome
    genome: dict = {}
    for rec in SeqIO.parse(genome_fasta, "fasta"):
        int_id = rec.id
        orig_header = id_2_name.get(int_id, f">{rec.id}").lstrip(">")
        orig_id = orig_header.split()[0]
        genome[int_id] = {"orig_id": orig_id, "seq": rec.seq}

    # Collect all GFF features from the pyMicrosporidiaAnnot output directory
    features = []   # list of _GFFRecord
    for gff_file in sorted(os.listdir(gff_dir)):
        if not gff_file.endswith(".gff"):
            continue
        with open(os.path.join(gff_dir, gff_file)) as fh:
            seqid = ""
            for line in fh:
                if line.startswith("#"):
                    seqid = line[1:].strip()
                    continue
                if not line.strip():
                    continue
                try:
                    rec = _GFFRecord(line)
                    rec.seqid = seqid
                    features.append(rec)
                except (ValueError, IndexError):
                    pass

    # Build gene/mRNA/CDS/exon hierarchy
    gene_counter = 0
    gff3_lines = ["##gff-version 3"]
    protein_records = []
    mrna_records = []
    cds_records = []

    for feat in features:
        if feat.ftype not in ("CDS", "gene"):
            continue
        if feat.ftype == "gene":
            # Frameshifted / truncated – output as a non-coding gene feature
            gene_counter += 1
            gid = f"{sp_code}{gene_counter:06d}"
            mid = f"{gid}-T1"
            strand = feat.strand
            seq_entry = genome.get(feat.seqid, {})
            orig_id = seq_entry.get("orig_id", feat.seqid)

            note = feat.attrs.get("Note", "").replace("%3B", ";")
            gff3_lines.append(
                f"{orig_id}\tMicroAnnot\tgene\t{feat.start}\t{feat.end}"
                f"\t.\t{strand}\t.\t"
                f"ID={gid};locus_tag={gid};Name={gid}"
            )
            gff3_lines.append(
                f"{orig_id}\tMicroAnnot\tmRNA\t{feat.start}\t{feat.end}"
                f"\t.\t{strand}\t.\t"
                f"ID={mid};Parent={gid};Note={note}"
            )
            gff3_lines.append(
                f"{orig_id}\tMicroAnnot\texon\t{feat.start}\t{feat.end}"
                f"\t.\t{strand}\t.\t"
                f"Parent={mid}"
            )
            continue

        # CDS feature
        gene_counter += 1
        gid = f"{sp_code}{gene_counter:06d}"
        mid = f"{gid}-T1"
        strand = feat.strand
        seq_entry = genome.get(feat.seqid, {})
        orig_id = seq_entry.get("orig_id", feat.seqid)
        seq = seq_entry.get("seq")

        note = feat.attrs.get("Note", "").replace("%3B", ";")
        hit = feat.attrs.get("hit", "")
        evalue = feat.attrs.get("evalue", "")

        # Extract CDS nucleotide sequence (strip < / > from partial coords)
        raw_start = feat.start
        raw_end = feat.end
        cds_len = raw_end - raw_start + 1

        if seq is not None:
            nt_seq = str(seq)[raw_start - 1: raw_end]
            if strand == "-":
                nt_seq = str(Seq(nt_seq).reverse_complement())
            prot_seq = str(Seq(nt_seq).translate()).rstrip("*")
        else:
            nt_seq = ""
            prot_seq = ""

        # GFF3 gene row
        extra_attrs = ""
        if hit:
            extra_attrs += f";hit={hit}"
        if evalue:
            extra_attrs += f";evalue={evalue}"

        gff3_lines.append(
            f"{orig_id}\tMicroAnnot\tgene\t{raw_start}\t{raw_end}"
            f"\t.\t{strand}\t.\t"
            f"ID={gid};locus_tag={gid};Name={gid}"
        )
        gff3_lines.append(
            f"{orig_id}\tMicroAnnot\tmRNA\t{raw_start}\t{raw_end}"
            f"\t.\t{strand}\t.\t"
            f"ID={mid};Parent={gid};Note={note}{extra_attrs}"
        )
        gff3_lines.append(
            f"{orig_id}\tMicroAnnot\texon\t{raw_start}\t{raw_end}"
            f"\t.\t{strand}\t.\t"
            f"Parent={mid}"
        )
        gff3_lines.append(
            f"{orig_id}\tMicroAnnot\tCDS\t{raw_start}\t{raw_end}"
            f"\t.\t{strand}\t0\t"
            f"ID={gid}_CDS;Parent={mid}"
        )

        if prot_seq:
            protein_records.append(
                SeqRecord(Seq(prot_seq), id=mid, description=gid)
            )
        if nt_seq:
            mrna_records.append(
                SeqRecord(Seq(nt_seq), id=mid, description=gid)
            )
            cds_records.append(
                SeqRecord(Seq(nt_seq), id=mid, description=gid)
            )

    # Write GFF3
    gff3_path = os.path.join(out_dir, "genome.gff3")
    with open(gff3_path, "w") as fh:
        fh.write("\n".join(gff3_lines) + "\n")

    # Copy genome FASTA (original headers)
    genome_out = os.path.join(out_dir, "genome.fasta")
    with open(genome_out, "w") as fh:
        for int_id, entry in sorted(genome.items(),
                                    key=lambda x: int(x[0])):
            fh.write(f">{entry['orig_id']}\n{entry['seq']}\n")

    # Protein FASTA
    with open(os.path.join(out_dir, "genome.proteins.fa"), "w") as fh:
        SeqIO.write(protein_records, fh, "fasta")

    # mRNA transcript FASTA
    with open(os.path.join(out_dir, "genome.mrna-transcripts.fa"), "w") as fh:
        SeqIO.write(mrna_records, fh, "fasta")

    # CDS transcript FASTA
    with open(os.path.join(out_dir, "genome.cds-transcripts.fa"), "w") as fh:
        SeqIO.write(cds_records, fh, "fasta")

    # Create empty training file so funannotate annotate doesn't complain
    open(os.path.join(out_dir, "funannotate_train.stringtie.gtf"), "w").close()
