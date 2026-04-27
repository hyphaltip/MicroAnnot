"""
Six-frame ORF extraction.

For every sequence in the input FASTA, all six reading frames are translated and
split at stop codons (``*``).  Candidate ORFs that meet the minimum size
threshold are written to ``orf.fa`` and registered in the annotation dictionary.

ORF identifiers follow the same convention as the original Perl code::

    {seq_id}_orf_{strand_frame}_{start}-{end}[_extr_5prim|_extr_3prim]

where ``strand_frame`` is e.g. ``+1``, ``-2`` etc., and start/end are 1-based
positions in the *translated* sequence (forward strand = + strand coordinates;
reverse strand = positions in the revcom sequence, later converted when building
the annotation features).

The ``annotation`` dict returned has the structure::

    {
        seq_id: {
            'bioseq': Bio.SeqRecord,
            orf_id: {
                'seqAA':        str,   # amino-acid sequence (no stop)
                'seqNT-add':    str,   # nt sequence + 20 nt upstream padding
                'seqNT-add-n':  int,   # actual upstream padding used (≤20)
            },
            ...
        },
        ...
    }
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq


def extract_orfs(
    input_fasta: str,
    orf_fasta: str,
    min_orf_size: int = 240,
) -> dict:
    """Extract ORFs from all sequences in *input_fasta*.

    Parameters
    ----------
    input_fasta:
        Path to the (preprocessed) FASTA file.
    orf_fasta:
        Path where the amino-acid ORF sequences will be written.
    min_orf_size:
        Minimum ORF length in *nucleotides* (default 240, minimum 200).

    Returns
    -------
    annotation : dict
        Nested dict keyed by sequence ID, then ORF ID (see module docstring).
    """
    min_orf_size = max(200, min_orf_size)
    annotation = {}

    with open(orf_fasta, "w") as orf_out:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            seq_id = rec.id
            seq_str = str(rec.seq).upper()
            seq_len = len(seq_str)
            annotation[seq_id] = {"bioseq": rec}

            revcom_str = str(rec.seq.reverse_complement()).upper()

            for i_cadre in range(6):
                n_cadre = i_cadre
                if i_cadre > 2:
                    n_cadre = i_cadre - 3
                    i_seq = revcom_str
                    cadre = f"-{n_cadre + 1}"
                else:
                    i_seq = seq_str
                    cadre = f"+{n_cadre + 1}"

                # Translate the whole sequence in this frame
                translated = str(Seq(i_seq[n_cadre:]).translate())
                split_orf = translated.split("*")

                bool_extrem_5prim = True
                start = 1 + n_cadre   # 1-based position in i_seq
                for idx, aa_seq in enumerate(split_orf):
                    nt_len = len(aa_seq) * 3
                    end = start + nt_len - 1   # 1-based inclusive

                    if (end - start + 1) >= min_orf_size:
                        extrem = ""
                        if bool_extrem_5prim:
                            extrem = "_extr_5prim"
                        elif idx == len(split_orf) - 1:
                            extrem = "_extr_3prim"

                        orf_id = (
                            f"{seq_id}_orf_{cadre}_{start}-{end}{extrem}"
                        )

                        orf_out.write(f">{orf_id}\n{aa_seq}\n")
                        annotation[seq_id][orf_id] = {
                            "seqAA": aa_seq,
                        }

                        # Store ≤ 20 nt upstream of the ORF (including
                        # padding for Met detection)
                        start_nt = start - 20
                        add_nt = 20
                        if start_nt <= 0:
                            add_nt = start - 1  # available upstream nt
                            start_nt = 1
                        annotation[seq_id][orf_id]["seqNT-add-n"] = add_nt
                        # The full nt run from (start-add_nt) to end
                        annotation[seq_id][orf_id]["seqNT-add"] = i_seq[
                            start_nt - 1: end
                        ]

                    bool_extrem_5prim = False
                    start += nt_len + 3   # skip the stop codon (3 nt)

    return annotation
