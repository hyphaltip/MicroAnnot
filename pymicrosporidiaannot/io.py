"""
FASTA pre-processing and related I/O helpers.

The original Perl script rewrites the input FASTA so that each sequence header
is replaced by a monotonically-increasing integer ID (``>1``, ``>2``, …).  The
original headers are preserved in a mapping dict.  This avoids BLAST or Glimmer
choking on headers that contain spaces or special characters.
"""

import os


def preprocess_fasta(input_file: str, out_dir: str) -> tuple:
    """Sanitise a FASTA file and renumber sequences.

    * Strips spaces, tabs, CR, LF from every line.
    * Upper-cases nucleotides.
    * Replaces each header with a sequential integer ID (``>1``, ``>2``, …).

    Parameters
    ----------
    input_file:
        Original FASTA file.
    out_dir:
        Directory where the sanitised copy will be written
        (``{out_dir}/{basename}``).

    Returns
    -------
    new_input_file : str
        Path to the sanitised FASTA file.
    id_2_name : dict
        Maps integer-string IDs (``"1"``, ``"2"``, …) to the original FASTA
        headers (including the ``>`` prefix).
    """
    os.makedirs(out_dir, exist_ok=True)
    basename = os.path.basename(input_file)
    new_input_file = os.path.join(out_dir, basename)

    id_2_name: dict = {}
    seq_id = 0

    with open(input_file) as fh, open(new_input_file, "w") as out:
        for line in fh:
            line = line.strip().replace(" ", "").replace("\t", "").upper()
            if not line:
                continue
            if line.startswith(">"):
                seq_id += 1
                id_2_name[str(seq_id)] = line
                out.write(f">{seq_id}\n")
            elif seq_id != 0:
                out.write(line + "\n")

    return new_input_file, id_2_name
