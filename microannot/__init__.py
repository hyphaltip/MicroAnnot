"""
MicroAnnot – Python port of the microannot.pl pipeline for microsporidia genome annotation.

The package provides:
  * A fully functional, importable Python API (``MicroAnnotPipeline``)
  * A stand-alone CLI entry-point (``microannot`` / ``python -m microannot``)
  * A funannotate shunt (``microannot.funannotate``) that redirects
    ``funannotate predict --organism microsporidia`` into this pipeline.
"""

from microannot.pipeline import MicroAnnotPipeline

__all__ = ["MicroAnnotPipeline"]
__version__ = "1.0.0"
