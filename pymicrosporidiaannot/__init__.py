"""
pyMicrosporidiaAnnot – Python port of the ``pymicrosporidiaannot``.pl pipeline for microsporidia genome annotation.

The package provides:
  * A fully functional, importable Python API (``MicrosporidiaAnnotPipeline``)
  * A stand-alone CLI entry-point (````pymicrosporidiaannot```` / ``python -m ``pymicrosporidiaannot````)
  * A funannotate shunt (````pymicrosporidiaannot``.funannotate``) that redirects
    ``funannotate predict --organism microsporidia`` into this pipeline.
"""

from pymicrosporidiaannot.pipeline import MicrosporidiaAnnotPipeline

__all__ = ["MicrosporidiaAnnotPipeline"]
__version__ = "1.0.0"
