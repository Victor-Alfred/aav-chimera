"""
aav-chimera
===========
Benchmark-validated chimeric read detection for recombinant AAV vectors
sequenced on Oxford Nanopore platforms.

Chimeric molecules — single reads spanning the vector genome, helper
plasmid, and host DNA — are detected from split-read alignment
signatures, inside a full vector-QC pipeline. A companion simulator
generates reads with known ground truth so detection accuracy can be
measured rather than asserted.

Modules
-------
aav_chimera.qc         : detection + QC pipeline (alignment, genome-structure
                         classification, chimeric/backbone/host detection,
                         HTML reporting).           CLI: `aav-chimera`
aav_chimera.simulator  : ground-truth Nanopore read simulator and the
                         benchmark harness.         CLI: `aav-chimera-sim`

This project is independent of Oxford Nanopore's `wf-aav-qc` workflow.
"""

__version__ = "1.0.0"
__all__ = ["__version__"]
