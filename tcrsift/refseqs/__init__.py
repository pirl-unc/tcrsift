"""Packaged reference data shipped with tcrsift.

Currently:

* ``canonical_constants.fasta`` — human TCR constant-region AA sequences
  for TRAC / TRBC1 / TRBC2 (generated from pyensembl GRCh38 release 110,
  verified against UniProt P01848 / P01850 / A0A5B9). Loaded by
  :mod:`tcrsift.assemble` at import time.

Named ``refseqs`` rather than ``data`` because :mod:`tcrsift.data` already
exists as a top-level module and a sibling package would shadow it.
"""
