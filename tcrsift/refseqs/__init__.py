"""Packaged reference data shipped with tcrsift.

Contents (see ``PROVENANCE.md`` in this directory for exact sources and build scripts):

* ``canonical_constants.fasta`` / ``canonical_constants_alleles.fasta`` — bare mature
  human ``TRAC``, ``TRBC1``, ``TRBC2`` canonical and allelic protein sequences.
* ``canonical_constants_mouse.fasta`` / ``canonical_constants_alleles_mouse.fasta`` —
  bare mature murine ``Trac``, ``Trbc1``, ``Trbc2`` sequences.
* ``vgene_leaders.csv.gz`` — human IMGT germline ``TRAV``/``TRBV`` signal peptides.
* ``kmer_pgen_{alpha,beta}.npz`` and ``kmer_ppost_{alpha,beta}.npz`` — bundled NumPy
  order-2 amino-acid Markov repertoire probability models.
* ``observed_pbmc_10x.csv.gz`` — offline 10x PBMC repertoire reference for model
  training/verification.
"""
