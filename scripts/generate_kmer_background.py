#!/usr/bin/env python
# Licensed under the Apache License, Version 2.0 (the "License").
"""Fit the shipped k-mer CDR3 background models — run OFFLINE at build time.

OLGA (GPL-3.0) is used here *only* to generate a synthetic reference
repertoire; the fitted :class:`tcrsift.seqprob.KmerProbabilityModel` params
(a plain numpy table) are committed to ``tcrsift/refseqs/``. tcrsift never
imports OLGA at runtime, so it stays Apache-2.0.

Usage::

    python scripts/generate_kmer_background.py --n 500000 --order 3 --seed 0

Regenerate both chains and overwrite the committed ``kmer_background_*.npz``.
"""

from __future__ import annotations

import argparse
import os
import time

import numpy as np


def _olga_seqgen(chain: str):
    import olga
    import olga.load_model as lm
    import olga.sequence_generation as sg

    name = "human_T_alpha" if chain == "alpha" else "human_T_beta"
    d = os.path.join(os.path.dirname(olga.__file__), "default_models", name)
    params = os.path.join(d, "model_params.txt")
    marg = os.path.join(d, "model_marginals.txt")
    vanc = os.path.join(d, "V_gene_CDR3_anchors.csv")
    janc = os.path.join(d, "J_gene_CDR3_anchors.csv")
    if chain == "alpha":
        g = lm.GenomicDataVJ()
        g.load_igor_genomic_data(params, vanc, janc)
        m = lm.GenerativeModelVJ()
        m.load_and_process_igor_model(marg)
        return sg.SequenceGenerationVJ(m, g), g
    g = lm.GenomicDataVDJ()
    g.load_igor_genomic_data(params, vanc, janc)
    m = lm.GenerativeModelVDJ()
    m.load_and_process_igor_model(marg)
    return sg.SequenceGenerationVDJ(m, g), g


def generate_cdr3(chain: str, n: int):
    """Return (cdr3s, v_genes, j_genes) — gene-aware OLGA synthetic Pgen set."""
    seqgen, gd = _olga_seqgen(chain)
    genV = [v[0] for v in gd.genV]
    genJ = [j[0] for j in gd.genJ]
    cdr3s, vs, js = [], [], []
    t = time.time()
    while len(cdr3s) < n:
        _nt, aa, vi, ji = seqgen.gen_rnd_prod_CDR3()
        if aa:
            cdr3s.append(aa)
            vs.append(genV[vi] if vi < len(genV) else "")
            js.append(genJ[ji] if ji < len(genJ) else "")
        if len(cdr3s) % 100000 == 0 and len(cdr3s) > 0:
            print(f"  [{chain}] {len(cdr3s)}/{n} ({len(cdr3s)/(time.time()-t):.0f}/s)")
    return cdr3s, vs, js


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=500_000,
                    help="synthetic CDR3s per chain (default 500k)")
    ap.add_argument("--order", type=int, default=2, help="Markov order")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", default=None,
                    help="default: tcrsift/refseqs next to the package")
    ap.add_argument("--observed", default=None,
                    help="paired observed repertoire CSV (cols: chain,seq,v,j) "
                         "for ppost; default: bundled refseqs/observed_pbmc_10x.csv.gz")
    args = ap.parse_args()

    import pandas as pd

    import tcrsift
    from tcrsift.seqprob import GeneAwareKmerModel

    out_dir = args.out_dir or os.path.join(
        os.path.dirname(tcrsift.__file__), "refseqs"
    )
    np.random.seed(args.seed)

    def _fit_save(seqs, vs, js, chain, role):
        print(f"Fitting gene-aware order-{args.order} k-mer {role} on "
              f"{len(seqs)} {chain} seqs...")
        model = GeneAwareKmerModel(order=args.order, chain=chain).fit(
            seqs, v_genes=vs, j_genes=js)
        path = os.path.join(out_dir, f"kmer_{role}_{chain}.npz")
        model.save(path)
        print(f"  wrote {path} ({os.path.getsize(path)/1024:.0f} KB, "
              f"n_train={model.n_train}, gene_aware={model.gene_aware})")

    # --- Pgen: OLGA-generated synthetic background, both chains (with V/J) ---
    for chain in ("beta", "alpha"):
        print(f"Generating {args.n} {chain} CDR3s via OLGA (pgen)...")
        cdr3s, vs, js = generate_cdr3(chain, args.n)
        _fit_save(cdr3s, vs, js, chain, "pgen")

    # --- Ppost: observed healthy-PBMC background, BOTH chains from the one
    # bundled, CC-licensed pooled 10x source (refseqs/observed_pbmc_10x.csv.gz;
    # see refseqs/PROVENANCE.md). No OLGA / no TCRpeg — a plain k-mer fit on
    # observed CDR3s. `--observed` overrides the bundled source. ---
    obs = args.observed or os.path.join(out_dir, "observed_pbmc_10x.csv.gz")
    if os.path.isfile(obs):
        src = pd.read_csv(obs)
        for chain_label, tag in (("TRA", "alpha"), ("TRB", "beta")):
            sub = src[src["chain"] == chain_label]
            print(f"Observed {tag} repertoire: {len(sub)} unique CDR3 from "
                  f"{os.path.basename(obs)}")
            _fit_save(sub["seq"].tolist(), sub["v"].tolist(), sub["j"].tolist(),
                      tag, "ppost")
    else:
        print(f"No bundled observed repertoire at {obs}; skipping ppost. "
              "Pass --observed a paired-chain CDR3 CSV (cols: chain,seq,v,j).")


if __name__ == "__main__":
    main()
