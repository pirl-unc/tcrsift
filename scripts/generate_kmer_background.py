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
    ap.add_argument("--observed-beta", default=None,
                    help="observed beta repertoire CSV (col 'seq') for ppost; "
                         "default: TCRpeg's bundled TCRs_train.csv")
    ap.add_argument("--observed-alpha", default=None,
                    help="observed alpha repertoire CSV (col 'seq') for ppost")
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

    # --- Ppost: observed-repertoire background (with V/J columns) ---
    def _read_obs(path):
        try:
            frame = pd.read_csv(path, compression="gzip")
        except Exception:
            frame = pd.read_csv(path)
        frame = frame.dropna(subset=["seq"])
        v = frame["v"].tolist() if "v" in frame.columns else None
        j = frame["j"].tolist() if "j" in frame.columns else None
        return frame["seq"].tolist(), v, j

    obs_beta = args.observed_beta
    if obs_beta is None:
        try:
            import tcrpeg
            obs_beta = os.path.join(os.path.dirname(tcrpeg.__file__),
                                    "data", "TCRs_train.csv")
        except Exception:
            obs_beta = None
    if obs_beta and os.path.isfile(obs_beta):
        seqs, vs, js = _read_obs(obs_beta)
        print(f"Observed beta repertoire: {len(seqs)} seqs from {obs_beta}")
        _fit_save(seqs, vs, js, "beta", "ppost")
    else:
        print("No observed beta repertoire; skipping beta ppost.")

    if args.observed_alpha and os.path.isfile(args.observed_alpha):
        seqs, vs, js = _read_obs(args.observed_alpha)
        print(f"Observed alpha repertoire: {len(seqs)} seqs")
        _fit_save(seqs, vs, js, "alpha", "ppost")
    else:
        print("No observed alpha repertoire; alpha ppost unavailable until a "
              "neutral reference is fetched (tcrsift pgen fetch --url ...).")


if __name__ == "__main__":
    main()
