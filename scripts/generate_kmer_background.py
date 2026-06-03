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
        return sg.SequenceGenerationVJ(m, g)
    g = lm.GenomicDataVDJ()
    g.load_igor_genomic_data(params, vanc, janc)
    m = lm.GenerativeModelVDJ()
    m.load_and_process_igor_model(marg)
    return sg.SequenceGenerationVDJ(m, g)


def generate_cdr3(chain: str, n: int) -> list[str]:
    seqgen = _olga_seqgen(chain)
    out = []
    t = time.time()
    while len(out) < n:
        aa = seqgen.gen_rnd_prod_CDR3()[1]  # (nt, aa, V, J)
        if aa:
            out.append(aa)
        if len(out) % 100000 == 0 and len(out) > 0:
            rate = len(out) / (time.time() - t)
            print(f"  [{chain}] {len(out)}/{n} ({rate:.0f}/s)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=500_000,
                    help="synthetic CDR3s per chain (default 500k)")
    ap.add_argument("--order", type=int, default=3, help="Markov order")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", default=None,
                    help="default: tcrsift/refseqs next to the package")
    ap.add_argument("--observed-beta", default=None,
                    help="observed beta repertoire CSV (col 'seq') for ppost; "
                         "default: TCRpeg's bundled TCRs_train.csv")
    ap.add_argument("--observed-alpha", default=None,
                    help="observed alpha repertoire CSV (col 'seq') for ppost")
    args = ap.parse_args()

    import tcrsift
    from tcrsift.seqprob import KmerProbabilityModel

    out_dir = args.out_dir or os.path.join(
        os.path.dirname(tcrsift.__file__), "refseqs"
    )
    np.random.seed(args.seed)

    def _fit_save(seqs, chain, role):
        print(f"Fitting order-{args.order} k-mer {role} on {len(seqs)} "
              f"{chain} seqs...")
        model = KmerProbabilityModel(order=args.order, chain=chain).fit(seqs)
        path = os.path.join(out_dir, f"kmer_{role}_{chain}.npz")
        model.save(path)
        print(f"  wrote {path} ({os.path.getsize(path)/1024:.0f} KB, "
              f"n_train={model.n_train})")

    # --- Pgen: OLGA-generated synthetic background, both chains ---
    for chain in ("beta", "alpha"):
        print(f"Generating {args.n} {chain} CDR3s via OLGA (pgen)...")
        _fit_save(generate_cdr3(chain, args.n), chain, "pgen")

    # --- Ppost: observed-repertoire background ---
    import pandas as pd

    def _read_seqs(path):
        try:
            return pd.read_csv(path, compression="gzip")["seq"].dropna().tolist()
        except Exception:
            return pd.read_csv(path)["seq"].dropna().tolist()

    obs_beta = args.observed_beta
    if obs_beta is None:
        try:
            import tcrpeg
            obs_beta = os.path.join(os.path.dirname(tcrpeg.__file__),
                                    "data", "TCRs_train.csv")
        except Exception:
            obs_beta = None
    if obs_beta and os.path.isfile(obs_beta):
        seqs = _read_seqs(obs_beta)
        print(f"Observed beta repertoire: {len(seqs)} seqs from {obs_beta}")
        _fit_save(seqs, "beta", "ppost")
    else:
        print("No observed beta repertoire; skipping beta ppost.")

    if args.observed_alpha and os.path.isfile(args.observed_alpha):
        seqs = _read_seqs(args.observed_alpha)
        print(f"Observed alpha repertoire: {len(seqs)} seqs")
        _fit_save(seqs, "alpha", "ppost")
    else:
        print("No observed alpha repertoire; alpha ppost will fall back "
              "to alpha pgen at load time.")


if __name__ == "__main__":
    main()
