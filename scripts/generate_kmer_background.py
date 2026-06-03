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
    args = ap.parse_args()

    import tcrsift
    from tcrsift.seqprob import KmerProbabilityModel

    out_dir = args.out_dir or os.path.join(
        os.path.dirname(tcrsift.__file__), "refseqs"
    )
    np.random.seed(args.seed)

    for chain in ("beta", "alpha"):
        print(f"Generating {args.n} {chain} CDR3s via OLGA...")
        seqs = generate_cdr3(chain, args.n)
        print(f"Fitting order-{args.order} k-mer model on {len(seqs)} {chain} seqs...")
        model = KmerProbabilityModel(order=args.order, chain=chain).fit(seqs)
        path = os.path.join(out_dir, f"kmer_background_{chain}.npz")
        model.save(path)
        size_kb = os.path.getsize(path) / 1024
        print(f"  wrote {path} ({size_kb:.0f} KB, n_train={model.n_train})")


if __name__ == "__main__":
    main()
