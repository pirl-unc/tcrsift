#!/usr/bin/env python
# Licensed under the Apache License, Version 2.0 (the "License").
"""Fit the shipped k-mer CDR3 **Ppost** background models — run OFFLINE.

Fits :class:`tcrsift.seqprob.GeneAwareKmerModel` on the bundled, CC-licensed
pooled 10x healthy-PBMC repertoire (``refseqs/observed_pbmc_10x.csv.gz``; see
``refseqs/PROVENANCE.md``) and overwrites the committed
``refseqs/kmer_ppost_{alpha,beta}.npz``. No OLGA / no TCRpeg — a plain k-mer
fit on observed CDR3s.

The **Pgen** models (``kmer_pgen_*.npz``) ship pre-built and are NOT regenerated
here: they were produced offline with OLGA, which is no longer a tcrsift
dependency.

Usage::

    python scripts/generate_kmer_background.py [--order 2] [--observed PATH]
"""

from __future__ import annotations

import argparse
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order", type=int, default=2, help="Markov order (default 2)")
    ap.add_argument("--out-dir", default=None,
                    help="default: tcrsift/refseqs next to the package")
    ap.add_argument("--observed", default=None,
                    help="paired observed repertoire CSV (cols: chain,seq,v,j); "
                         "default: bundled refseqs/observed_pbmc_10x.csv.gz")
    args = ap.parse_args()

    import pandas as pd

    import tcrsift
    from tcrsift.seqprob import GeneAwareKmerModel

    out_dir = args.out_dir or os.path.join(
        os.path.dirname(tcrsift.__file__), "refseqs"
    )
    obs = args.observed or os.path.join(out_dir, "observed_pbmc_10x.csv.gz")
    if not os.path.isfile(obs):
        raise SystemExit(
            f"No observed repertoire at {obs}; pass --observed a paired-chain "
            "CDR3 CSV (cols: chain,seq,v,j)."
        )
    src = pd.read_csv(obs)
    for chain_label, tag in (("TRA", "alpha"), ("TRB", "beta")):
        sub = src[src["chain"] == chain_label]
        print(f"Fitting order-{args.order} k-mer ppost on {len(sub)} {tag} "
              f"CDR3 from {os.path.basename(obs)}...")
        model = GeneAwareKmerModel(order=args.order, chain=tag).fit(
            sub["seq"].tolist(), v_genes=sub["v"].tolist(), j_genes=sub["j"].tolist())
        path = os.path.join(out_dir, f"kmer_ppost_{tag}.npz")
        model.save(path)
        print(f"  wrote {path} (n_train={model.n_train})")


if __name__ == "__main__":
    main()
