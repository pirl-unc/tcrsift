# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The shipped Ppost models are reproducible from the bundled source CSV.

Locks provenance: `observed_pbmc_10x.csv.gz` (pooled 10x healthy-PBMC VDJ-T,
see refseqs/PROVENANCE.md) must regenerate the committed kmer_ppost_*.npz
byte-for-byte, so the background can never silently drift from its documented,
single, both-chains-co-sourced origin.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import tcrsift
from tcrsift.seqprob import GeneAwareKmerModel

_REFDIR = Path(tcrsift.__file__).parent / "refseqs"


def test_bundled_source_exists_and_is_paired():
    src = pd.read_csv(_REFDIR / "observed_pbmc_10x.csv.gz")
    assert set(src.columns) == {"chain", "seq", "v", "j"}
    counts = src["chain"].value_counts().to_dict()
    # Both chains from the one bundled source.
    assert counts["TRA"] == 15948
    assert counts["TRB"] == 18110


def test_shipped_ppost_models_reproduce_from_bundled_source():
    src = pd.read_csv(_REFDIR / "observed_pbmc_10x.csv.gz")
    for chain, tag in (("TRA", "alpha"), ("TRB", "beta")):
        sub = src[src["chain"] == chain]
        model = GeneAwareKmerModel(order=2, chain=tag).fit(
            sub["seq"].tolist(), v_genes=sub["v"].tolist(), j_genes=sub["j"].tolist()
        )
        shipped = np.load(_REFDIR / f"kmer_ppost_{tag}.npz", allow_pickle=True)
        assert int(shipped["n_train"]) == model.n_train
        assert np.allclose(shipped["logp"], model._cdr3._logp, equal_nan=True), (
            f"{tag}: shipped Ppost model no longer matches the bundled source CSV"
        )
