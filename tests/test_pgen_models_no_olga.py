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

"""OLGA/SONIA are fully removed; Pgen ships pre-built and isn't runtime-trainable."""

from __future__ import annotations

import importlib

import pytest

import tcrsift


def test_pgen_training_raises_clearly():
    # Pgen models ship pre-built (OLGA, the synthetic-sequence generator, is no
    # longer a dependency), so a runtime train request must fail with a clear,
    # actionable error rather than ImportError-ing on a missing olga.
    from tcrsift.pgen_models import train_model

    with pytest.raises(ValueError, match="shipped pre-built"):
        train_model("beta", role="pgen")


def test_shipped_pgen_still_loads():
    # Degradation intact: the default k-mer Pgen background still resolves
    # without any OLGA — it's the committed numpy model.
    from tcrsift.pgen_models import ensure_model

    model = ensure_model("beta", backend="kmer", role="pgen", auto_train=False)
    assert model is not None


def test_olga_module_and_exports_gone():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("tcrsift.olga_ppost")
    for name in (
        "compute_pgen_ppost",
        "olga_sonia_available",
        "nearest_supported_allele",
        "supported_alleles",
        "normalize_gene_name",
        "flag_private_candidates",
    ):
        assert name not in tcrsift.__all__, f"{name} should be removed from the public API"
