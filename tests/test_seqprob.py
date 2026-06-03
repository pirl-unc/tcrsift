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

"""Tests for data-driven CDR3 sequence-probability models (seqprob)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift import seqprob as sp


def _toy_repertoire(n=2000, seed=0):
    """Synthetic CDR3-like strings: C…F with a biased middle composition."""
    rng = np.random.default_rng(seed)
    common = list("AGSTLEKVDR")  # over-represented residues
    seqs = []
    for _ in range(n):
        length = int(rng.integers(10, 18))
        mid = "".join(rng.choice(common, size=length - 2))
        seqs.append("C" + mid + "F")
    return seqs


class TestKmerModel:
    def test_fit_and_relative_probability(self):
        seqs = _toy_repertoire()
        m = sp.KmerProbabilityModel(order=3, chain="beta").fit(seqs)
        # an in-distribution sequence beats one full of rare residues
        common = m.log_prob([seqs[0]])[0]
        rare = m.log_prob(["CWWHWHWMPCWF"])[0]
        assert common > rare

    def test_unscorable_inputs_are_nan(self):
        m = sp.KmerProbabilityModel(order=2).fit(_toy_repertoire(500))
        lp = m.log_prob(["", "CASSX1F", None, 42])
        assert np.isnan(lp).all()

    def test_table_shape_matches_order(self):
        m = sp.KmerProbabilityModel(order=3).fit(_toy_repertoire(500))
        assert m._logp.shape == (sp.N_SYM**3, sp.N_SYM)

    def test_save_load_roundtrip(self, tmp_path):
        seqs = _toy_repertoire()
        m = sp.KmerProbabilityModel(order=3, chain="alpha").fit(seqs)
        p = tmp_path / "m.npz"
        m.save(p)
        m2 = sp.KmerProbabilityModel.load(p)
        assert m2.order == 3 and m2.chain == "alpha" and m2.n_train == len(seqs)
        a = m.log_prob(seqs[:20])
        b = m2.log_prob(seqs[:20])
        assert np.allclose(a, b, equal_nan=True)

    def test_order_validation(self):
        with pytest.raises(ValueError, match="order"):
            sp.KmerProbabilityModel(order=0)

    def test_unfitted_raises(self):
        with pytest.raises(RuntimeError, match="not fitted"):
            sp.KmerProbabilityModel().log_prob(["CASSF"])

    def test_fit_all_unscorable_raises(self):
        with pytest.raises(ValueError, match="no scorable"):
            sp.KmerProbabilityModel().fit(["", "X1Z", None])

    def test_probabilities_normalize_per_context(self):
        # Each context row is a proper distribution over the 22 symbols.
        m = sp.KmerProbabilityModel(order=2).fit(_toy_repertoire(800))
        row_sums = np.exp(m._logp).sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=1e-5)


class TestShippedDefaults:
    """Use the committed OLGA-generated k-mer backgrounds (if present)."""

    @pytest.mark.parametrize("chain", ["beta", "alpha"])
    def test_load_and_score(self, chain):
        try:
            model = sp.load_background_model(chain, "kmer")
        except FileNotFoundError:
            pytest.skip(f"shipped {chain} background not built yet")
        # A canonical CDR3 scores finite and beats AA gibberish.
        seq = "CASSLGTGELFF" if chain == "beta" else "CAVRDSNYQLIWF"
        good = model.log_prob([seq])[0]
        bad = model.log_prob(["CWWWWWWWWWWWF"])[0]
        assert np.isfinite(good)
        assert good > bad

    def test_score_log_pgen_dataframe(self):
        try:
            sp.load_background_model("beta", "kmer")
        except FileNotFoundError:
            pytest.skip("shipped beta background not built yet")
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLGTGELFF", "", "CASSPGQGAYEQYF"],
        })
        s = sp.score_log_pgen(df, chain="beta")
        assert s.name == "log_pgen"
        assert np.isfinite(s.iloc[0]) and np.isnan(s.iloc[1])

    def test_missing_column_raises(self):
        with pytest.raises(ValueError, match="missing"):
            sp.score_log_pgen(pd.DataFrame({"x": [1]}), chain="beta",
                              model=sp.KmerProbabilityModel(order=2).fit(
                                  _toy_repertoire(200)))


class TestBackendRegistry:
    def test_registry_has_both_backends(self):
        assert set(sp.BACKENDS) == {"kmer", "tcrpeg"}
        assert sp.BACKENDS["kmer"] is sp.KmerProbabilityModel

    def test_unknown_chain_for_default(self):
        with pytest.raises(ValueError, match="alpha.*beta"):
            sp.load_background_model("gamma", "kmer")


_HAS_TCRPEG = sp.TCRpegProbabilityModel.available()


@pytest.mark.skipif(_HAS_TCRPEG, reason="only when tcrpeg is NOT installed")
def test_tcrpeg_missing_deps_raises():
    with pytest.raises(ImportError, match=r"pip install tcrsift\[tcrpeg\]"):
        sp.TCRpegProbabilityModel().fit(["CASSF"])


@pytest.mark.skipif(not _HAS_TCRPEG, reason="needs the [tcrpeg] extra")
class TestTCRpegBackend:
    def test_fit_and_score_smoke(self):
        seqs = _toy_repertoire(400)
        m = sp.TCRpegProbabilityModel(epochs=1, batch_size=64, chain="beta")
        m.fit(seqs)
        lp = m.log_prob([seqs[0], "", "CASSX1F"])
        assert np.isfinite(lp[0])
        assert np.isnan(lp[1]) and np.isnan(lp[2])
