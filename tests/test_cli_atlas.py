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

"""CLI subcommands for the single-cell atlas path (#342).

`tcrsift atlas qc | embed | annotate` wrap cell_qc_funnel / embed_cells /
annotate_cells so the GEX atlas has CLI parity with the TCR-selection path.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift.cli import (
    _UNSET_CLI,
    cmd_atlas_annotate,
    cmd_atlas_qc,
    create_parser,
)


def _raw_adata(n=60):
    """Tiny raw-count AnnData: a CD8 T cluster and a B cluster over bulk filler
    (so CP10K normalization behaves), plus MT genes for the mito gate."""
    rng = np.random.default_rng(0)
    markers = ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMB", "NKG7",
               "MS4A1", "CD79A", "CD79B", "CD19", "BANK1"]
    # >=~150 genes so scanpy score_genes has enough control genes to bin (#312).
    filler = [f"G{i}" for i in range(200)]
    genes = markers + ["MT-CO1", "MT-ND1"] + filler
    gi = {g: i for i, g in enumerate(genes)}
    X = np.zeros((n, len(genes)))
    for f in filler:
        X[:, gi[f]] = rng.poisson(25, size=n)
    half = n // 2
    for g in ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMB", "NKG7"]:
        X[:half, gi[g]] = 40
    for g in ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"]:
        X[half:, gi[g]] = 40
    X[:, gi["MT-CO1"]] = 3
    X[:, gi["MT-ND1"]] = 3
    return ad.AnnData(
        X=X, var=pd.DataFrame(index=genes),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
    )


def _clustered_adata(n=60):
    adata = _raw_adata(n)
    half = n // 2
    adata.obs["leiden"] = pd.Categorical(["0"] * half + ["1"] * (n - half))
    return adata


class TestAtlasParser:
    def test_qc_wires_command_and_defaults(self):
        args = create_parser().parse_args(["atlas", "qc", "-i", "a.h5ad", "-o", "b.h5ad"])
        assert args.func is cmd_atlas_qc
        # Unpassed thresholds stay the sentinel → funnel uses its own defaults.
        assert args.min_genes is _UNSET_CLI
        assert args.max_genes is _UNSET_CLI
        assert args.pbmc_preset is False and args.solid_tumor is False

    def test_qc_threshold_number_and_none(self):
        args = create_parser().parse_args(
            ["atlas", "qc", "-i", "a", "-o", "b",
             "--min-genes", "300", "--max-genes", "none", "--max-mito", "12.5"]
        )
        assert args.min_genes == 300 and isinstance(args.min_genes, int)
        assert args.max_genes is None          # 'none' disables the bound
        assert args.max_mito == 12.5

    def test_qc_rejects_non_numeric_threshold(self):
        with pytest.raises(SystemExit):
            create_parser().parse_args(
                ["atlas", "qc", "-i", "a", "-o", "b", "--min-genes", "lots"]
            )

    def test_embed_wires_command(self):
        args = create_parser().parse_args(
            ["atlas", "embed", "-i", "a", "-o", "b", "--batch-key", "sample",
             "--resolution", "0.8"]
        )
        assert args.func.__name__ == "cmd_atlas_embed"
        assert args.batch_key == "sample" and args.resolution == 0.8

    def test_annotate_wires_command(self):
        args = create_parser().parse_args(
            ["atlas", "annotate", "-i", "a", "-o", "b", "--solid-tumor", "--no-suffix"]
        )
        assert args.func is cmd_atlas_annotate
        assert args.solid_tumor is True and args.no_suffix is True

    def test_atlas_requires_subcommand(self):
        with pytest.raises(SystemExit):
            create_parser().parse_args(["atlas"])


class TestAtlasQcCommand:
    def _run(self, tmp_path, extra=None):
        inp = tmp_path / "raw.h5ad"
        out = tmp_path / "qc.h5ad"
        _raw_adata().write_h5ad(inp)
        argv = ["atlas", "qc", "-i", str(inp), "-o", str(out),
                "--min-genes", "3", "--min-counts", "5"] + (extra or [])
        args = create_parser().parse_args(argv)
        args.func(args)
        return out

    def test_writes_filtered_output(self, tmp_path):
        out = self._run(tmp_path)
        assert out.exists()
        res = ad.read_h5ad(out)
        assert res.n_obs > 0
        assert "pct_counts_mt" in res.obs.columns

    def test_waterfall_csv_written(self, tmp_path):
        wf = tmp_path / "waterfall.csv"
        self._run(tmp_path, extra=["--waterfall-csv", str(wf)])
        assert wf.exists()
        df = pd.read_csv(wf)
        assert list(df.columns) == ["step", "reason", "removed", "remaining"]
        assert (df["step"] == "input").any()

    def test_none_disables_a_bound(self, tmp_path):
        # --max-genes none must parse and run (bound disabled, no crash).
        out = self._run(tmp_path, extra=["--max-genes", "none"])
        assert ad.read_h5ad(out).n_obs > 0


class TestAtlasAnnotateCommand:
    def test_annotates_and_writes_phenotype(self, tmp_path):
        inp = tmp_path / "clustered.h5ad"
        out = tmp_path / "annotated.h5ad"
        _clustered_adata().write_h5ad(inp)
        args = create_parser().parse_args(
            ["atlas", "annotate", "-i", str(inp), "-o", str(out)]
        )
        args.func(args)
        res = ad.read_h5ad(out)
        assert "phenotype" in res.obs.columns
        assert "phenotype_base" in res.obs.columns
        assert res.obs["phenotype"].notna().all()

    def test_missing_leiden_raises(self, tmp_path):
        inp = tmp_path / "noleiden.h5ad"
        out = tmp_path / "x.h5ad"
        _raw_adata().write_h5ad(inp)  # no leiden column
        args = create_parser().parse_args(
            ["atlas", "annotate", "-i", str(inp), "-o", str(out)]
        )
        with pytest.raises(ValueError, match="leiden"):
            args.func(args)


class TestAtlasEmbedCommand:
    def test_embed_produces_leiden(self, tmp_path):
        pytest.importorskip("igraph")
        pytest.importorskip("leidenalg")
        from tcrsift.cli import cmd_atlas_embed  # noqa: F401

        inp = tmp_path / "raw.h5ad"
        out = tmp_path / "embed.h5ad"
        _raw_adata(120).write_h5ad(inp)
        args = create_parser().parse_args(
            ["atlas", "embed", "-i", str(inp), "-o", str(out),
             "--n-pcs", "10", "--n-neighbors", "10", "--n-top-genes", "30"]
        )
        args.func(args)
        res = ad.read_h5ad(out)
        assert "leiden" in res.obs.columns
        assert "X_umap" in res.obsm
