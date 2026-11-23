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

"""Tests for v3.31.0 architectural streamlining and VDJdb 2026+ relational support."""

from __future__ import annotations

import gzip
import io
import json
import zipfile
from dataclasses import fields
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

import tcrsift
from tcrsift import (
    PipelineResult,
    PrioritizeConfig,
    TCRsiftConfig,
    run_pipeline,
    run_til_prioritize,
)
from tcrsift._schema import (
    clean_cdr3_series,
    clean_str_series,
    paired_cdr3ab_key,
    truthy_flag_series,
    valid_clone_mask,
)
from tcrsift.annotate import VIRAL_SPECIES_PATTERNS, load_vdjdb, match_clonotypes
from tcrsift.annotate_tcrs import naive_signature
from tcrsift.cohort import read_donor_clones
from tcrsift.config import _FLAT_TO_NESTED, _SECTION_CLASSES
from tcrsift.datacache import cached_path, download_database
from tcrsift.sample_sheet import Sample, SampleSheet
from tcrsift.signature_methods import SIGNATURES
from tcrsift.signatures import (
    ANTIGEN_RESPONSE_GENES_HGNC,
    CYTOLYTIC_GENES_HGNC,
    EFFECTOR_GENES_HGNC,
    NAIVE_STEM_GENES_HGNC,
)
from tcrsift.til_select import VIRAL_SPECIES_PATTERNS as TIL_SELECT_VIRAL_PATTERNS


def _write_synthetic_sample(sample_dir: Path, barcodes: list[str], cdr3s: list[tuple[str, str]]) -> None:
    """Create minimal 10x GEX Matrix Market + VDJ CSV sample directory."""
    mtx_dir = sample_dir / "filtered_feature_bc_matrix"
    mtx_dir.mkdir(parents=True, exist_ok=True)
    genes = [
        ("ENSG00000167286", "CD3D", "Gene Expression"),
        ("ENSG00000198851", "CD3E", "Gene Expression"),
        ("ENSG00000153563", "CD8A", "Gene Expression"),
        ("ENSG00000172116", "CD8B", "Gene Expression"),
        ("ENSG00000010610", "CD4", "Gene Expression"),
        ("ENSG00000111537", "IFNG", "Gene Expression"),
        ("ENSG00000100453", "GZMB", "Gene Expression"),
        ("ENSG00000180644", "PRF1", "Gene Expression"),
        ("ENSG00000115523", "GNLY", "Gene Expression"),
        ("ENSG00000105374", "NKG7", "Gene Expression"),
        ("ENSG00000081059", "TCF7", "Gene Expression"),
        ("ENSG00000138795", "LEF1", "Gene Expression"),
        ("ENSG00000126353", "CCR7", "Gene Expression"),
        ("ENSG00000188404", "SELL", "Gene Expression"),
        ("ENSG00000168685", "IL7R", "Gene Expression"),
        ("ENSG00000143297", "GZMK", "Gene Expression"),
        ("ENSG00000049249", "TNFRSF9", "Gene Expression"),
        ("ENSG00000186827", "TNFRSF4", "Gene Expression"),
        ("ENSG00000148773", "MKI67", "Gene Expression"),
        ("ENSG00000198888", "MT-ND1", "Gene Expression"),
    ]
    rng = np.random.default_rng(42)
    n_cells = len(barcodes)
    X = rng.integers(5, 40, size=(n_cells, len(genes)))
    X[:, 0] = 15  # CD3D
    X[:, 1] = 15  # CD3E
    X[:, 2] = 25  # CD8A
    X[:, 3] = 20  # CD8B
    X[:, 4] = 0   # CD4
    X[:, -1] = 15  # MT-ND1 (~4% mito, >250 counts)
    # Pad extra housekeeping genes so min_genes=10 passes easily
    extra_genes = [(f"ENSG99999{i:05d}", f"GENE{i}", "Gene Expression") for i in range(15)]
    X_extra = rng.integers(10, 30, size=(n_cells, len(extra_genes)))
    X_all = np.hstack([X, X_extra])
    all_genes = genes + extra_genes

    mtx_tmp = mtx_dir / "matrix.mtx"
    mmwrite(mtx_tmp, sp.coo_matrix(X_all.T))
    with open(mtx_tmp, "rb") as fh, gzip.open(mtx_dir / "matrix.mtx.gz", "wb") as gz:
        gz.write(fh.read())
    mtx_tmp.unlink()
    with gzip.open(mtx_dir / "barcodes.tsv.gz", "wt") as gz:
        gz.write("\n".join(barcodes) + "\n")
    with gzip.open(mtx_dir / "features.tsv.gz", "wt") as gz:
        for gid, gsym, gtype in all_genes:
            gz.write(f"{gid}\t{gsym}\t{gtype}\n")

    rows = []
    for i, (bc, (ca, cb)) in enumerate(zip(barcodes, cdr3s)):
        rows.append(
            {
                "barcode": bc,
                "is_cell": True,
                "contig_id": f"{bc}_contig_1",
                "high_confidence": True,
                "length": 500,
                "chain": "TRA",
                "v_gene": "TRAV12-2",
                "d_gene": "None",
                "j_gene": "TRAJ24",
                "c_gene": "TRAC",
                "full_length": True,
                "productive": True,
                "cdr3": ca,
                "cdr3_nt": "TGT" + "GCC" * (len(ca) - 2) + "TTT",
                "reads": 100,
                "umis": 10,
            }
        )
        rows.append(
            {
                "barcode": bc,
                "is_cell": True,
                "contig_id": f"{bc}_contig_2",
                "high_confidence": True,
                "length": 520,
                "chain": "TRB",
                "v_gene": "TRBV19",
                "d_gene": "TRBD1",
                "j_gene": "TRBJ2-7",
                "c_gene": "TRBC2",
                "full_length": True,
                "productive": True,
                "cdr3": cb,
                "cdr3_nt": "TGT" + "GCC" * (len(cb) - 2) + "TTT",
                "reads": 120,
                "umis": 12,
            }
        )
    pd.DataFrame(rows).to_csv(sample_dir / "filtered_contig_annotations.csv", index=False)


def test_schema_helpers_across_object_categorical_and_stringdtype():
    raw_vals = ["CASSIRSSYEQYF", "  CAVRDSNYQLIW  ", "", "_", "nan", "None", "<NA>", None]
    for dtype in (object, "category", "string"):
        s = pd.Series(raw_vals, dtype=dtype)
        cleaned = clean_cdr3_series(s)
        assert cleaned.dtype == object
        assert list(cleaned) == ["CASSIRSSYEQYF", "CAVRDSNYQLIW", "", "", "", "", "", ""]

    df = pd.DataFrame(
        {
            "CDR3_alpha": pd.Series(["CAVRDSNYQLIW", None, "CAVR", "_"], dtype="string"),
            "CDR3_beta": pd.Series(["CASSIRSSYEQYF", "CASSIRSSYEQYF", "nan", "CASSL"], dtype="category"),
        }
    )
    key, valid = paired_cdr3ab_key(df)
    assert list(valid) == [True, False, False, False]
    assert list(key) == ["CAVRDSNYQLIW_CASSIRSSYEQYF", "", "", ""]

    clones = pd.Series(
        ["CAVRDSNYQLIW_CASSIRSSYEQYF", "_CASSIRSSYEQYF", "CAVRDSNYQLIW_", "_", "None", None],
        dtype="string",
    )
    assert list(valid_clone_mask(clones, clone_col="CDR3ab")) == [True, False, False, False, False, False]
    assert list(valid_clone_mask(pd.Series(["CASSIRSSYEQYF", "_", None]), clone_col="CDR3_beta")) == [
        True,
        False,
        False,
    ]
    assert list(clean_str_series(pd.Series(["  foo ", "NaN", None], dtype="string"))) == ["foo", "", ""]

    flags = pd.Series(["True", "true", "1", "yes", "False", "None", "", None], dtype="category")
    assert list(truthy_flag_series(flags)) == [True, True, True, True, False, False, False, False]
    float_flags = pd.Series([1.0, 0.0, None, "1.0", "0.0"], dtype=object)
    assert list(truthy_flag_series(float_flags)) == [True, False, False, True, False]


def test_cohort_read_donor_clones_rejects_partial_single_chain_keys(tmp_path):
    donor_dir = tmp_path / "donorA" / "data"
    donor_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "CDR3ab": [
                "CAVRDSNYQLIW_CASSIRSSYEQYF",
                "_CASSIRSSYEQYF",
                "CAVRDSNYQLIW_",
                "_",
                "None",
            ],
            "cell_count": [5, 3, 2, 1, 1],
        }
    ).to_csv(donor_dir / "clonotypes.csv", index=False)

    out = read_donor_clones(tmp_path / "donorA")
    assert list(out["CDR3ab"]) == ["CAVRDSNYQLIW_CASSIRSSYEQYF"]
    assert list(out["CDR3_beta"]) == ["CASSIRSSYEQYF"]
    assert list(out["cells"]) == [5.0]


def test_vdjdb_2026_relational_release_download_and_load(tmp_path):
    """VDJdb 2026-11+ ships `vdjdb-YYYY-MM-DD/{chains.tsv.gz,records.tsv.gz}` without `vdjdb.txt`."""
    chains_df = pd.DataFrame(
        [
            {
                "record_id": "VDJDB0001",
                "gene": "TRA",
                "cdr3": "CAGAGSQGNLIF",
                "v_segm": "TRAV27*01",
                "j_segm": "TRAJ42*01",
            },
            {
                "record_id": "VDJDB0001",
                "gene": "TRB",
                "cdr3": "CASSIRSSYEQYF",
                "v_segm": "TRBV19*01",
                "j_segm": "TRBJ2-7*01",
            },
            {
                "record_id": "VDJDB0002",
                "gene": "TRB",
                "cdr3": "CASSLAPGATNEKLFF",
                "v_segm": "TRBV7-6*01",
                "j_segm": "TRBJ1-4*01",
            },
        ]
    )
    records_df = pd.DataFrame(
        [
            {
                "record_id": "VDJDB0001",
                "species": "HomoSapiens",
                "mhc_a": "HLA-A*02:01",
                "mhc_b": "B2M",
                "mhc_class": "MHCI",
                "antigen_epitope": "GILGFVFTL",
                "antigen_gene": "M",
                "antigen_species": "InfluenzaA",
                "reference_id": "PMID:28636592",
                "vdjdb_score": 3,
            },
            {
                "record_id": "VDJDB0002",
                "species": "HomoSapiens",
                "mhc_a": "HLA-B*07:02",
                "mhc_b": "B2M",
                "mhc_class": "MHCI",
                "antigen_epitope": "NLVPMVATV",
                "antigen_gene": "pp65",
                "antigen_species": "CMV",
                "reference_id": "PMID:12345678",
                "vdjdb_score": 2,
            },
        ]
    )

    def _gz_tsv(df: pd.DataFrame) -> bytes:
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            gz.write(df.to_csv(sep="\t", index=False).encode("utf-8"))
        return buf.getvalue()

    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w") as zf:
        zf.writestr("vdjdb-2026-11-06/chains.tsv.gz", _gz_tsv(chains_df))
        zf.writestr("vdjdb-2026-11-06/records.tsv.gz", _gz_tsv(records_df))

    api_payload = json.dumps(
        {
            "tag_name": "2026-11-06",
            "assets": [
                {
                    "name": "vdjdb-2026-11-06.zip",
                    "browser_download_url": "https://example.com/vdjdb-2026-11-06.zip",
                }
            ],
        }
    ).encode("utf-8")

    class _Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.close()
            return False

    payloads = {
        "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest": api_payload,
        "https://example.com/vdjdb-2026-11-06.zip": zip_buf.getvalue(),
    }
    resolved = download_database("vdjdb", data_dir=tmp_path, opener=lambda u: _Resp(payloads[u]))
    assert resolved == tmp_path / "vdjdb"
    assert cached_path("vdjdb", data_dir=tmp_path) == tmp_path / "vdjdb"

    db = load_vdjdb(resolved, verbose=False)
    assert len(db) == 2
    assert int(db["cdr3_alpha"].notna().sum()) == 1
    assert set(db["epitope"]) == {"GILGFVFTL", "NLVPMVATV"}
    assert bool(db["is_viral"].all())

    query = pd.DataFrame(
        {
            "CDR3ab": ["CAGAGSQGNLIF_CASSIRSSYEQYF", "CAVRDSNYQLIW_CASSLAPGATNEKLFF"],
            "CDR3_alpha": ["CAGAGSQGNLIF", "CAVRDSNYQLIW"],
            "CDR3_beta": ["CASSIRSSYEQYF", "CASSLAPGATNEKLFF"],
        }
    )
    strict = match_clonotypes(query, db, match_strictness="strict_ab", verbose=False, show_progress=False)
    assert list(strict["db_epitope"]) == ["GILGFVFTL", None]
    partial = match_clonotypes(
        query, db, match_strictness="ab_with_partial", verbose=False, show_progress=False
    )
    assert list(partial["db_epitope"]) == ["GILGFVFTL", "NLVPMVATV"]

    # Direct file path to chains.tsv.gz or records.tsv.gz also pairs sibling metadata automatically
    rel_subdir = resolved / "vdjdb-2026-11-06"
    for direct_file in (rel_subdir / "chains.tsv.gz", rel_subdir / "records.tsv.gz"):
        db_direct = load_vdjdb(direct_file, verbose=False)
        assert set(db_direct["epitope"]) == {"GILGFVFTL", "NLVPMVATV"}
        matched_direct = match_clonotypes(
            query, db_direct, match_strictness="strict_ab", verbose=False, show_progress=False
        )
        assert list(matched_direct["db_epitope"]) == ["GILGFVFTL", None]


def test_config_auto_derives_every_unambiguous_dataclass_field():
    shared = {"enabled", "min_cells", "match_by"}
    for sec_name, sec_cls in _SECTION_CLASSES.items():
        for f in fields(sec_cls):
            if f.name in shared:
                continue
            assert _FLAT_TO_NESTED.get(f.name) == (sec_name, f.name), (
                f"{sec_name}.{f.name} missing from _FLAT_TO_NESTED"
            )
    cfg = TCRsiftConfig._from_dict(
        {
            "enabled": False,
            "min_cells": 6,
            "match_by": "CDR3b_only",
            "til_match_by": "CDR3ab",
            "no_pgen_ppost": True,
            "no_assemble": True,
            "skip_plots": True,
        }
    )
    assert cfg.attribution.enabled is False
    assert cfg.filter.min_cells == 6
    assert cfg.annotate.match_by == "CDR3b_only"
    assert cfg.til.match_by == "CDR3ab"
    assert cfg.clonotype.add_pgen_ppost is False
    assert cfg.assemble.enabled is False
    assert cfg.output.generate_plots is False


def test_signature_and_viral_registries_share_single_source_of_truth():
    assert SIGNATURES["Cytolytic"].genes_up == CYTOLYTIC_GENES_HGNC
    assert SIGNATURES["Differentiated"].genes_up == EFFECTOR_GENES_HGNC
    assert SIGNATURES["Differentiated"].genes_down == NAIVE_STEM_GENES_HGNC[:4]
    assert SIGNATURES["DifferentiatedBroad"].genes_down == NAIVE_STEM_GENES_HGNC
    assert SIGNATURES["AntigenExperienced"].genes_up == EFFECTOR_GENES_HGNC
    assert SIGNATURES["AcuteActivation"].genes_up == ANTIGEN_RESPONSE_GENES_HGNC
    assert naive_signature().genes_up == NAIVE_STEM_GENES_HGNC
    assert TIL_SELECT_VIRAL_PATTERNS is VIRAL_SPECIES_PATTERNS
    assert tcrsift.__version__ == "3.31.0"


def test_run_pipeline_and_prioritize_config_python_apis(tmp_path):
    ctrl_dir = tmp_path / "ctrl"
    stim_dir = tmp_path / "stim"
    ctrl_bcs = [f"AAACCTGC{i:04d}-1" for i in range(6)]
    stim_bcs = [f"GGGTTTAC{i:04d}-1" for i in range(6)]
    ctrl_cdr3s = [
        ("CAGAGSQGNLIF", "CASSIRSSYEQYF"),
        ("CAVRDSNYQLIW", "CASSLAPGATNEKLFF"),
        ("CAVRDSNYQLIW", "CASSLAPGATNEKLFF"),
        ("CAVSDGQKLLF", "CASSQETQYF"),
        ("CAVNTGNQFYF", "CASSPGQGEGYEQYF"),
        ("CAVRPLYNQGGKLIF", "CASSYVGNTGELFF"),
    ]
    stim_cdr3s = [
        ("CAGAGSQGNLIF", "CASSIRSSYEQYF"),
        ("CAGAGSQGNLIF", "CASSIRSSYEQYF"),
        ("CAGAGSQGNLIF", "CASSIRSSYEQYF"),
        ("CAGAGSQGNLIF", "CASSIRSSYEQYF"),
        ("CAVRDSNYQLIW", "CASSLAPGATNEKLFF"),
        ("CAVRDSNYQLIW", "CASSLAPGATNEKLFF"),
    ]
    _write_synthetic_sample(ctrl_dir, ctrl_bcs, ctrl_cdr3s)
    _write_synthetic_sample(stim_dir, stim_bcs, stim_cdr3s)

    sheet = SampleSheet(
        samples=[
            Sample(
                sample="unstim",
                patient_id="P1",
                vdj_dir=str(ctrl_dir),
                gex_dir=str(ctrl_dir),
                source="control",
            ),
            Sample(
                sample="stim",
                patient_id="P1",
                vdj_dir=str(stim_dir),
                gex_dir=str(stim_dir),
                source="culture",
                control_sample="unstim",
            ),
        ]
    )
    cfg = TCRsiftConfig._from_dict(
        {
            "min_genes": 10,
            "min_counts": 50,
            "min_mito_pct": 0.0,
            "max_mito_pct": 50.0,
            "min_cells": 2,
            "min_frequency": 0.01,
            "skip_plots": True,
            "generate_report": False,
            "no_assemble": True,
        }
    )
    res = run_pipeline(
        sheet,
        tmp_path / "pipeline_out",
        config=cfg,
        cache_dir=tmp_path / "empty_cache",
        no_assemble=True,
    )
    assert isinstance(res, PipelineResult)
    assert len(res.clonotypes) == 2
    assert res["stimulation_response"] is not None
    assert "stim_log2_fc" in res.clonotypes.columns
    assert (tmp_path / "pipeline_out" / "data" / "stimulation_response.csv").is_file()

    # Also test PrioritizeConfig direct Python API without argparse.Namespace
    sheet_csv = tmp_path / "prio_sheet.csv"
    sheet_csv.write_text(
        f"sample,patient_id,tissue,vdj_dir,gex_dir\nstim,P1,blood,{stim_dir},{stim_dir}\n"
    )
    pcfg = PrioritizeConfig(
        sample_sheet=sheet_csv,
        output_dir=tmp_path / "prio_out",
        context="blood",
        min_genes=10,
        min_counts=50,
        max_mito_pct=50.0,
        min_cells=2,
        cache_dir=tmp_path / "empty_cache",
    )
    scored, candidates = run_til_prioritize(pcfg)
    assert len(scored) == 2
    assert isinstance(candidates, pd.DataFrame)
    assert (tmp_path / "prio_out" / "prioritization.json").is_file()

    # Empty string paths normalize to None and None list fields normalize to []
    pcfg_norm = PrioritizeConfig(
        sample_sheet=sheet_csv,
        output_dir=tmp_path / "prio_norm",
        vdjdb="",
        cache_dir="  ",
        min_expression=None,  # type: ignore[arg-type]
        include_v_gene=None,  # type: ignore[arg-type]
    )
    assert pcfg_norm.vdjdb is None
    assert pcfg_norm.cache_dir is None
    assert pcfg_norm.min_expression == []
    assert pcfg_norm.include_v_gene == []

    # Verify CLI `main(["run", ...])` exits cleanly (returns None, does not sys.exit(PipelineResult))
    from tcrsift.cli import main as cli_main

    cli_sheet = tmp_path / "cli_run_sheet.csv"
    cli_sheet.write_text(
        f"sample,patient_id,source,control_sample,vdj_dir,gex_dir\n"
        f"unstim,P1,control,,{ctrl_dir},{ctrl_dir}\n"
        f"stim,P1,culture,unstim,{stim_dir},{stim_dir}\n"
    )
    cfg_yaml = tmp_path / "cli_cfg.yaml"
    cfg.to_yaml(cfg_yaml)
    assert cfg.til.til_samples == []
    cli_rc = cli_main(
        [
            "run",
            "-s",
            str(cli_sheet),
            "-o",
            str(tmp_path / "cli_run_out"),
            "-c",
            str(cfg_yaml),
            "--cache-dir",
            str(tmp_path / "empty_cache"),
            "--no-assemble",
            "--skip-plots",
        ]
    )
    assert cli_rc is None
    assert (tmp_path / "cli_run_out" / "data" / "clonotypes.csv").is_file()
