"""Regression tests for v3.38.0 dead-code removal, helper deduplication, and til-select cleanup."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from tcrsift.assemble.codons import _revcomp
from tcrsift.assemble.germline import _v_func_rank
from tcrsift.assemble.qc import _restriction_hits
from tcrsift.assemble.vregion import _func_rank_key
from tcrsift.cli import create_parser
from tcrsift.mnemonic import tcr_name
from tcrsift.pgen import _j_usage_log
from tcrsift.til_prioritize import _ensure_clone_key


def test_unified_assemble_helpers():
    # _restriction_hits uses codons._revcomp to detect reverse-strand Type IIS sites
    bsai_fwd = "AAAGGTC" + "TCAAA"  # contains GGTCTC
    bsai_rev = _revcomp("GGTCTC")
    assert "BsaI" in _restriction_hits(bsai_fwd)
    assert "BsaI" in _restriction_hits(f"AAA{bsai_rev}AAA")

    # _v_func_rank is exact complement of _func_rank_key
    for func in ("F", "[F]", "(F)", "ORF", "P", ""):
        assert _v_func_rank(func) == 3 - _func_rank_key(func)
    assert _v_func_rank("F") > _v_func_rank("[F]") > _v_func_rank("ORF") > _v_func_rank("P")


def test_til_prioritize_ensure_clone_key_uses_paired_cdr3ab_key():
    adata = ad.AnnData(
        X=np.zeros((3, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "CDR3_alpha": ["CAVRDSNYQLIW", "", np.nan],
                "CDR3_beta": ["CASSLGGNTEVFF", "CASSF", "CASSG"],
            },
            index=["c0", "c1", "c2"],
        ),
        var=pd.DataFrame(index=["G1"]),
    )
    _ensure_clone_key(adata)
    assert adata.obs.loc["c0", "CDR3ab"] == "CAVRDSNYQLIW_CASSLGGNTEVFF"
    assert pd.isna(adata.obs.loc["c1", "CDR3ab"])
    assert pd.isna(adata.obs.loc["c2", "CDR3ab"])


def test_mnemonic_handles_non_string_and_mouse_j_usage_log():
    assert tcr_name(None) == "Anon"  # type: ignore[arg-type]
    assert tcr_name("   ") == "Anon"
    assert isinstance(tcr_name("CASSLGQAYEQYF"), str)
    # Mouse alpha J usage log works cleanly without dead TRAV alias lookup
    val = _j_usage_log("TRAJ18*01", chain="alpha", species="mouse")
    assert np.isfinite(val) and val < 0


def test_til_select_legacy_flags_still_accepted_but_suppressed_from_help():
    parser = create_parser()
    # Still accepted for backward compatibility
    args = parser.parse_args(
        [
            "til-select",
            "--data-dir",
            ".",
            "--trend-increase-ratio-min",
            "2.0",
            "--trend-decrease-ratio-max",
            "0.4",
            "--pyensembl-release",
            "110",
        ]
    )
    assert args.trend_increase_ratio_min == 2.0
    assert args.trend_decrease_ratio_max == 0.4
    assert args.pyensembl_release == 110

    # Suppressed from subcommand --help
    sub_actions = [a for a in parser._actions if hasattr(a, "choices") and isinstance(a.choices, dict)]
    til_select_parser = sub_actions[0].choices["til-select"]
    help_text = til_select_parser.format_help()
    assert "--trend-increase-ratio-min" not in help_text
    assert "--trend-decrease-ratio-max" not in help_text
    assert "--pyensembl-release" not in help_text
