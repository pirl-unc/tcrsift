"""
Tests for configuration handling.
"""

from tcrsift.config import TCRsiftConfig


def test_sct_flat_keys_are_mapped():
    """Flat SCT keys should map to the sct section."""
    config = TCRsiftConfig._from_dict(
        {
            "min_snr": 5.0,
            "min_reads_per_chain": 20,
            "require_mutation_match": False,
            "require_compact_match": True,
        }
    )

    assert config.sct.min_snr == 5.0
    assert config.sct.min_reads_per_chain == 20
    assert config.sct.require_mutation_match is False
    assert config.sct.require_compact_match is True


def test_attribution_split_flags_inherit_doublet_split():
    """dual_alpha_split / dual_beta_split default to None and inherit
    doublet_split via __post_init__, so old configs are unchanged (#181)."""
    from tcrsift.config import AttributionConfig

    assert AttributionConfig().dual_alpha_split is True
    assert AttributionConfig().dual_beta_split is True
    off = AttributionConfig(doublet_split=False)
    assert off.dual_alpha_split is False and off.dual_beta_split is False
    # Explicit value overrides the inherited default.
    mixed = AttributionConfig(doublet_split=False, dual_alpha_split=True)
    assert mixed.dual_alpha_split is True and mixed.dual_beta_split is False


def test_min_mito_default_retains_floor():
    """min_mito_pct keeps its 2.0 floor default (intentional empty/ambient
    removal); #168's fix is making the drop visible, not removing the floor."""
    config = TCRsiftConfig()
    assert config.load.min_mito_pct == 2.0


def test_min_cd3_reads_default_is_gentle_floor():
    """The now-functional CD3 gate defaults to a gentle floor of 3 (validated on
    pilot data); 0 disables it. Previously a dead parameter (#172)."""
    config = TCRsiftConfig()
    assert config.phenotype.min_cd3_reads == 3


def test_match_by_no_collision_on_merge():
    """annotate.match_by and til.match_by must not clobber each other through
    merge_with_args' flatten/rebuild (they share the field name match_by, F8)."""
    import argparse

    cfg = TCRsiftConfig._from_dict(
        {"annotate": {"match_by": "CDR3ab"}, "til": {"match_by": "CDR3b_only"}}
    )
    merged = cfg.merge_with_args(argparse.Namespace(func=None))
    assert merged.annotate.match_by == "CDR3ab"
    assert merged.til.match_by == "CDR3b_only"


def test_insilico_filter_survives_to_dict_and_merge():
    """insilico_filter must round-trip through to_dict/merge_with_args (F9) so it
    isn't silently dropped before the pipeline runs or when saving config.yaml."""
    import argparse

    spec = {"predicates": [{"score": "ppost_beta", "min_percentile": 50}]}
    cfg = TCRsiftConfig._from_dict({"insilico_filter": spec})
    assert cfg.to_dict()["insilico_filter"] == spec
    merged = cfg.merge_with_args(argparse.Namespace(func=None))
    assert merged.insilico_filter == spec


def test_clonotype_pgen_ppost_default_on():
    """Every run annotates clonotypes with k-mer Pgen/Ppost by default, so
    ppost_alpha/ppost_beta are available for PRISM without a separate step."""
    assert TCRsiftConfig().clonotype.add_pgen_ppost is True


def test_no_pgen_ppost_flag_parses_and_overrides():
    from tcrsift.cli import create_parser

    a = create_parser().parse_args(["run", "-s", "s.csv", "-o", "out", "--no-pgen-ppost"])
    assert a.no_pgen_ppost is True
    # default off (so the k-mer Pgen/Ppost annotation runs) when not passed
    b = create_parser().parse_args(["run", "-s", "s.csv", "-o", "out"])
    assert getattr(b, "no_pgen_ppost", False) is False


def test_baseline_markers_default_is_none():
    """Defaults to None so the format-module default (CTY) is used unchanged."""
    assert TCRsiftConfig().output.baseline_markers is None


def test_baseline_markers_nested_and_flat():
    """output.baseline_markers settable nested or via the flat alias (#208)."""
    nested = TCRsiftConfig._from_dict({"output": {"baseline_markers": ["CTY", "DMSO"]}})
    assert nested.output.baseline_markers == ["CTY", "DMSO"]
    flat = TCRsiftConfig._from_dict({"baseline_markers": ["DMSO"]})
    assert flat.output.baseline_markers == ["DMSO"]


def test_baseline_markers_round_trips_through_dict():
    cfg = TCRsiftConfig._from_dict({"output": {"baseline_markers": ["DMSO"]}})
    assert cfg.to_dict()["output"]["baseline_markers"] == ["DMSO"]


def test_cmd_run_applies_baseline_markers(tmp_path):
    """cmd_run pushes the configured markers into the format module (#208)."""
    import tcrsift.format as fmt
    from tcrsift.cli import create_parser

    original = fmt.BASELINE_MARKERS
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("output:\n  baseline_markers: [DMSO, CTY]\n")
    # An existing (but unusable) sheet passes arg validation, so execution
    # reaches the baseline block before the run bails on the empty sheet.
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("sample,vdj_dir\n")
    try:
        args = create_parser().parse_args(
            ["run", "-s", str(sheet), "-o", str(tmp_path / "out"), "-c", str(cfg)]
        )
        try:
            args.func(args)
        except Exception:
            pass
        assert fmt.BASELINE_MARKERS == ("DMSO", "CTY")
    finally:
        fmt.set_baseline_markers(*original)
