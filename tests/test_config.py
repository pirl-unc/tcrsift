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
