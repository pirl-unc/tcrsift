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


def test_min_mito_default_is_no_floor():
    """min_mito_pct is a FLOOR; the default must be 0 so out-of-the-box QC
    doesn't silently discard healthy low-mito cells (#168)."""
    config = TCRsiftConfig()
    assert config.load.min_mito_pct == 0.0


def test_min_cd3_reads_default_is_opt_in():
    """The CD3 gate must default to 0 (no gate) so it's a deliberate opt-in
    rather than a silent depth-correlated cull (#172)."""
    config = TCRsiftConfig()
    assert config.phenotype.min_cd3_reads == 0
