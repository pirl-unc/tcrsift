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
