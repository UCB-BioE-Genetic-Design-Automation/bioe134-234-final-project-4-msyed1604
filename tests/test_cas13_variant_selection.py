import pytest
from bio_functions import select_cas13_variant

@pytest.fixture
def test_select_cas13_variant_typical():
    """Test a typical case for selecting a Cas13 variant based on organism name."""
    result = select_cas13_variant("Drosophila melanogaster")
    expected_output = {"variant": "CasFX4", "efficiency": 90, "specificity": "high"}
    assert result == expected_output

def test_select_cas13_variant_invalid_organism():
    """Test edge case where organism is not recognized."""
    with pytest.raises(ValueError, match="No Cas13 variant data available for Unknown species."):
        select_cas13_variant("Unknown species")

def test_select_cas13_variant_case_insensitivity():
    """Test case insensitivity for organism name input."""
    result = select_cas13_variant("drosophila melanogaster")
    expected_output = {"variant": "CasFX4", "efficiency": 90, "specificity": "high"}
    assert result == expected_output

def test_select_cas13_variant_edge_case():
    """Test edge case with a valid but less common organism."""
    result = select_cas13_variant("Tribolium castaneum")
    expected_output = {"variant": "CasFB8", "efficiency": 70, "specificity": "moderate"}
    assert result == expected_output
