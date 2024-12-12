import pytest
from bio_functions import analyze_rna_secondary_structure

@pytest.fixture
def test_analyze_rna_secondary_structure_typical():
    """Test a typical RNA sequence for secondary structure analysis."""
    result = analyze_rna_secondary_structure("AUGCGCUAUGCUAGC")
    expected_output = {"hairpin_count": 2, "accessibility_score": 8}
    assert result == expected_output

def test_analyze_rna_secondary_structure_empty_sequence():
    """Test edge case where RNA sequence is empty."""
    with pytest.raises(ValueError, match="RNA sequence cannot be empty."):
        analyze_rna_secondary_structure("")

def test_analyze_rna_secondary_structure_complex_sequence():
    """Test a complex RNA sequence with many potential hairpins."""
    result = analyze_rna_secondary_structure("AUGCUUCGGAUCCGAUUGCUA")
    assert result["hairpin_count"] > 2

