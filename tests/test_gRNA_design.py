import pytest
from bio_functions import design_gRNA

@pytest.fixture
def test_design_gRNA_typical():
    """Test a typical case with valid PAM and gene sequence."""
    result = design_gRNA("NGG", "ATGCGTACGTAGCTAGCTAGNGG")
    expected_output = {"protospacer": "ATGCGTACGTAGCTA", "tracer_rna": "GTTTTAGAGCTAGAA"}
    assert result == expected_output

def test_design_gRNA_no_pam():
    """Test edge case where PAM sequence is not found in the gene sequence."""
    with pytest.raises(ValueError, match="PAM sequence not found in the gene sequence."):
        design_gRNA("NGG", "ATGCGTACGTAGCTAGCTAG")

def test_design_gRNA_short_gene():
    """Test edge case with a gene sequence too short to contain a valid protospacer."""
    with pytest.raises(ValueError, match="Gene sequence is too short for gRNA design."):
        design_gRNA("NGG", "ATGC")

def test_design_gRNA_multiple_pams():
    """Test case where multiple PAM sequences are present in the gene sequence."""
    result = design_gRNA("NGG", "ATGCGTACGTANGGCTAGNGG")
    expected_output = {"protospacer": "ATGCGTACGTANGG", "tracer_rna": "GTTTTAGAGCTAGAA"}
    assert result == expected_output

