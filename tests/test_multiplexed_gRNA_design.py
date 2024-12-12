import pytest
from bio_functions import design_multiplexed_gRNAs

@pytest.fixture
def test_design_multiplexed_gRNAs_typical():
    """Test a typical case with multiple PAMs and gene sequences."""
    result = design_multiplexed_gRNAs(
        ["NGG", "TTTV"], 
        ["ATGCGTACGTAGCTAGCTAGNGG", "TTTACGTAGCTTTTV"]
    )
    expected_output = [
        {"protospacer": "ATGCGTACGTAGCTA", "tracer_rna": "GTTTTAGAGCTAGAA"},
        {"protospacer": "TTTACGTAGCTTTT", "tracer_rna": "GTTTTAGAGCTAA"}
    ]
    assert result == expected_output

def test_design_multiplexed_gRNAs_mismatch_lengths():
    """Test edge case where PAM and gene lists have mismatched lengths."""
    with pytest.raises(ValueError, match="Number of PAM sequences must match number of gene sequences."):
        design_multiplexed_gRNAs(["NGG"], ["ATGCGTACGTANGGC", "TTTACGTANGGC"])
