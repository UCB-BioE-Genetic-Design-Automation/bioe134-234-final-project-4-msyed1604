import pytest
from bio_functions import create_cas13_construction_file

@pytest.fixture
def test_typical_case_tribolium():
    """Test a typical case for Tribolium castaneum with a valid RNA sequence."""
    result = create_cas13_construction_file(
        organism_name="Tribolium castaneum",
        rna_sequence="GCUAUGCUAGCGUACG"
    )
    expected_output = (
        "CRISPR/Cas13 Construction File\n"
        "Organism: Tribolium castaneum\n"
        "RNA Sequence: GCUAUGCUAGCGUACG\n"
        "Selected Cas13 Variant: CasFB8\n"
        "Efficiency: 70%\n"
        "Hairpin Count: 3\n"
        "Accessibility Score: 7\n"
        "Steps:\n"
        "1. Assemble crRNA with Cas13 (CasFB8).\n"
        "2. Deliver complex into insect tissues.\n"
        "3. Validate transcript knockdown via qPCR.\n"
    )
    assert result == expected_output, f"Expected output did not match. Got: {result}"

def test_edge_case_empty_rna_sequence():
    """Test edge case where RNA sequence is empty."""
    with pytest.raises(ValueError, match="RNA sequence cannot be empty"):
        create_cas13_construction_file(
            organism_name="Drosophila melanogaster",
            rna_sequence=""
        )

# Additional test cases for robustness
def test_edge_case_invalid_organism():
    """Test edge case where organism is not recognized."""
    with pytest.raises(ValueError, match="Organism not supported"):
        create_cas13_construction_file(
            organism_name="Unknown species",
            rna_sequence="AUGCGCUAUGCUAGC"
        )

def test_edge_case_long_rna_sequence():
    """Test edge case with an unusually long RNA sequence."""
    long_rna_sequence = "A" * 10000
    result = create_cas13_construction_file(
        organism_name="Drosophila melanogaster",
        rna_sequence=long_rna_sequence
    )
    assert isinstance(result, str), "Function should return a string even for long sequences."
