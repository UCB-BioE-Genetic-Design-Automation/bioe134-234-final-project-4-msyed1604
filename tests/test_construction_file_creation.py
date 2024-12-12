import pytest
from bio_functions import create_construction_file

@pytest.fixture
def test_create_construction_file_typical():
    """Test a typical case for generating a construction file."""
    result = create_construction_file("Drosophila melanogaster", "ATGCGTACGTAGCTAGCTAGNGG")
    expected_output = (
        "CRISPR Knockout Construction File\n"
        "Organism: Drosophila melanogaster\n"
        "Gene Sequence: ATGCGTACGTAGCTAGCTAGNGG\n"
        "Selected Toolkit: Cas9 with PAM NGG\n"
        "Protospacer: ATGCGTACGTAGCTA\n"
        "Tracer RNA: GTTTTAGAGCTAGAA\n"
        "Steps:\n"
        "1. Use the protospacer and tracer RNA to assemble gRNA.\n"
        "2. Deliver gRNA and Cas9 into insect tissues.\n"
        "3. Validate knockout efficiency through sequencing.\n"
    )
    assert result == expected_output

def test_create_construction_file_invalid_organism():
    """Test edge case where organism is not supported."""
    with pytest.raises(ValueError, match="No toolkit available for Unknown species."):
        create_construction_file("Unknown species", "ATGCGTACGTAGCTAGCTAGNGG")

def test_create_construction_file_empty_gene_sequence():
    """Test edge case where gene sequence is empty."""
    with pytest.raises(ValueError, match="Gene sequence cannot be empty."):
        create_construction_file("Drosophila melanogaster", "")
