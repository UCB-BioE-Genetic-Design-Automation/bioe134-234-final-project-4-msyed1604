import pytest
from bio_functions import select_toolkit

@pytest.fixture
def test_select_toolkit_typical():
    """Test a typical case for selecting a toolkit based on organism name."""
    result = select_toolkit("Drosophila melanogaster")
    expected_output = {"pam_sequence": "NGG", "system": "Cas9"}
    assert result == expected_output

def test_select_toolkit_invalid_organism():
    """Test edge case where organism is not in the database."""
    with pytest.raises(ValueError, match="No toolkit available for Unknown species."):
        select_toolkit("Unknown species")

def test_select_toolkit_case_insensitivity():
    """Test case insensitivity for organism name input."""
    result = select_toolkit("drosophila melanogaster")
    expected_output = {"pam_sequence": "NGG", "system": "Cas9"}
    assert result == expected_output
