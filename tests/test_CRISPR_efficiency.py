import pytest
from bio_functions import predict_crispr_efficiency

@pytest.fixture
def test_predict_crispr_efficiency_typical():
    """Test a typical case for predicting CRISPR efficiency."""
    result = predict_crispr_efficiency("NGG", "ATGCGTACGTAGCTA")
    assert 80 <= result <= 100

def test_predict_crispr_efficiency_no_pam_match():
    """Test edge case where PAM does not match target sequence."""
    result = predict_crispr_efficiency("TTTV", "ATGCGTACGTAGCTA")
    assert result < 60

def test_predict_crispr_efficiency_high_gc_content():
    """Test case with high GC content in the target sequence."""
    result = predict_crispr_efficiency("NGG", "GGGGGGGGGGGGGGGG")
    assert 90 <= result <= 100

