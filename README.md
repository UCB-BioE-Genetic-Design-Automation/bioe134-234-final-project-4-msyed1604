
# BioE 134 Final Project Submission

## Project Overview

This project focuses on developing a Python-based toolkit for CRISPR-mediated genome editing in insects. The toolkit is designed to facilitate gene knockouts in various insect species using different CRISPR/Cas systems. It involves designing guide RNAs (gRNAs), constructing insect-specific CRISPR vectors, and simulating transformation processes to deliver these vectors into insect tissues. The project aims to make the process of insect gene editing more accessible by providing a streamlined and reusable codebase.

---

## Scope of Work

As part of the final project for BioE 134, I developed several functions that are foundational for CRISPR-mediated genome editing:
1. gRNA Design: This function generates guide RNAs by taking a PAM sequence and a gene sequence as inputs. It outputs the protospacer and tracer RNA necessary for CRISPR targeting.
2. Toolkit Selection: This function selects the appropriate CRISPR toolkit based on the organism name, providing details such as PAM sequence and other properties.
3. Construction File Creation: This function generates a construction file outlining the steps for performing a gene knockout, integrating gRNA design and toolkit selection.
4. RNA Secondary Structure Analysis: This function evaluates RNA secondary structure to assess accessibility for gRNA binding.
5. CRISPR Efficiency Prediction: This function predicts the efficiency of a CRISPR gRNA design based on PAM sequence compatibility and target site properties.
6. Multiplexed gRNA Design: This function designs multiple gRNAs for simultaneous targeting of multiple genes.
7. Generate Construction File with Cas13 Integration: This function generates a construction file for RNA targeting using Cas13, integrating organism-specific Cas13 variant selection and crRNA design.

Each function includes input validation and error handling to ensure proper use, raising errors for invalid inputs or conditions that cannot be processed.

---

## Function Descriptions

### 1. gRNA Design (`design_gRNA`)

- **Description**: Generates guide RNAs by taking a PAM sequence and a gene sequence as inputs.
- **Input**:  PAM sequence (string), Gene sequence (string).
- **Output**: Dictionary containing protospacer and tracer RNA sequences.

**Example**:
```python
design_gRNA("NGG", "ATGCGTACGTAGCTAGCTAGNGG")
# Returns: {"protospacer": "ATGCGTACGTAGCTA", "tracer_rna": "GTTTTAGAGCTAGAA"}
```

### 2.  Toolkit Selection (`select_toolkit`)

- **Description**: Selects the appropriate CRISPR toolkit based on the organism name.
- **Input**: Organism name as represented by a string.
- **Output**: Dictionary containing PAM sequence and CRISPR system details.

**Example**:
```python
select_toolkit("Drosophila melanogaster")
# Returns: {"pam_sequence": "NGG", "system": "Cas9"}
```

### 3.  Construction File Creation (`create_construction_file`)

- **Description**: Generates a construction file outlining steps for gene knockout.
- **Input**: Organism name (string), Gene sequence (string).
- **Output**:  String detailing construction steps.

**Example**:
```python
create_construction_file("Drosophila melanogaster", "ATGCGTACGTAGCTAGCTAGNGG")
# Returns detailed construction file content
```

### 4.  RNA Secondary Structure Analysis (`analyze_rna_secondary_structure`)

- **Description**: Analyzes RNA secondary structure to assess accessibility for gRNA binding.
- **Input**: RNA sequence as represented by a string.
- **Output**: Dictionary containing hairpin count and accessibility score.

**Example**:
```python
analyze_rna_secondary_structure("AUGCGCUAUGCUAGC")
# Returns: {"hairpin_count": 2, "accessibility_score": 8}
```

### 5.  CRISPR Efficiency Prediction (`predict_crispr_efficiency`)

- **Description**: Predicts efficiency of a CRISPR gRNA design.
- **Input**: PAM sequence (string), Target sequence (string).
- **Output**: Efficiency score (float).

**Example**:
```python
predict_crispr_efficiency("NGG", "ATGCGTACGTAGCTA")
# Returns efficiency score
```

### 6.  Multiplexed gRNA Design (`design_multiplexed_gRNAs`)

- **Description**: Designs multiple gRNAs for simultaneous targeting.
- **Input**: List of PAM sequences, List of gene sequences.
- **Output**: List of dictionaries with protospacer and tracer RNA sequences.

**Example**:
```python
design_multiplexed_gRNAs(["NGG", "TTTV"], ["ATGCGTACGTAGCTAGCTAGNGG", "TTTACGTAGCTTTTV"])
# Returns list of gRNA designs
```

### 7.  Generate Construction File with Cas13 Integration (`create_cas13_construction_file`)

- **Description**: Generates a construction file for RNA targeting using Cas13.
- **Input**: Organism name (string), RNA sequence (string).
- **Output**: String detailing construction steps with Cas13 integration.

**Example**:
```python
create_cas13_construction_file("Drosophila melanogaster", "AUGCGCUAUGCUAGC")
# Returns detailed Cas13 construction file content
```
---

## Error Handling

### gRNA Design
- Raises `ValueError`  if the PAM sequence is not found in the gene sequence.
- Raises `ValueError` if the gene sequence is too short to contain a valid protospacer.

### Toolkit Selection
- Raises `ValueError` if the organism name is not recognized or supported.

### Construction File Creation
- Raises `ValueError` if the organism name is not recognized.
- Raises `ValueError` if the gene sequence is empty or invalid.

### RNA Secondary Structure Analysis
- Raises `ValueError` if the RNA sequence is empty or contains invalid characters.

### CRISPR Efficiency Prediction
- Returns lower efficiency scores if the PAM sequence does not match the target sequence.
- Raises `ValueError` if invalid characters are present in the target sequence.

### Multiplexed gRNA Design
- Raises `ValueError` if the number of PAM sequences does not match the number of gene sequences.
- Raises `ValueError` if any gene sequence is too short or invalid.


### Generate Construction File with Cas13 Integration
- Raises `ValueError` if the RNA sequence is empty.
- Raises `ValueError` if the organism name is not recognized or supported.

---

## Testing

All functions have been tested using **pytest**, covering both standard and edge cases to ensure robustness and accuracy. The tests include:
Valid inputs (e.g., typical PAM sequences, gene/RNA sequences).
Edge cases (e.g., empty inputs, mismatched lists for multiplexed designs).
Invalid inputs (e.g., unsupported organisms, non-nucleotide characters in sequences).

**Test Files**: 
- `tests/test_gRNA_design.py`
- `tests/test_toolkit_selection.py`
- `tests/test_construction_file.py`
- `tests/test_rna_secondary_structure.py`
- `tests/test_crispr_efficiency.py`
- `tests/test_multiplexed_gRNA.py`
- `tests/test_cas13_construction_file.py`

The tests validate:
1. Correct outputs for typical inputs.
2. Proper error handling for invalid inputs.
3. Functionality across diverse scenarios (e.g., long RNA sequences, multiple PAMs)

---

## Usage Instructions

Clone the repository and install the required dependencies listed in `requirements.txt`. The functions can be imported from the `bio_functions.py` module.

**Example**:

```bash
pip install -r requirements.txt
```

Once installed, you can use the functions as follows:

```python
from bio_functions import design_gRNA,
    select_toolkit,
    create_construction_file,
    analyze_rna_secondary_structure,
    predict_crispr_efficiency,
    design_multiplexed_gRNAs,
    create_cas13_construction_file

# Example: Design a gRNA
gRNA = design_gRNA("NGG", "ATGCGTACGTAGCTAGCTAGNGG")
print(gRNA)

# Example: Select a toolkit
toolkit = select_toolkit("Drosophila melanogaster")
print(toolkit)

# Example: Generate a construction file
construction_file = create_construction_file("Drosophila melanogaster", "ATGCGTACGTAGCTAGCTAGNGG")
print(construction_file)
```

---

## Conclusion

This project provides essential tools for CRISPR-mediated genome editing workflows in bioinformatics, supporting precise genetic manipulation across various insect species. The functions are robustly implemented, thoroughly tested, and well-documented to ensure reliability and ease of use in real-world applications.
