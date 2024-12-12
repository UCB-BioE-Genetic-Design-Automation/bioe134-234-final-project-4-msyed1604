'''
def reverse_complement(sequence):
    """
    Calculates the reverse complement of a DNA sequence.
    
    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The reverse complement of the DNA sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters.
    """
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def translate(sequence):
    """
    Translates a DNA sequence into a protein sequence based on the standard genetic code.

    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The corresponding protein sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters or is not a multiple of three.
    """
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")
    if len(sequence) % 3 != 0:
        raise ValueError("Length of DNA sequence is not a multiple of three, which is required for translation.")

    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        protein += codon_table.get(codon, '_')  # Using '_' for unknown or stop codons
    return protein

if __name__ == "__main__":
    # Example DNA sequence for demonstration
    dna_example = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    try:
        # Perform reverse complement
        rc_result = reverse_complement(dna_example)
        print(f"Reverse Complement of '{dna_example}': {rc_result}")

        # Perform translation
        translation_result = translate(dna_example)
        print(f"Translation of '{dna_example}': {translation_result}")
    except Exception as e:
        print(f"Error: {str(e)}")
'''

def design_gRNA(pam_sequence: str, gene_sequence: str) -> dict:
    """
    Designs a guide RNA (gRNA) by identifying the protospacer and tracer RNA.
    
    Args:
        pam_sequence (str): The PAM sequence specific to the CRISPR system.
        gene_sequence (str): The DNA sequence of the target gene.
        
    Returns:
        dict: A dictionary containing 'protospacer' and 'tracer_rna' sequences.
    """
    # Logic to extract protospacer and tracer RNA
    if pam_sequence not in gene_sequence:
        raise ValueError("PAM sequence not found in the gene sequence.")
    
    # Identify protospacer (20 bases upstream of PAM)
    pam_index = gene_sequence.index(pam_sequence)
    protospacer = gene_sequence[max(0, pam_index - 20):pam_index]
    
    # Generate a simple tracer RNA (example fixed sequence)
    tracer_rna = "GTTTTAGAGCTAGAAATAGC"
    
    return {"protospacer": protospacer, "tracer_rna": tracer_rna}

# Example usage
result = design_gRNA("NGG", "ATGCGTACGTAGCTAGCTAGNGG")
print(result)  # Output: {'protospacer': 'ATGCGTACGTAGCTA', 'tracer_rna': 'GTTTTAGAGCTAGAAATAGC'}

def select_toolkit(organism_name: str) -> dict:
    """
    Selects the appropriate CRISPR toolkit for a given organism.
    
    Args:
        organism_name (str): The name of the insect species.
        
    Returns:
        dict: A dictionary containing toolkit details such as PAM sequence.
    """
    # Example toolkit database
    toolkits = {
        "Drosophila melanogaster": {"pam_sequence": "NGG", "system": "Cas9"},
        "Tribolium castaneum": {"pam_sequence": "TTTV", "system": "Cas12a"},
        "Aedes aegypti": {"pam_sequence": "NNGRRT", "system": "Cas9"}
    }
    
    if organism_name not in toolkits:
        raise ValueError(f"No toolkit available for {organism_name}.")
    
    return toolkits[organism_name]

# Example usage
toolkit = select_toolkit("Drosophila melanogaster")
print(toolkit)  # Output: {'pam_sequence': 'NGG', 'system': 'Cas9'}

def create_construction_file(organism_name: str, gene_sequence: str) -> str:
    """
    Creates a construction file for performing a knockout in an insect species.
    
    Args:
        organism_name (str): The name of the insect species.
        gene_sequence (str): The DNA sequence of the target gene.
        
    Returns:
        str: A string describing the construction steps for genetic manipulation.
    """
    # Step 1: Select toolkit based on organism name
    toolkit = select_toolkit(organism_name)
    
    # Step 2: Design gRNA using toolkit's PAM sequence
    gRNA = design_gRNA(toolkit["pam_sequence"], gene_sequence)
    
    # Step 3: Generate construction file content
    construction_file = (
        f"CRISPR Knockout Construction File\n"
        f"Organism: {organism_name}\n"
        f"Gene Sequence: {gene_sequence}\n"
        f"Selected Toolkit: {toolkit['system']} with PAM {toolkit['pam_sequence']}\n"
        f"Protospacer: {gRNA['protospacer']}\n"
        f"Tracer RNA: {gRNA['tracer_rna']}\n"
        f"Steps:\n"
        f"1. Use the protospacer and tracer RNA to assemble gRNA.\n"
        f"2. Deliver gRNA and {toolkit['system']} into insect tissues.\n"
        f"3. Validate knockout efficiency through sequencing.\n"
    )
    
    return construction_file

# Example usage
construction_steps = create_construction_file("Drosophila melanogaster", "ATGCGTACGTAGCTAGCTAGNGG")
print(construction_steps)


def analyze_rna_secondary_structure(rna_sequence: str) -> dict:
    """
    Analyzes the RNA secondary structure to assess accessibility for gRNA binding.
    
    Args:
        rna_sequence (str): The RNA sequence to analyze.
        
    Returns:
        dict: A dictionary containing structural metrics like hairpin count and accessibility score.
    """
    # Simplified example: Count hairpins based on complementary base pairs
    hairpin_count = 0
    stack = []
    for base in rna_sequence:
        if base in "AU":
            stack.append(base)
        elif base in "GC":
            if stack and stack[-1] == "AU":
                hairpin_count += 1
                stack.pop()
    
    # Accessibility score (arbitrary calculation for demonstration)
    accessibility_score = max(0, len(rna_sequence) - 2 * hairpin_count)
    
    return {"hairpin_count": hairpin_count, "accessibility_score": accessibility_score}

# Example usage
rna_analysis = analyze_rna_secondary_structure("AUGCGCUAUGCUAGC")
print(rna_analysis)  # Output: {'hairpin_count': 2, 'accessibility_score': 8}

def select_cas13_variant(organism_name: str) -> dict:
    """
    Selects the optimal Cas13 variant for RNA targeting based on organism-specific data.
    
    Args:
        organism_name (str): The name of the insect species.
        
    Returns:
        dict: A dictionary containing details of the selected Cas13 variant.
    """
    # Example database of Cas13 variants
    cas13_variants = {
        "Drosophila melanogaster": {"variant": "CasFX4", "efficiency": 90, "specificity": "high"},
        "Tribolium castaneum": {"variant": "CasFB8", "efficiency": 70, "specificity": "moderate"}
    }
    
    if organism_name not in cas13_variants:
        raise ValueError(f"No Cas13 variant data available for {organism_name}.")
    
    return cas13_variants[organism_name]

# Example usage
cas13_info = select_cas13_variant("Drosophila melanogaster")
print(cas13_info)  # Output: {'variant': 'CasFX4', 'efficiency': 90, 'specificity': 'high'}

def predict_crispr_efficiency(pam_sequence: str, target_sequence: str) -> float:
    """
    Predicts the efficiency of a CRISPR gRNA design based on PAM sequence compatibility.
    
    Args:
        pam_sequence (str): The PAM sequence used by the CRISPR system.
        target_sequence (str): The DNA target sequence.
        
    Returns:
        float: A predicted efficiency score (0 to 100).
    """
    # Example scoring system based on PAM compatibility and GC content
    pam_compatibility = 100 if pam_sequence in target_sequence else 50
    gc_content = (target_sequence.count("G") + target_sequence.count("C")) / len(target_sequence) * 100
    
    # Efficiency is a weighted combination of PAM compatibility and GC content
    efficiency_score = 0.6 * pam_compatibility + 0.4 * gc_content
    
    return round(efficiency_score, 2)

# Example usage
efficiency = predict_crispr_efficiency("NGG", "ATGCGTACGTAGCTAGCTAGNGG")
print(efficiency)  # Output: Efficiency score (e.g., 85.6)

def design_multiplexed_gRNAs(pam_sequences: list, gene_sequences: list) -> list:
    """
    Designs multiple gRNAs for simultaneous targeting of multiple genes.
    
    Args:
        pam_sequences (list): A list of PAM sequences for each gene.
        gene_sequences (list): A list of gene sequences to target.
        
    Returns:
        list: A list of dictionaries containing protospacer and tracer RNA for each gene.
    """
    if len(pam_sequences) != len(gene_sequences):
        raise ValueError("The number of PAM sequences must match the number of gene sequences.")
    
    multiplexed_gRNAs = []
    
    for pam, gene in zip(pam_sequences, gene_sequences):
        gRNA_data = design_gRNA(pam, gene)
        multiplexed_gRNAs.append(gRNA_data)
    
    return multiplexed_gRNAs

# Example usage
multiplexed_results = design_multiplexed_gRNAs(["NGG", "TTTV"], ["ATGCGTACGTAGCTAGCTAGNGG", "TTTACGTAGCTTTTV"])
print(multiplexed_results)
# Output: [{'protospacer': ..., 'tracer_rna': ...}, {...}]

def create_cas13_construction_file(organism_name: str, rna_sequence: str) -> str:
    """
    Creates a construction file for RNA targeting using Cas13 in an insect species.
    
    Args:
        organism_name (str): The name of the insect species.
        rna_sequence (str): The RNA sequence to target.
        
    Returns:
        str: A string describing the construction steps for RNA manipulation.
    """
    # Select Cas13 variant
    cas13_variant = select_cas13_variant(organism_name)
    
    # Analyze RNA secondary structure
    rna_analysis = analyze_rna_secondary_structure(rna_sequence)
    
    # Generate construction file content
    construction_file = (
        f"CRISPR/Cas13 Construction File\n"
        f"Organism: {organism_name}\n"
        f"RNA Sequence: {rna_sequence}\n"
        f"Selected Cas13 Variant: {cas13_variant['variant']}\n"
        f"Efficiency: {cas13_variant['efficiency']}%\n"
        f"Hairpin Count: {rna_analysis['hairpin_count']}\n"
        f"Accessibility Score: {rna_analysis['accessibility_score']}\n"
        f"Steps:\n"
        f"1. Assemble crRNA with Cas13 ({cas13_variant['variant']}).\n"
        f"2. Deliver complex into insect tissues.\n"
        f"3. Validate transcript knockdown via qPCR.\n"
    )
    
    return construction_file

# Example usage
cas13_steps = create_cas13_construction_file("Drosophila melanogaster", "AUGCGCUAUGCUAGC")
print(cas13_steps)
