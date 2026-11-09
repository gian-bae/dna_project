# ============================================================================
# [Day1:Genomic(DNA/RNA) Sequence Matching and String Algorithm Visualization]
# ============================================================================

#src/matcher.py

from ast import pattern
from email.mime import text


def load_fasta(file_path):
    """
    Reads a FASTA file and parses headers and DNA sequences. Merges multi-line sequences into a single contiguouse string,
    """
    sequences = {}
    current_header = None

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue # Ignore empty lines

                if line.startswith('>'):
                    current_header = line[1:] # Remove the '>' symbol and set as key
                    sequences[current_header] = ""
                else:
                    if current_header:
                        sequences[current_header] += line.upper() # Convert to uppercase
                        
        return sequences
                    
    except FileNotFoundError:
        print(f"[Error] File not found: {file_path}")
        return {}

def compute_lps(pattern):
    """
    Computes the LPS (Longest Prefix which is also Suffix)table
    for the KMP string matching algorithm.
    """
    lps = [0] * len(pattern)
    length = 0
    i = 1

    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            lps [i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def kmp_search(text, pattern):
"""
Performs pattern matching using the Knuth-Morris-Pratt (KMP) algorithm.
Time Complexity: O(N + M)
"""
n = len(text)
m = len(pattern)
lps = compute_lps(pattern)

matches = []
i = 0 # index for text
j = 0 # index for pattern

while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == m:
            matches.append(i - j)
            j = lps[j - 1]
        elif i < n and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return matches

# Test execution code (runs only when this file is executed directly)
if __name__ == "__main__":
    # Specify the path to the data directory in the parent folder
    fasta_data = load_fasta ("data/sample.fasta")

    print("== FASTA Parsing Result ===")
    for header, seq in fasta_data.items():
        print(f"Header: {header}")
        print(f"Length: {len(seq)} bp")
        print(f"Sequence: {seq}\n")

# ============================================================================
# [Day2: Basic Sequence Statistics and GC Content Analysis]
# ============================================================================

def calculate_gc_content(sequence):
    """
    Calculates the GC content (percentage of G and C nucleotides)
    of a given DNA sequence.
    """
    seq_upper = sequence.upper()
    g_count = seq_upper.count('G')
    c_count = seq_upper.count('C')
    total_length = len(seq_upper)

    if total_length == 0:
        return 0.0

    gc_percentage = ((g_count + c_count) / total_length) * 100
    return round(gc_percentage, 2)

def get_base_counts(sequence):
    """
    Counts the occurrence of each nucleotide (A, T, G, C)
    in a given DNA sequence.
    """
    seq_upper = sequence.upper()
    return {
        'A': seq_upper.count('A'),
        'T': seq_upper.count('T'),
        'G': seq_upper.count('G'),
        'C': seq_upper.count('C'),
        'Other': len(seq_upper) - (seq_upper.count('A') + seq_upper.count('T') + seq_upper.count('G') + seq_upper.count('C'))
    }

# ==========================================
# 🧬 [Day 3 Features: Motif & Restriction Site Search]
# ==========================================

#Dictionary of common restriction enzymes and thier recognition motifs
RESTRICTION_ENZYMES = {
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC',
    'HindIII': 'AAGCTT',
    'NotI': 'GCGGCCGC',
}

def search_motif(sequence, motif):
    """
    Searches for a specific motif (short DNA pattern) in the sequence
    using the KMP search algorithm and returns all starting indices.
    """
    return kmp_search(sequence, motif)

def find_restriction_sites(sequence, enzyme_name):
    """
    Finds cleavage sites for a specific restriction enzyme in the DNA sequence.
    Returns a dictionary containing the enzyme name, motif, found poditions,
    and total count.
    """
    if enzyme_name not in RESTRICTION_ENZYMES:
        raise ValueError(f"Enzyme '{enzyme_name}' is not supported.")

    motif = RESTRICTION_ENZYMES[enzyme_name]
    sites = kmp_search(sequence, motif)

    return {
        'enzyme': enzyme_name,
        'motif': motif,
        'sites': sites,
        'count': len(sites)
    }

if __name__ == "__main__":
    # Test code for Day 3
    test_seq = "AGCTGAATTCCGATGAATTCGCTAAGCTT"
    print(f"Test Sequence: {test_seq}")

    # Test EcoRI motif ssearch
    ecori_result = find_restriction_sites(test_seq, 'EcoRI')
    print(f"EcoRI Analysis Result: {ecori_result}")