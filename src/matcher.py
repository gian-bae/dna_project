# ============================================================================
# [Day1:Genomic(DNA/RNA) Sequence Matching and String Algorithm Visualization]
# ============================================================================

#src/matcher.py

def load_fasta(file_path):
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

