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

