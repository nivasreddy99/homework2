def read_rtf_file(filename):
    sequence = ""
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            content = file.read()
            # Extract only capital letters (DNA/Protein sequence)
            sequence = ''.join(char for char in content if char.isupper())
        print(f"Read {len(sequence)} characters from {filename}")
        if len(sequence) == 0:
            print(f"Warning: No uppercase letters found in {filename}")
    except FileNotFoundError:
        print(f"Error: File {filename} not found")
    except Exception as e:
        print(f"Error reading file: {str(e)}")
    return sequence

class Node:
    def __init__(self, char, freq):
        self.char = char
        self.freq = freq
        self.left = None
        self.right = None

def calculate_frequency(text):
    """Calculate frequency of each character in the text."""
    if not text:
        return {}
    freq_dict = {}
    for char in text:
        freq_dict[char] = freq_dict.get(char, 0) + 1
    return freq_dict

def build_huffman_tree(freq_dict):
    """Build Huffman tree from frequency dictionary."""
    if not freq_dict:
        return None
        
    nodes = [Node(char, freq) for char, freq in freq_dict.items()]
    
    while len(nodes) > 1:
        nodes = sorted(nodes, key=lambda x: x.freq)
        left = nodes.pop(0)
        right = nodes.pop(0)
        
        internal = Node(None, left.freq + right.freq)
        internal.left = left
        internal.right = right
        
        nodes.append(internal)
    
    return nodes[0] if nodes else None

def generate_huffman_codes(root, code="", codes=None):
    """Generate Huffman codes for each character."""
    if codes is None:
        codes = {}
    
    if root:
        if root.char:
            codes[root.char] = code if code else "0"  # Handle single character case
        generate_huffman_codes(root.left, code + "0", codes)
        generate_huffman_codes(root.right, code + "1", codes)
    
    return codes

def binary_coding(text):
    """Implement binary coding for the sequence."""
    if not text:
        return {}
        
    unique_chars = sorted(set(text))
    n_chars = len(unique_chars)
    bits_needed = max(1, len(bin(n_chars - 1)[2:]))  # Minimum 1 bit
    
    binary_codes = {}
    for i, char in enumerate(unique_chars):
        binary_code = bin(i)[2:].zfill(bits_needed)
        binary_codes[char] = binary_code
    
    return binary_codes

def compress_sequence(text, coding_dict):
    """Compress sequence using provided coding dictionary."""
    if not text or not coding_dict:
        return ""
    return ''.join(coding_dict[char] for char in text)

def analyze_compression(text, name="Sequence"):
    """Analyze compression ratios for both binary and Huffman coding."""
    if not text:
        print(f"\nError: Empty {name} sequence")
        return 0, 0

    print(f"\nAnalysis for {name}:")
    print(f"Original sequence length: {len(text)} characters")
    print(f"Original size (ASCII): {len(text) * 8} bits")
    print(f"Unique characters: {len(set(text))}")
    
    # Binary coding analysis
    binary_codes = binary_coding(text)
    if not binary_codes:
        print("Error: Binary coding failed")
        return 0, 0
        
    binary_compressed = compress_sequence(text, binary_codes)
    binary_size = len(binary_compressed)
    print(f"\nBinary Coding:")
    print(f"Bits per character: {len(next(iter(binary_codes.values())))}")
    print(f"Compressed size: {binary_size} bits")
    if binary_size > 0:
        print(f"Compression ratio: {(len(text) * 8) / binary_size:.2f}x")
    print("Binary codes:", binary_codes)
    
    # Huffman coding analysis
    freq_dict = calculate_frequency(text)
    if not freq_dict:
        print("Error: Frequency calculation failed")
        return binary_size, 0
        
    huffman_tree = build_huffman_tree(freq_dict)
    if not huffman_tree:
        print("Error: Huffman tree construction failed")
        return binary_size, 0
        
    huffman_codes = generate_huffman_codes(huffman_tree)
    if not huffman_codes:
        print("Error: Huffman code generation failed")
        return binary_size, 0
        
    huffman_compressed = compress_sequence(text, huffman_codes)
    huffman_size = len(huffman_compressed)
    
    print(f"\nHuffman Coding:")
    print(f"Compressed size: {huffman_size} bits")
    if huffman_size > 0:
        print(f"Compression ratio: {(len(text) * 8) / huffman_size:.2f}x")
    print("Character frequencies:", freq_dict)
    print("Huffman codes:", huffman_codes)
    
    return binary_size, huffman_size

def main():
    # Read sequences from files
    print("Reading DNA sequence file...")
    dna_seq = read_rtf_file(r"C:\Users\src\homework2\DNA_SEQ-1.rtf")
    print("\nReading protein sequence file...")
    prot_seq = read_rtf_file(r"C:\Users\src\homework2\PROT_SEQ-1.rtf")
    
    if not dna_seq and not prot_seq:
        print("Error: Both sequence files are empty or couldn't be read")
        return
    
    if dna_seq:
        print("\n=== DNA Sequence Analysis ===")
        dna_binary_size, dna_huffman_size = analyze_compression(dna_seq, "DNA")
    else:
        dna_binary_size, dna_huffman_size = 0, 0
    
    if prot_seq:
        print("\n=== Protein Sequence Analysis ===")
        prot_binary_size, prot_huffman_size = analyze_compression(prot_seq, "Protein")
    else:
        prot_binary_size, prot_huffman_size = 0, 0
    
    # Print comparative analysis only if we have valid data
    if dna_seq or prot_seq:
        print("\n=== Comparative Analysis ===")
        
        if dna_seq and dna_binary_size > 0 and dna_huffman_size > 0:
            print("\nDNA Sequence:")
            print(f"Binary vs ASCII savings: {((len(dna_seq) * 8) - dna_binary_size) / (len(dna_seq) * 8) * 100:.2f}%")
            print(f"Huffman vs ASCII savings: {((len(dna_seq) * 8) - dna_huffman_size) / (len(dna_seq) * 8) * 100:.2f}%")
            print(f"Huffman vs Binary savings: {(dna_binary_size - dna_huffman_size) / dna_binary_size * 100:.2f}%")
        
        if prot_seq and prot_binary_size > 0 and prot_huffman_size > 0:
            print("\nProtein Sequence:")
            print(f"Binary vs ASCII savings: {((len(prot_seq) * 8) - prot_binary_size) / (len(prot_seq) * 8) * 100:.2f}%")
            print(f"Huffman vs ASCII savings: {((len(prot_seq) * 8) - prot_huffman_size) / (len(prot_seq) * 8) * 100:.2f}%")
            print(f"Huffman vs Binary savings: {(prot_binary_size - prot_huffman_size) / prot_binary_size * 100:.2f}%")

if __name__ == "__main__":
    main()