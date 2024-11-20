import re
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as scd
import matplotlib.pyplot as plt

def clean_rtf_text(rtf_text):
    """
    Clean RTF formatted text and extract plain text
    """
    # Remove RTF control words
    text = re.sub(r"\\[a-z0-9*\-]+", "", rtf_text)
    # Remove curly braces
    text = re.sub(r"[{}]", "", text)
    # Remove any remaining backslashes
    text = text.replace("\\", "")
    return text.strip()

def parse_fasta_from_rtf(rtf_text):
    """
    Parse FASTA format sequences from RTF text
    """
    # Clean RTF formatting
    clean_text = clean_rtf_text(rtf_text)
    
    sequences = {}
    current_name = None
    current_seq = []
    
    # Process each line
    for line in clean_text.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_name:
                sequences[current_name] = ''.join(current_seq)
            current_name = line[1:].split()[0]  # Get identifier before space
            current_seq = []
        elif line:
            # Remove any remaining non-amino acid characters
            cleaned_line = re.sub(r'[^A-Za-z]', '', line)
            if cleaned_line:
                current_seq.append(cleaned_line)
    
    # Add the last sequence
    if current_name:
        sequences[current_name] = ''.join(current_seq)
    
    return sequences

def lcs(seq1, seq2):
    """
    Compute the Longest Common Subsequence between two sequences
    """
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Fill dp table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    
    # Backtrack to find LCS
    lcs_seq = []
    i, j = m, n
    while i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]:
            lcs_seq.append(seq1[i-1])
            i -= 1
            j -= 1
        elif dp[i-1][j] > dp[i][j-1]:
            i -= 1
        else:
            j -= 1
            
    return dp[m][n], ''.join(reversed(lcs_seq))

def compute_all_pairs_lcs(sequences):
    """
    Compute LCS between all pairs of sequences
    """
    results = {}
    print("\nComputing LCS for all sequence pairs...")
    for name1 in sequences:
        for name2 in sequences:
            if name1 < name2:  # Avoid redundant comparisons
                length, sequence = lcs(sequences[name1], sequences[name2])
                similarity = length / min(len(sequences[name1]), len(sequences[name2]))
                results[(name1, name2)] = {
                    'length': length,
                    'sequence': sequence,
                    'similarity': similarity
                }
                print(f"{name1} vs {name2}: {similarity:.2%} similarity")
    return results

def save_results_to_file(results, sequences, filename="tulp3_analysis_results.txt"):
    """
    Save analysis results to a file
    """
    with open(filename, 'w') as f:
        # Write sequence information
        f.write("TULP3 Sequence Analysis Results\n")
        f.write("=" * 30 + "\n\n")
        
        f.write("Sequence Lengths:\n")
        f.write("-" * 20 + "\n")
        for name, seq in sequences.items():
            f.write(f"{name}: {len(seq)} amino acids\n")
        f.write("\n")
        
        # Write LCS results
        f.write("Pairwise LCS Results:\n")
        f.write("-" * 20 + "\n")
        for (seq1, seq2), result in results.items():
            f.write(f"\n{seq1} vs {seq2}:\n")
            f.write(f"LCS Length: {result['length']}\n")
            f.write(f"Similarity: {result['similarity']:.2%}\n")
            f.write(f"LCS Sequence (first 50 aa): {result['sequence'][:50]}...\n")
            f.write("-" * 50 + "\n")

def create_similarity_matrix(results, sequences):
    """
    Create a similarity matrix from pairwise LCS results
    """
    # Get unique sequence names
    sequence_names = list(sequences.keys())
    n = len(sequence_names)
    
    # Initialize similarity matrix
    similarity_matrix = np.zeros((n, n))
    np.fill_diagonal(similarity_matrix, 1.0)  # Diagonal is always 1
    
    # Fill similarity matrix
    for i in range(n):
        for j in range(i+1, n):
            seq1, seq2 = sequence_names[i], sequence_names[j]
            
            # Find the result in either direction
            if (seq1, seq2) in results:
                similarity = results[(seq1, seq2)]['similarity']
            elif (seq2, seq1) in results:
                similarity = results[(seq2, seq1)]['similarity']
            else:
                similarity = 0
            
            # Symmetrize the matrix
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    
    return similarity_matrix, sequence_names

def cluster_sequences(similarity_matrix, sequence_names, method='ward'):
    """
    Perform hierarchical clustering on sequences
    """
    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix
    
    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(scd.squareform(distance_matrix), method=method)
    
    # Determine optimal number of clusters using the elbow method
    last = linkage_matrix[-10:, 2]
    acceleration = np.diff(last, 2)
    acceleration_rev = acceleration[::-1]
    k = acceleration_rev.argmax() + 2
    
    # Cut the dendrogram to create clusters
    clusters = sch.fcluster(linkage_matrix, t=k, criterion='maxclust')
    
    # Organize sequences into groups
    clustering_groups = {}
    for i, cluster_id in enumerate(clusters):
        if cluster_id not in clustering_groups:
            clustering_groups[cluster_id] = []
        clustering_groups[cluster_id].append(sequence_names[i])
    
    return clustering_groups

def visualize_dendrogram(similarity_matrix, sequence_names):
    """
    Create and save a dendrogram visualization
    """
    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix
    
    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(scd.squareform(distance_matrix), method='ward')
    
    # Create dendrogram
    plt.figure(figsize=(10, 6))
    sch.dendrogram(
        linkage_matrix,
        labels=sequence_names,
        leaf_rotation=90,
        leaf_font_size=8
    )
    plt.title('Hierarchical Clustering of TULP3 Sequences')
    plt.xlabel('Sequence Names')
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.savefig('tulp3_dendrogram.png')
    plt.close()

def sequence_clustering_analysis(sequences, results):
    """
    Perform clustering analysis on sequences
    """
    # Create similarity matrix
    similarity_matrix, sequence_names = create_similarity_matrix(results, sequences)
    
    # Perform clustering
    clustering_groups = cluster_sequences(similarity_matrix, sequence_names)
    
    # Visualize dendrogram
    visualize_dendrogram(similarity_matrix, sequence_names)
    
    # Save clustering results
    with open('tulp3_clustering_results.txt', 'w') as f:
        f.write("TULP3 Sequence Clustering Analysis\n")
        f.write("=" * 40 + "\n\n")
        
        for cluster_id, group_sequences in clustering_groups.items():
            f.write(f"Cluster {cluster_id}:\n")
            for seq in group_sequences:
                f.write(f"  - {seq} (Length: {len(sequences[seq])} aa)\n")
            f.write("\n")
    
    # Print clustering summary
    print("\nClustering Results:")
    for cluster_id, group_sequences in clustering_groups.items():
        print(f"Cluster {cluster_id}: {len(group_sequences)} sequences")
        print("  Sequences:", ", ".join(group_sequences))
        print()

def main():
    # Read RTF file
    try:
        with open(r'C:\Users\src\homework2\tulp3_relatives-2-1.rtf', 'r', encoding='utf-8') as f:
            rtf_content = f.read()
    except UnicodeDecodeError:
        # Try with a different encoding if utf-8 fails
        with open('tulp3_relatives-2-1.rtf', 'r', encoding='latin-1') as f:
            rtf_content = f.read()
    
    # Parse sequences
    sequences = parse_fasta_from_rtf(rtf_content)
    
    # Print initial sequence info
    print("\nFound sequences:")
    for name, seq in sequences.items():
        print(f"{name}: {len(seq)} amino acids")
    
    # Compute LCS for all pairs
    results = compute_all_pairs_lcs(sequences)
    
    # Save initial results to file
    save_results_to_file(results, sequences)
    
    # Perform clustering analysis
    sequence_clustering_analysis(sequences, results)
    
    print("\nResults have been saved to:")
    print("1. tulp3_analysis_results.txt")
    print("2. tulp3_clustering_results.txt")
    print("3. tulp3_dendrogram.png")
    
    # Print summary of highly similar sequences (>70% similarity)
    print("\nHighly similar sequences (>70% similarity):")
    for (seq1, seq2), result in results.items():
        if result['similarity'] > 0.7:
            print(f"{seq1} - {seq2}: {result['similarity']:.2%}")

if __name__ == "__main__":
    main()