import pandas as pd
import os
import glob

# Paths
input_dir = "/Users/mateuslisboa/Desktop/Data/LRRK2/alignment_STAR"
output_file = "/Users/mateuslisboa/Desktop/Data/LRRK2/LRRK2_mRNA_counts_matrix.csv"

# Find all ReadsPerGene files
files = glob.glob(os.path.join(input_dir, "*_ReadsPerGene.out.tab"))
files.sort()

# Dictionary to hold dataframes
dfs = []

for f in files:
    # Extract sample name from filename (e.g., SC1014_mRNA_ReadsPerGene.out.tab -> SC1014_mRNA)
    sample_name = os.path.basename(f).replace("_ReadsPerGene.out.tab", "")
    
    # Read the file
    # STAR output: 1:geneID, 2:unstranded, 3:stranded_forward, 4:stranded_reverse
    # We use column 0 (geneID) and column 1 (unstranded counts)
    df = pd.read_csv(f, sep='\t', header=None, skiprows=4, names=['gene_id', sample_name, 'fwd', 'rev'])
    
    # Keep only gene_id and the unstranded count
    df = df[['gene_id', sample_name]]
    df.set_index('gene_id', inplace=True)
    dfs.append(df)

# Merge all dataframes on gene_id
matrix = pd.concat(dfs, axis=1)

# Save to CSV
matrix.to_csv(output_file)
print(f"Count matrix created with {matrix.shape[0]} genes and {matrix.shape[1]} samples.")
print(f"Saved to: {output_file}")
