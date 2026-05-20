import pandas as pd
import re

gtf_path = "/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/reference/gencode.v45.primary_assembly.annotation.gtf"
matrix_path = "/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/de/2_DESeq2_SNCA_Triplication_vs_Control_Basal.csv"
output_path = "/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/de/SNCA_Triplication_counts_with_symbols.csv"

# 1. Parse GTF to get gene_id to gene_name mapping
print("Parsing GTF for gene symbols...")
mapping = {}
with open(gtf_path, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        if '\tgene\t' in line:
            gene_id = re.search(r'gene_id "([^"]+)"', line).group(1)
            gene_name = re.search(r'gene_name "([^"]+)"', line).group(1)
            mapping[gene_id] = gene_name

# 2. Load matrix
print("Loading count matrix...")
df = pd.read_csv(matrix_path)

# 3. Add symbols
print("Adding gene symbols...")
df['gene_symbol'] = df['gene_id'].map(mapping)

# 4. Reorder columns to put symbol near ID
cols = df.columns.tolist()
cols = [cols[0], 'gene_symbol'] + cols[1:-1]
df = df[cols]

# 5. Save
df.to_csv(output_path, index=False)
print(f"Final matrix saved to: {output_path}")
print(f"Total genes mapped: {df['gene_symbol'].notna().sum()} / {len(df)}")
