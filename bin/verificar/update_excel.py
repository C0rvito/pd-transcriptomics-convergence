import pandas as pd
import os

# Paths
base_dir = "/Users/mateuslisboa/Desktop/Data/LRRK2/outputs"
combined_res = os.path.join(base_dir, "differential_expression", "DESeq2_Results_Combined.csv")
symbols_map = os.path.join(base_dir, "counts", "LRRK2_mRNA_counts_with_symbols.csv")
excel_out = os.path.join(base_dir, "differential_expression", "LRRK2_DEA_Results_Summary.xlsx")

# Load mapping
print("Loading gene symbols mapping...")
map_df = pd.read_csv(symbols_map, usecols=['gene_id', 'gene_symbol'])
mapping = dict(zip(map_df['gene_id'], map_df['gene_symbol']))

# 1. Annotate Combined Results
print("Annotating Combined Results...")
df = pd.read_csv(combined_res)
df.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
df['gene_symbol'] = df['gene_id'].map(mapping)
cols = ['gene_id', 'gene_symbol'] + [c for c in df.columns if c not in ['gene_id', 'gene_symbol']]
df = df[cols]
combined_annotated = combined_res.replace(".csv", "_Annotated.csv")
df.to_csv(combined_annotated, index=False)

# 2. Update Excel with all three sheets
print("Updating Excel summary...")
isogenic_annotated = os.path.join(base_dir, "differential_expression", "DESeq2_Results_Isogenic_Annotated.csv")
non_isogenic_annotated = os.path.join(base_dir, "differential_expression", "DESeq2_Results_Non_Isogenic_Annotated.csv")
metadata_csv = os.path.join(base_dir, "counts", "metadata.csv")

with pd.ExcelWriter(excel_out) as writer:
    pd.read_csv(combined_annotated).to_excel(writer, sheet_name='Combined_Comparison', index=False)
    pd.read_csv(isogenic_annotated).to_excel(writer, sheet_name='Isogenic_Comparison', index=False)
    pd.read_csv(non_isogenic_annotated).to_excel(writer, sheet_name='Non_Isogenic_Comparison', index=False)
    pd.read_csv(metadata_csv).to_excel(writer, sheet_name='Metadata', index=False)

print("Final Excel version created with 4 sheets.")
