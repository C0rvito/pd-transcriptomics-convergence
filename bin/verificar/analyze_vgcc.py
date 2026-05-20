import pandas as pd
import os

# Paths
base_dir = "/Users/mateuslisboa/Desktop/Data/LRRK2/outputs/differential_expression"
output_file = "/Users/mateuslisboa/Desktop/Data/LRRK2/outputs/VGCC_Analysis.csv"

# VGCC Gene families
# CACNA1* (Alpha subunits)
# CACNB* (Beta subunits)
# CACNA2D* (Alpha2-Delta subunits)
# CACNG* (Gamma subunits)

files = {
    'Isogenic': 'DESeq2_Results_Isogenic_Annotated.csv',
    'Combined': 'DESeq2_Results_Combined_Annotated.csv'
}

all_vgcc_results = []

for label, filename in files.items():
    path = os.path.join(base_dir, filename)
    if not os.path.exists(path): continue
    
    df = pd.read_csv(path)
    
    # Filter for VGCC genes
    # We want exact matches for the main families
    vgcc_mask = df['gene_symbol'].str.match(r'^(CACNA1[A-Z]|CACNA1[0-9]|CACNB[1-4]|CACNA2D[1-4]|CACNG[1-8])$', na=False)
    vgcc_df = df[vgcc_mask].copy()
    vgcc_df['Comparison'] = label
    all_vgcc_results.append(vgcc_df)

final_df = pd.concat(all_vgcc_results)
final_df = final_df.sort_values(['gene_symbol', 'Comparison'])

# Select relevant columns
cols = ['gene_symbol', 'Comparison', 'log2FoldChange', 'padj', 'baseMean']
final_df[cols].to_csv(output_file, index=False)

# Identify significant ones
sig_vgcc = final_df[final_df['padj'] < 0.05]
print("Significant VGCC Genes (padj < 0.05):")
if not sig_vgcc.empty:
    print(sig_vgcc[cols].to_string(index=False))
else:
    print("No significant VGCC genes found.")

# Summary of trends
print("\nGeneral Trends for VGCCs (Log2FC):")
summary = final_df.groupby('gene_symbol')['log2FoldChange'].mean()
print(summary)
