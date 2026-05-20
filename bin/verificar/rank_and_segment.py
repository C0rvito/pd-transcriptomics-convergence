import pandas as pd
import os

# Paths
base_dir = "/Users/mateuslisboa/Desktop/Data/LRRK2/outputs/differential_expression"
excel_out = os.path.join(base_dir, "LRRK2_DEA_Results_Summary_Ranked.xlsx")
metadata_path = "/Users/mateuslisboa/Desktop/Data/LRRK2/outputs/counts/metadata.csv"

# Files to process
files = {
    'Combined': 'DESeq2_Results_Combined_Annotated.csv',
    'Isogenic': 'DESeq2_Results_Isogenic_Annotated.csv',
    'Non_Isogenic': 'DESeq2_Results_Non_Isogenic_Annotated.csv'
}

print("Ranking and splitting results by regulation direction...")

with pd.ExcelWriter(excel_out) as writer:
    # 1. Metadata sheet
    pd.read_csv(metadata_path).to_excel(writer, sheet_name='Metadata', index=False)
    
    for label, filename in files.items():
        path = os.path.join(base_dir, filename)
        if not os.path.exists(path):
            print(f"Warning: {path} not found.")
            continue
            
        # Load and rank by p-value (lowest first)
        df = pd.read_csv(path)
        df_ranked = df.sort_values(by='pvalue', ascending=True)
        
        # Save the full ranked results
        df_ranked.to_excel(writer, sheet_name=f'{label}_Full_Ranked', index=False)
        
        # Define significance and regulation
        # padj < 0.05 is the standard threshold
        sig_mask = (df_ranked['padj'] < 0.05)
        
        # Up-regulated (log2FoldChange > 0)
        up_reg = df_ranked[sig_mask & (df_ranked['log2FoldChange'] > 0)]
        # Down-regulated (log2FoldChange < 0)
        down_reg = df_ranked[sig_mask & (df_ranked['log2FoldChange'] < 0)]
        
        # Save segmented sheets
        up_reg.to_excel(writer, sheet_name=f'{label}_UpRegulated', index=False)
        down_reg.to_excel(writer, sheet_name=f'{label}_DownRegulated', index=False)
        
        # Save individual CSVs for convenience
        up_reg.to_csv(os.path.join(base_dir, f"{label}_UpRegulated_padj05.csv"), index=False)
        down_reg.to_csv(os.path.join(base_dir, f"{label}_DownRegulated_padj05.csv"), index=False)
        
        print(f"{label}: {len(up_reg)} Up, {len(down_reg)} Down (padj < 0.05)")

print(f"\nFinal ranked and segmented Excel created: {excel_out}")
