import pandas as pd

# Files
iso_res = "/Users/mateuslisboa/Desktop/Data/LRRK2/DESeq2_Results_Isogenic.csv"
non_iso_res = "/Users/mateuslisboa/Desktop/Data/LRRK2/DESeq2_Results_Non_Isogenic.csv"
symbols_map = "/Users/mateuslisboa/Desktop/Data/LRRK2/LRRK2_mRNA_counts_with_symbols.csv"

# Load mapping
print("Loading gene symbols mapping...")
map_df = pd.read_csv(symbols_map, usecols=['gene_id', 'gene_symbol'])
mapping = dict(zip(map_df['gene_id'], map_df['gene_symbol']))

def process_results(path, name):
    print(f"Annotating {name}...")
    df = pd.read_csv(path)
    # The first column is unnamed and contains gene IDs
    df.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
    df['gene_symbol'] = df['gene_id'].map(mapping)
    
    # Reorder
    cols = ['gene_id', 'gene_symbol'] + [c for c in df.columns if c not in ['gene_id', 'gene_symbol']]
    df = df[cols]
    
    # Save annotated version
    out_path = path.replace(".csv", "_Annotated.csv")
    df.to_csv(out_path, index=False)
    
    # Get top 20 significant
    top_sig = df[df['padj'] < 0.05].sort_values('padj').head(20)
    print(f"\nTop 20 Significant Genes for {name}:")
    print(top_sig[['gene_symbol', 'log2FoldChange', 'padj']].to_string(index=False))
    print("-" * 30)

process_results(iso_res, "Isogenic Comparison")
process_results(non_iso_res, "Non-Isogenic Comparison")
