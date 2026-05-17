import pandas as pd
import os
import re

def annotate_de_results(res_path, symbols_map_path, output_path=None):
    """Annotates DESeq2 results with gene symbols."""
    map_df = pd.read_csv(symbols_map_path)
    mapping = dict(zip(map_df['gene_id'], map_df['gene_symbol']))
    
    # Check if it's tsv or csv
    sep = '\t' if res_path.endswith('.tsv') else ','
    df = pd.read_csv(res_path, sep=sep)
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
    
    df['gene_symbol'] = df['gene_id'].map(mapping)
    
    cols = ['gene_id', 'gene_symbol'] + [c for c in df.columns if c not in ['gene_id', 'gene_symbol']]
    df = df[cols]
    
    if output_path:
        out_sep = '\t' if output_path.endswith('.tsv') else ','
        df.to_csv(output_path, index=False, sep=out_sep)
    return df

def rank_and_segment_results(res_path, output_xlsx, metadata_path=None, padj_threshold=0.05):
    """Ranks by p-value and segments into Up/Down regulated sheets. Appends if file exists."""
    sep = '\t' if res_path.endswith('.tsv') else ','
    df = pd.read_csv(res_path, sep=sep)
    df_ranked = df.sort_values(by='pvalue', ascending=True)
    
    label = os.path.basename(res_path).replace(".csv", "").replace(".tsv", "")
    # Clean label to be used as sheet prefix
    prefix = label.replace("_DE_Mutant_vs_Control", "").replace("_Combined", "")

    sig_mask = (df_ranked['padj'] < padj_threshold)
    up_reg = df_ranked[sig_mask & (df_ranked['log2FoldChange'] > 0)]
    down_reg = df_ranked[sig_mask & (df_ranked['log2FoldChange'] < 0)]

    mode = 'a' if os.path.exists(output_xlsx) else 'w'
    if mode == 'a':
        # Use openpyxl engine for appending
        with pd.ExcelWriter(output_xlsx, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            df_ranked.to_excel(writer, sheet_name=f'{prefix}_Ranked', index=False)
            up_reg.to_excel(writer, sheet_name=f'{prefix}_Up', index=False)
            down_reg.to_excel(writer, sheet_name=f'{prefix}_Down', index=False)
    else:
        with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
            if metadata_path and os.path.exists(metadata_path):
                m_sep = '\t' if metadata_path.endswith('.tsv') else ','
                pd.read_csv(metadata_path, sep=m_sep).to_excel(writer, sheet_name='Metadata', index=False)
            
            df_ranked.to_excel(writer, sheet_name=f'{prefix}_Ranked', index=False)
            up_reg.to_excel(writer, sheet_name=f'{prefix}_Up', index=False)
            down_reg.to_excel(writer, sheet_name=f'{prefix}_Down', index=False)
    
    return up_reg, down_reg

def analyze_vgcc_family(res_path, output_csv):
    """Filters for Voltage-Gated Calcium Channel genes."""
    df = pd.read_csv(res_path)
    
    # CACNA1*, CACNB*, CACNA2D*, CACNG*
    vgcc_pattern = r'^(CACNA1[A-Z]|CACNA1[0-9]|CACNB[1-4]|CACNA2D[1-4]|CACNG[1-8])$'
    vgcc_mask = df['gene_symbol'].str.match(vgcc_pattern, na=False)
    vgcc_df = df[vgcc_mask].copy()
    
    cols = ['gene_symbol', 'log2FoldChange', 'padj', 'baseMean']
    vgcc_df[cols].to_csv(output_csv, index=False)
    return vgcc_df
