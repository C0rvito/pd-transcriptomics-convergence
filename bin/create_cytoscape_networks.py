#!/usr/bin/env python3
"""
create_cytoscape_networks.py
---------------------------
Orchestrates Cytoscape via py4cytoscape to create STRING networks 
and export PNG images for all DE comparisons.
"""

import os
import sys
import glob
import pandas as pd
import py4cytoscape as p4c
import time

def create_and_export_network(de_file, out_root, gene_name, comparison_name):
    print(f"\n>>> Processing {gene_name} | Comparison: {comparison_name}")
    
    try:
        df = pd.read_csv(de_file)
    except Exception as e:
        print(f"  [ERROR] Could not read {de_file}: {e}")
        return

    # Filter significant genes (padj < 0.05)
    sig_df = df[df['padj'] < 0.05].copy()
    if sig_df.empty:
        print(f"  [SKIP] No significant genes found.")
        return

    # Prepare sets
    tasks = [
        ("Up", sig_df[sig_df['log2FoldChange'] > 0].sort_values('padj').head(100)),
        ("Down", sig_df[sig_df['log2FoldChange'] < 0].sort_values('padj').head(100)),
        ("Combined", pd.concat([
            sig_df[sig_df['log2FoldChange'] > 0].sort_values('padj').head(100),
            sig_df[sig_df['log2FoldChange'] < 0].sort_values('padj').head(100)
        ]))
    ]

    for direction, selected_genes in tasks:
        if selected_genes.empty:
            print(f"  [SKIP] No {direction} genes found.")
            continue

        gene_list = ",".join(selected_genes['symbol'].dropna().astype(str).tolist())
        network_title = f"{gene_name}_{comparison_name}_{direction}"

        try:
            print(f"  Querying STRING for {len(selected_genes)} {direction} genes...")
            
            # ── Execute STRING Query ──
            string_cmd = (
                f'string protein query query="{gene_list}" '
                f'species="Homo sapiens" '
                f'limit=20 '
                f'cutoff=0.5 '
                f'networkType="full STRING network"'
            )
            p4c.commands.commands_run(string_cmd)
            time.sleep(3) # Give Cytoscape time to build
            
            p4c.networks.rename_network(network_title)

            # ── Apply Visual Styling ──
            # 1. Prepare data for mapping (Absolute Log2FC and Color)
            selected_genes['abs_log2FC'] = selected_genes['log2FoldChange'].abs()
            selected_genes['node_color'] = selected_genes['log2FoldChange'].apply(lambda x: '#FF0000' if x > 0 else '#0000FF')
            
            # 2. Upload metadata to Cytoscape node table
            # We use 'display name' or 'query term' usually created by STRING app as key
            p4c.tables.load_table_data(selected_genes[['symbol', 'abs_log2FC', 'node_color']], 
                                       data_key_column='symbol', 
                                       table_key_column='display name')

            # 3. Define the style
            style_name = f"Style_{direction}"
            if style_name not in p4c.styles.get_visual_style_names():
                p4c.styles.create_visual_style(style_name)
            
            # Node Color (Passthrough mapping)
            p4c.styles.set_node_color_mapping('node_color', mapping_type='p', style_name=style_name)
            
            # Node Size (Continuous mapping based on abs_log2FC)
            # Map abs_log2FC (typical range 1-10) to Size (typical range 30-100)
            p4c.styles.set_node_size_mapping('abs_log2FC', [1.0, 5.0], [30, 100], mapping_type='c', style_name=style_name)
            
            # Labels and Font
            p4c.styles.set_node_label_mapping('display name', style_name=style_name)
            p4c.styles.set_visual_property_default('NODE_LABEL_FONT_SIZE', 12, style_name=style_name)
            p4c.styles.set_visual_property_default('NODE_LABEL_COLOR', '#333333', style_name=style_name)
            
            # Edge Styling (Subtle grey lines)
            p4c.styles.set_visual_property_default('EDGE_WIDTH', 1.5, style_name=style_name)
            p4c.styles.set_visual_property_default('EDGE_STROKE_UNSELECTED_PAINT', '#CCCCCC', style_name=style_name)

            p4c.styles.set_visual_style(style_name)

            # ── Apply Layout ──
            try:
                p4c.layouts.layout_network(layout_name='yfiles-circular')
            except:
                p4c.layouts.layout_network(layout_name='circular')
            
            # ── Export Image ──
            target_dir = os.path.join(out_root, gene_name, comparison_name, direction)
            os.makedirs(target_dir, exist_ok=True)
            img_path = os.path.join(target_dir, f"network_{direction}_styled.png")
            
            p4c.export_image(img_path, type='PNG', overwrite_file=True)
            print(f"  [SUCCESS] {direction} PNG exported to: {img_path}")

        except Exception as e:
            print(f"  [ERROR] Failed to create/style {direction} network: {e}")

def main():
    root_dir = "/home/jess/Documents/estudos/pd-transcriptomics-convergence"
    out_root = os.path.join(root_dir, "data/cytoscape")
    
    # Check if Cytoscape is running
    try:
        p4c.cytoscape_ping()
        print("Connected to Cytoscape.")
    except Exception:
        print("CRITICAL ERROR: Cytoscape is not running or CyREST is not available.")
        print("Please open Cytoscape before running this script.")
        sys.exit(1)

    # Find all DE result CSVs
    de_files = glob.glob(os.path.join(root_dir, "data/de/**/results/DE_results_*.csv"), recursive=True)
    
    for f in sorted(de_files):
        # Infer gene and comparison from path
        parts = f.split(os.sep)
        gene = parts[-3]
        comp = parts[-1].replace("DE_results_", "").replace(".csv", "")
        
        create_and_export_network(f, out_root, gene, comp)

    print("\nAll tasks finished.")

if __name__ == "__main__":
    main()
