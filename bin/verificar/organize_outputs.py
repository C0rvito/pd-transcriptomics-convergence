import pandas as pd
import os
import shutil

# Paths
base_dir = "/Users/mateuslisboa/Desktop/Data/LRRK2"
output_dir = os.path.join(base_dir, "outputs")

# Files to move and convert
isogenic_csv = os.path.join(base_dir, "DESeq2_Results_Isogenic_Annotated.csv")
non_isogenic_csv = os.path.join(base_dir, "DESeq2_Results_Non_Isogenic_Annotated.csv")
counts_csv = os.path.join(base_dir, "LRRK2_mRNA_counts_with_symbols.csv")
metadata_csv = os.path.join(base_dir, "metadata.csv")

# 1. Create Excel Workbook for DE results
excel_out = os.path.join(output_dir, "differential_expression", "LRRK2_DEA_Results_Summary.xlsx")
print(f"Creating Excel summary at {excel_out}...")

with pd.ExcelWriter(excel_out) as writer:
    if os.path.exists(isogenic_csv):
        pd.read_csv(isogenic_csv).to_excel(writer, sheet_name='Isogenic_Comparison', index=False)
    if os.path.exists(non_isogenic_csv):
        pd.read_csv(non_isogenic_csv).to_excel(writer, sheet_name='Non_Isogenic_Comparison', index=False)
    if os.path.exists(metadata_csv):
        pd.read_csv(metadata_csv).to_excel(writer, sheet_name='Metadata', index=False)

# 2. Convert counts to Excel
counts_excel = os.path.join(output_dir, "counts", "LRRK2_mRNA_Counts_Annotated.xlsx")
if os.path.exists(counts_csv):
    print(f"Converting counts to {counts_excel}...")
    pd.read_csv(counts_csv).to_excel(counts_excel, index=False)

# 3. Organize other files
print("Moving files to subfolders...")

# Move Plots
for f in os.listdir(base_dir):
    if f.endswith(".png"):
        shutil.move(os.path.join(base_dir, f), os.path.join(output_dir, "plots", f))

# Move Scripts
scripts = ["align_sequential.py", "create_matrix.py", "add_symbols.py", "run_deseq2.R", "annotate_results.py"]
for s in scripts:
    path = os.path.join(base_dir, s)
    if os.path.exists(path):
        shutil.move(path, os.path.join(output_dir, "scripts", s))

# Move original CSVs to a 'raw' or keep them in DE folder
csv_files = ["DESeq2_Results_Isogenic.csv", "DESeq2_Results_Non_Isogenic.csv", 
             "DESeq2_Results_Isogenic_Annotated.csv", "DESeq2_Results_Non_Isogenic_Annotated.csv"]
for f in csv_files:
    path = os.path.join(base_dir, f)
    if os.path.exists(path):
        shutil.move(path, os.path.join(output_dir, "differential_expression", f))

# Move count matrix and metadata
if os.path.exists(os.path.join(base_dir, "LRRK2_mRNA_counts_matrix.csv")):
    shutil.move(os.path.join(base_dir, "LRRK2_mRNA_counts_matrix.csv"), os.path.join(output_dir, "counts", "LRRK2_mRNA_counts_matrix.csv"))
if os.path.exists(counts_csv):
    shutil.move(counts_csv, os.path.join(output_dir, "counts", "LRRK2_mRNA_counts_with_symbols.csv"))
if os.path.exists(metadata_csv):
    shutil.move(metadata_csv, os.path.join(output_dir, "counts", "metadata.csv"))
if os.path.exists(os.path.join(base_dir, "sequential_alignment.log")):
    shutil.move(os.path.join(base_dir, "sequential_alignment.log"), os.path.join(output_dir, "scripts", "sequential_alignment.log"))

print("Organization complete.")
