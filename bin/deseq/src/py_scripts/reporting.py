import pandas as pd
import os

def create_excel_summary(files_dict, output_xlsx, metadata_path=None):
    """Creates a multi-sheet Excel summary from a dictionary of CSV/TSV files."""
    with pd.ExcelWriter(output_xlsx) as writer:
        if metadata_path and os.path.exists(metadata_path):
            m_sep = '\t' if metadata_path.endswith('.tsv') else ','
            pd.read_csv(metadata_path, sep=m_sep).to_excel(writer, sheet_name='Metadata', index=False)
            
        for sheet_name, file_path in files_dict.items():
            if os.path.exists(file_path):
                f_sep = '\t' if file_path.endswith('.tsv') else ','
                pd.read_csv(file_path, sep=f_sep).to_excel(writer, sheet_name=sheet_name, index=False)
            else:
                print(f"Warning: {file_path} not found.")
    return output_xlsx
