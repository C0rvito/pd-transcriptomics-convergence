import pandas as pd
import re
import subprocess
import os
import glob

def parse_gtf_mapping(gtf_path):
    """Parses GTF to get gene_id to gene_name mapping."""
    mapping = {}
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            if '\tgene\t' in line:
                gene_id_match = re.search(r'gene_id "([^"]+)"', line)
                gene_name_match = re.search(r'gene_name "([^"]+)"', line)
                if gene_id_match and gene_name_match:
                    mapping[gene_id_match.group(1)] = gene_name_match.group(1)
    return mapping

def add_gene_symbols(matrix_path, gtf_path, output_path):
    """Adds gene symbols to a count matrix using a GTF file."""
    mapping = parse_gtf_mapping(gtf_path)
    df = pd.read_csv(matrix_path)
    df['gene_symbol'] = df['gene_id'].map(mapping)
    
    cols = df.columns.tolist()
    # Assume first column is gene_id
    cols = [cols[0], 'gene_symbol'] + [c for c in cols[1:] if c != 'gene_symbol']
    df = df[cols]
    
    df.to_csv(output_path, index=False)
    return df

def run_star_alignment(sample_list, genome_dir, trimmed_dir, out_dir, threads=8, star_path="STAR"):
    """Runs STAR alignment for a list of samples."""
    for sample in sample_list:
        fastq = os.path.join(trimmed_dir, f"{sample}_trimmed.fastq.gz")
        prefix = os.path.join(out_dir, f"{sample}_")
        
        command = [
            star_path,
            "--runThreadN", str(threads),
            "--genomeDir", genome_dir,
            "--readFilesIn", fastq,
            "--readFilesCommand", "gunzip -c",
            "--outFileNamePrefix", prefix,
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--quantMode", "GeneCounts"
        ]
        
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error aligning {sample}: {e}")
            break

def create_count_matrix(input_dir, output_file):
    """Merges STAR ReadsPerGene files into a single matrix."""
    files = glob.glob(os.path.join(input_dir, "*_ReadsPerGene.out.tab"))
    files.sort()
    
    dfs = []
    for f in files:
        sample_name = os.path.basename(f).replace("_ReadsPerGene.out.tab", "")
        df = pd.read_csv(f, sep='\t', header=None, skiprows=4, names=['gene_id', sample_name, 'fwd', 'rev'])
        df = df[['gene_id', sample_name]]
        df.set_index('gene_id', inplace=True)
        dfs.append(df)
    
    matrix = pd.concat(dfs, axis=1)
    matrix.to_csv(output_file)
    return matrix
