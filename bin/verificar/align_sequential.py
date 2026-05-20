import subprocess
import os

samples = [
    "SC1007_mRNA",
    "SC1015_mRNA",
    "SC1034_mRNA",
    "SC1041_mRNA",
    "SC1055_mRNA",
    "iso_GS_1_mRNA",
    "iso_GS_2_mRNA",
    "iso_MC_1_mRNA",
    "iso_MC_2_mRNA"
]

genome_dir = "/Users/mateuslisboa/Desktop/Data/Potapova_2016_RNAseq/Reference/STAR_Index"
trimmed_dir = "/Users/mateuslisboa/Desktop/Data/LRRK2/trimmed_mRNA"
out_dir = "/Users/mateuslisboa/Desktop/Data/LRRK2/alignment_STAR"

for sample in samples:
    print(f"Starting alignment for {sample}...")
    fastq = os.path.join(trimmed_dir, f"{sample}_trimmed.fastq.gz")
    prefix = os.path.join(out_dir, f"{sample}_")
    
    command = [
        "STAR",
        "--runThreadN", "8",
        "--genomeDir", genome_dir,
        "--readFilesIn", fastq,
        "--readFilesCommand", "gunzip -c",
        "--outFileNamePrefix", prefix,
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--quantMode", "GeneCounts"
    ]
    
    try:
        subprocess.run(command, check=True)
        print(f"Finished alignment for {sample}")
    except subprocess.CalledProcessError as e:
        print(f"Error aligning {sample}: {e}")
        # Stop if there's an error to avoid cascade failures
        break
