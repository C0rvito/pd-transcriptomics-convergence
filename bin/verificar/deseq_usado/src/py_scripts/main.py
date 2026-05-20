import argparse
import sys
import os
from preprocessing import run_star_alignment, create_count_matrix, add_gene_symbols
from postprocessing import annotate_de_results, rank_and_segment_results, analyze_vgcc_family
from reporting import create_excel_summary

def main():
    parser = argparse.ArgumentParser(description="RNA-Seq Pipeline Helper (Python Module)")
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands")

    # Preprocess
    preprocess_parser = subparsers.add_parser("preprocess", help="Preprocess steps")
    preprocess_parser.add_argument("--step", choices=["align", "matrix", "symbols"], required=True)
    preprocess_parser.add_argument("--samples", nargs="+", help="List of samples (for align)")
    preprocess_parser.add_argument("--genome-dir", help="STAR Genome directory")
    preprocess_parser.add_argument("--input-dir", help="Input directory")
    preprocess_parser.add_argument("--output", required=True, help="Output file or directory")
    preprocess_parser.add_argument("--gtf", help="GTF file path (for symbols)")
    preprocess_parser.add_argument("--threads", type=int, default=8)

    # Postprocess
    postprocess_parser = subparsers.add_parser("postprocess", help="Postprocess steps")
    postprocess_parser.add_argument("--step", choices=["annotate", "rank", "vgcc"], required=True)
    postprocess_parser.add_argument("--input", required=True, help="Input CSV file")
    postprocess_parser.add_argument("--map", help="Symbols mapping CSV (for annotate)")
    postprocess_parser.add_argument("--output", required=True, help="Output file")
    postprocess_parser.add_argument("--metadata", help="Metadata CSV (for rank)")
    postprocess_parser.add_argument("--padj", type=float, default=0.05)

    # Report
    report_parser = subparsers.add_parser("report", help="Reporting steps")
    report_parser.add_argument("--output", required=True, help="Output XLSX file")
    report_parser.add_argument("--inputs", nargs="+", help="List of SheetName:FilePath pairs (e.g. DE:results.csv)")
    report_parser.add_argument("--metadata", help="Metadata CSV")

    args = parser.parse_args()

    if args.command == "preprocess":
        if args.step == "align":
            run_star_alignment(args.samples, args.genome_dir, args.input_dir, args.output, args.threads)
        elif args.step == "matrix":
            create_count_matrix(args.input_dir, args.output)
        elif args.step == "symbols":
            add_gene_symbols(args.input_dir, args.gtf, args.output)

    elif args.command == "postprocess":
        if args.step == "annotate":
            annotate_de_results(args.input, args.map, args.output)
        elif args.step == "rank":
            rank_and_segment_results(args.input, args.output, args.metadata, args.padj)
        elif args.step == "vgcc":
            analyze_vgcc_family(args.input, args.output)

    elif args.command == "report":
        files_dict = {}
        for item in args.inputs:
            name, path = item.split(":")
            files_dict[name] = path
        create_excel_summary(files_dict, args.output, args.metadata)

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
