"""
Gets the context of each variant in the input TSV file.

Usage:
    python context.py --input_file <input_file> \
        --reference_fasta <reference_fasta> --output_file <output_file> \

"""

import sys
import os
import argparse
import pandas as pd
from tqdm import tqdm
from pyfaidx import Fasta

def get_context(reference, chrom, pos, length):
    start = pos - (length // 2)
    end = pos + (length // 2)
    return reference[chrom][start:end].seq

def main():
    """
    Main entry point
    """


    input_file = args.input_file
    reference_fasta = args.reference_fasta
    output_file = args.output_file
    assembly_report = args.assembly_report

    # Check if the input file, reference FASTA file, and assembly report file exist
    if not os.path.exists(input_file):
        sys.exit(f"Input file {input_file} does not exist.")
    if not os.path.exists(reference_fasta):
        sys.exit(f"Reference FASTA file {reference_fasta} does not exist.")
    if not os.path.exists(assembly_report):
        sys.exit(f"Assembly report file {assembly_report} does not exist.")

    # Check if output file already exists
    if os.path.exists(output_file):
        sys.exit(f"Output file {output_file} already exists.")


    # Load the reference genome
    reference = Fasta(reference_fasta)
    
    # Context lengths
    context_lengths = [3, 5, 7, 9, 12, 15]
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Header for the output TSV
        header = ["chr", "pos", "variant_id"] + [f"context_{length}bp" for length in context_lengths]
        outfile.write("\t".join(header) + "\n")
        
        for line in tqdm(infile.readlines()):
            chrom, pos, variant = line.strip().split()
            pos = int(pos)
            
            # Obtain the context for each length and append to contexts list
            contexts = [get_context(reference, chrom, pos, length) for length in context_lengths]
            
            # Write to the output TSV
            outfile.write("\t".join([chrom, str(pos), variant] + contexts) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gets the context of each variant in the input TSV file.")
    parser.add_argument("input_file", type=str, help="Input TSV file with variants.")
    parser.add_argument("reference_fasta", type=str, help="Reference genome FASTA file.")
    parser.add_argument("assembly_report", type=str, help="Assembly report file.")
    parser.add_argument("output_file",type=str, help="Output TSV file with variant contexts.")    
    main(parser.parse_args())
