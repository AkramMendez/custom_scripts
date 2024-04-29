#!/usr/bin/env python
import sys
import os
import argparse
import pysam
import gzip
sys.path.append("/home/akram/gitrepos/bioinfopipes/python")
import bam_parsing.bam_parsing_functions as bp


parser = argparse.ArgumentParser(description='Script for filtering unmapped reads from a BAM file')

# Add arguments for input and output files
parser.add_argument('input_file', type=str, help='Path to the input file.')
parser.add_argument('output_file', type=str, help='Name for the output file (it is saved to the same path as the input BAM file by default).')

# Parse the command-line arguments
args = parser.parse_args()

input_filepath = args.input_file #"/home/akram/Documents/data/bam/tNCC_Day7_Input_S31_R1_001.ht2.srt.bam"
output_filename = args.output_file
print("Input:",input_filepath)
print("Output:", output_filename)

input_bam = pysam.AlignmentFile(input_filepath,"rb") #sys.argv[0] #"/home/akram/Documents/data/bam/tNCC_Day7_Input_S31_R1_001.ht2.srt.bam"
output_filepath= os.path.join(os.path.dirname(input_filepath),output_filename) #"/home/akram/Documents/data/bam/output_unmapped.bam" sys.argv[1]

#input_bam = pysam.AlignmentFile(sys.argv[0],"rb")
#output_bam = pysam.AlignmentFile(output_filepath,"wb",template=input_bam)

#bp.extractUnmappedReads(input_bam,output_filepath)

#output_bam = pysam.AlignmentFile(output_filepath,"rb")
bp.bamToFastQ(output_filepath)



""" bam = pysam.AlignmentFile(output_filepath,"rb")
output_filename = output_filepath.replace(".bam",".fastq.gz") """


""" print("Converting BAM to FASTQ.gz")
with gzip.open(output_filename,'wt') as output_fastq:
            for read in bam.fetch(until_eof=True):
                output_fastq.write(pysam.fastq(output_filepath))
        
print("Done.")
 """