#!/usr/bin/bash

inputdir=$(realpath $1) #Input directory path containing bascalled fastq.gz files
multi_reads_dir=$(realpath $2) # Output directory path

# Exit the shell script with a status of 1 using exit 1 command.
[ $# -eq 0 ] && { echo "Usage: $0 inputdir multi_reads_dir"; exit 1; }

# This script concatenate all basecalled FASTQ files (generated after basecalling dRNA Oxford Nanopore data using guppy for example) 
# Concatenated fastq files are then converted RNA to DNA and compressed into a single FASTQ file.

echo "Concatenating and converting FASTQ files (RNA to DNA):"

cat ${inputdir}/*.fastq.gz | seqkit seq --rna2dna | gzip >> ${filename_base}.fastq.gz

if [[ -d ${multi_reads_dir} ]];then
    echo "Directory exists"
else
    echo "Creating ${multi_reads_dir} directory"
    
    mkdir -p ${multi_reads_dir}
    
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory."
    fi

fi

echo "Done."
