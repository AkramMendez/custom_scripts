#!/bin/bash -l
# Script for fetching TCGA data using manifest files, the script searches for available genomic alignments in a TCGA manifest file and extracts the raw reads into Fastq files
# Extracted raw reads (detected as paired- or single-end reads) are further remapped to the T2T CHM13 human genome assembly using HISAT2.

manifest=$(realpath $1)
n=8 # Number of threads
nSamples=30 # Number of samples per batch to process
token=$(realpath $2) # Token txt file path
hisat2_index_basename=$(realpath $3) #HISAT2 index basename
raw_data=$(realpath $4) # Output directory to deposit extracted reads



cd ${raw_data}

echo "Downloading data from ${manifest} file:"

gdc-client download -m <(awk 'BEGIN{FS=OFS="\t"}{if(NR==1 || $2 ~ /genomic/){print $0}}' ${manifest} | head -${nSamples}) -t ${token}


isPaired(){
# Check if BAM file is paired-end
    local paired=$(samtools view -c -f 1 $1 -@ 2)
    echo "${paired}"
}

extractReads(){
    if [[ $is_paired -gt 0 ]]
    then 

    echo "Sorting BAM file"
    samtools sort -o ${bam%.bam}.sortedByName.bam -n ${bam} -@ ${n}
    echo "Getting paired-end reads from BAM: ${bam}"
    bamToFastq -i ${bam%.bam}.sortedByName.bam -fq ${bam%.bam}_R1.fq -fq2 ${bam%.bam}_R2.fq

    echo "Compressing FASTQ files"
    gzip ${bam%.bam}_R1.fq
    gzip ${bam%.bam}_R2.fq
    
    else 

    echo "Sorting BAM file"
    samtools sort -o ${bam%.bam}.sortedByName.bam -n ${bam} -@ ${n}
    echo "Getting single-end reads from BAM: ${bam}"
    bamToFastq -i ${bam%.bam}.sortedByName.bam -fq ${bam%.bam}_R1.fq

    echo "Compressing FASTQ files"
    gzip ${bam%.bam}_R1.fq

    fi

    echo "Cleaning intermediate files"
    rm ${bam%.bam}.sortedByName.bam
    rm ${bam}

    echo "Extracting reads finished."
}

mapReadsPairedEndToT2T(){
    sample_name=$(basename $1)
    R1=$1
	R2=$2
	
	echo "Mapping paired-end reads to index genome using HISAT2"
	
	hisat2 -p ${n} -x ${hisat2_index_basename} --no-discordant --no-mixed -I 1 -X 1000 -1 ${R1} -2 ${R2} --summary-file hisat2_summary_${sample_name%.fq.gz}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name%.fq.gz}.hisat2.sorted.bam -@ ${n}

	samtools index ${sample_name%.fq.gz}.hisat2.sorted.bam -@ ${n}
}

mapReadsSingleEndToT2T(){
    sample_name=$(basename $1)
    R1=$1

	echo "Mapping single-end reads to index genome using HISAT2"
	hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --summary-file hisat2_summary_${sample_name%.fq.gz}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name%.fq.gz}.hisat2.sorted.bam -@ ${n}
    
}


for bam in $(ls *.bam)
do
    echo "Processing BAM file: ${bam}"
    isPaired ${bam} && extractReads ${bam}

    if [[ $is_paired -gt 0 ]]
    then
    mapReadsPairedEndToT2T ${bam%.bam}_R1.fq.gz ${bam%.bam}_R2.fqz

    else
    mapReadsSingleEndToT2T ${bam%.bam}_R1.fq.gz

    fi
done