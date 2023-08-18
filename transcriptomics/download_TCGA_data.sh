#!/bin/bash -l
# Script for fetching TCGA data

manifest=$(realpath $1)
n=8 # Number of threads
nSamples=30
token="/home/akram/gdc-user-token.2023-03-15T09_50_54.309Z.txt"
raw_data=$(realpath /data3/akram/tcga_data)
chrEndCoords="/data3/akram/ref/chm13v2/chrEndCoords_chm13v2.bed"
hisat2_index_basename="/data3/akram/ref/chm13v2/hisat2idx_t2tv2.0/chm13v2.0"
outdir_path=$(realpath /data3/akram/tcga_filtered_mappings)

outname=$(basename ${manifest} | sed -E 's/(gdc_manifest_TCGA_|.txt)//g')

outdir="${outdir_path}/${outname}"


cd ${raw_data}



    #echo "[1] Downloading data from ${manifest} file:"
    #echo "-------------------------------------------"
    #gdc-client download ${id} -t ${token} --n-processes 8
    #gdc-client download -m <(awk 'BEGIN{FS=OFS="\t"}{if(NR==1 || $2 ~ /genomic/){print $0}}' ${manifest} | head -${nSamples}) -t ${token} --n-processes 8


isPaired(){
# Check if BAM file is paired-end
    is_paired=$(samtools view -c -f 1 $1 -@ 2)
    echo "${is_paired}"
}

extractReads(){
    if [[ $is_paired -gt 0 ]]
    then

    bam=$(basename ${bamfile} | sed 's/.rna_seq.genomic.gdc_realn//g')

    echo "Sorting BAM file"
    samtools sort -o ${id}.sortedByName.bam -n ${bamfile} -@ ${n}
    echo "[2] Extracting paired-end reads from BAM: ${bamfile}"
    echo "-------------------------------------------"
    bamToFastq -i ${id}.sortedByName.bam -fq ${id}_R1.fq -fq2 ${id}_R2.fq

    echo "Compressing FASTQ files"
    gzip ${id}_R1.fq
    gzip ${id}_R2.fq
    
    else 

    echo "Sorting BAM file"
    samtools sort -o ${id}.sortedByName.bam -n ${bamfile} -@ ${n}
    echo "Extracting single-end reads from BAM: ${bam}"
    echo "-------------------------------------------"
    bamToFastq -i ${id}.sortedByName.bam -fq ${id}_R1.fq

    echo "Compressing FASTQ files"
    gzip ${id}_R1.fq

    fi

    echo "Cleaning intermediate files"
    rm ${id}.sortedByName.bam
    rm ${bamfile}

    echo "Extracting reads finished."
}

mapReadsSingleEndToT2T(){
    sample_name=$(basename $1)
    R1=$1

	echo "[3] Mapping single-end reads to index genome using HISAT2"
    echo "-------------------------------------------"
	hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --summary-file hisat2_summary_${sample_name%_R1.fq.gz}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name%_R1.fq.gz}.hisat2.sorted.bam -@ ${n}

    echo "Indexing BAM file:"
	samtools index ${sample_name%_R1.fq.gz}.hisat2.sorted.bam -@ ${n}

    echo "Cleaning FASTQ files"
    rm ${R1}
    
}

mapReadsPairedEndToT2T(){
    sample_name=$(basename $1)
    R1=$1
	R2=$2
	
	echo "[3] Mapping paired-end reads to index genome using HISAT2"
    echo "-------------------------------------------"
	
	hisat2 -p ${n} -x ${hisat2_index_basename} --no-discordant --no-mixed -I 1 -X 1000 -1 ${R1} -2 ${R2} --summary-file hisat2_summary_${sample_name%_R1.fq.gz}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name%_R1.fq.gz}.hisat2.sorted.bam -@ ${n}

    echo "Indexing BAM file:"
	samtools index ${sample_name%_R1.fq.gz}.hisat2.sorted.bam -@ ${n}

    echo "Cleaning FASTQ files"
    rm ${R1}
    rm ${R2}
}

filterMappingsSingle(){
sorted_bam=$(realpath $1)
sorted_name=$(basename ${sorted_bam})
echo "[4] Filtering BAM files to retrieve chromosome end coordinates (50kb)"
echo "-------------------------------------------"
samtools view -h -b -L ${chrEndCoords} ${sorted_bam} -@ ${n} | samtools sort -o ${outdir}/${sorted_name%.bam}_chrEnds50kb.bam -@ ${n}
echo "[4] Filtering BAM files to retrieve unmmaped reads"
echo "-------------------------------------------"
samtools view -h -b -f 0x4 ${sorted_bam} -@ ${n} > ${outdir}/${sorted_name%.bam}.unmapped.bam 

echo "Cleaning intermediate BAM files"
rm ${sorted_bam}
rm ${sorted_bam}.bai
}


filterMappingsPaired(){
sorted_bam=$(realpath $1)
sorted_name=$(basename ${sorted_bam})
echo "[4] Filtering BAM files to retrieve chromosome end coordinates (50kb)"
echo "-------------------------------------------"
samtools view -h -b -f 2 -L ${chrEndCoords} ${sorted_bam} -@ ${n} | samtools sort -o ${outdir}/${sorted_name%.bam}_chrEnds50kb.bam -@ ${n}
echo "[4] Filtering BAM files to retrieve unmmaped reads"
echo "-------------------------------------------"
samtools view -h -b -f 0x4 ${sorted_bam} -@ ${n} > ${outdir}/${sorted_name%.bam}.unmapped.bam 

echo "Cleaning intermediate BAM files"
rm ${sorted_bam}
rm ${sorted_bam}.bai

}

#####

#for bamfile in $(find . -iname "*.bam")
#do
while read -r line; do

    columns=(${line}) 
    
    id=${columns[0]}
    bamfile="${id}/${columns[1]}"

    echo "Processing UUID file: ${id}"
    
    #echo "[1] Downloading data from ${manifest} file:"
    #echo "-------------------------------------------"
    gdc-client download ${id} -t ${token} --n-processes 8 && isPaired ${bamfile} && extractReads ${bamfile}

    if [[ $is_paired -gt 0 ]]
    then
    mapReadsPairedEndToT2T ${id}_R1.fq.gz ${id}_R2.fq.gz && filterMappingsPaired ${sample_name%_R1.fq.gz}.hisat2.sorted.bam

    else
    mapReadsSingleEndToT2T ${id}_R1.fq.gz && filterMappingsSingle ${sample_name%_R1.fq.gz}.hisat2.sorted.bam

    fi

#done
done < <(awk 'BEGIN{FS=OFS="\t"}{if($2 ~ /genomic/){print $0}}' ${manifest})

echo "-------------------------------------------"
echo "Done."