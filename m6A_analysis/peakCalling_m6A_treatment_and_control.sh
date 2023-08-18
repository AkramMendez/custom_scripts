#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J macs_peakcalling

module load bioinfo-tools
module load MACS/2.2.6

inputdir_fwd=$(realpath $1) #Path of BAM files mapped to the forward strand
inputdir_rev=$(realpath $2) #Path of BAM files mapped to the reverse strand
outdir=$(realpath $3) # Output directory path
gsize=$4 # Effective Genome Size
extsize=$5 # Extension size

# Script for performing m6A peak calling using MACS2 from stranded sequencing short-reads data corresponding to 'Control' and 'Treated' conditions.

#Transcriptome size Hg38: UCSC Table Browser > GENCODE v3 > knownGene > summarystatistics > blocl total:
#374732024
#gsize=2.7e9
#gsize=374732024
#pvalue=0.01
#extsize=75

cd ${SNIC_TMP}

cp ${inputdir_fwd}/*.bam .
cp ${inputdir_rev}/*.bam .

cp ${inputdir_fwd}/*.bai .
cp ${inputdir_rev}/*.bai .


callPeaks(){
	IP=$1
	input=$2
	name=$3
	gsize=$4
	extsize=$5
	#pvalue=$5
	outdir=$6
	

	echo "Calling peaks: ${IP} and ${input}:"

	#macs2 callpeak -t ${IP} -c ${input} --name ${name}_p${pvalue}_ext${extsize} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --extsize ${extsize} --keep-dup auto --call-summits --outdir ${outdir}/${name} 2> macs2.${name}.out
	macs2 callpeak -t ${IP} -c ${input} --name ${name}_ext${extsize} --format="BAM" --gsize=${gsize} --nomodel --bdg --cutoff-analysis --extsize ${extsize} --keep-dup auto --call-summits --outdir ${outdir}/${name} 2> macs2.${name}.out
}

#### Alignments Fwd:
ctrl_input_fwd="ctrl_input_fwd.bam"
ctrl_m6AIP_fwd="ctrl_m6AIP_fwd.bam"

treated_input_fwd="treated_input_fwd.bam"
treated_m6AIP_fwd="treated_m6AIP_fwd.bam"

#### Alignments Rev:
ctrl_input_rev="ctrl_input_rev.bam"
ctrl_m6AIP_rev="ctrl_m6AIP_rev.bam"

treated_input_rev="treated_input_rev.bam"
treated_m6AIP_rev="treated_m6AIP_rev.bam"

#### Output directory names:

name_ctrl_fwd="control_fwd"
name_treated_fwd="treated_fwd"

name_ctrl_rev="control_rev"
name_treated_rev="treated_rev"

echo "Performing peak calling:"

callPeaks ${ctrl_m6AIP_fwd} ${ctrl_input_fwd} ${name_ctrl_fwd} ${gsize} ${extsize} ${outdir}
callPeaks ${ctrl_m6AIP_rev} ${ctrl_input_rev} ${name_ctrl_rev} ${gsize} ${extsize} ${outdir}

callPeaks ${treated_m6AIP_fwd} ${treated_input_fwd} ${name_treated_fwd} ${gsize} ${extsize} ${outdir}
callPeaks ${treated_m6AIP_rev} ${treated_input_rev} ${name_treated_rev} ${gsize} ${extsize} ${outdir}

cp *.out ${outdir}

echo "Done."

