#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 08:00:00
#SBATCH -J macs_peakcalling

# Script to perform peak calling using MACS2 from m6A-RIP and input samples.
module load bioinfo-tools
module load MACS/2.2.6

ip_sample=$(realpath $1) # Path to m6A-RIP BAM file
input_sample=$(realpath $2) # Path to input BAM file
output_name=$3 # Output name
outdir=$(realpath $4) # Output directory
gsize=$5 # Genome size
extsize=$6 # Extension size

cd ${SNIC_TMP}

cp ${ip_sample} .
cp ${input_sample} .

callPeaks(){
	IP=$1
	input=$2
	name=$3
	gsize=$4
	extsize=$5
	outdir=$6
	
	echo "Calling peaks: ${IP} and ${input}:"
	macs2 callpeak -t ${IP} -c ${input} --name ${name}_ext${extsize} --format="BAM" --gsize=${gsize} --nomodel --bdg --cutoff-analysis --extsize ${extsize} --keep-dup auto --call-summits --outdir ${outdir}/${name} 2> ${outdir}/macs2.${name}.out
}

echo "Performing peak calling"
callPeaks ${ip_sample} ${input_sample} ${output_name} ${gsize} ${extsize} ${outdir}
echo "Done."

