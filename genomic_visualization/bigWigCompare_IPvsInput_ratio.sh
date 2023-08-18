#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J bamCoverage


module load bioinfo-tools
module load deepTools/3.3.2

inputdir=$(realpath $1) # Path to input directory containing the strand-separated bigwig tracks (fwd and rev)
outdir=$(realpath $2) # Output directory

# Script for calculating the ratio of m6A-IP vs Input coverage tracks using deepTools bigwigCompare fucntion.

cd ${SNIC_TMP}

cp ${inputdir}/*_rev.bw .
cp ${inputdir}/*_fwd.bw .

# Input and m6A-IP track bigwig files
cond1_input_rev="cond1_input_rev.bw"
cond1_m6aIP_rev="cond1_m6AIP_rev.bw"

cond1_input_fwd="cond1_input_fwd.bw"
cond1_m6aIP_fwd="cond1_m6AIP_fwd.bw"

# Output file names
name_cond1_rev="cond1_rev" 
name_cond1_fwd="cond1_fwd"

	
getBigwigRatio(){
	
	b1=$1
	b2=$2
	name=$3

	bigwigCompare --bigwig1 ${b1} --bigwig2 ${b2} --binSize 10 --outFileName ${name}_ratio_bs10.bw --operation ratio --outFileFormat bigwig --numberOfProcessors ${SLURM_NTASKS} --verbose 
}


echo "Making BigWig files"

getBigwigRatio ${cond1_m6aIP_rev} ${cond1_input_rev} ${name_cond1_rev}
getBigwigRatio ${cond1_m6aIP_fwd} ${cond1_input_fwd} ${name_cond1_fwd}

echo "Copying files to ${outdir}"

cp *_ratio_bs10.bw ${outdir}

echo "Done."