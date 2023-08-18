#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 01:30:00
#SBATCH -J macs_peakcalling


module load bioinfo-tools
module load MACS/2.2.6

#inputdir=$(realpath $1) 
#outdir=$(realpath $2)
inputdir="/mappings/infected_cells_mapped_to_sarscov2/infected_nodup_alns/infected_nodup_cov2only"
outdir="peaksInfectedCells_mappedTo_ViralGenome"
#/mappings/infected_cells_mapped_to_sarscov2/infected_nodup_alns

#------ Peak calling Sars-cov2 infected cell mapped to viral genome ---#

gsize=30000
pvalue=0.05
#qvalue=0.01

#Input files Fwd:
input_B1_fwd=${inputdir}/fwd/B1_Infected_cell_Input_1_S7_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam
input_B1351_fwd=${inputdir}/fwd/B1351_Infected_cell_Input_1_S9_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam
input_B117_fwd=${inputdir}/fwd/B117_Infected_cell_Input_1_S8_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam
input_vero_fwd=${inputdir}/fwd/Vero_Input_1_S13_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam

#IP files Fwd:
IP_B1_fwd=${inputdir}/fwd/B1_Infected_cell_m6A_S10_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam
IP_B1351_fwd=${inputdir}/fwd/B1351_Infected_cell_m6A_S12_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam
IP_B117_fwd=${inputdir}/fwd/B117_Infected_cell_m6A_renamed.hisat2.nodup.uniq.sorted.sarscov2Only.fwd.bam
IP_vero_fwd=${inputdir}/fwd/Vero_m6A_S14_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam

#Names Fwd output file sufix:
name_B1_fwd="peaks_nodup_uniq_B1_infectedMappedToCov2_fwd"
name_B1351_fwd="peaks_nodup_uniq_B1351_infectedMappedToCov2_fwd"
name_B117_fwd="peaks_nodup_uniq_B117_infectedMappedToCov2_fwd"
name_vero_fwd="peaks_nodup_uniq_vero_infectedMappedToCov2_fwd"

input_B1_rev=${inputdir}/rev/B1_Infected_cell_Input_1_S7_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam
input_B1351_rev=${inputdir}/rev/B1351_Infected_cell_Input_1_S9_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam
input_B117_rev=${inputdir}/rev/B117_Infected_cell_Input_1_S8_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam
input_vero_rev=${inputdir}/rev/Vero_Input_1_S13_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam

IP_B1_rev=${inputdir}/rev/B1_Infected_cell_m6A_S10_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam
IP_B1351_rev=${inputdir}/rev/B1351_Infected_cell_m6A_S12_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam
IP_B117_rev=${inputdir}/rev/B117_Infected_cell_m6A_renamed.hisat2.nodup.uniq.sorted.sarscov2Only.rev.bam
IP_vero_rev=${inputdir}/rev/Vero_m6A_S14_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam

name_B1_rev="peaks_nodup_uniq_B1_infectedMappedToCov2_rev"
name_B1351_rev="peaks_nodup_uniq_B1351_infectedMappedToCov2_rev"
name_B117_rev="peaks_nodup_uniq_B117_infectedMappedToCov2_rev"
name_vero_rev="peaks_nodup_uniq_vero_infectedMappedToCov2_rev"

mkdir -p ${outdir}/${name_B1_fwd}
mkdir -p ${outdir}/${name_B1351_fwd}
mkdir -p ${outdir}/${name_B117_fwd}
mkdir -p ${outdir}/${name_vero_fwd}

mkdir -p ${outdir}/${name_B1_rev}
mkdir -p ${outdir}/${name_B1351_rev}
mkdir -p ${outdir}/${name_B117_rev}
mkdir -p ${outdir}/${name_vero_rev}

echo "Calling peaks for Fwd samples"
macs2 callpeak -t ${IP_B1_fwd} -c ${input_B1_fwd} --name ${name_B1_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --extsize 75 --outdir ${outdir}/${name_B1_fwd} 2> ${outdir}/${name_B1_fwd}/macs2.${name_B1_fwd}.out
macs2 callpeak -t ${IP_B1351_fwd} -c ${input_B1351_fwd} --name ${name_B1351_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --extsize 75 --outdir ${outdir}/${name_B1351_fwd} 2> ${outdir}/${name_B1351_fwd}/macs2.${name_B1351_fwd}.out
macs2 callpeak -t ${IP_B117_fwd} -c ${input_B117_fwd} --name ${name_B117_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --extsize 75 --outdir ${outdir}/${name_B117_fwd} 2> ${outdir}/${name_B117_fwd}/macs2.${name_B117_fwd}.out
macs2 callpeak -t ${IP_vero_fwd} -c ${input_vero_fwd} --name ${name_vero_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_vero_fwd} 2> ${outdir}/${name_vero_fwd}/macs2.${name_vero_fwd}.out

echo "Calling peaks for Rev samples"
macs2 callpeak -t ${IP_B1_rev} -c ${input_B1_rev} --name ${name_B1_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --extsize 75 --outdir ${outdir}/${name_B1_rev} 2> ${outdir}/${name_B1_rev}/macs2.${name_B1_rev}.out
macs2 callpeak -t ${IP_B1351_rev} -c ${input_B1351_rev} --name ${name_B1351_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --extsize 75 --outdir ${outdir}/${name_B1351_rev} 2> ${outdir}/${name_B1351_rev}/macs2.${name_B1351_rev}.out
macs2 callpeak -t ${IP_B117_rev} -c ${input_B117_rev} --name ${name_B117_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --extsize 75 --outdir ${outdir}/${name_B117_rev} 2> ${outdir}/${name_B117_rev}/macs2.${name_B117_rev}.out
macs2 callpeak -t ${IP_vero_rev} -c ${input_vero_rev} --name ${name_vero_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -p ${pvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_vero_rev} 2> ${outdir}/${name_vero_rev}/macs2.${name_vero_rev}.out
echo "Done."