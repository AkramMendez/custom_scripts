#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 04:00:00
#SBATCH -J macs_peakcalling


module load bioinfo-tools
module load MACS/2.2.6

#inputdir=$(realpath $1) 
#outdir=$(realpath $2)
inputdir="/mappings/covid19_infected_cells/nodup_uniq_alns"
outdir="/peakCalling/macs2_peaks/peaks_infectedCells"


#------ Peak calling Sars-cov2 viral data samples ---#

gsize=2756670577
#pvalue=0.05
qvalue=0.01

#Input files Fwd:
input_WU_fwd=${inputdir}/fwd/COVID_19_Infected_cell_Input_1_S7_R1_001.chroms.nodup.uniq.sorted.fwd.bam
input_SA_fwd=${inputdir}/fwd/COVID_SA_Infected_cell_Input_1_S9_R1_001.chroms.nodup.uniq.sorted.fwd.bam
input_UK_fwd=${inputdir}/fwd/COVID_UK_Infected_cell_Input_1_S8_R1_001.chroms.nodup.uniq.sorted.fwd.bam
input_vero_fwd=${inputdir}/fwd/Vero_Input_1_S13_R1_001.chroms.nodup.uniq.sorted.fwd.bam

#IP files Fwd:
IP_WU_fwd=${inputdir}/fwd/COVID_19_Infected_cell_m6A_S10_R1_001.chroms.nodup.uniq.sorted.fwd.bam
IP_SA_fwd=${inputdir}/fwd/COVID_SA_Infected_cell_m6A_S12_R1_001.chroms.nodup.uniq.sorted.fwd.bam
IP_UK_fwd=${inputdir}/fwd/COVID_UK_Infected_cell_m6A_S11_R1_001_renamed.chroms.nodup.uniq.sorted.fwd.bam
IP_vero_fwd=${inputdir}/fwd/Vero_m6A_S14_R1_001.chroms.nodup.uniq.sorted.fwd.bam

#Names Fwd output file sufix:
name_WU_fwd="peaks_nodup_uniq_WU_infectedAll_fwd"
name_SA_fwd="peaks_nodup_uniq_SA_infectedAll_fwd"
name_UK_fwd="peaks_nodup_uniq_UK_infectedAll_fwd"
name_vero_fwd="peaks_nodup_uniq_vero_infectedAll_fwd"

input_WU_rev=${inputdir}/rev/COVID_19_Infected_cell_Input_1_S7_R1_001.chroms.nodup.uniq.sorted.rev.bam
input_SA_rev=${inputdir}/rev/COVID_SA_Infected_cell_Input_1_S9_R1_001.chroms.nodup.uniq.sorted.rev.bam
input_UK_rev=${inputdir}/rev/COVID_UK_Infected_cell_Input_1_S8_R1_001.chroms.nodup.uniq.sorted.rev.bam
input_vero_rev=${inputdir}/rev/Vero_Input_1_S13_R1_001.chroms.nodup.uniq.sorted.rev.bam

IP_WU_rev=${inputdir}/rev/COVID_19_Infected_cell_m6A_S10_R1_001.chroms.nodup.uniq.sorted.rev.bam
IP_SA_rev=${inputdir}/rev/COVID_SA_Infected_cell_m6A_S12_R1_001.chroms.nodup.uniq.sorted.rev.bam
IP_UK_rev=${inputdir}/rev/COVID_UK_Infected_cell_m6A_S11_R1_001_renamed.chroms.nodup.uniq.sorted.rev.bam
IP_vero_rev=${inputdir}/rev/Vero_m6A_S14_R1_001.chroms.nodup.uniq.sorted.rev.bam

name_WU_rev="peaks_nodup_uniq_WU_infectedAll_rev"
name_SA_rev="peaks_nodup_uniq_SA_infectedAll_rev"
name_UK_rev="peaks_nodup_uniq_UK_infectedAll_rev"
name_vero_rev="peaks_nodup_uniq_vero_infectedAll_rev"

mkdir -p ${outdir}/${name_WU_fwd}
mkdir -p ${outdir}/${name_SA_fwd}
mkdir -p ${outdir}/${name_UK_fwd}
mkdir -p ${outdir}/${name_vero_fwd}

mkdir -p ${outdir}/${name_WU_rev}
mkdir -p ${outdir}/${name_SA_rev}
mkdir -p ${outdir}/${name_UK_rev}
mkdir -p ${outdir}/${name_vero_rev}

echo "Calling peaks for Fwd samples"
macs2 callpeak -t ${IP_WU_fwd} -c ${input_WU_fwd} --name ${name_WU_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_WU_fwd} 2> ${outdir}/${name_WU_fwd}/macs2.${name_WU_fwd}.out
macs2 callpeak -t ${IP_SA_fwd} -c ${input_SA_fwd} --name ${name_SA_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_SA_fwd} 2> ${outdir}/${name_SA_fwd}/macs2.${name_SA_fwd}.out
macs2 callpeak -t ${IP_UK_fwd} -c ${input_UK_fwd} --name ${name_UK_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_UK_fwd} 2> ${outdir}/${name_UK_fwd}/macs2.${name_UK_fwd}.out
macs2 callpeak -t ${IP_vero_fwd} -c ${input_vero_fwd} --name ${name_vero_fwd} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_vero_fwd} 2> ${outdir}/${name_vero_fwd}/macs2.${name_vero_fwd}.out

echo "Calling peaks for Rev samples"
macs2 callpeak -t ${IP_WU_rev} -c ${input_WU_rev} --name ${name_WU_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_WU_rev} 2> ${outdir}/${name_WU_rev}/macs2.${name_WU_rev}.out
macs2 callpeak -t ${IP_SA_rev} -c ${input_SA_rev} --name ${name_SA_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_SA_rev} 2> ${outdir}/${name_SA_rev}/macs2.${name_SA_rev}.out
macs2 callpeak -t ${IP_UK_rev} -c ${input_UK_rev} --name ${name_UK_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_UK_rev} 2> ${outdir}/${name_UK_rev}/macs2.${name_UK_rev}.out
macs2 callpeak -t ${IP_vero_rev} -c ${input_vero_rev} --name ${name_vero_rev} --format="BAM" --gsize=${gsize} --nomodel --bdg -q ${qvalue} --keep-dup auto --call-summits --outdir ${outdir}/${name_vero_rev} 2> ${outdir}/${name_vero_rev}/macs2.${name_vero_rev}.out
echo "Done."