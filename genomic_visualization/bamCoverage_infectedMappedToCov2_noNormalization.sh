#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 02:00:00
#SBATCH -J bamCoverage


module load bioinfo-tools
module load deepTools/3.3.2

rev="mappings/infected_cells_mapped_to_sarscov2/infected_nodup_alns/infected_nodup_cov2only/rev"

fwd="mappings/infected_cells_mapped_to_sarscov2/infected_nodup_alns/infected_nodup_cov2only/fwd"

outdir="bigwig/infected_mappedTocov2_noNormalization"

sf_wu_inf_in1=1
sf_wu_inf_in2=1
sf_wu_inf_m6a=1

sf_sa_inf_in1=1
sf_sa_inf_in2=1
sf_sa_inf_m6a=1

sf_uk_inf_in1=1
sf_uk_inf_in2=1
sf_uk_inf_m6a=1

sf_vero_in1=1
sf_vero_in2=1
sf_vero_m6a=1

wu_inf_input1_rev="B1_Infected_cell_Input_1_S7_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
wu_inf_input2_rev="B1_Infected_cell_Input_2_S15_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
wu_inf_m6a_rev="B1_Infected_cell_m6A_S10_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"

sa_inf_input1_rev="B1351_Infected_cell_Input_1_S9_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
sa_inf_input2_rev="B1351_Infected_cell_Input_2_S17_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
sa_inf_m6a_rev="B1351_Infected_cell_m6A_S12_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"

uk_inf_input1_rev="B117_Infected_cell_Input_1_S8_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
uk_inf_input2_rev="B117_Infected_cell_Input_2_S16_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
uk_inf_m6a_rev="B117_Infected_cell_m6A_renamed.hisat2.nodup.uniq.sorted.sarscov2Only.rev.bam"

vero_input1_rev="Vero_Input_1_S13_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
vero_input2_rev="Vero_Input_2_S18_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"
vero_m6a_rev="Vero_m6A_S14_R1_001.hisat2.uniq.sorted.rev.nodup.sorted.sarscov2only.bam"

wu_inf_input1_fwd="B1_Infected_cell_Input_1_S7_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
wu_inf_input2_fwd="B1_Infected_cell_Input_2_S15_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
wu_inf_m6a_fwd="B1_Infected_cell_m6A_S10_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"

sa_inf_input1_fwd="B1351_Infected_cell_Input_1_S9_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
sa_inf_input2_fwd="B1351_Infected_cell_Input_2_S17_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
sa_inf_m6a_fwd="B1351_Infected_cell_m6A_S12_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"

uk_inf_input1_fwd="B117_Infected_cell_Input_1_S8_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
uk_inf_input2_fwd="B117_Infected_cell_Input_2_S16_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
uk_inf_m6a_fwd="B117_Infected_cell_m6A_renamed.hisat2.nodup.uniq.sorted.sarscov2Only.fwd.bam"

vero_input1_fwd="Vero_Input_1_S13_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
vero_input2_fwd="Vero_Input_2_S18_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"
vero_m6a_fwd="Vero_m6A_S14_R1_001.hisat2.uniq.sorted.fwd.nodup.sorted.sarscov2only.bam"

echo "Making coverage files:"

#WU rev
 bamCoverage -b ${rev}/${wu_inf_input1_rev} -o ${outdir}/B1_infected_input1_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${wu_inf_input2_rev} -o ${outdir}/B1_infeced_input2_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${wu_inf_m6a_rev} -o ${outdir}/B1_infected_m6a_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 #SA rev
 bamCoverage -b ${rev}/${sa_inf_input1_rev} -o ${outdir}/B1351_infected_input1_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${sa_inf_input2_rev} -o ${outdir}/B1351_infected_input2_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${sa_inf_m6a_rev} -o ${outdir}/B1351_infected_m6a_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 #UK rev
 bamCoverage -b ${rev}/${uk_inf_input1_rev} -o ${outdir}/B117_infected_input1_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${uk_inf_input2_rev} -o ${outdir}/B117_infected_input2_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${uk_inf_m6a_rev} -o ${outdir}/B117_infected_m6a_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 #Vero rev
 bamCoverage -b ${rev}/${vero_input1_rev} -o ${outdir}/Vero_input1_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${vero_input2_rev} -o ${outdir}/Vero_input2_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${vero_m6a_rev} -o ${outdir}/Vero_m6a_normalizedToSNratio_rev.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}

#bamCoverage Fwd:
#WU fwd
bamCoverage -b ${fwd}/${wu_inf_input1_fwd} -o ${outdir}/B1_infected_input1_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${wu_inf_input2_fwd} -o ${outdir}/B1_infected_input2_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${wu_inf_m6a_fwd} -o ${outdir}/B1_infected_m6a_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
#SA fwd
bamCoverage -b ${fwd}/${sa_inf_input1_fwd} -o ${outdir}/B1351_infected_input1_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${sa_inf_input2_fwd} -o ${outdir}/B1351_infected_input2_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${sa_inf_m6a_fwd} -o ${outdir}/B1351_infected_m6a_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
#UK fwd
bamCoverage -b ${fwd}/${uk_inf_input1_fwd} -o ${outdir}/B117_infected_input1_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${uk_inf_input2_fwd} -o ${outdir}/B117_infected_input2_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${uk_inf_m6a_fwd} -o ${outdir}/B117_infected_m6a_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
#Vero fwd
bamCoverage -b ${fwd}/${vero_input1_fwd} -o ${outdir}/Vero_input1_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${vero_input2_fwd} -o ${outdir}/Vero_input2_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${vero_m6a_fwd} -o ${outdir}/Vero_m6a_normalizedToSNratio_fwd.bedgraph --outFileFormat bedgraph --binSize 1 --verbose -p ${SLURM_NTASKS}

echo "Done."
