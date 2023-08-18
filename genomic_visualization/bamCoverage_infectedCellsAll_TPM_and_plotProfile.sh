#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 12:00:00
#SBATCH -J bamCoverage


module load bioinfo-tools
module load deepTools/3.3.2

ref="/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"

all="/mappings/covid19_infected_cells/nodup_uniq_alns/all"

rev="/mappings/covid19_infected_cells/nodup_uniq_alns/rev"

fwd="/mappings/covid19_infected_cells/nodup_uniq_alns/fwd"

outdir="/bigwig/infectedCellsAll_TPM_all"

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

#All (forward and reverse)
wu_inf_input1_all="COVID_19_Infected_cell_Input_1_S7_R1_001.chroms.nodup.uniq.sorted.bam"
wu_inf_input2_all="COVID_19_Infected_cell_Input_2_S15_R1_001.chroms.nodup.uniq.sorted.bam"
wu_inf_m6a_all="COVID_19_Infected_cell_m6A_S10_R1_001.chroms.nodup.uniq.sorted.bam"

sa_inf_input1_all="COVID_SA_Infected_cell_Input_1_S9_R1_001.chroms.nodup.uniq.sorted.bam"
sa_inf_input2_all="COVID_SA_Infected_cell_Input_2_S17_R1_001.chroms.nodup.uniq.sorted.bam"
sa_inf_m6a_all="COVID_SA_Infected_cell_m6A_S12_R1_001.chroms.nodup.uniq.sorted.bam"

uk_inf_input1_all="COVID_UK_Infected_cell_Input_1_S8_R1_001.chroms.nodup.uniq.sorted.bam"
uk_inf_input2_all="COVID_UK_Infected_cell_Input_2_S16_R1_001.chroms.nodup.uniq.sorted.bam"
uk_inf_m6a_all="COVID_UK_Infected_cell_m6A_S11_R1_001_renamed.chroms.nodup.uniq.sorted.bam"

vero_input1_all="Vero_Input_1_S13_R1_001.chroms.nodup.uniq.sorted.bam"
vero_input2_all="Vero_Input_2_S18_R1_001.chroms.nodup.uniq.sorted.bam"
vero_m6a_all="Vero_m6A_S14_R1_001.chroms.nodup.uniq.sorted.bam"

#Reverse
wu_inf_input1_rev="COVID_19_Infected_cell_Input_1_S7_R1_001.chroms.nodup.uniq.sorted.rev.bam"
wu_inf_input2_rev="COVID_19_Infected_cell_Input_2_S15_R1_001.chroms.nodup.uniq.sorted.rev.bam"
wu_inf_m6a_rev="COVID_19_Infected_cell_m6A_S10_R1_001.chroms.nodup.uniq.sorted.rev.bam"

sa_inf_input1_rev="COVID_SA_Infected_cell_Input_1_S9_R1_001.chroms.nodup.uniq.sorted.rev.bam"
sa_inf_input2_rev="COVID_SA_Infected_cell_Input_2_S17_R1_001.chroms.nodup.uniq.sorted.rev.bam"
sa_inf_m6a_rev="COVID_SA_Infected_cell_m6A_S12_R1_001.chroms.nodup.uniq.sorted.rev.bam"

uk_inf_input1_rev="COVID_UK_Infected_cell_Input_1_S8_R1_001.chroms.nodup.uniq.sorted.rev.bam"
uk_inf_input2_rev="COVID_UK_Infected_cell_Input_2_S16_R1_001.chroms.nodup.uniq.sorted.rev.bam"
uk_inf_m6a_rev="COVID_UK_Infected_cell_m6A_S11_R1_001_renamed.chroms.nodup.uniq.sorted.rev.bam"

vero_input1_rev="Vero_Input_1_S13_R1_001.chroms.nodup.uniq.sorted.rev.bam"
vero_input2_rev="Vero_Input_2_S18_R1_001.chroms.nodup.uniq.sorted.rev.bam"
vero_m6a_rev="Vero_m6A_S14_R1_001.chroms.nodup.uniq.sorted.rev.bam"

#Forward
wu_inf_input1_fwd="COVID_19_Infected_cell_Input_1_S7_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
wu_inf_input2_fwd="COVID_19_Infected_cell_Input_2_S15_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
wu_inf_m6a_fwd="COVID_19_Infected_cell_m6A_S10_R1_001.chroms.nodup.uniq.sorted.fwd.bam"

sa_inf_input1_fwd="COVID_SA_Infected_cell_Input_1_S9_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
sa_inf_input2_fwd="COVID_SA_Infected_cell_Input_2_S17_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
sa_inf_m6a_fwd="COVID_SA_Infected_cell_m6A_S12_R1_001.chroms.nodup.uniq.sorted.fwd.bam"

uk_inf_input1_fwd="COVID_UK_Infected_cell_Input_1_S8_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
uk_inf_input2_fwd="COVID_UK_Infected_cell_Input_2_S16_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
uk_inf_m6a_fwd="COVID_UK_Infected_cell_m6A_S11_R1_001_renamed.chroms.nodup.uniq.sorted.fwd.bam"

vero_input1_fwd="Vero_Input_1_S13_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
vero_input2_fwd="Vero_Input_2_S18_R1_001.chroms.nodup.uniq.sorted.fwd.bam"
vero_m6a_fwd="Vero_m6A_S14_R1_001.chroms.nodup.uniq.sorted.fwd.bam"

echo "Making coverage files:"
#All (Fwd and Rev):
#WU all
 bamCoverage -b ${all}/${wu_inf_input1_all} -o ${outdir}/WU_infected_input1_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${wu_inf_input2_all} -o ${outdir}/WU_infected_input2_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${wu_inf_m6a_all} -o ${outdir}/WU_infected_m6a_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 #SA all
 bamCoverage -b ${all}/${sa_inf_input1_all} -o ${outdir}/SA_infected_input1_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${sa_inf_input2_all} -o ${outdir}/SA_infected_input2_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${sa_inf_m6a_all} -o ${outdir}/SA_infected_m6a_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 #UK all
 bamCoverage -b ${all}/${uk_inf_input1_all} -o ${outdir}/UK_infected_input1_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${uk_inf_input2_all} -o ${outdir}/UK_infected_input2_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${uk_inf_m6a_all} -o ${outdir}/UK_infected_m6a_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 #Vero all
 bamCoverage -b ${all}/${vero_input1_all} -o ${outdir}/Vero_input1_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${vero_input2_all} -o ${outdir}/Vero_input2_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${all}/${vero_m6a_all} -o ${outdir}/Vero_m6a_TPM_all.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}


#WU rev
 bamCoverage -b ${rev}/${wu_inf_input1_rev} -o ${outdir}/WU_infected_input1_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${wu_inf_input2_rev} -o ${outdir}/WU_infected_input2_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${wu_inf_m6a_rev} -o ${outdir}/WU_infected_m6a_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 #SA rev
 bamCoverage -b ${rev}/${sa_inf_input1_rev} -o ${outdir}/SA_infected_input1_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${sa_inf_input2_rev} -o ${outdir}/SA_infected_input2_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${sa_inf_m6a_rev} -o ${outdir}/SA_infected_m6a_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 #UK rev
 bamCoverage -b ${rev}/${uk_inf_input1_rev} -o ${outdir}/UK_infected_input1_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${uk_inf_input2_rev} -o ${outdir}/UK_infected_input2_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${uk_inf_m6a_rev} -o ${outdir}/UK_infected_m6a_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 #Vero rev
 bamCoverage -b ${rev}/${vero_input1_rev} -o ${outdir}/Vero_input1_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${vero_input2_rev} -o ${outdir}/Vero_input2_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
 bamCoverage -b ${rev}/${vero_m6a_rev} -o ${outdir}/Vero_m6a_TPM_rev.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}

# #bamCoverage Fwd:
#WU fwd
bamCoverage -b ${fwd}/${wu_inf_input1_fwd} -o ${outdir}/WU_infected_input1_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${wu_inf_input2_fwd} -o ${outdir}/WU_infected_input2_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${wu_inf_m6a_fwd} -o ${outdir}/WU_infected_m6a_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
#SA fwd
bamCoverage -b ${fwd}/${sa_inf_input1_fwd} -o ${outdir}/SA_infected_input1_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${sa_inf_input2_fwd} -o ${outdir}/SA_infected_input2_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${sa_inf_m6a_fwd} -o ${outdir}/SA_infected_m6a_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
#UK fwd
bamCoverage -b ${fwd}/${uk_inf_input1_fwd} -o ${outdir}/UK_infected_input1_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${uk_inf_input2_fwd} -o ${outdir}/UK_infected_input2_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${uk_inf_m6a_fwd} -o ${outdir}/UK_infected_m6a_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
#Vero fwd
bamCoverage -b ${fwd}/${vero_input1_fwd} -o ${outdir}/Vero_input1_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${vero_input2_fwd} -o ${outdir}/Vero_input2_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}
bamCoverage -b ${fwd}/${vero_m6a_fwd} -o ${outdir}/Vero_m6a_TPM_fwd.bw --outFileFormat bigwig --binSize 10 --normalizeUsing BPM --verbose -p ${SLURM_NTASKS}

cd $SNIC_TMP

computeMatrix scale-regions -S $(ls ${outdir}/*_all.bw | tr "\n" " ") -R ${ref} -a 3000 -b 3000 --regionBodyLength 3000 --outFileName scale_regions_infectedCells_TPM_All.mat.gz --outFileNameMatrix scale_regions_infectedCells_TPM_All.tabular --samplesLabel $(ls ${outdir}/*_all.bw | sed -E 's/(_infected_|_TPM_|.bw)//g' | tr "\n" " ") --verbose -p ${SLURM_NTASKS} --skipZeros

plotProfile --matrixFile scale_regions_infectedCells_TPM_All.mat.gz --perGroup --plotFileFormat pdf --outFileName ${outdir}/profile_scaleRegions_TPM_infectedCells_all.pdf --plotTitle "" --regionsLabel "" --outFileNameData scale_regions_TPM_infectedCells_All.tsv

# #Rev
# computeMatrix scale-regions -S $(ls ${outdir}/*_rev.bw | tr "\n" " ") -R ${ref} -a 3000 -b 3000 --regionBodyLength 3000 --outFileName scale_regions_infectedCells_TPM__Rev.mat.gz --outFileNameMatrix scale_regions_infectedCells_TPM_Rev.tabular --samplesLabel $(ls ${outdir}/*_rev.bw | sed -E 's/(_m6avsInput1_TPM|.bw)//g' | tr "\n" " ") --verbose -p 8 --skipZeros

# plotProfile --matrixFile scale_regions_infectedCells_TPM__Rev.mat.gz --perGroup --plotFileFormat pdf --outFileName ${outdir}/profile_scaleRegions_TPM_infectedCells_rev.pdf --plotTitle "" --regionsLabel "" --outFileNameData scale_regions_TPM_infectedCells_Rev.tsv

# #Fwd
# computeMatrix scale-regions -S $(ls ${outdir}/*_fwd.bw | tr "\n" " ") -R ${ref} -a 3000 -b 3000 --regionBodyLength 3000 --outFileName scale_regions_infectedCells_TPM_Fwd.mat.gz --outFileNameMatrix scale_regions_infectedCells_TPM_Fwd.tabular --samplesLabel $(ls ${outdir}/*_fwd.bw | sed -E 's/(_m6avsInput1_TPM|.bw)//g' | tr "\n" " ") --verbose -p 8 --skipZeros

# plotProfile --matrixFile scale_regions_infectedCells_TPM_Fwd.mat.gz --perGroup --plotFileFormat pdf --outFileName ${outdir}/profile_scaleRegions_TPM_infectedCells_Fwd.pdf --plotTitle "" --regionsLabel "" --outFileNameData scale_regions_infectedCells_TPM_Fwd.tsv

echo "Copying result files to ${outdir}:"

cp *.gz ${outdir}
cp *.tabular ${outdir}


echo "Done."
