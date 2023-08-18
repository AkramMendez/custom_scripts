#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 1-00:00:00
#SBATCH -J compMat


module load bioinfo-tools
module load deepTools/3.3.2

ref="/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
inputdir="/aln/aln_vero_batch1/bigwig_CPMcoverage_infectedVero/all/bigwig_CPM_inputs"
inputdir2="/aln/aln_vero_batch1/bigwig_CPMcoverage_infectedVero/all/bigwig_spikeInNorm_IP"
outdir="/profilePlots"

cd $SNIC_TMP

cp ${inputdir}/*.bw .
cp ${inputdir2}/*.bw .

echo "ComputeMatrix calculation"
computeMatrix scale-regions -S $(ls *.bw | tr "\n" " ") -R ${ref} \
-a 3000 \
-b 3000 \
--regionBodyLength 3000 \
--outFileName scale_regions_infNonInf_vero_All.mat.gz \
--outFileNameMatrix scale_regions_infNonInf_vero_All.tabular \
--samplesLabel $(ls *.bw | sed -E 's/(_bs10|.bw)//g' | tr "\n" " ") \
--verbose \
-p ${SLURM_NTASKS} \
--skipZeros

echo "Computing plotProfile"
plotProfile --matrixFile scale_regions_infNonInf_vero_All.mat.gz \
--perGroup \
--plotFileFormat pdf \
--outFileName ${outdir}/profile_scaleRegions_infNonInf_vero_infectedCells_all.pdf \
--plotTitle "" \
--regionsLabel "" \
--outFileNameData ${outdir}/scale_regions_infNonInf_vero_infectedCells_All.tsv

echo "Copying result files to ${outdir}:"

cp *.gz ${outdir}
cp *.tabular ${outdir}

echo "Done."
