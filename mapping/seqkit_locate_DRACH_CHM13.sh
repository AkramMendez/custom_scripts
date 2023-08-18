#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J seqkit_locate

module load bioinfo-tools
module load SeqKit/0.15.0


echo "Searching DRACH motif:"

seqkit locate --degenerate -p DRACH --bed CHM13_T2T.fa > seqkit_locate_DRACHmotif_CHM13_hg38chrY_ecoliK12.bed

echo "Done."
