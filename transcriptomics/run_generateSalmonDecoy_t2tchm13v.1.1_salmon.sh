#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J salmonDecoy


ml Salmon BEDTools

ref=$(realpath $1)
gtf=$(realpath $2)

echo "Starting generateDecoy"

/crex/proj/nb_storage/private/terra_project/scripts/generateDecoyTranscriptome.sh -m /domus/h1/amendez/private/MashMap/mashmap -a ${gtf} -g ${ref} -t ../../ref/CHM13/chm13.draft_v1.1.transcriptome.gffread.fasta -o .

echo "Done."
