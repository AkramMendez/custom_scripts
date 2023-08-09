#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 32
#SBATCH -t 56:00:00
#SBATCH -J xPorePrep

# This script cocantenates the huge tables (>300G) resulting from Nanopolish eventalign analysis
# The tables can be fetched from aws s3 bucket or copied locally.

#inputdir=$(realpath $1)
#outdir=$(realpath $2)
module load bioinfo-tools
module load awscli

name=$1
outdir=$(realpath $2)

gtf="/crex/proj/nb_storage/private/terra_project/ref/CHM13/CHM13_hg38chrY_AND_chromEndCoordinates_50kb_notUCSCformat_custom.gtf"

#### Transcriptome including sequences from 50kb chromosome ends (Reference transcriptome + whole 50kb chr. end sequences concatenated) to simulate a big transcript for each chromosome end since those regio$
transcript_fasta="/proj/nb_storage/private/terra_project/ref/chm13.draft_v1.1_transcriptome_gffread_concatenatedTo_chm13_chromEndCoord50kb_sequences.fa"



cd ${SNIC_TMP}

echo "Copying eventalign tables"

aws s3 cp s3://momo23/nanopolish/mettl3kd_customTxm_withSecondary_001_eventalign.txt .

cp /crex/proj/nb_storage/private/terra_ont/nanopolish/eventalign_summary_mettl3kd_customTxm_withSecondary_002.txt .

#rsync -avP ${inputdir}/*_eventalign.txt .

echo "Concatenating tables"

cat *_eventalign.txt > ${name}_eventalign_concat.txt

eventalign="${name}_eventalign_concat.txt"

echo "Activating xPore conda environment"

conda activate /proj/nb_project/private/conda_envs/nanopolish/

xpore dataprep --eventalign ${eventalign} \
--out_dir=${outdir} \
--gtf_or_gff=/crex/proj/nb_storage/private/terra_project/ref/CHM13/chm13.draft_v1.1.gene_annotation.v4.gffread.noUCSCformat.gtf \
--transcript_fasta=/crex/proj/nb_storage/private/terra_project/ref/CHM13/chm13.draft_v1.1.transcriptome.gffread.wnocds.fasta \
--chunk_size 10000 \
--genome \
--n_processes=${SLURM_NTASKS}

echo "Done."

echo "Copying file to ${outdir}"

cp ${name}_eventalign_concat.txt ${outdir}

echo "Done."

