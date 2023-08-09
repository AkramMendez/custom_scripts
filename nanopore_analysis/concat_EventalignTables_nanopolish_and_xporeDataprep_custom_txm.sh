#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p node -n 32
#SBATCH -t 56:00:00
#SBATCH -J xPorePrep

# xPore data preparation for m6A detection after event align of dRNA Oxford Nanopore data:
# This script cocantenates the huge tables (>300G) resulting from Nanopolish eventalign analysis
# The tables can be fetched from aws s3 bucket (if previously stored after nanopolish eventalign analysis ) or copied locally.
# After the eventalign files are ready, they are processed using the xPore data preparation pipeline before m6A detection.
# The resulting prepared files are deposited in the output directory, whereas the whole eventalign table is then stored into the S3 bucket.

module load bioinfo-tools
module load awscli

name=$1
inputdir=$(realpath $2)
outdir=$(realpath $3)

gtf="CHM13_hg38chrY.gtf"

#### T2T CHM13 Transcriptome
transcript_fasta="chm13.draft_v1.1_transcriptome.fa"

cd ${SNIC_TMP}

echo "Copying eventalign tables"

aws s3 cp s3://momo23/nanopolish/nanopolish_customTxm_001_eventalign.txt .

cp ${inputdir}/nanopolish_customTxm_002_eventalign.txt .

#rsync -avP ${inputdir}/*_eventalign.txt .

echo "Concatenating tables"

cat *_eventalign.txt > ${name}_eventalign_concat.txt

echo "Size of concatenated table"

du -sh ${name}_eventalign_concat.txt

eventalign="${name}_eventalign_concat.txt"

echo "Activating xPore conda environment"

conda activate /proj/nb_project/private/conda_envs/nanopolish/

xpore dataprep --eventalign ${eventalign} \
--out_dir=${outdir} \
--gtf_or_gff=${gtf} \
--transcript_fasta=${transcript_fasta} \
--chunk_size 10000 \
--genome \
--n_processes=${SLURM_NTASKS}

echo "xPore dataprep done."

conda deactivate

echo "Copying file to s3 bucket"

aws s3 cp ${name}_eventalign_concat.txt s3://momo23/nanopolish/

echo "Done."

