#!/bin/bash -l
#SBATCH -A naiss2023-22-153
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH -J filterBAM

module load bioinfo-tools
module load sambamba/0.7.1
module load samtools/1.12

inputdir=$(realpath $1)
outdir=$(realpath $2)

cd $SNIC_TMP

cp ${inputdir}/*.markdup.bam* .

splitBAM(){
    bam_file=$1
echo "Splitting aln file: $1"

	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "[XS]=='+'" ${bam_file} | samtools sort -o ${bam_file%.bam}.fwd.bam -@ ${SLURM_NTASKS}
	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "[XS]=='-'" ${bam_file} | samtools sort -o ${bam_file%.bam}.rev.bam -@ ${SLURM_NTASKS}

	echo "Indexing filtered alignment"
	samtools index ${bam_file%.bam}.fwd.bam -@ ${SLURM_NTASKS}
	samtools index ${bam_file%.bam}.rev.bam -@ ${SLURM_NTASKS}
}

echo "Filtering markdup BAM file to keep human chromosomes (only chr, no patches), plus E.coli genomes"

for bam in $(ls *.markdup.bam)
do 
    echo "Filtering file: ${bam}"
	sambamba view -q -h -f bam -F "[NH]==1 and ref_name=~/^(chr|U00096.2)/ and not unmapped" ${bam}  | samtools sort -o ${bam%.markdup.bam}.chrms.mapped.uniq.bam -@ ${SLURM_NTASKS}

    echo "Indexing BAM file"
    samtools index ${bam%.markdup.bam}.chrms.mapped.uniq.bam -@ ${SLURM_NTASKS}

    echo "Copying BAM and index to ${outdir}"
    cp ${bam%.markdup.bam}.chrms.mapped.uniq.bam ${outdir}
    cp ${bam%.markdup.bam}.chrms.mapped.uniq.bam.bai ${outdir}

    echo "Splitting BAM file"

    splitBAM ${bam%.markdup.bam}.chrms.mapped.uniq.bam
done

echo "Copying BAM files to ${outdir}"
cp *.fwd.bam* ${outdir}
cp *.rev.bam* ${outdir}

echo "Done."


