#!/bin/bash -l
#SBATCH -A naiss2023-22-153
#SBATCH -p core -n 8
#SBATCH -t 04:00:00
#SBATCH -J bam_stats

module load bioinfo-tools
module load sambamba/0.7.1
module load samtools/1.12

inputdir=$(realpath $1) # Path to BAM files to calculate statistics from
outdir=$(realpath $2) # Output directory
project=$3 # Project name for labeling the output statistics file

# Script for calculating BAM statistics  using Sambamba for reads mapped to the human reference genome and E.coli genome (used as spiked-in control).

cd $SNIC_TMP

cp ${inputdir}/*.markdup.bam .

echo -e "Sample\tTotal_reads\tMapped\tUnmapped\tDuplicate\tmapped_uniq_withdup_human\tmapped_uniq_withdup_ecoli\tnodup_uniq_human\tnodup_uniq_human_chr\tnodup_uniq_ecoli\tnodup_uniq_human_fwd\tnodup_uniq_human_rev\tmapped_uniq_withdup_human_fwd\tmapped_uniq_withdup_human_rev\tmapped_uniq_withdup_ecoli_fwd\tmapped_uniq_withdup_ecoli_rev" > bamstats_sambamba_${project}.tsv

for bam in $(ls *.markdup.bam)
do 
	echo "Getting stats for aln file: ${bam}"
	sample=${bam%.markdup.bam}
	total_reads=$(sambamba view -q -c ${bam} --nthreads ${SLURM_NTASKS})
	mapped=$(sambamba view -q -c -F "not unmapped" ${bam} --nthreads ${SLURM_NTASKS})
	unmapped=$(sambamba view -q -c -F "unmapped" ${bam} --nthreads ${SLURM_NTASKS})
	duplicate=$(sambamba view -q -c -F "duplicate" ${bam} --nthreads ${SLURM_NTASKS})
    mapped_uniq_withdup_human=$(sambamba view -q -c -F "[NH]==1 and not (ref_name=~/^(U00096.2)/) and not unmapped" ${bam} --nthreads ${SLURM_NTASKS})
    mapped_uniq_withdup_ecoli=$(sambamba view -q -c -F "[NH]==1 and ref_name=='U00096.2' and not unmapped" ${bam} --nthreads ${SLURM_NTASKS})
	nodup_uniq_human=$(sambamba view -q -c -F "[NH]==1 and not (ref_name=~/^(U00096.2)/) and not unmapped and not duplicate" ${bam} --nthreads ${SLURM_NTASKS})
	nodup_uniq_human_chr=$(sambamba view -q -c -F "[NH]==1 and ref_name=~/^(chr)/ and not unmapped and not duplicate" ${bam} --nthreads ${SLURM_NTASKS})
	nodup_uniq_ecoli=$(sambamba view -q -c -F "[NH]==1 and ref_name=='U00096.2' and not unmapped and not duplicate" ${bam} --nthreads ${SLURM_NTASKS})
	nodup_uniq_human_fwd=$(sambamba view -q -c -F "ref_name=~/^chr/ and [NH]==1 and [XS]=='+' and not unmapped and not duplicate" ${bam} --nthreads ${SLURM_NTASKS})
	nodup_uniq_human_rev=$(sambamba view -q -c -F "ref_name=~/^chr/ and [NH]==1 and [XS]=='-' and not unmapped and not duplicate" ${bam} --nthreads ${SLURM_NTASKS})
    mapped_uniq_withdup_human_fwd=$(sambamba view -q -c -F "[NH]==1 and [XS]=='+' and not (ref_name=~/^(U00096.2)/) and not unmapped" ${bam} --nthreads ${SLURM_NTASKS})
    mapped_uniq_withdup_human_rev=$(sambamba view -q -c -F "[NH]==1 and [XS]=='-' and not (ref_name=~/^(U00096.2)/) and not unmapped" ${bam} --nthreads ${SLURM_NTASKS})
    mapped_uniq_withdup_ecoli_fwd=$(sambamba view -q -c -F "[NH]==1 and [XS]=='+' and ref_name=='U00096.2' and not unmapped" ${bam} --nthreads ${SLURM_NTASKS})
    mapped_uniq_withdup_ecoli_rev=$(sambamba view -q -c -F "[NH]==1 and [XS]=='-' and ref_name=='U00096.2' and not unmapped" ${bam} --nthreads ${SLURM_NTASKS})
echo -e "${sample}\t${total_reads}\t${mapped}\t${unmapped}\t${duplicate}\t${mapped_uniq_withdup_human}\t${mapped_uniq_withdup_ecoli}\t${nodup_uniq_human}\t${nodup_uniq_human_chr}\t${nodup_uniq_ecoli}\t${nodup_uniq_human_fwd}\t${nodup_uniq_human_rev}\t${mapped_uniq_withdup_human_fwd}\t${mapped_uniq_withdup_human_rev}\t${mapped_uniq_withdup_ecoli_fwd}\t${mapped_uniq_withdup_ecoli_rev}" >> bamstats_sambamba_${project}.tsv
done

echo "Saving alignment stats to ${outdir}"
	cp bamstats_sambamba_${project}.tsv ${outdir}

echo "Done."