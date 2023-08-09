#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J downsampPicard

module load bioinfo-tools
module load picard/2.23.4
module load samtools/1.12

inputdir="/crex/proj/nb_storage/private/batchMar2022/aln/tumor/nodup"
outdir="/crex/proj/nb_storage/private/batchMar2022/aln/tumor/normSpikeIn"

cd $SNIC_TMP

cp ${inputdir}/*.bam .
cp ${inputdir}/*.bai .

downsampleBAM(){
        bam=$1
	scalefactor=$2
        outdir=$3
        echo "Downsampling file: ${bam}"

        java -jar $PICARD_ROOT/picard.jar DownsampleSam -I ${bam} \
        -O ${outdir}/${bam%.bam}.spikeInNorm.bam \
        -STRATEGY HighAccuracy \
        -P ${scalefactor}
}


#BAM files
SKNF1_Csh_input="SKNF1_Csh_Input_S27_R1_001_trimmed.fq.gz.hisat2.sorted.chrms.nodup.bam"
SKNF1_Csh_m6a="SKNF1_Csh_m6A_S29_R1_001_trimmed.fq.gz.hisat2.sorted.chrms.nodup.bam"
SKNF1_sh1_input="SKNF1_Sh1_1_Input_S28_R1_001_trimmed.fq.gz.hisat2.sorted.chrms.nodup.bam"
SKNF1_sh1_m6a="SKNF1_Sh1_1_m6A_S30_R1_001_trimmed.fq.gz.hisat2.sorted.chrms.nodup.bam"

#Scaling factors:

sf_SKNF1_Csh_input=1
sf_SKNF1_Csh_m6a=1
sf_SKNF1_sh1_input=0.9606281796
sf_SKNF1_sh1_m6a=0.6752390373

#Downsample BAM files for batch1:

downsampleBAM ${SKNF1_Csh_input} ${sf_SKNF1_Csh_input} ${outdir}
downsampleBAM ${SKNF1_Csh_m6a} ${sf_SKNF1_Csh_m6a} ${outdir}
downsampleBAM ${SKNF1_sh1_input} ${sf_SKNF1_sh1_input} ${outdir}
downsampleBAM ${SKNF1_sh1_m6a} ${sf_SKNF1_sh1_m6a} ${outdir}

echo "Done."