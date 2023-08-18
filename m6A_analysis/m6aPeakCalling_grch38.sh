#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J macs_peakcalling


module load bioinfo-tools
module load bowtie/1.2.3
module load MACS/2.2.6
module load MEMEsuite/5.1.1
module load BEDTools/2.29.2
module load samtools




gtf="gencode.v36.chr_patch_hapl_scaff.basic.annotation.gtf"
hg38="GRCh38.primary_assembly.genome.fa"
IP=$1
input=$2
n=$(echo "$(($SLURM_NTASKS-1))")


#NOTE: The bowtie index needs to be located in the same folder of the analysis pipeline.
#echo "Starting Bowtie mapping for IP"

#bowtie -m 5 -a --sam --best --strata GRCh38_noalt_as/GRCh38_noalt_as ${IP} --threads 8 --chunkmbs 1024 | samtools view -bS -@ 8 > IP.bam

#echo "Starting Bowtie mapping for input.bample"

#bowtie -m 5 -a --sam --best --strata GRCh38_noalt_as/GRCh38_noalt_as ${input} --threads 8 --chunkmbs 1024 | samtools view -bS -@ 8 > input.bam

#echo "Finished Bowtie mapping"

echo "Starting MACS peak calling"

echo "Calling peaks (IP vs input)"

macs2 callpeak -t ${IP} -c ${input} --name=m6a --format="BAM" --gsize=2747877777 --tsize=75 --nomodel --extsize=100 --scale-to small --bdg 2> macs.out

echo "Calling negative peaks (flipping input vs IP)"

macs2 callpeak -t ${input} -c ${IP} --name=m6a_negativePeaks --format="BAM" --gsize=2747877777 --tsize=75 --nomodel --extsize=100 --scale-to small --bdg 2> macs_negativePeaks.out

echo "MACS2 peak calling finished"

echo "Sort peaks and select those with a FRD < 5%"

awk '{if($9 <= 5) print }' m6a_peaks.xls > m6a_sig_peaks.xls

awk '{if($9 <= 5) print }' m6a_peaks.xls > m6a_negativePeaks_sig_peaks.xls

echo "Sort peaks with FDR < 5% by fold change and retrieve peak-summit coordinates (50bp flanking)"

# Note: In macs14, the 5th column contained the summit position relative to the start position of the peak, 
# this changed in MACS2 where it now contains the absolute summit position.

sort -k8,8 -n -r m6a_sig_peaks.xls | head -1000 | awk '{summit=$5; print $1 "\t" summit-51 "\t" summit+50}' > bestPeaks.location
sort -k8,8 -n -r m6a_negativePeaks_sig_peaks.xls | head -1000 | awk '{summit=$5; print $1 "\t" summit-51 "\t" summit+50}' > bestPeaks_negativePeaks.location

echo "Map summits to gene annotation and retrieve the sequences from the sense strand"

intersectBed -wo -a bestPeaks.location -b ${gtf} | awk -v OFS="\t" '{print $1, $2,$3,"*","*",$10}' | uniq > bestPeaks.bed
intersectBed -wo -a bestPeaks_negativePeaks.location -b ${gtf} | awk -v OFS="\t" '{print $1, $2,$3,"*","*",$10}' | uniq > bestPeaks_negativePeaks.bed

fastaFromBed -s -fi ${hg38} -bed bestPeaks.bed -fo bestPeaks.fa
fastaFromBed -s -fi ${hg38} -bed bestPeaks_negativePeaks.bed -fo bestPeaks_negativePeaks.fa

echo "Removing duplicate sequences from fasta"

awk '/^>/{f=!d[$1];d[$1]=1}f' bestPeaks.fa > bestPeaks_uniq.fa
awk '/^>/{f=!d[$1];d[$1]=1}f' bestPeaks_negativePeaks.fa > bestPeaks_negativePeaks_uniq.fa

echo "De novo motif finding with MEME"

meme bestPeaks_uniq.fa -dna -nmotifs 25 -maxsize 1000000 -maxw 7 -o bestPeaks.top25.maxw7_meme -p ${n}

meme bestPeaks_negativePeaks_uniq.fa -dna -nmotifs 25 -maxsize 1000000 -maxw 7 -o bestPeaks.top25.maxw7_negativePeaks_meme -p ${n}

echo "Calculate summit-to-motif distance with CentriMo"

echo "Generate location file with peak-summit sequences (150bp) for CentriMo"

awk '{if($1~/[^#]/){summit=$5; print $1 "\t" summit-151 "\t" summit+150 }}' m6a_sig_peaks.xls > m6a_sig_peaks_summit.location

awk '{if($1~/[^#]/){summit=$5; print $1 "\t" summit-151 "\t" summit+150 }}' m6a_negativePeaks_sig_peaks.xls > m6a_sig_peaks_negativePeaks_summit.location

echo "Map summits to gene annotation"

intersectBed -wo -a m6a_sig_peaks_summit.location -b ${gtf} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | uniq > m6a_sig_peaks_summit.bed
intersectBed -wo -a m6a_sig_peaks_negativePeaks_summit.location -b ${gtf} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | uniq > m6a_sig_peaks_negativePeaks_summit.bed

echo "Retrieve sequences of peaks"
fastaFromBed -s -fi ${hg38} -bed m6a_sig_peaks_summit.bed -fo m6a_sig_peaks_summit.fa
fastaFromBed -s -fi ${hg38} -bed m6a_sig_peaks_negativePeaks_summit.bed -fo m6a_sig_peaks_negativePeaks_summit.fa


echo "Running CentriMo"
centrimo --motif 1 --o peaks_motif_centrimo --norc m6a_sig_peaks_summit.fa bestPeaks_meme/meme.txt --inc '*'
centrimo --motif 1 --o peaks_motif_centrimo_negativePeaks --norc m6a_sig_peaks_negativePeaks_summit.fa bestPeaks_negativePeaks_meme/meme.txt --inc '*'

echo "Done."
