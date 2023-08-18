#!/bin/bash -l

hg38="hg38.ucsc.chromosomes.fa"

# Script for concatenating the human reference genome (UCSC hg38) with the humanSTF500 subtelomer assembly.

#awk 'BEGIN{FS=OFS="\t"}{print $1,0,500000}' chromosome.labels.humanSTF500.igv.shortlabels.txt | grep -P "chr(\d+|X|Y)(q)" > chromosome.labels.humanSTF500.igv.shortlabels.chrQ.bed
#awk 'BEGIN{FS=OFS="\t"}{print $1,0,500000}' chromosome.labels.humanSTF500.igv.shortlabels.txt | grep -P "chr(\d+|X|Y)(p)" > chromosome.labels.humanSTF500.igv.shortlabels.chrP.bed
#seqkit replace humanSTF500.igv.shortlabels.fa -p "chrXpYp" -r "chrXp" > humanSTF500.igv.shortlabels.relabel.fa
#seqkit subseq humanSTF500.igv.shortlabels.relabel.fa --bed chromosome.labels.humanSTF500.igv.shortlabels.chrP.bed   > humanSTF500.igv.shortlabels.chrP.seqkit.subseq.fa
#seqkit subseq humanSTF500.igv.shortlabels.relabel.fa --bed chromosome.labels.humanSTF500.igv.shortlabels.chrQ.bed | seqkit seq --seq-type DNA --reverse --complement > humanSTF500.igv.shortlabels.chrQ.seqkit.subseq.revcomp.fa

#cat humanSTF500.igv.shortlabels.chrP.seqkit.subseq.fa | seqkit replace -p "p_.*" -r "" > humanSTF500.igv.shortlabels.chrP.seqkit.subseq.relabel.fa
#cat humanSTF500.igv.shortlabels.chrQ.seqkit.subseq.revcomp.fa | seqkit replace -p "q_.*" -r "" > humanSTF500.igv.shortlabels.chrQ.seqkit.subseq.revcomp.relabel.fa


#bedtools complement -i subtel.coordinates.formasking.hg38.bed -g chromsizes.ucsc.hg38.sorted.tsv > complement.bedtools.coordinates.notsubtel.hg38.bed

echo "sorting bed"
bedtools sort -i starting.chrP.arms.hg38.before.subtel.bed > starting.chrP.arms.hg38.before.subtel.sorted.bed
bedtools sort -i middle.intrachrom.hg38.before.chrQ.arms.bed > middle.intrachrom.hg38.before.chrQ.arms.sorted.bed
bedtools sort -i ending.chrQ.arms.hg38.after.subtel.bed > ending.chrQ.arms.hg38.after.subtel.sorted.bed

cat humanSTF500.igv.shortlabels.chrQ.seqkit.subseq.revcomp.relabel.fa.fai | awk 'BEGIN{FS=OFS="\t"}{if($1=="chr13"){print $1,0,466457}else{print $1,0,$2}}' | bedtools sort -i - > subtel.coordinates.forextractingfasta.chrQarms.sorted.bed
cat humanSTF500.igv.shortlabels.chrP.seqkit.subseq.relabel.fa.fai | awk 'BEGIN{FS=OFS="\t"}{if($1=="chr7"){print $1,29405,$2}else{print $1,0,$2}}' | bedtools sort -i - > subtel.coordinates.forextractingfasta.chrP.arms.sorted.bed

echo "extracting chr P and Q arms subtelomeric sequences:"
seqkit subseq --bed subtel.coordinates.forextractingfasta.chrQarms.sorted.bed humanSTF500.igv.shortlabels.chrQ.seqkit.subseq.revcomp.relabel.fa | seqkit rmdup | seqkit replace -p "_.*:." -r "" - > humanSTF500.igv.shortlabels.chrQ.seqkit.subseq.revcomp.relabel.limitscheked.fa
seqkit subseq --bed subtel.coordinates.forextractingfasta.chrP.arms.sorted.bed humanSTF500.igv.shortlabels.chrP.seqkit.subseq.relabel.fa | seqkit rmdup | seqkit replace -p "_.*:." -r "" - > humanSTF500.igv.shortlabels.chrP.seqkit.subseq.relabel.limitschecked.fa

cat humanSTF500.igv.shortlabels.chrP.seqkit.subseq.relabel.limitschecked.fa <(echo -e ">chr13\n>chr14\n>chr15\n>chr21\n>chr22\n>chrY") > humanSTF500.igv.shortlabels.chrP.seqkit.subseq.relabel.limitschecked.addedlabels.fa

echo "extracting starting sequences from hg38"
seqkit subseq --bed starting.chrP.arms.hg38.before.subtel.sorted.bed ${hg38} | seqkit rmdup | seqkit replace -p "_1-\d+:." -r "" > subseq.seqkit.hg38.starting.chrP.before.subtel.relabel.fa
cat subseq.seqkit.hg38.starting.chrP.before.subtel.relabel.fa <(echo -e ">chr7") > subseq.seqkit.hg38.starting.chrP.before.subtel.relabel.addedlabel.fa

echo "extracting middle sequences from hg38"
seqkit subseq --bed middle.intrachrom.hg38.before.chrQ.arms.sorted.bed ${hg38} | seqkit rmdup | seqkit replace -p "_.*:." -r "" > subseq.seqkit.hg38.intrachrom.between.subtels.relabel.fa

echo "extracting ending sequences from hg38"
seqkit subseq --bed ending.chrQ.arms.hg38.after.subtel.sorted.bed ${hg38} | seqkit rmdup | seqkit replace -p "_.*:." -r "" > subseq.seqkit.hg38.endings.chrQ.after.subtel.relabel.fa
cat subseq.seqkit.hg38.endings.chrQ.after.subtel.relabel.fa <(echo -e ">chr13") > subseq.seqkit.hg38.endings.chrQ.after.subtel.relabel.addedlabel.fa



echo "concatenating starting seqs and subtel P arms"
seqkit concat subseq.seqkit.hg38.starting.chrP.before.subtel.relabel.addedlabel.fa humanSTF500.igv.shortlabels.chrP.seqkit.subseq.relabel.limitschecked.addedlabels.fa > concat.start.hg38.subtel.chrP.fa

echo "concatenating starting-subtelP arms with middle sequences"
seqkit concat concat.start.hg38.subtel.chrP.fa subseq.seqkit.hg38.intrachrom.between.subtels.relabel.fa > concat.starting.hg38.chrP.subtel.intrachrom.hg38.fa

echo "concatenating starting-subtelP arms-middle with chrQ subtel arms"
seqkit concat concat.starting.hg38.chrP.subtel.intrachrom.hg38.fa humanSTF500.igv.shortlabels.chrQ.seqkit.subseq.revcomp.relabel.limitscheked.fa > concat.starting.hg38.chrP.subtel.intrachrom.hg38.chrQ.subtel.fa

echo "concatenating starting-subtelP arms-middle-subtelQ arms with ending sequences"
seqkit concat concat.starting.hg38.chrP.subtel.intrachrom.hg38.chrQ.subtel.fa subseq.seqkit.hg38.endings.chrQ.after.subtel.relabel.addedlabel.fa | seqkit sort -N - > concat.hg38.subtel.hybgenome.v4.sorted.fa

echo "merged genome done."

# Calculate chromosome lengths 
 seqkit fx2tab --length --name --header-line subseq.seqkit.hg38.starting.chrP.before.subtel.relabel.addedlabel.fa | sort -k1,1 -k2,2n

