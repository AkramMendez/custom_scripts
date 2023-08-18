#!/usr/bin/bash


inputdir=$(realpath $1)
regions=$(realpath $2)
total_reads=$(realpath $3)


tot_csh_pooled_in=
tot_mettl3kd_pooled_in=
tot_mf1_in=
tot_mf2_in=
tot_sknf1_csh_in=
tot_sknf1_sh1_in=

echo "Converting BAM counts to CPM"

bedtools coverage -a ${bam_csh_pooled_in} -b ${regions} -s | awk '{ print $0"\t"($5/(${tot_csh_pooled_in}/1000000))}'