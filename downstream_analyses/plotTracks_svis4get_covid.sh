#!/bin/bash -l

gtf="reference_genomes/wuhCor1/wuhCor1.ncbiGenes.gtf"
ref="wuhCor1_ecoliK12_concat.fa"
outfile_fwd="svis4get_peaks_sarscov2_fwd.pdf"
outfile_rev="svis4get_peaks_sarscov2_rev.pdf"
configFile="svist4get_configFile_covid.txt"

wu_fwd="peaks_nodup_uniq_WU_infectedMappedToCov2_fwd_peaks_bedtoolsMerged.bdg"
uk_fwd="peaks_nodup_uniq_UK_infectedMappedToCov2_fwd_peaks_bedtoolsMerged.bdg"
sa_fwd="peaks_nodup_uniq_SA_infectedMappedToCov2_fwd_peaks_bedtoolsMerged.bdg"

wu_rev="peaks_nodup_uniq_WU_infectedMappedToCov2_rev_peaks_bedtoolsMerged.bdg"
uk_rev="peaks_nodup_uniq_UK_infectedMappedToCov2_rev_peaks_bedtoolsMerged.bdg"
sa_rev="peaks_nodup_uniq_SA_infectedMappedToCov2_rev_peaks_bedtoolsMerged.bdg"

peaksLi2021="peaks_CellReports_Li2021_star_bedtoolsMerged.bdg"
#Peaks from Li2021 are in Tables S1, sheet: "star_no_dup_peaks"
peaksZhang2021="peaks_mBio_Zhang2021_bedtoolsMerged.bdg"
peaksLiu2021_56h_fwd="peaks_vero56h_luiGenomeResearch_fwd.bdg"
peaksLiu2021_56h_rev="peaks_vero56h_luiGenomeResearch_rev.bdg"

wu_in_rev="WU_infected_input1_rev.bedgraph"
uk_in_rev="UK_infected_input1_rev.bedgraph"
sa_in_rev="SA_infected_input1_rev.bedgraph"

wu_m6a_rev="WU_infected_m6a_rev.bedgraph"
uk_m6a_rev="UK_infected_m6a_rev.bedgraph"
sa_m6a_rev="SA_infected_m6a_rev.bedgraph"
# getBedgraph(){
# 	narrowPeak=$1
# 	cat ${narrowPeak} | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,0.25}'
# }

# peaksLi2021_bdg=$(echo "$(getBedgraph ${peaksLi2021})")
# peaksLiu2021_bdg=$(echo "$(getBedgraph ${peaksLiu2021})")
# peaksZhang2021_bdg=$(echo "$(getBedgraph ${peaksZhang2021})")

# wu_fwd_bdg=$(getBedgraph ${wu_fwd})
# uk_fwd_bdg=$(getBedgraph ${uk_fwd})
# sa_fwd_bdg=$(getBedgraph ${sa_fwd})

# wu_rev_bdg=$(getBedgraph ${wu_rev})
# uk_rev_bdg=$(getBedgraph ${uk_rev})
# sa_rev_bdg=$(getBedgraph ${sa_rev})

# echo "${peaksLi2021_bdg}"

echo "Plotting tracks Positive Strand"
svist4get -c ${configFile} -gtf ${gtf} -fa ${ref} -bg ${peaksLi2021} ${peaksLiu2021_56h_fwd} ${peaksZhang2021} ${wu_fwd} ${uk_fwd} ${sa_fwd} -bl 'Li et al.,2021' 'Liu et al. 2021' 'Zhang et al., 2021' 'B.1' 'B.1.1.7' 'B.1.351' -blp center -bgb max -w NC_045512v2 1 30000 -it "Peaks Genomic Positive Strand" -o ${outfile_fwd}


echo "Plotting tracks Negative Strand"
svist4get -c ${configFile} -gtf ${gtf} -fa ${ref} -bg ${peaksLiu2021_56h_rev} ${wu_rev} ${uk_rev} ${sa_rev} -bl 'Liu et al. 2021' 'B.1' 'B.1.1.7' 'B.1.351' -blp center -bgb max -w NC_045512v2 1 30000 -it "Peaks Genomic Negative Strand" -hf 28000 30000 N -o ${outfile_rev}

echo "Plotting Zoom-in tracks Negative Strand"
svist4get -gtf ${gtf} -fa ${ref} -bg ${peaksLiu2021_56h_rev} ${wu_m6a_rev} ${wu_in_rev} -bl 'Liu et al. 2021' 'B.1 m6A' "B.1 input" -blp center -bgb none -w NC_045512v2 28000 30000 -it "Peaks Genomic Negative Strand - B.1 Zoom-in" -bul 300 300 -o ${outfile_rev%.pdf}_zoomIn_B.1.pdf

svist4get -gtf ${gtf} -fa ${ref} -bg ${peaksLiu2021_56h_rev} ${uk_m6a_rev} ${uk_in_rev} -bl 'Liu et al. 2021' 'B.1.1.7 m6A' 'B.1.1.7 input' -blp center -blp center -bgb none -w NC_045512v2 28000 30000 -it "Peaks Genomic Negative Strand - B.1.1.7 Zoom-in" -bul 300 300 -o ${outfile_rev%.pdf}_zoomIn_B.1.1.7.pdf

svist4get -gtf ${gtf} -fa ${ref} -bg ${peaksLiu2021_56h_rev} ${sa_m6a_rev} ${sa_in_rev} -bl 'Liu et al. 2021' 'B.1.351 m6A' 'B.1.351 input' -blp center -bgb none -w NC_045512v2 28000 30000 -it "Peaks Genomic Negative Strand - B.1.351 Zoom-in" -bul 300 300 -o ${outfile_rev%.pdf}_zoomIn_B.1.351.pdf

echo "Done."
#-hf 28000 30000 N
