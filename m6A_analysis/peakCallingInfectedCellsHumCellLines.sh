#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 16:00:00
#SBATCH -J macs_peakcalling

module load bioinfo-tools
module load MACS/2.2.6

#inputdir=$(realpath $1) 
#outdir=$(realpath $2)
inputdir="/crex/proj/nb_storage/private/BEA21P074_Roshan_Vaid/aln/batch202203/BEA2216Pa/normSpikeIn"
outdir="/crex/proj/nb_storage/private/BEA21P074_Roshan_Vaid/peakCalling/peaks_LCsamples"

cd ${SNIC_TMP}

#------ Peak calling Sars-cov2 infected LC samples ---#
#vero
#gsize=2756670577
#human
gsize=2.7e9
pvalue=0.01
#extsize=75
#qvalue=0.01

#Names Fwd output file sufix:
name_LCcov2d4_fwd="peaks_LC_cov2_d4_fwd"
name_LCcov2d7_fwd="peaks_LC_cov2_d7_fwd"
name_LCmockd4_fwd="peaks_LC_mock_d4_fwd"
name_LCmockd7_fwd="peaks_LC_mock_d7_fwd"
name_LC27cov2_fwd="peaks_LCN27_cov2_fwd"
name_LC27mock_fwd="peaks_LCN27_mock_fwd"

#Names Rev output file sufix:
name_LCcov2d4_rev="peaks_LC_cov2_d4_rev"
name_LCcov2d7_rev="peaks_LC_cov2_d7_rev"
name_LCmockd4_rev="peaks_LC_mock_d4_rev"
name_LCmockd7_rev="peaks_LC_mock_d7_rev"
name_LC27cov2_rev="peaks_LCN27_cov2_rev"
name_LC27mock_rev="peaks_LCN27_mock_rev"

echo "Making output directories"

# mkdir -p ${outdir}/${name_LCcov2d4_fwd}
# mkdir -p ${outdir}/${name_LCcov2d7_fwd}
# mkdir -p ${outdir}/${name_LCmockd4_fwd}
# mkdir -p ${outdir}/${name_LCmockd7_fwd}
# mkdir -p ${outdir}/${name_LC27cov2_fwd}
# mkdir -p ${outdir}/${name_LC27mock_fwd}

# mkdir -p ${outdir}/${name_LCcov2d4_rev}
# mkdir -p ${outdir}/${name_LCcov2d7_rev}
# mkdir -p ${outdir}/${name_LCmockd4_rev}
# mkdir -p ${outdir}/${name_LCmockd7_rev}
# mkdir -p ${outdir}/${name_LC27cov2_rev}
# mkdir -p ${outdir}/${name_LC27mock_rev}

#Input files Fwd:
LC_cov2_d4_i2_fwd=${inputdir}/fwd/LC_CoV2_D4_Input_2_S8_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LC_cov2_d4_m6A_fwd=${inputdir}/fwd/LC_CoV2_D4_m6A_RIP_S9_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LC_cov2_d7_i1_fwd=${inputdir}/fwd/LC_CoV2_D7_Input_1_S10_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LC_cov2_d7_m6a_fwd=${inputdir}/fwd/LC_CoV2_D7_m6A_RIP_S12_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LC_mock_d4_i1_fwd=${inputdir}/fwd/LC_Mock_D4_Input_1_S1_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LC_mock_d4_m6a_fwd=${inputdir}/fwd/LC_Mock_D4_m6A_RIP_S3_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LC_mock_d7_i1_fwd=${inputdir}/fwd/LC_Mock_D7_Input_1_S4_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LC_mock_d7_m6a_fwd=${inputdir}/fwd/LC_Mock_D7_m6A_RIP_S6_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LCN27_cov2_i1_fwd=${inputdir}/fwd/LC_N27_CoV2_Input_S15_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LCN27_cov2_m6a_fwd=${inputdir}/fwd/LC_N27_CoV2_m6A_S16_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LCN27_mock_i1_fwd=${inputdir}/fwd/LC_N27_Mock_Input_S13_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam
LCN27_mock_m6a_fwd=${inputdir}/fwd/LC_N27_Mock_m6A_S14_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.fwd.bam

#Input files Rev:
LC_cov2_d4_i2_rev=${inputdir}/rev/LC_CoV2_D4_Input_2_S8_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LC_cov2_d4_m6A_rev=${inputdir}/rev/LC_CoV2_D4_m6A_RIP_S9_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LC_cov2_d7_i1_rev=${inputdir}/rev/LC_CoV2_D7_Input_1_S10_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LC_cov2_d7_m6a_rev=${inputdir}/rev/LC_CoV2_D7_m6A_RIP_S12_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LC_mock_d4_i1_rev=${inputdir}/rev/LC_Mock_D4_Input_1_S1_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LC_mock_d4_m6a_rev=${inputdir}/rev/LC_Mock_D4_m6A_RIP_S3_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LC_mock_d7_i1_rev=${inputdir}/rev/LC_Mock_D7_Input_1_S4_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LC_mock_d7_m6a_rev=${inputdir}/rev/LC_Mock_D7_m6A_RIP_S6_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LCN27_cov2_i1_rev=${inputdir}/rev/LC_N27_CoV2_Input_S15_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LCN27_cov2_m6a_rev=${inputdir}/rev/LC_N27_CoV2_m6A_S16_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LCN27_mock_i1_rev=${inputdir}/rev/LC_N27_Mock_Input_S13_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam
LCN27_mock_m6a_rev=${inputdir}/rev/LC_N27_Mock_m6A_S14_R1_001.ht2.srt.chrms.nodup.uniq.spikeInNorm.rev.bam

callPeaks(){
    IP=$1
    input=$2
    name=$3
    gsize=$4
    pvalue=$5
    outdir=$6
    echo "Calling peaks ${IP} vs ${input}"

    macs2 callpeak -t ${IP} -c ${input} \
    --name ${name} --format="BAM" --gsize=${gsize} \
    --nomodel --bdg -p ${pvalue} --keep-dup auto \
    --call-summits \
    --outdir ${outdir}/${name} 2> macs2.${name}.out

}

echo "Calling peaks for Fwd samples"

callPeaks ${LC_cov2_d4_m6A_fwd} ${LC_cov2_d4_i2_fwd} ${name_LCcov2d4_fwd} ${gsize} ${pvalue} ${outdir}

callPeaks ${LC_cov2_d7_m6a_fwd} ${LC_cov2_d7_i1_fwd} ${name_LCcov2d7_fwd} ${gsize} ${pvalue} ${outdir}

callPeaks ${LC_mock_d4_m6a_fwd} ${LC_mock_d4_i1_fwd} ${name_LCmockd4_fwd} ${gsize} ${pvalue} ${outdir}

callPeaks ${LC_mock_d7_m6a_fwd} ${LC_mock_d7_i1_fwd} ${name_LCmockd7_fwd} ${gsize} ${pvalue} ${outdir}

callPeaks ${LCN27_cov2_m6a_fwd} ${LCN27_cov2_i1_fwd} ${name_LC27cov2_fwd} ${gsize} ${pvalue} ${outdir}

callPeaks ${LCN27_mock_m6a_fwd} ${LCN27_mock_i1_fwd} ${name_LC27mock_fwd} ${gsize} ${pvalue} ${outdir}


echo "Calling peaks for Rev samples"

callPeaks ${LC_cov2_d4_m6A_rev} ${LC_cov2_d4_i2_rev} ${name_LCcov2d4_rev} ${gsize} ${pvalue} ${outdir}

callPeaks ${LC_cov2_d7_m6a_rev} ${LC_cov2_d7_i1_rev} ${name_LCcov2d7_rev} ${gsize} ${pvalue} ${outdir}

callPeaks ${LC_mock_d4_m6a_rev} ${LC_mock_d4_i1_rev} ${name_LCmockd4_rev} ${gsize} ${pvalue} ${outdir}

callPeaks ${LC_mock_d7_m6a_rev} ${LC_mock_d7_i1_rev} ${name_LCmockd7_rev} ${gsize} ${pvalue} ${outdir}

callPeaks ${LCN27_cov2_m6a_rev} ${LCN27_cov2_i1_rev} ${name_LC27cov2_rev} ${gsize} ${pvalue} ${outdir}

callPeaks ${LCN27_mock_m6a_rev} ${LCN27_mock_i1_rev} ${name_LC27mock_rev} ${gsize} ${pvalue} ${outdir}

echo "Copying log files to ${outdir}"
cp *.out ${outdir}

echo "Done."