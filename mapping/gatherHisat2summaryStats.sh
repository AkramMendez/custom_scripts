#!/bin/bash -l

#Gather summary files within the same folder.
#Based on: https://gist.github.com/slavailn/cafd7c50276a37b3aad61756ba994353

inputdir=$(realpath $1)

echo -e "sampleID\ttotal_reads\tunmapped\taligned_one_time\taligned_multiple\talignment_rate"
for file in ${inputdir}/*summary*.txt
do
    filename=`basename $file`
    samplename=${filename%.summary}
    
    total_reads=""
    unmapped=""
    aligned_one_time=""
    aligned_multiple=""
    alignment_rate=""
   
    while IFS= read -r line
    do
        if [[ $line == *reads* ]];
        then
            total_reads=`echo $line | awk '{print $1}'`;
        fi
        
        if [[ $line == *aligned[[:space:]]*0[[:space:]]*times* ]];
        then 
            unmapped_reads=`echo $line | awk '{print $1}'`;
        fi

        if [[ $line == *aligned[[:space:]]*exactly[[:space:]]*1[[:space:]]*time* ]]
        then
           aligned_one_time=`echo $line | awk '{print $1}'`;
        fi

        if [[ $line == *aligned[[:space:]]*\>1[[:space:]]*times*  ]]
        then
           aligned_multiple=`echo $line | awk '{print $1}'`;
        fi

        if [[ $line == *overall[[:space:]]*alignment[[:space:]]*rate*  ]]
        then
           alignment_rate=`echo $line | awk '{print $1}'`;
           alignment_rate=`echo $alignment_rate | sed s/%//`;
        fi
      
    done < "$file"
    echo -e "$samplename\t$total_reads\t$unmapped_reads\t$aligned_one_time\t$aligned_multiple\t$alignment_rate"
done