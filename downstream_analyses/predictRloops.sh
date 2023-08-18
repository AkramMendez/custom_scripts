#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 04:00:00
#SBATCH -J RloopFind

#R-loop finder
#QmRLFS-finder
#http://rloop.bii.a-star.edu.sg/?pg=qmrlfs
#python QmRLFS-finder.py --model m1 -bed -i snrpn.fasta -o snrpn
fasta=$(realpath $1) # Reference genome FASTA file

conda activate /proj/nb_project/private/conda_envs/py2.7

# Script for predicting genome-wide R-loops in a reference genome using RloopFinder.

RloopFinder="/domus/h1/amendez/private/Rloop_finder.py"

python2.7 ${RloopFinder} --model m1 -bed -i ${fasta} -o ${name}_RloppFinder_QmRLFS.bed