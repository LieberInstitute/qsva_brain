#!/bin/bash 
#$ -cwd
# Specify log file names
#$ -o logs/quantify_top1000.txt
#$ -e logs/quantify_top1000.txt
#$ -pe local 8
#$ -m e
#$ -l bluejay,mem_free=15G,h_vmem=15G,h_fsize=100G
#$ -N quantify_top1000

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

date
module load wiggletools/default
module load ucsctools
Rscript quantify_top1000.R
date
