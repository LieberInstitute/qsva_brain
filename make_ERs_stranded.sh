#!/bin/bash 
#$ -cwd
# Specify log file names
#$ -o logs/make_ERs_stranded.txt
#$ -e logs/make_ERs_stranded.txt
#$ -pe local 4
#$ -m e
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N make_ERs_stranded

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

date
module load wiggletools/default
module load ucsctools
Rscript make_ERs_stranded.R
date
