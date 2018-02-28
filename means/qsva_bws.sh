#!/bin/bash 
#$ -cwd
# Specify log file names
#$ -o logs/qsva_bws.txt
#$ -e logs/qsva_bws.txt
#$ -m e
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

date
module load wiggletools/default
module load ucsctools
Rscript qsva_bws.R
date