#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=200G,h_fsize=200G
#$ -N merge_qsva_data
#$ -o ./logs/merge_qsva_data.txt
#$ -e ./logs/merge_qsva_data.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

#module load conda_R/3.4.x
Rscript merge_data.R

echo "**** Job ends ****"
date
