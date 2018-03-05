## Based on parts of /dcl01/lieber/ajaffe/lab/brainseq_phase2/get_degradation_regions.R

library('SummarizedExperiment')
library('rtracklayer')
library('recount.bwtool')
library('BiocParallel')
library('devtools')

## For file paths to BrainSeq phase 2 bigwigs
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata')

## Degradation regions
regs <- import('/dcl01/ajaffe/data/lab/qsva_brain/ERs/bed/DLPFC_Plus_HIPPO_RiboZero_degradation_top1000.bed')


region <- rep(rse_gene$Region, elementNROWS(rse_gene$SAMPLE_ID))
region_dir <- ifelse(region == 'DLPFC', 'DLPFC', 'Hippo')
region_path <- paste0('/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/', region_dir, '_RiboZero/Coverage/')

forwardBw = paste0(region_path, unlist(rse_gene$SAMPLE_ID), ".Forward.bw")
reverseBw = paste0(region_path, unlist(rse_gene$SAMPLE_ID), ".Reverse.bw")

stopifnot(all(file.exists(c(forwardBw,reverseBw))))

names(forwardBw) = names(reverseBw) = unlist(rse_gene$SAMPLE_ID)


covForward = coverage_bwtool(forwardBw, regs, strand = "+", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covForward$bigwig_path = NULL
covForward$bigwig_file = NULL

covReverse = coverage_bwtool(reverseBw, regs, strand = "-", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## combine
cov_rse = rbind(covForward, covReverse)	
rownames(cov_rse) = rowData(cov_rse)$name
cov_rse = cov_rse[regs$name,]

## divide by number of reads
assays(cov_rse)$counts = assays(cov_rse)$counts/100 # divide by read length

## make positive
assays(cov_rse)$counts = abs(assays(cov_rse)$counts) 

dir.create('rdas', showWarnings = FALSE)
save(cov_rse, file = 'rdas/degradation_rse_phase2_usingJoint.rda')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
