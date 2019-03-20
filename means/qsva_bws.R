library('recount')
library('recount.bwtool')
library('devtools')
library('SummarizedExperiment')
load('/dcl01/ajaffe/data/lab/qsva_brain/expr_data/rda/rse_gene.Rdata')

## Forward
bw_f_dlpfc <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'DLPFC'], '.Forward.bw')

bw_f_hippo <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'HIPPO'], '.Forward.bw')

bw_f<- c(bw_f_dlpfc, bw_f_hippo)

compute_mean(bws = bw_f, outfile ='mean_forward',tempdir = 'forward_tmp')



## Reverse
bw_r_dlpfc <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'DLPFC'], '.Reverse.bw')

bw_r_hippo <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'HIPPO'], '.Reverse.bw')

bw_r<- c(bw_r_dlpfc, bw_r_hippo)

compute_mean(bws = bw_r, outfile ='mean_reverse',tempdir = 'reverse_tmp')

Sys.time()
proc.time()
options(width = 120)
session_info()
