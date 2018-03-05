library('SummarizedExperiment')
library('ggplot2')
library('devtools')

## For getting the sample groups
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata')

## Load degradation quantifications
load('rdas/degradation_rse_phase2_usingJoint.rda')

n_samples <- elementNROWS(rse_gene$SAMPLE_ID)

errdf <- do.call(rbind, lapply(2:4, function(n) {
    n_i <- which(n_samples == n)
    do.call(rbind, lapply(n_i, function(i) {
        dat <- assays(cov_rse[, unlist(rse_gene$SAMPLE_ID[i]) ])$counts
        mean <- rowMeans(dat)
        err <- sweep(dat, 1, mean)
        data.frame(rmse = sqrt(colMeans(err)^2), sample = colnames(err), n = n, i = i)
    }))
}))

dim(errdf)

dir.create('pdf', showWarnings = FALSE)
pdf('pdf/RMSE.pdf')
ggplot(errdf, aes(x = as.factor(n), y = rmse)) + geom_boxplot() + theme_light(base_size = 18) + xlab('Number of replicates') + ylab('RMSE across all 1000 degradation ERs')
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
