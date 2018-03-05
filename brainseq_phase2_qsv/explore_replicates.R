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
dir.create('rdas', showWarnings = FALSE)
save(errdf, file = 'rdas/errdf.Rdata')

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

# > session_info()
# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.3 Patched (2018-01-20 r74142)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-03-05
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  base                 * 3.4.3     2018-01-20 local
#  Biobase              * 2.38.0    2017-11-07 Bioconductor
#  BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
#  colorout             * 1.2-0     2018-02-19 Github (jalvesaq/colorout@2f01173)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
#  compiler               3.4.3     2018-01-20 local
#  datasets             * 3.4.3     2018-01-20 local
#  DelayedArray         * 0.4.1     2017-11-07 Bioconductor
#  devtools             * 1.13.4    2017-11-09 CRAN (R 3.4.2)
#  digest                 0.6.15    2018-01-28 cran (@0.6.15)
#  GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
#  GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
#  GenomicRanges        * 1.30.2    2018-02-17 Bioconductor
#  ggplot2              * 2.2.1     2016-12-30 CRAN (R 3.4.1)
#  graphics             * 3.4.3     2018-01-20 local
#  grDevices            * 3.4.3     2018-01-20 local
#  grid                   3.4.3     2018-01-20 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
#  htmlwidgets            0.9       2017-07-10 CRAN (R 3.4.1)
#  httpuv                 1.3.5     2017-07-04 CRAN (R 3.4.1)
#  IRanges              * 2.12.0    2017-11-29 Bioconductor
#  labeling               0.3       2014-08-23 CRAN (R 3.4.1)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
#  Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
#  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.3     2018-01-20 local
#  munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
#  parallel             * 3.4.3     2018-01-20 local
#  pillar                 1.1.0     2018-01-14 CRAN (R 3.4.2)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
#  png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
#  Rcpp                   0.12.14   2017-11-23 CRAN (R 3.4.2)
#  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
#  rlang                  0.1.6     2017-12-21 CRAN (R 3.4.2)
#  rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
#  S4Vectors            * 0.16.0    2017-11-29 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
#  servr                  0.8       2017-11-06 CRAN (R 3.4.3)
#  stats                * 3.4.3     2018-01-20 local
#  stats4               * 3.4.3     2018-01-20 local
#  SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
#  tibble                 1.4.1     2017-12-25 CRAN (R 3.4.2)
#  tools                  3.4.3     2018-01-20 local
#  utils                * 3.4.3     2018-01-20 local
#  withr                  2.1.1     2017-12-19 CRAN (R 3.4.2)
#  XVector                0.18.0    2017-11-29 Bioconductor
#  zlibbioc               1.24.0    2017-11-07 Bioconductor
