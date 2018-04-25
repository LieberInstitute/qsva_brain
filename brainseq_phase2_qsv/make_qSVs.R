

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
# [1] 57.800 15.700  3.360  2.700  2.360  1.160  1.040  0.848  0.750  0.661
# [11]  0.587  0.511  0.461  0.419  0.356  0.312  0.291  0.284  0.271  0.259
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/brainseq_phase2_qsvs.Rdata')
pdf('pdf/qsvs_var_explained.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()




### Re-make qSVs without age>17 and no HIPPO gold samples
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop age<=17 and HIPPO Gold samples
hipHMR <- as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
hipHMR[rse_gene$Region == 'DLPFC'] <- TRUE

keepIndex <- which(rse_gene$Age > 17 & hipHMR)
length(keepIndex)
# [1] 712
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 70.200  5.540  3.070  2.600  1.180  1.070  0.954  0.770  0.646  0.577
# [11]  0.512  0.422  0.400  0.334  0.316  0.315  0.278  0.276  0.243  0.235
# [21]  0.208  0.204
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'rdas/brainseq_phase2_qsvs_age17_noHGold.Rdata')
pdf('pdf/qsvs_var_explained_age17_noHGold.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()







### Re-make qSVs without age>17 and no HIPPO gold samples
## Only HIPPO
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop age<=17 and HIPPO Gold samples
hipHMR <- as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
hipHMR[rse_gene$Region == 'DLPFC'] <- TRUE

keepIndex <- which(rse_gene$Age > 17 & hipHMR & rse_gene$Region == 'HIPPO')
length(keepIndex)
# [1] 333
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 60.200  7.850  4.630  2.460  2.040  1.640  1.240  1.060  0.974  0.770
# [11]  0.665  0.611  0.580  0.489  0.458  0.432
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata')
pdf('pdf/qsvs_var_explained_age17_noHGold_HIPPO.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()




### Re-make qSVs without age>17 and no HIPPO gold samples
## Only DLPFC
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop age<=17 and HIPPO Gold samples
hipHMR <- as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
hipHMR[rse_gene$Region == 'DLPFC'] <- TRUE

keepIndex <- which(rse_gene$Age > 17 & hipHMR & rse_gene$Region == 'DLPFC')
length(keepIndex)
# [1] 379
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 69.600  5.620  2.730  1.900  1.420  1.130  0.875  0.869  0.845  0.633
# [11]  0.576  0.498  0.469  0.424  0.401
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata')
pdf('pdf/qsvs_var_explained_age17_noHGold_DLPFC.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()




## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.3 Patched (2018-01-20 r74142)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-04-25
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  annotate               1.56.2    2018-04-18 Bioconductor
#  AnnotationDbi          1.40.0    2017-11-29 Bioconductor
#  base                 * 3.4.3     2018-01-20 local
#  Biobase              * 2.38.0    2017-11-07 Bioconductor
#  BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
#  BiocParallel         * 1.12.0    2017-11-29 Bioconductor
#  bit                    1.1-12    2014-04-09 CRAN (R 3.4.1)
#  bit64                  0.9-7     2017-05-08 CRAN (R 3.4.1)
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
#  blob                   1.1.1     2018-03-25 CRAN (R 3.4.3)
#  cellranger             1.1.0     2016-07-27 CRAN (R 3.4.1)
#  colorout             * 1.2-0     2018-02-19 Github (jalvesaq/colorout@2f01173)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
#  compiler               3.4.3     2018-01-20 local
#  datasets             * 3.4.3     2018-01-20 local
#  DBI                    0.8       2018-03-02 CRAN (R 3.4.3)
#  DelayedArray         * 0.4.1     2017-11-07 Bioconductor
#  devtools             * 1.13.5    2018-02-18 CRAN (R 3.4.3)
#  digest                 0.6.15    2018-01-28 cran (@0.6.15)
#  genefilter           * 1.60.0    2017-11-29 Bioconductor
#  GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
#  GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
#  GenomicRanges        * 1.30.3    2018-04-18 Bioconductor
#  ggplot2                2.2.1     2016-12-30 CRAN (R 3.4.1)
#  graphics             * 3.4.3     2018-01-20 local
#  grDevices            * 3.4.3     2018-01-20 local
#  grid                   3.4.3     2018-01-20 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.4.3)
#  httpuv                 1.3.6.2   2018-03-02 CRAN (R 3.4.3)
#  IRanges              * 2.12.0    2017-11-29 Bioconductor
#  jaffelab             * 0.99.20   2018-04-19 Github (LieberInstitute/jaffelab@04c470a)
#  later                  0.7.1     2018-03-07 CRAN (R 3.4.3)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
#  limma                  3.34.9    2018-04-18 Bioconductor
#  Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
#  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.3     2018-01-20 local
#  mgcv                 * 1.8-22    2017-09-24 CRAN (R 3.4.3)
#  munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
#  nlme                 * 3.1-131   2017-02-06 CRAN (R 3.4.3)
#  parallel             * 3.4.3     2018-01-20 local
#  pillar                 1.2.1     2018-02-27 CRAN (R 3.4.3)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
#  png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.4.1)
#  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.4.1)
#  Rcpp                   0.12.16   2018-03-13 CRAN (R 3.4.3)
#  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
#  readxl               * 1.1.0     2018-04-20 CRAN (R 3.4.3)
#  rlang                  0.2.0     2018-02-20 CRAN (R 3.4.3)
#  rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
#  RSQLite                2.1.0     2018-03-29 CRAN (R 3.4.3)
#  S4Vectors            * 0.16.0    2017-11-29 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.4.2)
#  servr                  0.9       2018-03-25 CRAN (R 3.4.3)
#  splines                3.4.3     2018-01-20 local
#  stats                * 3.4.3     2018-01-20 local
#  stats4               * 3.4.3     2018-01-20 local
#  SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
#  survival               2.41-3    2017-04-04 CRAN (R 3.4.3)
#  sva                  * 3.26.0    2017-11-07 Bioconductor
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.4.3)
#  tools                  3.4.3     2018-01-20 local
#  utils                * 3.4.3     2018-01-20 local
#  withr                  2.1.2     2018-03-15 CRAN (R 3.4.3)
#  xfun                   0.1       2018-01-22 CRAN (R 3.4.3)
#  XML                    3.98-1.10 2018-02-19 CRAN (R 3.4.3)
#  xtable                 1.8-2     2016-02-05 CRAN (R 3.4.1)
#  XVector                0.18.0    2017-11-29 Bioconductor
#  zlibbioc               1.24.0    2017-11-07 Bioconductor
#
