library('SummarizedExperiment')
library('sva')
library('jaffelab')
library('sessioninfo')

## Create directories
dir.create('brainseq_phase4_5_qsv', showWarnings = FALSE)
dir.create('brainseq_phase4_5_qsv/pdf', showWarnings = FALSE)
dir.create('brainseq_phase4_5_qsv/rdas', showWarnings = FALSE)

## BSP4 and 5

### mPFC
load('/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase4and5/count_data/degradation_rse_phase4and5_mPFC.rda', verbose = TRUE)
load('/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase4and5/count_data/brainseq_phases4and5_hg38_rseGene_merged_n492.rda', verbose = TRUE)

## Match by sample id
m <- match(colData(cov_rse)$sample, colData(rse_gene)$SAMPLE_ID)

## Check that all samples matched
stopifnot(all(!is.na(m)))

## Replace the degradation RSE colData() with the gene-level one
colData(cov_rse) <- colData(rse_gene)[m, ]

## Add genotype data
mds <- read.csv('/dcl01/lieber/ajaffe/lab/brainseq_phase4and5/genotype_data/BrainSeq_Phase4and5_RiboZero_Genotypes_n492_snpPCs.csv')

## Match brains to genotype
m_geno <- match(colData(cov_rse)$SAMPLE_ID, mds$X)
stopifnot(all(!is.na(m_geno)))
colData(cov_rse) <- cbind(colData(cov_rse), mds[m_geno, 2:11])

## Fix mitoRate, rRNA_rate and totalAssignedGene
# colData(cov_rse)$mitoRate = sapply(colData(cov_rse)$mitoRate,mean)
# colData(cov_rse)$rRNA_rate = sapply(colData(cov_rse)$rRNA_rate,mean)
# colData(cov_rse)$totalAssignedGene = sapply(colData(cov_rse)$totalAssignedGene,mean)

## Drop age 17 or less
keepIndex <- which(colData(cov_rse)$Age > 17)
length(keepIndex)
# [1] 299

cov_rse <- cov_rse[, keepIndex]
with(colData(cov_rse), table(Dx))
# Dx
# Bipolar Control     MDD  Schizo
#      38     107      60      94

## Find qsvs
mod <- model.matrix(
    ~ Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(cov_rse)
)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
k
# [1] 14
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 54.800  7.720  4.300  4.040  1.770  1.490  1.360  1.040  0.910  0.769
# [11]  0.719  0.637  0.572  0.536
modQsva = cbind(mod, qSVs)

pdf('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/pdf/qsvs_var_explained_mPFC_age17.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/rdas/brainseq_phase4_5_mPFC_qsvs_age17.Rdata')

## Save subsetted cov_rse with the snpPCs already added
save(cov_rse, file = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/rdas/brainseq_phase4_5_mPFC_cov_rse_age17.Rdata')




### HIPPO
rm(list = ls())
load('/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase4and5/count_data/degradation_rse_phase4and5_HIPPO.rda', verbose = TRUE)
load('/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase4and5/count_data/brainseq_phases4and5_hg38_rseGene_merged_n492.rda', verbose = TRUE)

## Match by sample id
m <- match(colData(cov_rse)$sample, colData(rse_gene)$SAMPLE_ID)

## Check that all samples matched
stopifnot(all(!is.na(m)))

## Replace the degradation RSE colData() with the gene-level one
colData(cov_rse) <- colData(rse_gene)[m, ]

## Add genotype data
mds <- read.csv('/dcl01/lieber/ajaffe/lab/brainseq_phase4and5/genotype_data/BrainSeq_Phase4and5_RiboZero_Genotypes_n492_snpPCs.csv')

## Match brains to genotype
m_geno <- match(colData(cov_rse)$SAMPLE_ID, mds$X)
stopifnot(all(!is.na(m_geno)))
colData(cov_rse) <- cbind(colData(cov_rse), mds[m_geno, 2:11])

## Fix mitoRate, rRNA_rate and totalAssignedGene
# colData(cov_rse)$mitoRate = sapply(colData(cov_rse)$mitoRate,mean)
# colData(cov_rse)$rRNA_rate = sapply(colData(cov_rse)$rRNA_rate,mean)
# colData(cov_rse)$totalAssignedGene = sapply(colData(cov_rse)$totalAssignedGene,mean)

## Drop age 17 or less
keepIndex <- which(colData(cov_rse)$Age > 17)
length(keepIndex)
# [1] 92

cov_rse <- cov_rse[, keepIndex]
with(colData(cov_rse), table(Dx))
# Dx
# Bipolar     MDD
#      35      57

## Find qsvs
mod <- model.matrix(
    ~ Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(cov_rse)
)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
k
# [1] 8
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
# [1] 46.90 12.70  5.73  3.24  2.55  2.19  1.77  1.36
modQsva = cbind(mod, qSVs)

pdf('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/pdf/qsvs_var_explained_HIPPO_age17.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/rdas/brainseq_phase4_5_HIPPO_qsvs_age17.Rdata')

## Save subsetted cov_rse with the snpPCs already added
save(cov_rse, file = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/rdas/brainseq_phase4_5_HIPPO_cov_rse_age17.Rdata')






### DLPFC
rm(list = ls())
load('/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase4and5/count_data/degradation_rse_phase4and5_DLPFC.rda', verbose = TRUE)
load('/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase4and5/count_data/brainseq_phases4and5_hg38_rseGene_merged_n492.rda', verbose = TRUE)

## Match by sample id
m <- match(colData(cov_rse)$sample, colData(rse_gene)$SAMPLE_ID)

## Check that all samples matched
stopifnot(all(!is.na(m)))

## Replace the degradation RSE colData() with the gene-level one
colData(cov_rse) <- colData(rse_gene)[m, ]

## Add genotype data
mds <- read.csv('/dcl01/lieber/ajaffe/lab/brainseq_phase4and5/genotype_data/BrainSeq_Phase4and5_RiboZero_Genotypes_n492_snpPCs.csv')

## Match brains to genotype
m_geno <- match(colData(cov_rse)$SAMPLE_ID, mds$X)
stopifnot(all(!is.na(m_geno)))
colData(cov_rse) <- cbind(colData(cov_rse), mds[m_geno, 2:11])

## Fix mitoRate, rRNA_rate and totalAssignedGene
# colData(cov_rse)$mitoRate = sapply(colData(cov_rse)$mitoRate,mean)
# colData(cov_rse)$rRNA_rate = sapply(colData(cov_rse)$rRNA_rate,mean)
# colData(cov_rse)$totalAssignedGene = sapply(colData(cov_rse)$totalAssignedGene,mean)

## Drop age 17 or less
keepIndex <- which(colData(cov_rse)$Age > 17)
length(keepIndex)
# [1] 101

cov_rse <- cov_rse[, keepIndex]
with(colData(cov_rse), table(Dx))
# Dx
# Bipolar     MDD
#      39      62

## Find qsvs
mod <- model.matrix(
    ~ Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(cov_rse)
)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
k
# [1] 8
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
# [1] 49.40  9.51  6.31  4.96  2.98  1.67  1.29  1.10
modQsva = cbind(mod, qSVs)

pdf('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/pdf/qsvs_var_explained_DLPFC_age17.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/rdas/brainseq_phase4_5_DLPFC_qsvs_age17.Rdata')

## Save subsetted cov_rse with the snpPCs already added
save(cov_rse, file = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase4_5_qsv/rdas/brainseq_phase4_5_DLPFC_cov_rse_age17.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.3 Patched (2019-03-11 r76311)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-06-05
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  annotate               1.60.1    2019-03-07 [1] Bioconductor
#  AnnotationDbi          1.44.0    2018-10-30 [1] Bioconductor
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bit                    1.1-14    2018-05-29 [2] CRAN (R 3.5.1)
#  bit64                  0.9-7     2017-05-08 [2] CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 [2] CRAN (R 3.5.0)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.5.3)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DBI                    1.0.0     2018-05-02 [2] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  genefilter           * 1.64.0    2018-10-30 [1] Bioconductor
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.5.1)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.5.0     2019-03-15 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.3)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.5.1)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.3)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  mgcv                 * 1.8-27    2019-02-06 [3] CRAN (R 3.5.3)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  nlme                 * 3.1-137   2018-04-07 [3] CRAN (R 3.5.3)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.3.2     2019-03-15 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.1     2019-03-17 [1] CRAN (R 3.5.3)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  rlang                  0.3.3     2019-03-29 [1] CRAN (R 3.5.3)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  RSQLite                2.1.1     2018-05-06 [2] CRAN (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  survival               2.43-3    2018-11-26 [3] CRAN (R 3.5.3)
#  sva                  * 3.30.1    2019-01-04 [2] Bioconductor
#  tibble                 2.1.1     2019-03-16 [1] CRAN (R 3.5.3)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.6       2019-04-02 [1] CRAN (R 3.5.3)
#  XML                    3.98-1.19 2019-03-06 [2] CRAN (R 3.5.1)
#  xtable                 1.8-3     2018-08-29 [2] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
