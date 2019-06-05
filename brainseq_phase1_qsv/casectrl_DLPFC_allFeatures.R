# Adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/casectrl_DLPFC_allFeatures.R
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')


### Run with the qSVs made without the age>17 samples
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda", verbose = TRUE)
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseExon_merged_n732.rda", verbose = TRUE)
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseJxn_merged_n732.rda", verbose = TRUE)
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseTx_merged_n732.rda", verbose = TRUE)
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase1_qsv/rdas/brainseq_phase1_qsvs_age17.Rdata', verbose = TRUE)

## Check https://github.com/LieberInstitute/qsva_brain/blob/master/process_bsp1_and_bsp3.R#L36-L39
stopifnot(all(table(rse_gene[, keepIndex]$Dx) - c(67, 221, 147, 179) == 0))

## Drop samples absent in mod and modQsVA
rse_gene <- rse_gene[, keepIndex]
rse_jxn <- rse_jxn[, keepIndex]
rse_exon <- rse_exon[, keepIndex]
rse_tx <- rse_tx[, keepIndex]

## Keep diagnosis-specific samples
keepIndex = which(rse_gene$Dx %in% c('Control', 'Schizo'))
table(rse_gene[, keepIndex]$Dx)
# Control  Schizo
#     221     179
rse_gene <- rse_gene[, keepIndex]
rse_jxn <- rse_jxn[, keepIndex]
rse_exon <- rse_exon[, keepIndex]
rse_tx <- rse_tx[, keepIndex]

## Fix the mod/modQsva
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase1_qsv/rdas/brainseq_phase1_cov_rse_age17.Rdata', verbose = TRUE)
cov_rse <- cov_rse[, keepIndex]
stopifnot(identical(table(cov_rse$Dx), table(rse_gene$Dx)))
mod <- model.matrix(
    ~ Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(cov_rse)
)
modQsva = cbind(mod, qSVs[keepIndex, ])

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
pdf('pdf/voom_qsva_BSP1_DLPFC_gene.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]



##### Exon ######
dee = DGEList(counts = assays(rse_exon)$counts,
	genes = rowData(rse_exon))
dee = calcNormFactors(dee)
pdf('pdf/voom_qsva_BSP1_DLPFC_exon.pdf', useDingbats = FALSE)
vExon = voom(dee,modQsva, plot=TRUE)
dev.off()
fitExon = lmFit(vExon)
eBExon = eBayes(fitExon)
sigExon = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon))
outExon = sigExon[rownames(rse_exon),]

##### Junction ######
dje = DGEList(counts = assays(rse_jxn)$counts,
	genes = rowData(rse_jxn))
dje = calcNormFactors(dje)
pdf('pdf/voom_qsva_BSP1_DLPFC_jxn.pdf', useDingbats = FALSE)
vJxn = voom(dje,modQsva, plot=TRUE)
dev.off()
fitJxn = lmFit(vJxn)
eBJxn = eBayes(fitJxn)
sigJxn = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jxn))
outJxn = sigJxn[rownames(rse_jxn),]

##### Transcript ######
fitTx = lmFit(log2(assays(rse_tx)$tpm + 0.5), modQsva)
eBTx = eBayes(fitTx)
sigTx = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx))
outTx = sigTx[rownames(rse_tx),]
outTx <- cbind(outTx, rowData(rse_tx))

save(outGene, outExon, outJxn,outTx,
	file = "rdas/dxStats_dlpfc_filtered_qSVA_BSP1_DLPFC.rda")


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
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.5.1)
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  callr                  3.2.0     2019-03-15 [2] CRAN (R 3.5.1)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.5.3)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  desc                   1.2.0     2018-05-01 [2] CRAN (R 3.5.1)
#  devtools             * 2.0.1     2018-10-26 [1] CRAN (R 3.5.1)
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  edgeR                * 3.24.3    2019-01-02 [1] Bioconductor
#  fs                     1.2.7     2019-03-19 [2] CRAN (R 3.5.1)
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
#  limma                * 3.38.3    2018-12-02 [1] Bioconductor
#  locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.3)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgbuild               1.0.3     2019-03-20 [2] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  pkgload                1.0.2     2018-10-29 [2] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  processx               3.3.0     2019-03-10 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  ps                     1.3.0     2018-12-21 [2] CRAN (R 3.5.1)
#  purrr                  0.3.2     2019-03-15 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.1     2019-03-17 [1] CRAN (R 3.5.3)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  remotes                2.0.2     2018-10-30 [1] CRAN (R 3.5.1)
#  rlang                  0.3.3     2019-03-29 [1] CRAN (R 3.5.3)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo            1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  testthat               2.0.1     2018-10-13 [1] CRAN (R 3.5.1)
#  tibble                 2.1.1     2019-03-16 [1] CRAN (R 3.5.3)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  usethis              * 1.4.0     2018-08-14 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.6       2019-04-02 [1] CRAN (R 3.5.3)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
