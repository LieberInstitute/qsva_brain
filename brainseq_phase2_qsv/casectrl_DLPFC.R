
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

#load qsvBonf, qSVs, mod, modQsva object from brainseq data
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs.Rdata", verbose = TRUE)
#filter for age and DLPFC; no age in rse_gene
keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('pdf/dlpfc_voom_qsva.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('pdf/dlpfc_voom_noqsva.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]


## no adjustment vars
pdf('pdf/dlpfc_voom_noadj.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]

stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "rdas/dxStats_dlpfc_filtered_qSVA_geneLevel.rda")

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
#  date     2018-04-18
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
#  devtools             * 1.13.5    2018-02-18 CRAN (R 3.4.3)
#  digest                 0.6.15    2018-01-28 cran (@0.6.15)
#  edgeR                * 3.20.9    2018-04-18 Bioconductor
#  GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
#  GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
#  GenomicRanges        * 1.30.3    2018-04-18 Bioconductor
#  ggplot2                2.2.1     2016-12-30 CRAN (R 3.4.1)
#  graphics             * 3.4.3     2018-01-20 local
#  grDevices            * 3.4.3     2018-01-20 local
#  grid                   3.4.3     2018-01-20 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
#  htmlwidgets            1.0       2018-01-20 CRAN (R 3.4.3)
#  httpuv                 1.3.6.2   2018-03-02 CRAN (R 3.4.3)
#  IRanges              * 2.12.0    2017-11-29 Bioconductor
#  jaffelab             * 0.99.18   2018-02-22 Github (LieberInstitute/jaffelab@a8e6430)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
#  limma                * 3.34.9    2018-04-18 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.4.1)
#  Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
#  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.3     2018-01-20 local
#  munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
#  parallel             * 3.4.3     2018-01-20 local
#  pillar                 1.2.1     2018-02-27 CRAN (R 3.4.3)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
#  png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.4.1)
#  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.4.1)
#  Rcpp                   0.12.16   2018-03-13 CRAN (R 3.4.3)
#  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
#  rlang                  0.2.0     2018-02-20 CRAN (R 3.4.3)
#  rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
#  S4Vectors            * 0.16.0    2017-11-29 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.4.2)
#  servr                  0.9       2018-03-25 CRAN (R 3.4.3)
#  stats                * 3.4.3     2018-01-20 local
#  stats4               * 3.4.3     2018-01-20 local
#  SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.4.3)
#  tools                  3.4.3     2018-01-20 local
#  utils                * 3.4.3     2018-01-20 local
#  withr                  2.1.2     2018-03-15 CRAN (R 3.4.3)
#  xfun                   0.1       2018-01-22 CRAN (R 3.4.3)
#  XVector                0.18.0    2017-11-29 Bioconductor
#  zlibbioc               1.24.0    2017-11-07 Bioconductor


