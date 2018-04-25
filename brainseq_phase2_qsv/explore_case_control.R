library('clusterProfiler')
library('gplots')
library('GenomicRanges')
library('devtools')

## Load case-control results
files <- c(
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_matchHIPPO.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda'

)
outGene <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene)
})

outGene0 <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene0)
})

outGeneNoAdj <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGeneNoAdj)
})


names(outGeneNoAdj) <- names(outGene0) <- names(outGene) <- c('HIPPO_allQSV', 'DLPFC_allQSV', 'HIPPO_noHGoldQSV', 'DLPFC_noHGoldQSV', 'HIPPO_matchQSV', 'DLPFC_matchQSV')

## Load BrainSeq Phase 1 and Common Mind results
load('/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda', verbose = TRUE)

prev <- list(
    'DLPFC_P1' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = outStatsExprs$Gene$fdr_qsva
        ),
    'CMC' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = p.adjust(outStatsExprs$Gene$CMC_pval_qsva, method = 'fdr')
        )
)

outGene <- c(outGene, prev)


## Get number of DE genes at different FDR cutoffs
n_de <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
    xx <- sapply(outGene, function(x) {
        table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')))
    })
    cbind(xx, cutoff = cut)
}))
n_de

data.frame(colnames(n_de))


## Explore number of DE genes across models
make_venn <- function(i, txt = 'HIPPO_', cut = 0.1) {
    vinfo <- lapply(outGene[i], function(x) {
        x$ensemblID[x$adj.P.Val < cut]
    })
    names(vinfo) <- gsub(txt, '', names(vinfo))
    venn(vinfo) + title(paste('FDR cutoff:', cut))
}

make_venn2 <- function(i, txt = 'QSV|IPPO|LPFC') {
    make_venn(i, txt = txt, cut = 0.05)
    make_venn(i, txt = txt)
}

pdf('pdf/venn_across_models.pdf', useDingbats = FALSE)
make_venn2(c(1, 3, 5))
make_venn2(c(2, 4, 6))
make_venn2(3:6)
make_venn2(c(4, 6, 7, 8))
make_venn2(5:6)
make_venn2(c(5, 6, 7, 8))
dev.off()


## Degradation results
load("/dcl01/ajaffe/data/lab/qsva_brain/ERs/rdas/DLPFC_Plus_HIPPO_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda", verbose = TRUE)
load("/dcl01/ajaffe/data/lab/qsva_brain/ERs/rdas/DLPFC_HIPPO_degradationStats_hg38.rda", verbose = TRUE)

## For DE_qual plots
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}

plot_dequal <- function(out_input, degrade_input, var = 't', xlabtxt = 'case-control', main, ylabtxt = '') {
	both <- intersect(rownames(out_input), rownames(degrade_input))
	degrade <- degrade_input[both, ]
	interest <- out_input[both, ]

	stopifnot(identical(rownames(degrade), rownames(interest)))
	corr = signif(cor( degrade[, var],  interest[, var]), 3)
	plot(y = degrade[, var], x = interest[, var], xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlabtxt), ylab = paste(ylabtxt, var, 'degradation'), main = main, col = add.alpha('black', 1/10), pch = 16)
	legend('topleft', legend = paste('r =', corr))
}


## Make de_qual plots
plot_dequal2 <- function(i) {
    geneinfo <- list(outGeneNoAdj[[i]], outGeneNoAdj[[i]], outGene0[[i]], outGene0[[i]], outGene[[i]], outGene[[i]])
    makedequal <- function(v, gene, main, xlab) {
        plot_dequal(gene, degradeStats, var = v, xlabtxt = xlab, main = main, ylabtxt = 'Combined')
        plot_dequal(gene, degradeStats_HIPPO, var = v, xlabtxt = xlab, main = main, ylabtxt = 'HIPPO')
        plot_dequal(gene, degradeStats_DLPFC, var = v, xlabtxt = xlab, main = main, ylabtxt = 'DLPFC')
        plot_dequal(gene, degradeStatsInt, var = v, xlabtxt = xlab, main = main, ylabtxt = 'Combined adj interaction')
    }
    par(mfcol = c(4, 2))
    mapply(makedequal,
           v = rep(c('t', 'logFC'), 3),
           gene = geneinfo,
           main = rep(names(outGene)[i], 6),
           xlab = rep(c('case-control (Dx only)', 'case-control (without qSVs)', 'case-control (with qSVs)'), each = 2))
}


pdf('pdf/dequal_plots.pdf', useDingbats = FALSE, width = 8, height = 16)
for(i in 1:6) plot_dequal2(i)
dev.off()


## Gene ontology analysis
de_genes <- lapply(outGene[5:6], function(x) {
    x$ensemblID[x$adj.P.Val < 0.05]
})


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
