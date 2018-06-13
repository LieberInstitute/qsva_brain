library('clusterProfiler')
library('gplots')
library('GenomicRanges')
library('devtools')
library('VennDiagram')
library('RColorBrewer')
library('ggplot2')
library('jaffelab')
library('limma')
library('edgeR')
library('SummarizedExperiment')

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
    'BSP1' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = outStatsExprs$Gene$fdr_qsva,
        logFC = outStatsExprs$Gene$log2FC_qsva,
        t = outStatsExprs$Gene$tstat_qsva
        ),
    'CMC' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = p.adjust(outStatsExprs$Gene$CMC_pval_qsva, method = 'fdr'),
        logFC = outStatsExprs$Gene$CMC_log2FC_qsva,
        t = outStatsExprs$Gene$CMC_tstat_qsva
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

options(width = 160)
n_de
#       HIPPO_allQSV DLPFC_allQSV HIPPO_noHGoldQSV DLPFC_noHGoldQSV HIPPO_matchQSV DLPFC_matchQSV  BSP1   CMC cutoff
# FALSE        24652        24307            24635            24190          24604          24407 23939 23746   0.05
# TRUE             0          345               17              462             48            245   183   376   0.05
# FALSE        24648        23923            24592            23800          24551          24020 23616 23220   0.10
# TRUE             4          729               60              852            101            632   506   902   0.10
# FALSE        24601        23543            24502            23360          24489          23618 23205 22661   0.15
# TRUE            51         1109              150             1292            163           1034   917  1461   0.15
# FALSE        24524        23112            24369            22857          24320          23103 22736 22002   0.20
# TRUE           128         1540              283             1795            332           1549  1386  2120   0.20

data.frame(colnames(n_de))
#     colnames.n_de.
# 1     HIPPO_allQSV
# 2     DLPFC_allQSV
# 3 HIPPO_noHGoldQSV
# 4 DLPFC_noHGoldQSV
# 5   HIPPO_matchQSV
# 6   DLPFC_matchQSV
# 7             BSP1
# 8              CMC
# 9           cutoff

n_de_sign <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
    xx <- lapply(outGene[5:8], function(x) {
        y <- table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')), factor(sign(x$logFC), levels = c(-1, 0, 1)))
        data.frame(de_status = rownames(y), n = as.vector(y), sign = rep(colnames(y), each = 2), cutoff = cut)
    })
    for(i in seq_len(length(xx))) { xx[[i]]$model = names(xx)[i] }
    names(xx) <- NULL
    do.call(rbind, xx)
}))
n_de_sign$group <- ifelse(n_de_sign$sign == 0, 'none', ifelse(n_de_sign$sign == -1, 'Control', 'Schizo'))
n_de_sign
#    de_status     n sign cutoff          model   group
# 1      FALSE 12825   -1   0.05 HIPPO_matchQSV Control
# 2       TRUE    27   -1   0.05 HIPPO_matchQSV Control
# 3      FALSE     0    0   0.05 HIPPO_matchQSV    none
# 4       TRUE     0    0   0.05 HIPPO_matchQSV    none
# 5      FALSE 11779    1   0.05 HIPPO_matchQSV  Schizo
# 6       TRUE    21    1   0.05 HIPPO_matchQSV  Schizo
# 7      FALSE 13065   -1   0.05 DLPFC_matchQSV Control
# 8       TRUE   142   -1   0.05 DLPFC_matchQSV Control
# 9      FALSE     0    0   0.05 DLPFC_matchQSV    none
# 10      TRUE     0    0   0.05 DLPFC_matchQSV    none
# 11     FALSE 11342    1   0.05 DLPFC_matchQSV  Schizo
# 12      TRUE   103    1   0.05 DLPFC_matchQSV  Schizo
# 13     FALSE 11415   -1   0.05           BSP1 Control
# 14      TRUE    50   -1   0.05           BSP1 Control
# 15     FALSE     0    0   0.05           BSP1    none
# 16      TRUE     0    0   0.05           BSP1    none
# 17     FALSE 12524    1   0.05           BSP1  Schizo
# 18      TRUE   133    1   0.05           BSP1  Schizo
# 19     FALSE 12101   -1   0.05            CMC Control
# 20      TRUE   243   -1   0.05            CMC Control
# 21     FALSE    21    0   0.05            CMC    none
# 22      TRUE     0    0   0.05            CMC    none
# 23     FALSE 11624    1   0.05            CMC  Schizo
# 24      TRUE   133    1   0.05            CMC  Schizo
# 25     FALSE 12800   -1   0.10 HIPPO_matchQSV Control
# 26      TRUE    52   -1   0.10 HIPPO_matchQSV Control
# 27     FALSE     0    0   0.10 HIPPO_matchQSV    none
# 28      TRUE     0    0   0.10 HIPPO_matchQSV    none
# 29     FALSE 11751    1   0.10 HIPPO_matchQSV  Schizo
# 30      TRUE    49    1   0.10 HIPPO_matchQSV  Schizo
# 31     FALSE 12828   -1   0.10 DLPFC_matchQSV Control
# 32      TRUE   379   -1   0.10 DLPFC_matchQSV Control
# 33     FALSE     0    0   0.10 DLPFC_matchQSV    none
# 34      TRUE     0    0   0.10 DLPFC_matchQSV    none
# 35     FALSE 11192    1   0.10 DLPFC_matchQSV  Schizo
# 36      TRUE   253    1   0.10 DLPFC_matchQSV  Schizo
# 37     FALSE 11288   -1   0.10           BSP1 Control
# 38      TRUE   177   -1   0.10           BSP1 Control
# 39     FALSE     0    0   0.10           BSP1    none
# 40      TRUE     0    0   0.10           BSP1    none
# 41     FALSE 12328    1   0.10           BSP1  Schizo
# 42      TRUE   329    1   0.10           BSP1  Schizo
# 43     FALSE 11796   -1   0.10            CMC Control
# 44      TRUE   548   -1   0.10            CMC Control
# 45     FALSE    21    0   0.10            CMC    none
# 46      TRUE     0    0   0.10            CMC    none
# 47     FALSE 11403    1   0.10            CMC  Schizo
# 48      TRUE   354    1   0.10            CMC  Schizo
# 49     FALSE 12768   -1   0.15 HIPPO_matchQSV Control
# 50      TRUE    84   -1   0.15 HIPPO_matchQSV Control
# 51     FALSE     0    0   0.15 HIPPO_matchQSV    none
# 52      TRUE     0    0   0.15 HIPPO_matchQSV    none
# 53     FALSE 11721    1   0.15 HIPPO_matchQSV  Schizo
# 54      TRUE    79    1   0.15 HIPPO_matchQSV  Schizo
# 55     FALSE 12583   -1   0.15 DLPFC_matchQSV Control
# 56      TRUE   624   -1   0.15 DLPFC_matchQSV Control
# 57     FALSE     0    0   0.15 DLPFC_matchQSV    none
# 58      TRUE     0    0   0.15 DLPFC_matchQSV    none
# 59     FALSE 11035    1   0.15 DLPFC_matchQSV  Schizo
# 60      TRUE   410    1   0.15 DLPFC_matchQSV  Schizo
# 61     FALSE 11126   -1   0.15           BSP1 Control
# 62      TRUE   339   -1   0.15           BSP1 Control
# 63     FALSE     0    0   0.15           BSP1    none
# 64      TRUE     0    0   0.15           BSP1    none
# 65     FALSE 12079    1   0.15           BSP1  Schizo
# 66      TRUE   578    1   0.15           BSP1  Schizo
# 67     FALSE 11491   -1   0.15            CMC Control
# 68      TRUE   853   -1   0.15            CMC Control
# 69     FALSE    21    0   0.15            CMC    none
# 70      TRUE     0    0   0.15            CMC    none
# 71     FALSE 11149    1   0.15            CMC  Schizo
# 72      TRUE   608    1   0.15            CMC  Schizo
# 73     FALSE 12681   -1   0.20 HIPPO_matchQSV Control
# 74      TRUE   171   -1   0.20 HIPPO_matchQSV Control
# 75     FALSE     0    0   0.20 HIPPO_matchQSV    none
# 76      TRUE     0    0   0.20 HIPPO_matchQSV    none
# 77     FALSE 11639    1   0.20 HIPPO_matchQSV  Schizo
# 78      TRUE   161    1   0.20 HIPPO_matchQSV  Schizo
# 79     FALSE 12286   -1   0.20 DLPFC_matchQSV Control
# 80      TRUE   921   -1   0.20 DLPFC_matchQSV Control
# 81     FALSE     0    0   0.20 DLPFC_matchQSV    none
# 82      TRUE     0    0   0.20 DLPFC_matchQSV    none
# 83     FALSE 10817    1   0.20 DLPFC_matchQSV  Schizo
# 84      TRUE   628    1   0.20 DLPFC_matchQSV  Schizo
# 85     FALSE 10949   -1   0.20           BSP1 Control
# 86      TRUE   516   -1   0.20           BSP1 Control
# 87     FALSE     0    0   0.20           BSP1    none
# 88      TRUE     0    0   0.20           BSP1    none
# 89     FALSE 11787    1   0.20           BSP1  Schizo
# 90      TRUE   870    1   0.20           BSP1  Schizo
# 91     FALSE 11132   -1   0.20            CMC Control
# 92      TRUE  1212   -1   0.20            CMC Control
# 93     FALSE    21    0   0.20            CMC    none
# 94      TRUE     0    0   0.20            CMC    none
# 95     FALSE 10849    1   0.20            CMC  Schizo
# 96      TRUE   908    1   0.20            CMC  Schizo

pdf('pdf/n_de_sign.pdf', useDingbats = FALSE, width = 10, height = 10)
ggplot(n_de_sign, aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
ggplot(subset(n_de_sign, de_status == 'TRUE'), aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
ggplot(subset(n_de_sign, de_status == 'TRUE' & sign != 0), aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
dev.off()

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

## DE genes by sign
de_genes_sign <- mapply(function(x, cut, sign) {
    x$ensemblID[x$adj.P.Val < cut & sign(x$logFC) == sign]
}, outGene[c(5, 5, 6, 6)], rep(c(0.2, 0.1), each = 2), sign = rep(c(-1, 1), 2))
## Sign -1 corresponds to control, +1 to schizo
# load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
# > colnames(modQsva)[grep('Dx', colnames(modQsva))]
# [1] "DxSchizo"
names(de_genes_sign) <- paste0(rep(c('HIPPO', 'DLPFC'), each = 2), '_', c('control', 'schizo'))

sapply(de_genes_sign, length)
# HIPPO_control  HIPPO_schizo DLPFC_control  DLPFC_schizo
#           171           161           379           253


## Top DE genes by sign
de_genes_sign_top <- mapply(function(x, sign) {
    x <- x[sign(x$logFC) == sign, ]
    head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], 400)
}, outGene[c(5, 5, 6, 6)], sign = rep(c(-1, 1), 2), SIMPLIFY = FALSE)
names(de_genes_sign_top) <- names(de_genes_sign)
sapply(de_genes_sign_top, length)

## Pretty venn code
venn_cols <- brewer.pal('Set1', n = 4)
names(venn_cols) <- names(de_genes_sign)
make_venn_pretty <- function(genes, title = 'DE features grouped by gene id') {
    v <- venn.diagram(genes, filename = NULL,
        main = title,
        col = 'transparent', fill = venn_cols[seq_len(length(genes))],
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = venn_cols[seq_len(length(genes))])
    grid.newpage()
    grid.draw(v)
}

pdf('pdf/venn_de_genes_by_sign.pdf', useDingbats = FALSE)
make_venn_pretty(de_genes_sign, 'DLPFC FDR10%, HIPPO FDR20%')
make_venn_pretty(list('HIPPO' = de_genes_sign[[1]], 'DLPFC' = de_genes_sign[[3]]), 'Control (DLPFC FDR10%, HIPPO FDR20%)')
make_venn_pretty(list('HIPPO' = de_genes_sign[[2]], 'DLPFC' = de_genes_sign[[4]]), 'Schizo (DLPFC FDR10%, HIPPO FDR20%)')
make_venn_pretty(de_genes_sign_top, 'top 400 in each group')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[1]], 'DLPFC' = de_genes_sign_top[[3]]), 'Control (top 400)')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[2]], 'DLPFC' = de_genes_sign_top[[4]]), 'Schizo (top 400)')
make_venn_pretty(lapply(de_genes_sign_top, head, n = 200), 'top 200 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 200), 'DLPFC' = head(de_genes_sign_top[[3]], 200)), 'Control (top 200)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 200), 'DLPFC' = head(de_genes_sign_top[[4]], 200)), 'Schizo (top 200)')
make_venn_pretty(lapply(de_genes_sign_top, head, n = 150), 'top 150 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 150), 'DLPFC' = head(de_genes_sign_top[[3]], 150)), 'Control (top 150)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 150), 'DLPFC' = head(de_genes_sign_top[[4]], 150)), 'Schizo (top 150)')
make_venn_pretty(lapply(de_genes_sign_top, head, n = 100), 'top 100 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 100), 'DLPFC' = head(de_genes_sign_top[[3]], 100)), 'Control (top 100)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 100), 'DLPFC' = head(de_genes_sign_top[[4]], 100)), 'Schizo (top 100)')
make_venn_pretty(lapply(de_genes_sign_top, head, n = 50), 'top 50 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 50), 'DLPFC' = head(de_genes_sign_top[[3]], 50)), 'Control (top 50)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 50), 'DLPFC' = head(de_genes_sign_top[[4]], 50)), 'Schizo (top 50)')
dev.off()
system('rm VennDiagram*.log')


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
uni <- outGene[[3]]$ensemblID


run_go <- function(genes, ont = c('BP', 'MF', 'CC')) {
    ## Change to ENSEMBL ids
    genes_ens <- lapply(genes, function(x) { gsub('\\..*', '', x) })

    #genes_venn <- venn(genes_ens, show.plot = FALSE)

    ## Run GO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running GO analysis for', bp))
        tryCatch(compareCluster(genes_ens, fun = "enrichGO",
            universe = uni, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
            readable = TRUE, keyType = 'ENSEMBL'),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    genes_ncbi <- lapply(lapply(genes_ens, bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(uni, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        universe = uni_ncbi, organism = 'hsa', pAdjustMethod = 'BH',
        pvalueCutoff = 0.1, qvalueCutoff = 0.05, keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })
    
    return(go_cluster)
}


## Development DE genes
if(!file.exists('rdas/go_de_genes.Rdata')) {
    system.time( go_de_genes <- run_go(de_genes_sign) )
    message(paste(Sys.time(), 'saving rda/go_de_genes.Rdata'))
    save(go_de_genes, file = 'rdas/go_de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rdas/go_de_genes.Rdata'))
    load('rdas/go_de_genes.Rdata', verbose = TRUE)
}
sapply(go_de_genes, class)

if(!file.exists('rdas/go_de_genes_top.Rdata')) {
    go_de_genes_top <- lapply(c(50, 100, 150, 200), function(n) {
        run_go(lapply(de_genes_sign_top, head, n = n))
    })
    names(go_de_genes_top) <- c(50, 100, 150, 200)
    message(paste(Sys.time(), 'saving rda/go_de_genes_top.Rdata'))
    save(go_de_genes_top, file = 'rdas/go_de_genes_top.Rdata')
} else {
    message(paste(Sys.time(), 'loading rdas/go_de_genes_top.Rdata'))
    load('rdas/go_de_genes_top.Rdata', verbose = TRUE)
}
lapply(go_de_genes_top, function(x) sapply(x, class))


simplify_go <- function(x) {
    #gsub('QSV|IPPO|LPFC', '', x)
    #gsub('_matchQSV', '', x)
    gsub('IPPO|LPFC|ontrol|chizo', '', x)
}

plot_go <- function(go_cluster, cat = 10) {
    lapply(names(go_cluster), function(bp) {
        go <- go_cluster[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        ## Simplify names
        go@compareClusterResult$Cluster <- simplify_go(go@compareClusterResult$Cluster)
        names(go@geneClusters) <- simplify_go(names(go@geneClusters))

        print(plot(go, title = paste('ontology:', bp), font.size = 18, showCategory = cat, includeAll = TRUE))
        return(NULL)
    })
}


pdf('pdf/go_de_genes.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_go(go_de_genes)
dev.off()

pdf('pdf/go_all_de_genes.pdf', width = 16, height = 70, useDingbats = FALSE)
plot_go(go_de_genes, cat = NULL)
dev.off()

for(i in names(go_de_genes_top)) {
    pdf(paste0('pdf/go_de_genes_top', i, '.pdf'), width = 14, height = 9, useDingbats = FALSE)
    plot_go(go_de_genes_top[[i]])
    dev.off()

    pdf(paste0('pdf/go_all_de_genes_top', i, '.pdf'), width = 16, height = 70, useDingbats = FALSE)
    plot_go(go_de_genes_top[[i]], cat = NULL)
    dev.off()
}

run_gse <-  function(region, ont = c('BP', 'MF', 'CC')) {

    genes <- outGene[[paste0(region, '_matchQSV')]]$t
    names(genes) <- outGene[[paste0(region, '_matchQSV')]]$ensemblID
    genes <- genes[order(genes, decreasing = TRUE)]

    ## Run gseGO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running gseGO analysis for', bp))
        tryCatch(gseGO(genes, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            keyType = 'ENSEMBL', verbose = TRUE),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    
    
    genes_tab <- bitr(names(genes), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
    genes_tab <- genes_tab[!duplicated(genes_tab$ENTREZID), ]
        
    genes_m <- match(genes_tab$ENSEMBL, names(genes))
    genes_ncbi <- genes[genes_m]
    names(genes_ncbi) <- genes_tab$ENTREZID

    go_cluster$KEGG <- tryCatch(gseKEGG(genes_ncbi, organism = 'hsa',
        pAdjustMethod = "BH",
        keyType = 'ncbi-geneid', verbose = TRUE),
        error = function(e) { return(NULL) })
        
    return(go_cluster)
}


system.time( gse_hippo <- run_gse('HIPPO') )
system.time( gse_dlpfc <- run_gse('DLPFC') )
save(gse_hippo, gse_dlpfc, file = 'rdas/gse.Rdata')


plot_gse <- function(gse) {
    lapply(names(gse), function(bp) {
        go <- gse[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        print(dotplot(go, title = paste('ontology:', bp), font.size = 18))
        return(NULL)
    })
}


pdf('pdf/gse_hippo.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_gse(gse_hippo)
dev.off()

pdf('pdf/gse_dlpfc.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_gse(gse_dlpfc)
dev.off()



## Scatter of logFC
comp_log <- function(x, y, xlab, ylab, var = 'logFC', de = FALSE, n = 150, onlyx = FALSE) {
    if(de) {
        if(!onlyx) {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n),
                head(y$ensemblID[order(y$adj.P.Val, decreasing = FALSE)], n)
            ))
        } else {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n)
            ))
        }
        
    } else {
        common <- intersect(x$ensemblID, y$ensemblID)
    }
    x <- x[match(common, x$ensemblID), ]
    y <- y[match(common, y$ensemblID), ]
    corr = signif(cor(x[, var], y[, var], use = 'pairwise.complete.obs'), 3)

    plot(x = x[, var], y = y[, var],
         xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlab),
         ylab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), ylab),
         col = add.alpha('black', ifelse(de, 1/2, 1/10)), pch = 16)
    legend('topleft', legend = paste('r =', corr))
    lines(loess.smooth(y = y[, var], x = x[, var]), col = 'red')
    abline(lm(y[, var] ~ x[, var]), col = 'blue')
    abline(h = 0, col = 'grey20')
    abline(v = 0, col = 'grey20')
}


pdf('pdf/scatter_models.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC')
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1')
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC')
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1')
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC')

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't')
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't')
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't')
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't')
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't')
dev.off()

pdf('pdf/scatter_models_top150de.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', de = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', de = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', de = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', de = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', de = TRUE)

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't', de = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't', de = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't', de = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't', de = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't', de = TRUE)
dev.off()

pdf('pdf/scatter_models_top400de.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', de = TRUE, n = 400)

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't', de = TRUE, n = 400)
dev.off()


pdf('pdf/scatter_models_top400de_onlyX.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', de = TRUE, n = 400, onlyx = TRUE)

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
dev.off()




## Load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
rse_gene_HIPPO <- rse_gene[, keepIndex]
mod_HIPPO <- mod
modQsva_HIPPO <- modQsva
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)
rse_gene_DLPFC <- rse_gene[, keepIndex]
mod_DLPFC <- mod
modQsva_DLPFC <- modQsva

## Get normalized expression and clean it
cleaned <- mapply(function(expr, model) {
    dge = DGEList(counts = assays(expr)$counts,
	genes = rowData(expr))
    #calculate library-size adjustment
    dge = calcNormFactors(dge)
    vGene = voom(dge, model, plot=FALSE)
    list('norm' = vGene$E, 'cleaned' = cleaningY(vGene$E, model, P = 2))
},
    list('HIPPO' = rse_gene_HIPPO, 'DLPFC' = rse_gene_DLPFC),
    list('HIPPO' = modQsva_HIPPO, 'DLPFC' = modQsva_DLPFC),
    SIMPLIFY = FALSE
)

## Adapted some functions from
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/region_specific/explore_reg_specific.R
# and
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/region_specific/explore_reg_specific_top.R
get_ylab <- function(type, cleaned = FALSE) {
    if(type %in% c('gene', 'exon', 'jxn')) {
        res <- 'log2(CPM + 0.5)'
    } else {
        res <- 'log2(TPM + 0.5)'
    }

    if(cleaned) res <- paste(res, '- covariate effects removed')
    return(res)
}

get_main <- function(i, region, group = 'schizo', type = 'gene') {
    if(type %in% c('gene', 'exon')) {
        var <- 'gencodeID'
        vars <- 'Symbol'
    } else if (type == 'jxn') {
        var <- 'gencodeGeneID'
        vars <- 'Symbol'
    } else {
        var <- 'gene_id'
        vars <- 'gene_name'
    }

    rse <- if(region == 'HIPPO') rse_gene_HIPPO else rse_gene_DLPFC

    topnow <- outGene[[paste0(region, '_matchQSV')]]

    j <- which(topnow$ensemblID == de_genes_sign_top[[paste0(region, '_', group)]][i])
    k <- which(names(rowRanges(rse)) == rownames(topnow)[j])

    paste(if(type != 'jxn') mcols(rowRanges(rse))[, var][k] else rownames(topnow)[j],
          if(is.na(mcols(rowRanges(rse))[, vars][k])) '' else mcols(rowRanges(rse))[, vars][k],
          'FDR',
          signif(topnow$adj.P.Val[j], 3),
          group)
}

get_ylim_mult <- function(rang) {
    c(
        ifelse(sign(rang[1]) == 1, 0.95, 1.05),
        ifelse(sign(rang[2]) == 1, 1.05, 0.95)
    )
}


plot_top_sign <- function(i, region, group = 'schizo', normtype = 'norm') {
    g <- de_genes_sign_top[[paste0(region, '_', group)]][i]
    j <- which(gsub('\\..*', '', rownames(cleaned[[region]][[normtype]])) == g)

    set.seed(20180426)
    dx <- if(region == 'HIPPO') colData(rse_gene_HIPPO)$Dx else colData(rse_gene_DLPFC)$Dx
    dx <- factor(dx, levels = c('Control', 'Schizo'))

    y <- cleaned[[region]][[normtype]][j, ]
    boxplot(y ~ dx,
            ylab = get_ylab('gene', cleaned = normtype == 'cleaned'),
            main = get_main(i, region, group),
            col = c('orchid1', 'aquamarine1'),
            ylim = abs(range(y)) * get_ylim_mult(range(y)) * sign(range(y)),
            outline = FALSE)
    points(y ~ jitter(as.integer(dx), 1), pch = 21,
           bg = c('orchid4', 'aquamarine4')[as.integer(dx)])
}

for(reg in c('HIPPO', 'DLPFC')) {
    pdf(paste0('pdf/top_50each_', reg, '_norm.pdf'))
    for(i in 1:50) {
        plot_top_sign(i, reg)
        plot_top_sign(i, reg, 'control')
    }
    dev.off()
    pdf(paste0('pdf/top_50each_', reg, '_cleaned.pdf'))
    for(i in 1:50) {
        plot_top_sign(i, reg, normtype = 'cleaned')
        plot_top_sign(i, reg, 'control', normtype = 'cleaned')
    }
    dev.off()
}


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
