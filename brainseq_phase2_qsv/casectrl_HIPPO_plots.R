
## Load degradeStats, degradeStatsInt
## degradeStats from main model, mod = model.matrix(~pd$DegradationTime + pd$Region + factor(pd$BrNum)
## degradeInt from interaction model, modInt = model.matrix(~pd$DegradationTime*pd$Region + factor(pd$BrNum))
load("../ERs/rdas/DLPFC_Plus_HIPPO_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda", verbose = TRUE)
## Load outGene, outGene0
load("rdas/dxStats_hippo_filtered_qSVA_geneLevel.rda", verbose = TRUE)

plot_dequal <- function(out_input, degrade_input, var = 't', xlabtxt = 'case-control') {
	both <- intersect(rownames(out_input), rownames(degrade_input))
	degrade <- degrade_input[both, ]
	interest <- out_input[both, ]
	
	stopifnot(identical(rownames(degrade), rownames(interest)))
	plot(y = degrade[, var], x = interest[, var], xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlabtxt), ylab = paste(var, 'degradation'), main = signif(cor( degrade[, var],  interest[, var]), 3))
}

pdf('pdf/HIPPO_dequal_mainmodel.pdf', useDingbats = FALSE)
plot_dequal(out_input = outGene0, degrade_input = degradeStats, var = 't', xlabtxt = 'case-control (without qSVs)')
plot_dequal(out_input = outGene, degrade_input = degradeStats, var = 't', xlabtxt = 'case-control (with qSVs)')

plot_dequal(out_input = outGene0, degrade_input = degradeStats, var = 'logFC', xlabtxt = 'case-control (without qSVs)')
plot_dequal(out_input = outGene, degrade_input = degradeStats, var = 'logFC', xlabtxt = 'case-control (with qSVs)')
dev.off()


pdf('pdf/HIPPO_dequal_intmodel.pdf', useDingbats = FALSE)
plot_dequal(out_input = outGene0, degrade_input = degradeStatsInt, var = 't', xlabtxt = 'case-control (without qSVs)')
plot_dequal(out_input = outGene, degrade_input = degradeStatsInt, var = 't', xlabtxt = 'case-control (with qSVs)')

plot_dequal(out_input = outGene0, degrade_input = degradeStatsInt, var = 'logFC', xlabtxt = 'case-control (without qSVs)')
plot_dequal(out_input = outGene, degrade_input = degradeStatsInt, var = 'logFC', xlabtxt = 'case-control (with qSVs)')
dev.off()




