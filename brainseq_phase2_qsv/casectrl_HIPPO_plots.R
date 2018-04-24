
## Load degradeStats, degradeStatsInt
## degradeStats from main model, mod = model.matrix(~pd$DegradationTime + pd$Region + factor(pd$BrNum)
## degradeInt from interaction model, modInt = model.matrix(~pd$DegradationTime*pd$Region + factor(pd$BrNum))
load("../ERs/rdas/DLPFC_Plus_HIPPO_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda", verbose = TRUE)
## Load outGene, outGene0, outGeneNoAdj
load("rdas/dxStats_hippo_filtered_qSVA_geneLevel.rda", verbose = TRUE)
## load degradeStats_DLPFC, degradeStats_HIPPO from model.matrix( ~ DegradationTime + factor(BrNum)))
load("/dcl01/ajaffe/data/lab/qsva_brain/ERs/rdas/DLPFC_HIPPO_degradationStats_hg38.rda", verbose = TRUE)

plot_dequal <- function(out_input, degrade_input, var = 't', xlabtxt = 'case-control', main) {
	both <- intersect(rownames(out_input), rownames(degrade_input))
	degrade <- degrade_input[both, ]
	interest <- out_input[both, ]
	
	stopifnot(identical(rownames(degrade), rownames(interest)))
	corr = signif(cor( degrade[, var],  interest[, var]), 3)
	plot(y = degrade[, var], x = interest[, var], xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlabtxt), ylab = paste(var, 'degradation'), main = main)
	legend('topleft', legend = paste('r =', corr))
}

pdf('pdf/HIPPO_dequal_mainmodel.pdf', useDingbats = FALSE)
plot_dequal(out_input = outGene0, degrade_input = degradeStats, var = 't', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Main Model")
plot_dequal(out_input = outGene, degrade_input = degradeStats, var = 't', xlabtxt = 'case-control (with qSVs)', main = "HIPPO Main Model")

plot_dequal(out_input = outGene0, degrade_input = degradeStats, var = 'logFC', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Main Model")
plot_dequal(out_input = outGene, degrade_input = degradeStats, var = 'logFC', xlabtxt = 'case-control (with qSVs)', main = "HIPPO Main Model")
dev.off()


pdf('pdf/HIPPO_dequal_region.pdf', useDingbats = FALSE)
plot_dequal(out_input = outGene0, degrade_input = degradeStats_HIPPO, var = 't', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Region-specific Main Model")
plot_dequal(out_input = outGene, degrade_input = degradeStats_HIPPO, var = 't', xlabtxt = 'case-control (with qSVs)', main = "HIPPO Region-specific Main Model")

plot_dequal(out_input = outGene0, degrade_input = degradeStats_HIPPO, var = 'logFC', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Region-specific Main Model")
plot_dequal(out_input = outGene, degrade_input = degradeStats_HIPPO, var = 'logFC', xlabtxt = 'case-control (with qSVs)', main = "HIPPO Region-specific Main Model")
dev.off()


pdf('pdf/HIPPO_dequal_dx.pdf', useDingbats = FALSE)
plot_dequal(out_input = outGeneNoAdj, degrade_input = degradeStats_HIPPO, var = 't', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Dx Model")
plot_dequal(out_input = outGeneNoAdj, degrade_input = degradeStats_HIPPO, var = 'logFC', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Dx Model")
dev.off()


pdf('pdf/HIPPO_dequal_intmodel.pdf', useDingbats = FALSE)
plot_dequal(out_input = outGene0, degrade_input = degradeStatsInt, var = 't', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Interaction Model")
plot_dequal(out_input = outGene, degrade_input = degradeStatsInt, var = 't', xlabtxt = 'case-control (with qSVs)', main = "HIPPO Interaction Model")

plot_dequal(out_input = outGene0, degrade_input = degradeStatsInt, var = 'logFC', xlabtxt = 'case-control (without qSVs)', main = "HIPPO Interaction Model")
plot_dequal(out_input = outGene, degrade_input = degradeStatsInt, var = 'logFC', xlabtxt = 'case-control (with qSVs)', main = "HIPPO Interaction Model")
dev.off()



