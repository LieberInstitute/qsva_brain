

library(jaffelab)
library(SummarizedExperiment)
library(sva)

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
pdf('pdf/qsvs_var_explained.pdf')
plot(getPcaVars(qsvBonf), pch=20)
dev.off()
