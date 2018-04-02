
library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)

load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda")
getRPKM = recount::getRPKM
geneRpkm = getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)
pd$Dx = factor(pd$Dx)
pd$Dx = relevel(pd$Dx, "Control")

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))
pcaVars = getPcaVars(pca)

## check output
plot(pca$x[,1] ~ pd$totalAssignedGene)
plot(pca$x[,2] ~ pd$RIN)
plot(pca$x[,3] ~ pd$mitoRate)
plot(pca$x[,3] ~ factor(pd$BatchLab))


## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse[rowData(cov_rse)$bonfSig,])$counts+1)))
getPcaVars(qsvBonf)[1:5]

pdf("qcChecks/qSVs_bonfRegions.pdf")
mypar(2,2)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=19,bg="grey",
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,2] ~ factor(pd$BatchLab),
     xlab = "Batch",	ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
plot(qsvBonf$x[,3] ~ factor(pd$BatchLab),
     xlab = "Batch",	ylab=paste0("qSV3: ",getPcaVars(qsvBonf)[3],"% Var Expl"))
plot(qsvBonf$x[,4] ~ factor(pd$BatchLab),
     xlab = "Batch",	ylab=paste0("qSV4: ",getPcaVars(qsvBonf)[4],"% Var Expl"))
dev.off()

summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene))
summary(lm(qsvBonf$x[,1] ~ pd$Dx))
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$totalAssignedGene))


## get qSVs for top 1000
qsvTop = prcomp(t(log2(assays(cov_rse)$counts+1)))
getPcaVars(qsvTop)[1:5]
plot(qsvTop$x[,1] ~ pd$totalAssignedGene)
plot(qsvTop$x[,1] ~ factor(pd$Dx))
plot(qsvTop$x[,2] ~ factor(pd$Batch))
plot(qsvTop$x[,3] ~ factor(pd$Batch))
plot(qsvTop$x[,4] ~ factor(pd$Batch))
plot(qsvTop$x[,4] ~ pd$overallMapRate)

cc = cor(qsvBonf$x[,1:10], qsvTop$x[,1:10])
signif(cc,2)
