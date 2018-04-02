
library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)

# gene-level expression for degradation data
load("/dcl01/ajaffe/data/lab/qsva_brain/expr_data/rda/rse_gene.Rdata", verbose = TRUE)
rse_gene <- rse_gene[, colData(rse_gene)$Region %in% c('HIPPO', 'DLPFC')]


getRPKM = recount::getRPKM
geneRpkm = getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))
# percent of variance explained by pcas
pcaVars = getPcaVars(pca)

## check output
#region
dir.create('pdf', showWarnings = FALSE)
pdf('pdf/name.pdf')
plot(pca$x[,1] ~ pd$totalAssignedGene)
plot(pca$x[,2] ~ pd$RIN)
plot(pca$x[,3] ~ pd$mitoRate)
plot(pca$x[,3] ~ factor(pd$BatchLab))
dev.off()



load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)


colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)
colData(rse_gene)$Kit = ifelse(colData(rse_gene)$mitoRate < 0.05, "Gold", "HMR")

geneRpkm = getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)
pd$Dx = factor(pd$Dx)
pd$Dx = relevel(pd$Dx, "Control")



# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))
# percent of variance explained by pcas
pcaVars = getPcaVars(pca)

## check output
#region
pdf('pdf/brainseqplots.pdf')
plot(pca$x[,1] ~ pd$totalAssignedGene)
plot(pca$x[,2] ~ pd$RIN)
plot(pca$x[,3] ~ pd$mitoRate)
plot(pca$x[,3] ~ factor(pd$BatchLab))
plot(pca$x[,1] ~ pd$Dx)
plot(pca$x[,1] ~ factor(pd$Kit))
plot(pca$x[,1] ~ pd$totalAssignedGene, col = c('orange', 'light blue')[factor(pd$Kit)], pch = 21)
dev.off()


# top 1000 expressed regions associated with degradation
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))
getPcaVars(qsvBonf)[1:5]

pdf("pdf/qSVs_brainseq.pdf")
mypar(2,2)
#1st qsv across all 900 samples
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=19,bg="grey",
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))

plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=19,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'light blue')[c('Gold' = 1, 'HMR' = 2)[pd$Kit]])


plot(qsvBonf$x[,1] ~ pd$Dx,
     xlab = "Batch",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))


plot(qsvBonf$x[,2] ~ factor(pd$BatchLab),
     xlab = "Batch",	ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))

plot(qsvBonf$x[,2] ~ pd$Dx,
     xlab = "Batch",	ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))

plot(qsvBonf$x[,3] ~ factor(pd$BatchLab),
     xlab = "Batch",	ylab=paste0("qSV3: ",getPcaVars(qsvBonf)[3],"% Var Expl"))

plot(qsvBonf$x[,3] ~ pd$Dx,
     xlab = "Batch",	ylab=paste0("qSV3: ",getPcaVars(qsvBonf)[3],"% Var Expl"))

plot(qsvBonf$x[,4] ~ factor(pd$BatchLab),
     xlab = "Batch",	ylab=paste0("qSV4: ",getPcaVars(qsvBonf)[4],"% Var Expl"))
dev.off()

print('Regression: qsv1 vs totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene))
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene + pd$Kit))
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Kit))

# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx)
# 
# Residuals:
#    Min      1Q  Median      3Q     Max
# -56.723 -10.200   1.171  12.296  39.468
# 
# Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
# (Intercept)   2.2592     0.6896   3.276  0.00109 **
# pd$DxSchizo  -7.1092     1.2233  -5.812 8.59e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 17.09 on 898 degrees of freedom
# Multiple R-squared:  0.03625,	Adjusted R-squared:  0.03518
# F-statistic: 33.78 on 1 and 898 DF,  p-value: 8.586e-09



summary(lm(qsvBonf$x[,1] ~ pd$Dx))
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$totalAssignedGene))


