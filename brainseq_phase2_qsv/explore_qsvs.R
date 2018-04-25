
library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)
library('readxl')
library('devtools')

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
getPcaVars(pca)[1:5]
# [1] 26.60 12.90  6.10  5.81  3.82

## degradation plots
dir.create('pdf', showWarnings = FALSE)
pdf('pdf/degradation.pdf', useDingbats = FALSE)
plot(pca$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,1] ~ factor(pd$Region),
     xlab = "Region",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))
plot(pca$x[,2] ~ factor(pd$Region),
     xlab = "Region",
     ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"))
plot(pca$x[,1] ~ pd$RIN,
     xlab = "RIN", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
dev.off()


load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)
# Kit definition previously used only applicable to HIPPO not DLPFC
# colData(rse_gene)$Kit = ifelse(colData(rse_gene)$mitoRate < 0.05, "Gold", "HMR")


## Properly add kit
hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')
hipHMR <- as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
colData(rse_gene)$kit <- ifelse(hipHMR, 'HMR', 'Gold')
colData(rse_gene)$regkit <- with(colData(rse_gene), paste0(Region, '_', kit))

with(colData(rse_gene), addmargins(table(kit, Region)))



geneRpkm = getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)
pd$Dx = factor(pd$Dx)
pd$Dx = relevel(pd$Dx, "Control")
pd$Sex = factor(pd$Sex)
pd$Sex = relevel(pd$Sex, "M")

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))
# percent of variance explained by pcas
pcaVars = getPcaVars(pca)
getPcaVars(pca)[1:5]
# [1] 33.50 17.90  6.17  4.65  2.75

## brainseq plots
pdf('pdf/brainseqplots.pdf', useDingbats = FALSE)
plot(pca$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,1] ~ pd$RIN,
     xlab = "RIN", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(pca$x[,1] ~ factor(pd$Region),
     xlab = "Region",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))
plot(pca$x[,1] ~ pd$Dx,
     xlab = "Dx",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))
plot(pca$x[,1] ~ pd$Sex,
     xlab = "Sex",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))

plot(pca$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('topleft', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('dark orange', 'blue', 'skyblue3'), bty = 'n')

plot(pca$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('topleft', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('orange', 'blue', 'skyblue3'), bty = 'n')

plot(pca$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('bottomright', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('orange', 'blue', 'skyblue3'), bty = 'n')

plot(pca$x[,1] ~ factor(pd$regkit),
     xlab = "Region and kit",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))

plot(pca$x[,1] ~ pd$Age,
     xlab = "Age", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('topleft', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('dark orange', 'blue', 'skyblue3'), bty = 'n')

plot(pca$x[,2] ~ pd$Age,
     xlab = "Age", pch=20,
     ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('bottomright', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('dark orange', 'blue', 'skyblue3'), bty = 'n')

plot(pca$x[,2] ~ factor(ifelse(pd$Age < 0, 'Prenatal', 'Postnatal')),
     xlab = 'Age group',
     ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"))


dev.off()


# top 1000 expressed regions associated with degradation
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))
getPcaVars(qsvBonf)[1:5]
# [1] 57.80 15.70  3.36  2.70  2.36

## brainseq qsv plots
pdf("pdf/qSVs_brainseq.pdf", useDingbats = FALSE)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Dx,
     xlab = "Diagnosis",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ factor(pd$Region),
     xlab = "Region",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$RIN,
     xlab = "RIN",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
# exploring differences in prenatal samples
boxplot(qsvBonf$x[,2] ~ pd$Age<0,
	xlab = "Age<0",
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
boxplot(qsvBonf$x[,2] ~ factor(ifelse(pd$Age<0, 'Prenatal', 'Postnatal')) + pd$Region,
	xlab = "Age and Region", cex.axis=0.9,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ factor(pd$Sex),
     xlab = "Sex",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
boxplot(qsvBonf$x[,1] ~ pd$Region + pd$Dx,
	xlab = "Region and Diagnosis",
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))

plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('topleft', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('orange', 'blue', 'skyblue3'), bty = 'n')

plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('topleft', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('orange', 'blue', 'skyblue3'), bty = 'n')

plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('bottomright', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('orange', 'blue', 'skyblue3'), bty = 'n')

plot(qsvBonf$x[,1] ~ factor(pd$regkit),
     xlab = "Region and kit",
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))

plot(qsvBonf$x[,1] ~ pd$Age,
     xlab = "Age", pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('topleft', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('dark orange', 'blue', 'skyblue3'), bty = 'n')

plot(qsvBonf$x[,2] ~ pd$Age,
     xlab = "Age", pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('orange', 'blue', 'skyblue3')[factor(pd$regkit)])
legend('topright', c('DLPFC_Gold', 'HIPPO_Gold', 'HIPPO_HMR'), lwd = 2, col = c('dark orange', 'blue', 'skyblue3'), bty = 'n')

plot(qsvBonf$x[,2] ~ factor(ifelse(pd$Age < 0, 'Prenatal', 'Postnatal')),
     xlab = 'Age group',
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))

dev.off()


print('Regression: qsv1 vs totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene)

# Residuals:
#    Min      1Q  Median      3Q     Max
# -53.972 -10.494   2.319  12.299  44.369

# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)             7.011      3.581   1.958   0.0506 .
# pd$totalAssignedGene  -15.464      7.796  -1.984   0.0476 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 17.37 on 898 degrees of freedom
# Multiple R-squared:  0.004363,	Adjusted R-squared:  0.003254
# F-statistic: 3.935 on 1 and 898 DF,  p-value: 0.0476

print('Regression: qsv1 vs Region')
summary(lm(qsvBonf$x[,1] ~ pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Region)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -49.716  -7.472   2.994  10.634  32.709
#
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)      9.0177     0.6974   12.93   <2e-16 ***
# pd$RegionHIPPO -18.1564     0.9896  -18.35   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.84 on 898 degrees of freedom
# Multiple R-squared:  0.2726,	Adjusted R-squared:  0.2718
# F-statistic: 336.6 on 1 and 898 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssignedGene + Region')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene + pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene + pd$Region)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -59.545  -6.375   2.574   8.577  33.261
#
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)           -29.994      3.255  -9.214   <2e-16 ***
# pd$totalAssignedGene   95.604      7.819  12.227   <2e-16 ***
# pd$RegionHIPPO        -26.873      1.161 -23.140   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 13.75 on 897 degrees of freedom
# Multiple R-squared:  0.3765,	Adjusted R-squared:  0.3752
# F-statistic: 270.9 on 2 and 897 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssignedGene*Region')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Region)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -46.857  -5.155   1.868   7.036  39.507
#
# Coefficients:
#                                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)                          -79.199      4.367  -18.14   <2e-16 ***
# pd$totalAssignedGene                 216.190     10.608   20.38   <2e-16 ***
# pd$RegionHIPPO                        68.409      6.393   10.70   <2e-16 ***
# pd$totalAssignedGene:pd$RegionHIPPO -212.882     14.094  -15.10   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 12.28 on 896 degrees of freedom
# Multiple R-squared:  0.5031,	Adjusted R-squared:  0.5014
# F-statistic: 302.4 on 3 and 896 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx')
summary(lm(qsvBonf$x[,1] ~ pd$Dx))
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

print('Regression: qsv1 vs Region')
summary(lm(qsvBonf$x[,1] ~ pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Region)
#
# Residuals:
#    Min      1Q  Median      3Q     Max
# -49.716  -7.472   2.994  10.634  32.709
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)      9.0177     0.6974   12.93   <2e-16 ***
# pd$RegionHIPPO -18.1564     0.9896  -18.35   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.84 on 898 degrees of freedom
# Multiple R-squared:  0.2726,	Adjusted R-squared:  0.2718
# F-statistic: 336.6 on 1 and 898 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs rRNA_rate')
summary(lm(qsvBonf$x[,1] ~ pd$rRNA_rate))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$rRNA_rate)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -54.675 -10.304   2.222  12.500  42.473
#
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -1.4561     0.9547  -1.525   0.1276
# pd$rRNA_rate 36589.7405 19076.4582   1.918   0.0554 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 17.37 on 898 degrees of freedom
# Multiple R-squared:  0.00408,	Adjusted R-squared:  0.002971
# F-statistic: 3.679 on 1 and 898 DF,  p-value: 0.05542

print('Regression: qsv1 vs RIN')
summary(lm(qsvBonf$x[,1] ~ pd$RIN))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$RIN)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -49.443 -10.055  -1.155  10.203  36.338
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -80.7783     3.3846  -23.87   <2e-16 ***
# pd$RIN       10.4206     0.4327   24.08   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 13.57 on 898 degrees of freedom
# Multiple R-squared:  0.3924,    Adjusted R-squared:  0.3917
# F-statistic:   580 on 1 and 898 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$totalAssignedGene)
#
# Residuals:
#    Min      1Q  Median      3Q     Max
# -54.303 -10.640   1.557  11.930  42.312
#
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)             9.898      3.550   2.788  0.00541 **
# pd$DxSchizo            -7.188      1.221  -5.886 5.58e-09 ***
# pd$totalAssignedGene  -16.795      7.657  -2.193  0.02854 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 17.05 on 897 degrees of freedom
# Multiple R-squared:  0.04139,	Adjusted R-squared:  0.03925
# F-statistic: 19.36 on 2 and 897 DF,  p-value: 5.84e-09

print('Regression: qsv1 vs totalAssignedGene*Region*Dx')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Region * pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Region *
#     pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -48.162  -5.243   1.759   6.700  45.791
#
# Coefficients:
#                                                 Estimate Std. Error t value
# (Intercept)                                      -72.657      4.969 -14.621
# pd$totalAssignedGene                             201.743     11.908  16.941
# pd$RegionHIPPO                                    68.190      6.852   9.952
# pd$DxSchizo                                      -19.485      9.128  -2.135
# pd$totalAssignedGene:pd$RegionHIPPO             -204.575     15.203 -13.456
# pd$totalAssignedGene:pd$DxSchizo                  44.212     22.528   1.963
# pd$RegionHIPPO:pd$DxSchizo                       -54.842     16.247  -3.376
# pd$totalAssignedGene:pd$RegionHIPPO:pd$DxSchizo   79.862     34.640   2.305
#                                                 Pr(>|t|)
# (Intercept)                                      < 2e-16 ***
# pd$totalAssignedGene                             < 2e-16 ***
# pd$RegionHIPPO                                   < 2e-16 ***
# pd$DxSchizo                                     0.033060 *
# pd$totalAssignedGene:pd$RegionHIPPO              < 2e-16 ***
# pd$totalAssignedGene:pd$DxSchizo                0.050012 .
# pd$RegionHIPPO:pd$DxSchizo                      0.000768 ***
# pd$totalAssignedGene:pd$RegionHIPPO:pd$DxSchizo 0.021370 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 11.59 on 892 degrees of freedom
# Multiple R-squared:  0.5596,	Adjusted R-squared:  0.5561
# F-statistic: 161.9 on 7 and 892 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$Region + pd$rRNA_rate + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$Region + pd$rRNA_rate + pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -37.703  -4.870   0.716   5.248  47.000
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -7.520e+01  2.958e+00 -25.422  < 2e-16 ***
# pd$DxSchizo          -2.256e+00  6.949e-01  -3.246 0.001213 **
# pd$Age               -8.578e-02  1.500e-02  -5.719 1.47e-08 ***
# pd$SexF              -2.317e+00  6.340e-01  -3.654 0.000273 ***
# pd$mitoRate          -1.920e+02  6.006e+00 -31.970  < 2e-16 ***
# pd$RegionHIPPO       -2.681e+00  1.086e+00  -2.470 0.013710 *
# pd$rRNA_rate         -4.073e+04  1.066e+04  -3.819 0.000143 ***
# pd$totalAssignedGene  2.334e+02  6.859e+00  34.033  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 8.891 on 892 degrees of freedom
# Multiple R-squared:  0.7408,	Adjusted R-squared:  0.7388
# F-statistic: 364.3 on 7 and 892 DF,  p-value: < 2.2e-16


print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$Region + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$Region + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -32.658  -4.111  -0.153   4.211  29.812
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -1.114e+02  2.795e+00 -39.864  < 2e-16 ***
# pd$DxSchizo          -2.649e+00  5.473e-01  -4.840 1.53e-06 ***
# pd$Age                9.667e-03  1.249e-02   0.774   0.4391
# pd$SexF              -1.637e+00  4.999e-01  -3.275   0.0011 **
# pd$mitoRate          -1.292e+02  5.435e+00 -23.773  < 2e-16 ***
# pd$RegionHIPPO       -8.811e+00  8.936e-01  -9.859  < 2e-16 ***
# pd$rRNA_rate          2.209e+03  8.591e+03   0.257   0.7972
# pd$totalAssignedGene  1.791e+02  5.876e+00  30.474  < 2e-16 ***
# pd$RIN                6.555e+00  2.798e-01  23.425  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 6.998 on 891 degrees of freedom
# Multiple R-squared:  0.8396,    Adjusted R-squared:  0.8382
# F-statistic:   583 on 8 and 891 DF,  p-value: < 2.2e-16


## Save for later
pd_all <- pd


## Explore qSVs with no Hippo Gold samples
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold.Rdata', verbose = TRUE)
pd <- pd_all[keepIndex, ]

getPcaVars(qsvBonf)[1:5]
# [1] 70.20  5.54  3.07  2.60  1.18

## brainseq qsv plots
pdf("pdf/qSVs_brainseq_age17_noHGold.pdf", useDingbats = FALSE)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Dx,
     xlab = "Diagnosis",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ factor(pd$Region),
     xlab = "Region",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$RIN,
     xlab = "RIN",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
# exploring differences in prenatal samples
boxplot(qsvBonf$x[,2] ~ pd$Age<0,
	xlab = "Age<0",
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
boxplot(qsvBonf$x[,2] ~ factor(ifelse(pd$Age<0, 'Prenatal', 'Postnatal')) + pd$Region,
	xlab = "Age and Region", cex.axis=0.9,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('orange', 'skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ factor(pd$Sex),
     xlab = "Sex",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
boxplot(qsvBonf$x[,1] ~ pd$Region + pd$Dx,
	xlab = "Region and Diagnosis",
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
dev.off()

print('Regression: qsv1 vs totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -53.361 -10.384   2.836  12.973  45.066
#
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)            14.127      4.361   3.239  0.00125 **
# pd$totalAssignedGene  -30.559      9.321  -3.279  0.00109 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 17.92 on 710 degrees of freedom
# Multiple R-squared:  0.01491,   Adjusted R-squared:  0.01353
# F-statistic: 10.75 on 1 and 710 DF,  p-value: 0.001094

print('Regression: qsv1 vs Region')
summary(lm(qsvBonf$x[,1] ~ pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Region)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -49.476  -7.359   3.155  10.979  34.232
#
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)      9.5688     0.7647   12.51   <2e-16 ***
# pd$RegionHIPPO -20.4594     1.1182  -18.30   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.89 on 710 degrees of freedom
# Multiple R-squared:  0.3204,    Adjusted R-squared:  0.3195
# F-statistic: 334.8 on 1 and 710 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssignedGene + Region')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene + pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene + pd$Region)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -67.041  -4.424   2.791   6.923  34.204
#
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)           -67.301      4.052  -16.61   <2e-16 ***
# pd$totalAssignedGene  187.197      9.750   19.20   <2e-16 ***
# pd$RegionHIPPO        -41.134      1.408  -29.21   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 12.08 on 709 degrees of freedom
# Multiple R-squared:  0.5529,    Adjusted R-squared:  0.5516
# F-statistic: 438.4 on 2 and 709 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssignedGene*Region')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Region)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -55.648  -4.409   1.792   6.431  45.016
#
# Coefficients:
#                                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)                          -91.050      4.971 -18.314  < 2e-16 ***
# pd$totalAssignedGene                 245.030     12.019  20.387  < 2e-16 ***
# pd$RegionHIPPO                        29.407      9.275   3.171  0.00159 **
# pd$totalAssignedGene:pd$RegionHIPPO -147.633     19.203  -7.688 5.01e-14 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 11.62 on 708 degrees of freedom
# Multiple R-squared:  0.5873,    Adjusted R-squared:  0.5856
# F-statistic: 335.9 on 3 and 708 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssignedGene*Region*Dx')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Region * pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Region *
#     pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -53.450  -4.477   1.651   6.439  45.227
#
# Coefficients:
#                                                  Estimate Std. Error t value
# (Intercept)                                      -90.4169     6.6124 -13.674
# pd$totalAssignedGene                             243.7794    15.6559  15.571
# pd$RegionHIPPO                                    54.9894    12.0246   4.573
# pd$DxSchizo                                       -0.9119     9.9771  -0.091
# pd$totalAssignedGene:pd$RegionHIPPO             -190.2316    24.5915  -7.736
# pd$totalAssignedGene:pd$DxSchizo                   1.5528    24.3581   0.064
# pd$RegionHIPPO:pd$DxSchizo                       -41.6296    18.7407  -2.221
# pd$totalAssignedGene:pd$RegionHIPPO:pd$DxSchizo   65.2321    39.0797   1.669
#                                                 Pr(>|t|)
# (Intercept)                                      < 2e-16 ***
# pd$totalAssignedGene                             < 2e-16 ***
# pd$RegionHIPPO                                  5.68e-06 ***
# pd$DxSchizo                                       0.9272
# pd$totalAssignedGene:pd$RegionHIPPO             3.57e-14 ***
# pd$totalAssignedGene:pd$DxSchizo                  0.9492
# pd$RegionHIPPO:pd$DxSchizo                        0.0266 *
# pd$totalAssignedGene:pd$RegionHIPPO:pd$DxSchizo   0.0955 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 11.31 on 704 degrees of freedom
# Multiple R-squared:  0.6112,    Adjusted R-squared:  0.6073
# F-statistic: 158.1 on 7 and 704 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx')
summary(lm(qsvBonf$x[,1] ~ pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -57.120 -11.073   1.061  13.243  38.469
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   2.9401     0.8575   3.429 0.000641 ***
# pd$DxSchizo  -7.3194     1.3529  -5.410 8.62e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 17.7 on 710 degrees of freedom
# Multiple R-squared:  0.03959,   Adjusted R-squared:  0.03824
# F-statistic: 29.27 on 1 and 710 DF,  p-value: 8.621e-08

print('Regression: qsv1 vs Region')
summary(lm(qsvBonf$x[,1] ~ pd$Region))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Region)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -49.476  -7.359   3.155  10.979  34.232
#
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)      9.5688     0.7647   12.51   <2e-16 ***
# pd$RegionHIPPO -20.4594     1.1182  -18.30   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.89 on 710 degrees of freedom
# Multiple R-squared:  0.3204,    Adjusted R-squared:  0.3195
# F-statistic: 334.8 on 1 and 710 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs rRNA_rate')
summary(lm(qsvBonf$x[,1] ~ pd$rRNA_rate))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$rRNA_rate)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -54.492 -11.050   2.379  13.591  41.741
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -1.545      1.092  -1.415   0.1576
# pd$rRNA_rate 38085.535  21162.668   1.800   0.0723 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 18.02 on 710 degrees of freedom
# Multiple R-squared:  0.004541,  Adjusted R-squared:  0.003139
# F-statistic: 3.239 on 1 and 710 DF,  p-value: 0.07234

print('Regression: qsv1 vs Dx + totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -53.394 -10.524   1.658  12.417  42.726
#
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)            20.901      4.405   4.745 2.52e-06 ***
# pd$DxSchizo            -8.096      1.351  -5.994 3.26e-09 ***
# pd$totalAssignedGene  -38.177      9.188  -4.155 3.65e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 17.5 on 709 degrees of freedom
# Multiple R-squared:  0.06242,   Adjusted R-squared:  0.05978
# F-statistic:  23.6 on 2 and 709 DF,  p-value: 1.193e-10

print('Regression: qsv1 vs totalAssignedGene*Region*Dx')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Region * pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Region *
#     pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -53.450  -4.477   1.651   6.439  45.227
#
# Coefficients:
#                                                  Estimate Std. Error t value
# (Intercept)                                      -90.4169     6.6124 -13.674
# pd$totalAssignedGene                             243.7794    15.6559  15.571
# pd$RegionHIPPO                                    54.9894    12.0246   4.573
# pd$DxSchizo                                       -0.9119     9.9771  -0.091
# pd$totalAssignedGene:pd$RegionHIPPO             -190.2316    24.5915  -7.736
# pd$totalAssignedGene:pd$DxSchizo                   1.5528    24.3581   0.064
# pd$RegionHIPPO:pd$DxSchizo                       -41.6296    18.7407  -2.221
# pd$totalAssignedGene:pd$RegionHIPPO:pd$DxSchizo   65.2321    39.0797   1.669
#                                                 Pr(>|t|)
# (Intercept)                                      < 2e-16 ***
# pd$totalAssignedGene                             < 2e-16 ***
# pd$RegionHIPPO                                  5.68e-06 ***
# pd$DxSchizo                                       0.9272
# pd$totalAssignedGene:pd$RegionHIPPO             3.57e-14 ***
# pd$totalAssignedGene:pd$DxSchizo                  0.9492
# pd$RegionHIPPO:pd$DxSchizo                        0.0266 *
# pd$totalAssignedGene:pd$RegionHIPPO:pd$DxSchizo   0.0955 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 11.31 on 704 degrees of freedom
# Multiple R-squared:  0.6112,    Adjusted R-squared:  0.6073
# F-statistic: 158.1 on 7 and 704 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$Region + pd$rRNA_rate + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$Region + pd$rRNA_rate + pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -30.716  -4.728   0.301   4.898  49.398
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -8.045e+01  3.786e+00 -21.252  < 2e-16 ***
# pd$DxSchizo          -1.952e+00  6.917e-01  -2.822 0.004900 **
# pd$Age               -1.062e-01  2.248e-02  -4.725 2.78e-06 ***
# pd$SexF              -2.772e+00  7.079e-01  -3.916 9.88e-05 ***
# pd$mitoRate          -2.025e+02  8.560e+00 -23.658  < 2e-16 ***
# pd$RegionHIPPO       -2.618e+00  1.862e+00  -1.406 0.160263
# pd$rRNA_rate         -4.050e+04  1.101e+04  -3.679 0.000252 ***
# pd$totalAssignedGene  2.510e+02  8.170e+00  30.725  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 8.583 on 704 degrees of freedom
# Multiple R-squared:  0.776,     Adjusted R-squared:  0.7738
# F-statistic: 348.5 on 7 and 704 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$Region + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$Region + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -24.1189  -3.9843   0.2086   3.7032  27.9970
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -1.167e+02  3.229e+00 -36.130  < 2e-16 ***
# pd$DxSchizo          -1.930e+00  5.183e-01  -3.723 0.000213 ***
# pd$Age               -3.836e-02  1.709e-02  -2.245 0.025089 *
# pd$SexF              -1.569e+00  5.329e-01  -2.944 0.003346 **
# pd$mitoRate          -1.194e+02  7.327e+00 -16.299  < 2e-16 ***
# pd$RegionHIPPO       -1.223e+01  1.454e+00  -8.410 2.28e-16 ***
# pd$rRNA_rate          1.166e+04  8.543e+03   1.365 0.172633
# pd$totalAssignedGene  1.827e+02  6.779e+00  26.950  < 2e-16 ***
# pd$RIN                7.327e+00  3.122e-01  23.470  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 6.431 on 703 degrees of freedom
# Multiple R-squared:  0.8744,    Adjusted R-squared:  0.873
# F-statistic: 611.9 on 8 and 703 DF,  p-value: < 2.2e-16







## Explore qSVs with no Hippo Gold samples: only HIPPO
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
pd <- pd_all[keepIndex, ]

getPcaVars(qsvBonf)[1:5]
# [1] 60.20  7.85  4.63  2.46  2.04

## brainseq qsv plots
pdf("pdf/qSVs_brainseq_age17_noHGold_HIPPO.pdf", useDingbats = FALSE)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Dx,
     xlab = "Diagnosis",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ factor(pd$Region),
     xlab = "Region",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$RIN,
     xlab = "RIN",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('skyblue3')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
# exploring differences in prenatal samples
boxplot(qsvBonf$x[,2] ~ pd$Age<0,
	xlab = "Age<0",
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
boxplot(qsvBonf$x[,2] ~ factor(ifelse(pd$Age<0, 'Prenatal', 'Postnatal')) + pd$Region,
	xlab = "Age and Region", cex.axis=0.9,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('skyblue3')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ factor(pd$Sex),
     xlab = "Sex",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
boxplot(qsvBonf$x[,1] ~ pd$Region + pd$Dx,
	xlab = "Region and Diagnosis",
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
dev.off()


print('Regression: qsv1 vs totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -55.945  -7.603   3.437  10.194  24.320
#
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)           -51.619      9.778  -5.279 2.36e-07 ***
# pd$totalAssignedGene   99.061     18.703   5.297 2.16e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.51 on 331 degrees of freedom
# Multiple R-squared:  0.07813,   Adjusted R-squared:  0.07535
# F-statistic: 28.05 on 1 and 331 DF,  p-value: 2.157e-07

print('Regression: qsv1 vs Dx')
summary(lm(qsvBonf$x[,1] ~ pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -47.160  -8.292   2.477  10.051  40.485
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    3.784      1.016   3.722 0.000232 ***
# pd$DxSchizo   -9.474      1.608  -5.890 9.48e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.38 on 331 degrees of freedom
# Multiple R-squared:  0.09487,   Adjusted R-squared:  0.09214
# F-statistic: 34.69 on 1 and 331 DF,  p-value: 9.483e-09

print('Regression: qsv1 vs rRNA_rate')
summary(lm(qsvBonf$x[,1] ~ pd$rRNA_rate))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$rRNA_rate)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -44.716  -8.827   2.795  11.687  35.272
#
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -0.5813     1.1917  -0.488    0.626
# pd$rRNA_rate 20022.3735 29539.1256   0.678    0.498
#
# Residual standard error: 15.1 on 331 degrees of freedom
# Multiple R-squared:  0.001386,  Adjusted R-squared:  -0.001631
# F-statistic: 0.4594 on 1 and 331 DF,  p-value: 0.4984

print('Regression: qsv1 vs RIN')
summary(lm(qsvBonf$x[,1] ~ pd$RIN))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$RIN)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -42.105  -6.089  -0.071   5.438  35.164
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -86.5136     3.6578  -23.65   <2e-16 ***
# pd$RIN       11.4109     0.4779   23.88   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 9.158 on 331 degrees of freedom
# Multiple R-squared:  0.6327,    Adjusted R-squared:  0.6316
# F-statistic: 570.1 on 1 and 331 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -56.847  -7.629   3.474   9.561  31.124
#
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)           -38.812      9.765  -3.974 8.67e-05 ***
# pd$DxSchizo            -8.094      1.597  -5.068 6.72e-07 ***
# pd$totalAssignedGene   80.687     18.403   4.384 1.57e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14 on 330 degrees of freedom
# Multiple R-squared:  0.1447,    Adjusted R-squared:  0.1395
# F-statistic: 27.91 on 2 and 330 DF,  p-value: 6.309e-12

print('Regression: qsv1 vs totalAssignedGene*Dx')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -53.725  -7.689   3.415   9.720  26.306
#
# Coefficients:
#                                  Estimate Std. Error t value Pr(>|t|)
# (Intercept)                        -25.08      12.39  -2.025   0.0437 *
# pd$totalAssignedGene                54.68      23.39   2.338   0.0200 *
# pd$DxSchizo                        -43.03      19.57  -2.199   0.0286 *
# pd$totalAssignedGene:pd$DxSchizo    67.53      37.69   1.792   0.0741 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 13.95 on 329 degrees of freedom
# Multiple R-squared:  0.153,     Adjusted R-squared:  0.1452
# F-statistic:  19.8 on 3 and 329 DF,  p-value: 7.965e-12

print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$rRNA_rate + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$rRNA_rate + pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -20.830  -5.306  -0.077   4.648  23.595
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -1.217e+02  7.150e+00 -17.021  < 2e-16 ***
# pd$DxSchizo          -2.244e+00  9.964e-01  -2.253  0.02495 *
# pd$Age               -1.377e-01  3.192e-02  -4.314 2.12e-05 ***
# pd$SexF              -2.521e+00  9.590e-01  -2.629  0.00897 **
# pd$mitoRate          -2.661e+02  1.091e+01 -24.404  < 2e-16 ***
# pd$rRNA_rate          1.004e+04  1.652e+04   0.608  0.54381
# pd$totalAssignedGene  3.768e+02  1.618e+01  23.286  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 7.982 on 326 degrees of freedom
# Multiple R-squared:  0.7251,    Adjusted R-squared:  0.7201
# F-statistic: 143.3 on 6 and 326 DF,  p-value: < 2.2e-16


print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$rRNA_rate + pd$totalAssignedGene + pd$RIN)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -17.6581  -4.2996  -0.3485   3.2335  26.2595
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -1.390e+02  5.890e+00 -23.597  < 2e-16 ***
# pd$DxSchizo          -2.025e+00  8.011e-01  -2.528  0.01195 *
# pd$Age               -6.499e-02  2.622e-02  -2.478  0.01371 *
# pd$SexF              -1.560e+00  7.742e-01  -2.015  0.04472 *
# pd$mitoRate          -1.657e+02  1.153e+01 -14.365  < 2e-16 ***
# pd$rRNA_rate          3.777e+04  1.344e+04   2.810  0.00526 **
# pd$totalAssignedGene  2.613e+02  1.560e+01  16.747  < 2e-16 ***
# pd$RIN                6.292e+00  4.696e-01  13.399  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 6.416 on 325 degrees of freedom
# Multiple R-squared:  0.823,     Adjusted R-squared:  0.8191
# F-statistic: 215.8 on 7 and 325 DF,  p-value: < 2.2e-16









## Explore qSVs with no Hippo Gold samples: only DLPFC
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)
pd <- pd_all[keepIndex, ]

getPcaVars(qsvBonf)[1:5]
# [1] 69.60  5.62  2.73  1.90  1.42

## brainseq qsv plots
pdf("pdf/qSVs_brainseq_age17_noHGold_DLPFC.pdf", useDingbats = FALSE)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('dark orange')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Dx,
     xlab = "Diagnosis",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ factor(pd$Region),
     xlab = "Region",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$RIN,
     xlab = "RIN",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('dark orange')[factor(pd$Region)])
legend('topleft', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('dark orange')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('dark orange')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('dark orange')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,2] ~ pd$Age,
     xlab = "Age",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c('dark orange')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
# exploring differences in prenatal samples
boxplot(qsvBonf$x[,2] ~ pd$Age<0,
	xlab = "Age<0",
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
boxplot(qsvBonf$x[,2] ~ factor(ifelse(pd$Age<0, 'Prenatal', 'Postnatal')) + pd$Region,
	xlab = "Age and Region", cex.axis=0.9,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c('dark orange')[factor(pd$Region)])
legend('topright', c('DLPFC', 'HIPPO'), lwd = 2, col = c('dark orange', 'skyblue3'), bty = 'n')
plot(qsvBonf$x[,1] ~ factor(pd$Sex),
     xlab = "Sex",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
boxplot(qsvBonf$x[,1] ~ pd$Region + pd$Dx,
	xlab = "Region and Diagnosis",
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
dev.off()


print('Regression: qsv1 vs totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -31.860  -3.247   1.152   4.043  47.056
#
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -101.846      3.711  -27.45   <2e-16 ***
# pd$totalAssignedGene  248.021      8.972   27.64   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 8.672 on 377 degrees of freedom
# Multiple R-squared:  0.6697,    Adjusted R-squared:  0.6688
# F-statistic: 764.2 on 1 and 377 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx')
summary(lm(qsvBonf$x[,1] ~ pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -51.244  -6.031   2.782  10.044  30.291
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   2.3296     0.9857   2.363  0.01862 *
# pd$DxSchizo  -5.7707     1.5514  -3.720  0.00023 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.82 on 377 degrees of freedom
# Multiple R-squared:  0.0354,    Adjusted R-squared:  0.03284
# F-statistic: 13.84 on 1 and 377 DF,  p-value: 0.0002298

print('Regression: qsv1 vs rRNA_rate')
summary(lm(qsvBonf$x[,1] ~ pd$rRNA_rate))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$rRNA_rate)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -44.616  -6.732   2.499   9.493  34.317
#
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# (Intercept)   7.575e+00  1.388e+00   5.456 8.84e-08 ***
# pd$rRNA_rate -1.494e+05  2.322e+04  -6.433 3.79e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.32 on 377 degrees of freedom
# Multiple R-squared:  0.09892,   Adjusted R-squared:  0.09653
# F-statistic: 41.39 on 1 and 377 DF,  p-value: 3.79e-10

print('Regression: qsv1 vs RIN')
summary(lm(qsvBonf$x[,1] ~ pd$RIN))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$RIN)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -34.175  -5.544  -0.066   6.460  27.971
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -95.1783     4.3682  -21.79   <2e-16 ***
# pd$RIN       12.3931     0.5648   21.94   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 9.999 on 377 degrees of freedom
# Multiple R-squared:  0.5608,    Adjusted R-squared:  0.5597
# F-statistic: 481.4 on 1 and 377 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + totalAssigned')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -31.992  -3.366   1.201   4.000  47.106
#
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -101.5031     3.9068 -25.981   <2e-16 ***
# pd$DxSchizo            -0.2650     0.9318  -0.284    0.776
# pd$totalAssignedGene  247.4453     9.2081  26.873   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 8.682 on 376 degrees of freedom
# Multiple R-squared:  0.6697,    Adjusted R-squared:  0.668
# F-statistic: 381.2 on 2 and 376 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssignedGene*Dx')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Dx)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -31.995  -3.368   1.202   4.000  47.120
#
# Coefficients:
#                                   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                      -101.4808     5.0834 -19.963   <2e-16 ***
# pd$totalAssignedGene              247.3922    12.0357  20.555   <2e-16 ***
# pd$DxSchizo                        -0.3173     7.6701  -0.041    0.967
# pd$totalAssignedGene:pd$DxSchizo    0.1286    18.7257   0.007    0.995
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 8.694 on 375 degrees of freedom
# Multiple R-squared:  0.6697,    Adjusted R-squared:  0.6671
# F-statistic: 253.5 on 3 and 375 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$rRNA_rate + pd$totalAssignedGene))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$rRNA_rate + pd$totalAssignedGene)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -30.067  -3.727   1.014   4.565  38.528
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -9.264e+01  5.006e+00 -18.506  < 2e-16 ***
# pd$DxSchizo           1.007e-01  8.735e-01   0.115  0.90827
# pd$Age               -3.404e-02  2.937e-02  -1.159  0.24718
# pd$SexF              -2.499e+00  9.239e-01  -2.704  0.00716 **
# pd$mitoRate           1.154e+02  5.863e+01   1.967  0.04987 *
# pd$rRNA_rate         -8.762e+04  1.351e+04  -6.484 2.84e-10 ***
# pd$totalAssignedGene  2.361e+02  9.188e+00  25.696  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 8.059 on 372 degrees of freedom
# Multiple R-squared:  0.7185,    Adjusted R-squared:  0.7139
# F-statistic: 158.2 on 6 and 372 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN')
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Dx + pd$Age + pd$Sex + pd$mitoRate +
#     pd$rRNA_rate + pd$totalAssignedGene + pd$RIN)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -20.3239  -3.4711   0.7331   3.9102  23.6099
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -1.281e+02  4.381e+00 -29.229   <2e-16 ***
# pd$DxSchizo          -9.517e-01  6.686e-01  -1.423   0.1554
# pd$Age               -7.302e-03  2.244e-02  -0.325   0.7450
# pd$SexF              -1.444e+00  7.068e-01  -2.042   0.0418 *
# pd$mitoRate           5.455e+01  4.482e+01   1.217   0.2244
# pd$rRNA_rate         -1.669e+04  1.116e+04  -1.495   0.1358
# pd$totalAssignedGene  1.789e+02  7.819e+00  22.882   <2e-16 ***
# pd$RIN                7.219e+00  4.395e-01  16.427   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 6.14 on 371 degrees of freedom
# Multiple R-squared:  0.837,     Adjusted R-squared:  0.8339
# F-statistic: 272.2 on 7 and 371 DF,  p-value: < 2.2e-16


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
#  date     2018-04-25
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  acepack                1.4.1     2016-10-29 CRAN (R 3.4.1)
#  AnnotationDbi          1.40.0    2017-11-29 Bioconductor
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.4.1)
#  backports              1.1.2     2017-12-13 CRAN (R 3.4.2)
#  base                 * 3.4.3     2018-01-20 local
#  base64enc              0.1-3     2015-07-28 CRAN (R 3.4.1)
#  bindr                  0.1.1     2018-03-13 CRAN (R 3.4.3)
#  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.4.3)
#  Biobase              * 2.38.0    2017-11-07 Bioconductor
#  BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
#  BiocParallel         * 1.12.0    2017-11-29 Bioconductor
#  biomaRt                2.34.2    2018-02-17 Bioconductor
#  Biostrings             2.46.0    2017-11-29 Bioconductor
#  bit                    1.1-12    2014-04-09 CRAN (R 3.4.1)
#  bit64                  0.9-7     2017-05-08 CRAN (R 3.4.1)
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
#  blob                   1.1.1     2018-03-25 CRAN (R 3.4.3)
#  BSgenome               1.46.0    2017-11-29 Bioconductor
#  bumphunter             1.20.0    2017-11-29 Bioconductor
#  cellranger             1.1.0     2016-07-27 CRAN (R 3.4.1)
#  checkmate              1.8.5     2017-10-24 CRAN (R 3.4.2)
#  cluster                2.0.6     2017-03-10 CRAN (R 3.4.3)
#  codetools              0.2-15    2016-10-05 CRAN (R 3.4.3)
#  colorout             * 1.2-0     2018-02-19 Github (jalvesaq/colorout@2f01173)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
#  compiler               3.4.3     2018-01-20 local
#  data.table             1.10.4-3  2017-10-27 CRAN (R 3.4.2)
#  datasets             * 3.4.3     2018-01-20 local
#  DBI                    0.8       2018-03-02 CRAN (R 3.4.3)
#  DelayedArray         * 0.4.1     2017-11-07 Bioconductor
#  derfinder              1.12.6    2018-04-25 Bioconductor
#  derfinderHelper        1.12.0    2017-11-29 Bioconductor
#  devtools             * 1.13.5    2018-02-18 CRAN (R 3.4.3)
#  digest                 0.6.15    2018-01-28 cran (@0.6.15)
#  doRNG                  1.6.6     2017-04-10 CRAN (R 3.4.1)
#  downloader             0.4       2015-07-09 CRAN (R 3.4.1)
#  dplyr                  0.7.4     2017-09-28 CRAN (R 3.4.1)
#  foreach                1.4.4     2017-12-12 CRAN (R 3.4.2)
#  foreign                0.8-69    2017-06-22 CRAN (R 3.4.3)
#  Formula                1.2-2     2017-07-10 CRAN (R 3.4.1)
#  GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
#  GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
#  GenomicAlignments      1.14.2    2018-04-18 Bioconductor
#  GenomicFeatures        1.30.3    2018-02-17 Bioconductor
#  GenomicFiles           1.14.0    2017-11-29 Bioconductor
#  GenomicRanges        * 1.30.3    2018-04-18 Bioconductor
#  GEOquery               2.46.15   2018-03-06 Bioconductor
#  ggplot2                2.2.1     2016-12-30 CRAN (R 3.4.1)
#  glue                   1.2.0     2017-10-29 CRAN (R 3.4.2)
#  graphics             * 3.4.3     2018-01-20 local
#  grDevices            * 3.4.3     2018-01-20 local
#  grid                   3.4.3     2018-01-20 local
#  gridExtra              2.3       2017-09-09 CRAN (R 3.4.1)
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
#  Hmisc                  4.1-1     2018-01-03 CRAN (R 3.4.2)
#  hms                    0.4.2     2018-03-10 CRAN (R 3.4.3)
#  htmlTable              1.11.2    2018-01-20 CRAN (R 3.4.3)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.4.3)
#  httpuv                 1.3.6.2   2018-03-02 CRAN (R 3.4.3)
#  httr                   1.3.1     2017-08-20 CRAN (R 3.4.1)
#  IRanges              * 2.12.0    2017-11-29 Bioconductor
#  iterators              1.0.9     2017-12-12 CRAN (R 3.4.2)
#  jaffelab             * 0.99.20   2018-04-19 Github (LieberInstitute/jaffelab@04c470a)
#  jsonlite               1.5       2017-06-01 CRAN (R 3.4.1)
#  knitr                  1.20      2018-02-20 CRAN (R 3.4.3)
#  later                  0.7.1     2018-03-07 CRAN (R 3.4.3)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
#  latticeExtra           0.6-28    2016-02-09 CRAN (R 3.4.1)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
#  limma                  3.34.9    2018-04-18 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.4.1)
#  magrittr               1.5       2014-11-22 CRAN (R 3.4.1)
#  Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
#  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.3     2018-01-20 local
#  mime                   0.5       2016-07-07 CRAN (R 3.4.1)
#  munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
#  nnet                   7.3-12    2016-02-02 CRAN (R 3.4.3)
#  parallel             * 3.4.3     2018-01-20 local
#  pillar                 1.2.1     2018-02-27 CRAN (R 3.4.3)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.4.1)
#  pkgmaker               0.22      2014-05-14 CRAN (R 3.4.1)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
#  png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
#  prettyunits            1.0.2     2015-07-13 CRAN (R 3.4.1)
#  progress               1.1.2     2016-12-14 CRAN (R 3.4.1)
#  purrr                  0.2.4     2017-10-18 CRAN (R 3.4.2)
#  qvalue                 2.10.0    2017-11-29 Bioconductor
#  R.methodsS3            1.7.1     2016-02-16 CRAN (R 3.4.1)
#  R.oo                   1.22.0    2018-04-22 CRAN (R 3.4.3)
#  R.utils                2.6.0     2017-11-05 CRAN (R 3.4.2)
#  R6                     2.2.2     2017-06-17 CRAN (R 3.4.1)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.4.1)
#  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.4.1)
#  Rcpp                   0.12.16   2018-03-13 CRAN (R 3.4.3)
#  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
#  readr                  1.1.1     2017-05-16 CRAN (R 3.4.1)
#  readxl               * 1.1.0     2018-04-20 CRAN (R 3.4.3)
#  recount              * 1.4.6     2018-04-18 Bioconductor
#  recount.bwtool       * 0.99.20   2017-08-10 Github (LieberInstitute/recount.bwtool@db4155a)
#  registry               0.5       2017-12-03 CRAN (R 3.4.2)
#  rentrez                1.2.1     2018-03-05 CRAN (R 3.4.3)
#  reshape2               1.4.3     2017-12-11 CRAN (R 3.4.2)
#  rlang                  0.2.0     2018-02-20 CRAN (R 3.4.3)
#  rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
#  RMySQL                 0.10.14   2018-02-26 CRAN (R 3.4.3)
#  rngtools               1.2.4     2014-03-06 CRAN (R 3.4.1)
#  rpart                  4.1-12    2018-01-12 CRAN (R 3.4.3)
#  Rsamtools              1.30.0    2017-11-29 Bioconductor
#  RSQLite                2.1.0     2018-03-29 CRAN (R 3.4.3)
#  rstudioapi             0.7       2017-09-07 CRAN (R 3.4.1)
#  rtracklayer          * 1.38.3    2018-02-17 Bioconductor
#  S4Vectors            * 0.16.0    2017-11-29 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.4.2)
#  servr                  0.9       2018-03-25 CRAN (R 3.4.3)
#  splines                3.4.3     2018-01-20 local
#  stats                * 3.4.3     2018-01-20 local
#  stats4               * 3.4.3     2018-01-20 local
#  stringi                1.1.7     2018-03-12 CRAN (R 3.4.3)
#  stringr                1.3.0     2018-02-19 CRAN (R 3.4.3)
#  SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
#  survival               2.41-3    2017-04-04 CRAN (R 3.4.3)
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.4.3)
#  tidyr                  0.8.0     2018-01-29 CRAN (R 3.4.3)
#  tools                  3.4.3     2018-01-20 local
#  utils                * 3.4.3     2018-01-20 local
#  VariantAnnotation      1.24.5    2018-01-16 Bioconductor
#  withr                  2.1.2     2018-03-15 CRAN (R 3.4.3)
#  xfun                   0.1       2018-01-22 CRAN (R 3.4.3)
#  XML                    3.98-1.10 2018-02-19 CRAN (R 3.4.3)
#  xml2                   1.2.0     2018-01-24 CRAN (R 3.4.3)
#  xtable                 1.8-2     2016-02-05 CRAN (R 3.4.1)
#  XVector                0.18.0    2017-11-29 Bioconductor
#  zlibbioc               1.24.0    2017-11-07 Bioconductor
