
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
getPcaVars(pca)[1:5]

## degradation plots
dir.create('pdf', showWarnings = FALSE)
pdf('pdf/degradation.pdf')
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
getPcaVars(pca)[1:5]

## brainseq plots
pdf('pdf/brainseqplots.pdf')
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
plot(pca$x[,1] ~ factor(pd$Region),
     xlab = "Region", 
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))
plot(pca$x[,1] ~ pd$Dx,
     xlab = "Dx", 
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))
dev.off()


# top 1000 expressed regions associated with degradation
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))
getPcaVars(qsvBonf)[1:5]

## brainseq qsv plots
pdf("pdf/qSVs_brainseq.pdf")
mypar(2,2)
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

