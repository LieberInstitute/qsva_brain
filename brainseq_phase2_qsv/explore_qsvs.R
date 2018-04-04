
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
dir.create('pdf', showWarnings = FALSE)
pdf('pdf/degradation.pdf')
plot(pca$x[,1] ~ pd$totalAssignedGene)
plot(pca$x[,1] ~ factor(pd$Region))
plot(pca$x[,2] ~ factor(pd$Region))
plot(pca$x[,1] ~ pd$RIN)
plot(pca$x[,1] ~ pd$mitoRate)
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
pdf('pdf/brainseqplots.pdf')
plot(pca$x[,1] ~ pd$totalAssignedGene)
plot(pca$x[,1] ~ pd$RIN)
plot(pca$x[,1] ~ pd$mitoRate)
plot(pca$x[,2] ~ pd$mitoRate)
plot(pca$x[,1] ~ factor(pd$Region))
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
     xlab = "Gene Assignment Rate",pch=19,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c('orange', 'light blue')[c('Gold' = 1, 'HMR' = 2)[pd$Kit]])

plot(qsvBonf$x[,1] ~ pd$Dx,
     xlab = "Dx",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ factor(pd$Region),
     xlab = "Region",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$RIN,
     xlab = "RIN",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate",	ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
     
plot(qsvBonf$x[,1] ~ pd$Age,
     xlab = "Age",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,2] ~ pd$Age,
     xlab = "Age",	ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
     
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

print('Regression: qsv1 vs Kit')
summary(lm(qsvBonf$x[,1] ~ pd$Kit))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$Kit)
#
# Residuals:
#    Min      1Q  Median      3Q     Max 
# -49.085  -8.018   2.924  10.951  33.340 
#
# Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   8.3864     0.6668   12.58   <2e-16 ***
# pd$KitHMR   -18.4543     0.9892  -18.66   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 14.78 on 898 degrees of freedom
# Multiple R-squared:  0.2793,	Adjusted R-squared:  0.2785 
# F-statistic: 348.1 on 1 and 898 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssigned + Kit')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene + pd$Kit))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene + pd$Kit)
#
# Residuals:
#    Min      1Q  Median      3Q     Max 
# -63.119  -5.670   2.514   7.340  36.036 
#
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -46.659      3.311  -14.09   <2e-16 ***
# pd$totalAssignedGene  135.466      8.023   16.89   <2e-16 ***
# pd$KitHMR             -32.463      1.196  -27.13   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 12.88 on 897 degrees of freedom
# Multiple R-squared:  0.4531,	Adjusted R-squared:  0.4519 
# F-statistic: 371.6 on 2 and 897 DF,  p-value: < 2.2e-16

print('Regression: qsv1 vs totalAssigned*Kit')
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

print('Regression: qsv1 vs totalAssigned*Kit*Dx')
summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene * pd$Kit * pd$Dx))
# Call:
# lm(formula = qsvBonf$x[, 1] ~ pd$totalAssignedGene * pd$Kit * 
#    pd$Dx)
#
# Residuals:
#    Min      1Q  Median      3Q     Max 
# -50.785  -5.273   1.711   6.804  45.280 
#
# Coefficients:
#                                           Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                 -76.957      4.794 -16.054  < 2e-16
# pd$totalAssignedGene                        210.873     11.587  18.200  < 2e-16
# pd$KitHMR                                    58.966      7.635   7.723 3.05e-14
# pd$DxSchizo                                 -14.179      8.964  -1.582   0.1141
# pd$totalAssignedGene:pd$KitHMR             -188.785     16.380 -11.525  < 2e-16
# pd$totalAssignedGene:pd$DxSchizo             32.870     22.174   1.482   0.1386
# pd$KitHMR:pd$DxSchizo                       -42.810     15.946  -2.685   0.0074
# pd$totalAssignedGene:pd$KitHMR:pd$DxSchizo   58.970     33.995   1.735   0.0831
                                              
# (Intercept)                                ***
# pd$totalAssignedGene                       ***
# pd$KitHMR                                  ***
# pd$DxSchizo                                   
# pd$totalAssignedGene:pd$KitHMR             ***
# pd$totalAssignedGene:pd$DxSchizo              
# pd$KitHMR:pd$DxSchizo                      ** 
# pd$totalAssignedGene:pd$KitHMR:pd$DxSchizo .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 11.43 on 892 degrees of freedom
# Multiple R-squared:  0.5716,	Adjusted R-squared:  0.5682 
# F-statistic:   170 on 7 and 892 DF,  p-value: < 2.2e-16





