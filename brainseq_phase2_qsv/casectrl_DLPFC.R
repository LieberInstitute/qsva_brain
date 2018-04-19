
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

#load qsvBonf, qSVs, mod, modQsva object from brainseq data 
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs.Rdata", verbose = TRUE)
#filter for age and DLPFC; no age in rse_gene
keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('pdf/dlpfc_voom_qsva.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('pdf/dlpfc_voom_noqsva.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]
save(outGene, outGene0,file = "rdas/dxStats_dlpfc_filtered_qSVA_geneLevel.rda")



