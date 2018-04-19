
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

#load qsvBonf, qSVs, mod, modQsva object from brainseq data
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs.Rdata", verbose = TRUE)
#filter for age and hippocampus; no age in rse_gene
keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('pdf/voom_hippo_qsva.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('pdf/hippo_voom_noqsva.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('pdf/hippo_voom_noadj.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]

stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "rdas/dxStats_hippo_filtered_qSVA_geneLevel.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

