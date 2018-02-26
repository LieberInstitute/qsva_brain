## Based on /dcl01/lieber/ajaffe/lab/degradation_experiments/Joint/bipseq_sACC_Amygdala_RiboZero/make_ERs_stranded.R

###
library(rtracklayer)
library(derfinder)
library(jaffelab)
library(SummarizedExperiment)
library(recount.bwtool)
library(readxl)
library(BiocParallel)
library(limma)
library(edgeR)
library(biomaRt)
library(GenomicRanges)


dir.create("bed")
dir.create("rdas")
pdf('make_ERs_stranded.pdf')

## read in phenotype data
## subset by brain regions of interest, DLPFC and HIPPO
load('/dcl01/ajaffe/data/lab/qsva_brain/expr_data/rda/rse_gene.Rdata')
rse_gene <- rse_gene[, rse_gene$Region %in% c('DLPFC', 'HIPPO')]
pd = colData(rse_gene)

## chromosome info
chrInfo <- read.table('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode',
    header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
chrInfo <- subset(chrInfo, chr %in% paste0('chr', c(1:22, 'X', 'Y', 'M')))

## load stranded mean bigwigs
meanBWs = paste0("mean_",c("forward","reverse"), ".bw")
names(meanBWs) = c("Forward","Reverse")
meanList = lapply(meanBWs, import)
strand(meanList$Forward) = Rle("+")
strand(meanList$Reverse) = Rle("-")

## make ERs
meanListFilter = GRangesList(lapply(meanList, function(x) x[abs(x$score) >= 5]))
reduceList = endoapply(meanListFilter, reduce,min.gapwidth=2)

## at least 50 BP and on main chromosomes
erList = endoapply(reduceList, function(x) 
	x[width(x) >= 50 & seqnames(x) %in% chrInfo$chr])

##########################3


## keep both regions
pd$RNum = ss(pd$SAMPLE_ID, "_")

## use recount bw tool	

bw_f_dlpfc <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'DLPFC'], '.Forward.bw')

bw_f_hippo <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'HIPPO'], '.Forward.bw')

bw_f<- c(bw_f_dlpfc, bw_f_hippo)

bwFilesForward = c(bw_f_dlpfc, bw_f_hippo)


bw_r_dlpfc <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'DLPFC'], '.Reverse.bw')

bw_r_hippo <- paste0('/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_RiboZero/Coverage/', colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'HIPPO'], '.Reverse.bw')

bwFilesReverse = c(bw_r_dlpfc, bw_r_hippo)
names(bwFilesForward) =names(bwFilesReverse) = c(colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'DLPFC'], colData(rse_gene)$SAMPLE_ID[colData(rse_gene)$Region == 'HIPPO'])
all(file.exists(c(bwFilesForward, bwFilesReverse))) # TRUE

## get sums
covForward = coverage_bwtool(bwFilesForward, erList$Forward, 
	sumsdir = "ers", bpparam = MulticoreParam(8) ,strand = "+")
covReverse = coverage_bwtool(bwFilesReverse, erList$Reverse, 
	sumsdir = "ers", bpparam = MulticoreParam(8), strand="-")

covForward$bigwig_path = NULL
covForward$bigwig_file = NULL
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## divide by number of reads
assays(covForward)$counts = assays(covForward)$counts/100 # divide by read length
assays(covForward)$counts = abs(assays(covForward)$counts) 
assays(covReverse)$counts = assays(covReverse)$counts/100 # divide by read length
assays(covReverse)$counts = abs(assays(covReverse)$counts) 

## combine
covComb = SummarizedExperiment(
	assays = list(counts = rbind(assays(covForward)$counts,
		assays(covReverse)$counts)),
	colData = pd, rowData = c(rowRanges(covForward), rowRanges(covReverse)))
rownames(covComb) = paste0(seqnames(covComb), ":",start(covComb), "-",
	end(covComb), "(", strand(covComb), ")")

###############
## annotate ###
###############



## load genomic state
load('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/gs/gs_gencode_v25_hg38.Rdata')
gs = gs_gencode_v25_hg38$fullGenome

all <- GRanges(names(seqlengths(gs)), IRanges(start = 1, end = seqlengths(gs)), strand = '*')
new_intergenic <- lapply(c('+', '-'), function(str) {
    onestrand <- gs[strand(gs) == str | strand(gs) == '*']
    onestrand_simple <- onestrand
    mcols(onestrand_simple) <- NULL
    pieces <- disjoin(c(all, onestrand_simple), ignore.strand = TRUE)
    ov <- countOverlaps(pieces, onestrand_simple, ignore.strand = TRUE)
    new_intergenic <- pieces[ov == 0]
    strand(new_intergenic) <- str
    new_intergenic$theRegion <- paste0('intergenic', str)
    new_intergenic$tx_id <- IntegerList(NA)
    new_intergenic$tx_name <- CharacterList(NA)
    new_intergenic$gene <- IntegerList(NA)
    return(new_intergenic)
})
gs_stranded <- c(gs, unlist(GRangesList(new_intergenic)))

ensemblAnno = annotateRegions(rowRanges(covComb),gs_stranded,ignore.strand=FALSE,minoverlap=1)
countTable = ensemblAnno$countTable 
countTable[,2] = rowSums(countTable[,2:4])
countTable = countTable[,c(1,2,5)]

## gene annotation overlaps
geneMapGR = rowRanges(rse_gene) 
dA = distanceToNearest(rowRanges(covComb), geneMapGR)
rowRanges(covComb)$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
rowRanges(covComb)$nearestID = names(geneMapGR)[subjectHits(dA)]
rowRanges(covComb)$distToGene = mcols(dA)$distance
mcols(rowRanges(covComb)) = cbind(mcols(rowRanges(covComb)), countTable)

## add additional annotation
rowRanges(covComb)$annoClass = NA
rowRanges(covComb)$annoClass[rowRanges(covComb)$exon > 0 & 
	rowRanges(covComb)$intron == 0 &
	rowRanges(covComb)$intergenic == 0] = "strictExonic"
rowRanges(covComb)$annoClass[rowRanges(covComb)$exon == 0 & 
	rowRanges(covComb)$intron > 0 &
	rowRanges(covComb)$intergenic == 0] = "strictIntronic"
rowRanges(covComb)$annoClass[rowRanges(covComb)$exon == 0 & 
	rowRanges(covComb)$intron == 0 &
	rowRanges(covComb)$intergenic > 0] = "strictIntergenic"
rowRanges(covComb)$annoClass[rowRanges(covComb)$exon > 0 & 
	rowRanges(covComb)$intron > 0 &
	rowRanges(covComb)$intergenic == 0] = "exonIntron"
rowRanges(covComb)$annoClass[rowRanges(covComb)$exon > 0 & 
	rowRanges(covComb)$intergenic > 0] = "extendUTR"
	
save(covComb, file = "rdas/expressedRegions_DLPFC_Plus_HIPPO_RiboZero_degradation_cut5_hg38_n40.rda")

## main model
mod = model.matrix(~pd$DegradationTime +
	pd$Region + factor(pd$BrNum))
fit = lmFit(log2(assays(covComb)$counts+1), mod)
eb = ebayes(fit)
out = topTable(eBayes(fit),coef=2,n = nrow(covComb))

## interaction model
modInt = model.matrix(~pd$DegradationTime*pd$Region + 
	factor(pd$BrNum))
fitInt = lmFit(log2(assays(covComb)$counts+1), modInt)
ebInt = ebayes(fitInt)
outInt = topTable(eBayes(fitInt),coef=c(2,8), n = nrow(covComb))
outInt = outInt[rownames(out),]

sum(p.adjust(ebInt$p[,8],"fdr") < 0.05) # good

plot(outInt$F[out$, out$t)

## write out
dir.create("bed")
digBed =rowRanges(covComb)[rownames(out)[1:1000]]
export(digBed, con = "bed/DLPFC_Plus_HIPPO_RiboZero_degradation_top1000.bed")

#########################################
## make degradation features for genes
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=FALSE)
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)
degradeStats = topTable(eBayes(fitGene),coef=2,
	p.value = 1,number=nrow(rse_gene))
degradeStats$gencodeTx = NULL
degradeStats$bonf = NA
degradeStats$bonf[degradeStats$AveExpr > -2] = p.adjust(degradeStats$P.Value[degradeStats$AveExpr > -2] , "bonf")
degradeStats = degradeStats[rownames(rse_gene),]

## add interaction p-value
vGeneInt = voom(dge,mod,plot=FALSE)
fitGeneInt = lmFit(vGeneInt)
degradeStatsInt = topTable(eBayes(fitGeneInt),coef=2,
	p.value = 1,number=nrow(rse_gene))

save(degradeStats, file = "rdas/DLPFC_Plus_HIPPO_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda")

dev.off()

