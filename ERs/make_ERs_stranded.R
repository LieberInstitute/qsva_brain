## Based on /dcl01/lieber/ajaffe/lab/degradation_experiments/Joint/bipseq_sACC_Amygdala_RiboZero/make_ERs_stranded.R

###
library('rtracklayer')
library('derfinder')
library('jaffelab')
library('SummarizedExperiment')
library('recount.bwtool')
library('readxl')
library('BiocParallel')
library('limma')
library('edgeR')
library('biomaRt')
library('GenomicRanges')
library('devtools')


dir.create("bed", showWarnings = FALSE)
dir.create("rdas", showWarnings = FALSE)

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
meanBWs = paste0("/dcl01/ajaffe/data/lab/qsva_brain/means/mean_",c("forward","reverse"), ".bw")
names(meanBWs) = c("Forward","Reverse")
stopifnot(all(file.exists(meanBWs)))
meanList = lapply(meanBWs, function(x) {
    message(paste(Sys.time(), 'importing', x))
    import(x)
})
strand(meanList$Forward) = Rle("+")
strand(meanList$Reverse) = Rle("-")

## make ERs
message(paste(Sys.time(), 'filtering mean BWs'))
meanListFilter = GRangesList(lapply(meanList, function(x) x[abs(x$score) >= 5]))
reduceList = endoapply(meanListFilter, reduce,min.gapwidth=2)

## at least 50 BP and on main chromosomes
erList = endoapply(reduceList, function(x)
	x[width(x) >= 50 & seqnames(x) %in% chrInfo$chr])

## Clean up (to reduce mem burden)
rm(meanList, meanListFilter, reduceList)

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
stopifnot(all(file.exists(c(bwFilesForward, bwFilesReverse))))

## get sums
covForward = coverage_bwtool(bwFilesForward, erList$Forward,
	sumsdir = "ers", bpparam = MulticoreParam(1) ,strand = "+")
covReverse = coverage_bwtool(bwFilesReverse, erList$Reverse,
	sumsdir = "ers", bpparam = MulticoreParam(1), strand="-")

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
## Drop dummy var for last BrNum, otherwise the model is not well specified
fit = lmFit(log2(assays(covComb)$counts+1), mod[, -ncol(mod)])
eb = ebayes(fit)
out = topTable(eBayes(fit),coef=2,n = nrow(covComb))

## interaction model
modInt = model.matrix(~pd$DegradationTime*pd$Region +
	factor(pd$BrNum))
fitInt = lmFit(log2(assays(covComb)$counts+1), modInt[, - (ncol(modInt) - 1)])
ebInt = ebayes(fitInt)
outInt = topTable(eBayes(fitInt),coef=c(2, ncol(modInt) -1 ), n = nrow(covComb))
outInt = outInt[rownames(out),]


save(out, outInt, file = "rdas/DLPFC_Plus_HIPPO_RiboZero_ERlevel_degradationStats_forDEqual_hg38.rda")

sum(p.adjust(ebInt$p[, ncol(ebInt$p)],"fdr") < 0.05) # good
sum(p.adjust(eb$p[, 2],"fdr") < 0.05)
sum(p.adjust(ebInt$p[, 2],"fdr") < 0.05)

pdf('make_ERs_stranded.pdf')
plot(outInt$F, out$t)
plot(-log10(ebInt$p[, 2]), -log10(ebInt$p[, ncol(ebInt$p)]))
dev.off()

table(p.adjust(ebInt$p[head(rownames(out), n = 1000), ncol(ebInt$p)],"fdr") < 0.05)


## write out
dir.create("bed", showWarnings = FALSE)
digBed =rowRanges(covComb)[head(rownames(out), n = 1000)]
export(digBed, con = "bed/DLPFC_Plus_HIPPO_RiboZero_degradation_top1000.bed")

#########################################
## make degradation features for genes
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod[, -ncol(mod)],plot=FALSE)
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)
degradeStats = topTable(eBayes(fitGene),coef=2,
	p.value = 1,number=nrow(rse_gene))
degradeStats$gencodeTx = NULL
degradeStats$bonf = NA
degradeStats$bonf[degradeStats$AveExpr > -2] = p.adjust(degradeStats$P.Value[degradeStats$AveExpr > -2] , "bonf")
degradeStats = degradeStats[rownames(rse_gene),]

## add interaction p-value
vGeneInt = voom(dge,modInt[, -(ncol(modInt) - 1)],plot=FALSE)
fitGeneInt = lmFit(vGeneInt)
degradeStatsInt = topTable(eBayes(fitGeneInt),coef = (ncol(modInt) - 1),
	p.value = 1,number=nrow(rse_gene))

degradeStatsInt$bonf = NA
degradeStatsInt$bonf[degradeStatsInt$AveExpr > -2] = p.adjust(degradeStatsInt$P.Value[degradeStatsInt$AveExpr > -2] , "bonf")
degradeStatsInt = degradeStatsInt[rownames(rse_gene),]

save(degradeStats, degradeStatsInt, file = "rdas/DLPFC_Plus_HIPPO_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda")


d_stats <- function(region) {
    d <- dge[, rse_gene$Region == region]
    m <- with(pd[rse_gene$Region == region, ], model.matrix( ~ DegradationTime + factor(BrNum)))

    v = voom(d, m,plot=FALSE)
    f = lmFit(v)
    eb = eBayes(f)
    dS = topTable(eBayes(f),coef=2,
    	p.value = 1,number=nrow(rse_gene))
    dS$gencodeTx = NULL
    dS$bonf = NA
    dS$bonf[dS$AveExpr > -2] = p.adjust(dS$P.Value[dS$AveExpr > -2] , "bonf")
    dS = dS[rownames(rse_gene),]
    return(dS)
}

degradeStats_DLPFC <- d_stats('DLPFC')
degradeStats_HIPPO <- d_stats('HIPPO')

save(degradeStats_DLPFC, degradeStats_HIPPO, file = 'rdas/DLPFC_HIPPO_degradationStats_hg38.rda')

stopifnot(identical(rownames(degradeStats), rownames(degradeStatsInt)))
stopifnot(identical(rownames(degradeStats), rownames(degradeStats_DLPFC)))
stopifnot(identical(rownames(degradeStats), rownames(degradeStats_HIPPO)))

Sys.time()
proc.time()
options(width = 120)
session_info()
