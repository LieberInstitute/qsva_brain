## Based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/expr_cutoff/expr_cutoff.R
library('SummarizedExperiment')
library('recount')
library('jaffelab')
library('devtools')
library('readxl')

dirs <- dir('/dcl01/lieber/ajaffe/lab/degradation_experiments', pattern = 'RiboZero', full.names = TRUE)
names(dirs) <- dir('/dcl01/lieber/ajaffe/lab/degradation_experiments', pattern = 'RiboZero')

file_list <- lapply(dirs, function(d) {
    files <- dir(d, pattern = '^rse', full.names = TRUE)
    names(files) <- dir(d, pattern = '^rse', full.names = TRUE)
    return(files)    
})

## Remove pre-qc Hippo data
file_list[['Hippo_RiboZero']] <- file_list[['Hippo_RiboZero']][grep('postQC', file_list[['Hippo_RiboZero']] )]


## Get files, feature types and brain regions
files <- unlist(file_list)
types <- gsub('_.*', '', gsub('.*rse_', '', tolower(files)))
regions <- toupper(gsub('_.*', '', names(files)))

## Load the raw data and calculate RPKMs & RP10M when necessary
all <- mapply(function(f, type, region) {
    message(paste(Sys.time(), 'processing brain region', region, 'for feature type', type))
    load(f, verbose = TRUE)
    if(type == 'gene') {
        assays(rse_gene)$rpkm <- recount::getRPKM(rse_gene, 'Length')
        rowRanges(rse_gene)$meanExprs <- NA
        colData(rse_gene)$Region <- region
        rse <- rse_gene
    } else if (type == 'exon') {
        rowRanges(rse_exon)$meanExprs <- NA
        assays(rse_exon)$rpkm <- recount::getRPKM(rse_exon, 'Length')
        colData(rse_exon)$Region <- region
        rse <- rse_exon
    } else if (type == 'jx') {
        ## Try getRPKM based on https://github.com/LieberInstitute/brainseq_phase2/blob/54c73b2b4cd65af93a254ff8f38eed6a8d5c362a/caseControl_analysis_hippo.R#L42
        rowRanges(rse_jx)$Length <- 100
        #assays(rse_jx)$rp10m <- recount::getRPKM(rse_jx, 'Length')
        rowRanges(rse_jx)$meanExprs <- NA
        colData(rse_jx)$Region <- region
        rse <- rse_jx
    } else if (type == 'tx') {
        rowRanges(rse_tx)$meanExprs <- NA
        colData(rse_tx)$Region <- region
        rse <- rse_tx
    }
    
    if(region %in% c('AMYGDALA', 'SACC')) {
        ## Based on /dcl01/lieber/ajaffe/lab/degradation_experiments/prep_limbic_data.R
        
        pd = as.data.frame(read_excel(
        	"/dcl01/lieber/ajaffe/lab/degradation_experiments//RIN_Amyg_sACC_degradation_20171003.xlsx",skip=3))
        colnames(pd)[8:9] = c("Region","DegradationTime")
        pd$pH <- NA
        
        colData(rse) <- cbind(colData(rse), pd[match(ss(rse$SAMPLE_ID, '_'), pd$RNum), ])
        
    } else if (region == 'CAUDATE') {
        ## Adapted from /dcl01/lieber/ajaffe/lab/degradation_experiments/Caudate_RiboZero/make_ERs_stranded.R
        pd =  as.data.frame(read_excel("/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_mPFC_Caudate_degradation_05_12_2017.xlsx"))
        colnames(pd)[c(1,2,7,8,14)] = c("Position","ID","Region", "DegradationTime", "Dx")
        pd = pd[pd$Region == "Caudate",]
        pd$SampleID = paste0(pd$ID, "_", pd$Flowcell)
        rownames(pd) = pd$SampleID
        pd = pd[rownames(colData(rse)),]
        pd$DegradationTime = as.numeric(ss(pd$DegradationTime,"mi"))
        colData(rse) <- cbind(colData(rse), pd[,c(1,3:13)])
    } else if (region == 'DLPFC') {
       ## Based on /dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/make_ERs_stranded.R
       load("/dcl01/lieber/ajaffe/lab/degradation_experiments/overall_degradation_pheno.rda")
        pd = pd[pd$Dataset == "DLPFC_RiboZero",]
        colData(rse) <- cbind(colData(rse), pd[rownames(colData(rse)), ])
        
    } else if (region == 'HIPPO') {
        ## Based on /dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_RiboZero/make_ERs_stranded.R
        pheno = as.data.frame(read_excel("/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_mPFC_Caudate_degradation_05_12_2017.xlsx"))
        colnames(pheno)[c(1,2,7,8,14)] = c("Position","ID","Region", "DegradationTime", "Dx")

        ## update hippo IDs
        pheno$ID[pheno$Region == "Hippo"] = pheno$BrNum[pheno$Region == "Hippo"]
        pheno$BrNum[pheno$Region == "Hippo"] = ss(pheno$ID[pheno$Region == "Hippo"], "-")

        pheno$SampleID = pheno$ID
        pheno$SampleID[pheno$Region == "Hippo"] = paste0(pheno$BrNum[pheno$Region == "Hippo"], 
        	"-H", rep(1:5, each=4), rep(c("a","b","c","d"), times=4))
        rownames(pheno) = pheno$SampleID
        
        pheno$DegradationTime = as.numeric(ss(pheno$DegradationTime,"mi"))
        
        colData(rse) <- cbind(colData(rse), pheno[rownames(colData(rse)),])
    } else if (region == 'MPFC') {
        ## Based on /dcl01/lieber/ajaffe/lab/degradation_experiments/mPFC_RiboZero/make_ERs_stranded.R
        pd =  as.data.frame(read_excel("/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_mPFC_Caudate_degradation_05_12_2017.xlsx"))
        colnames(pd)[c(1,2,7,8,14)] = c("Position","ID","Region", "DegradationTime", "Dx")
        pd = pd[pd$Region == "mPFC",]
        pd$SampleID = paste0(pd$ID, "_", pd$Flowcell)
        rownames(pd) = pd$SampleID
        pd = pd[rownames(colData(rse)),]
        pd$DegradationTime = as.numeric(ss(pd$DegradationTime,"mi"))
        colData(rse) <- cbind(colData(rse), pd[,c(1,3:13)])
    }
    
    return(rse)
}, files, types, regions)


rse_merge_pair <- function(rse1, rse2) {
    if(!identical(mcols(rse1)$Symbol, mcols(rse2)$Symbol)) {
        stopifnot(nrow(rse1) == nrow(rse2))
        message(paste(Sys.time(), 'using the "Symbol" information from the first region'))
        mcols(rse2)$Symbol <- mcols(rse1)$Symbol
    }
    rses <- list(rse1, rse2)
    cols <- sapply(rses, function(x) colnames(colData(x)))
    common <- intersect(cols[[1]], cols[[2]])
    do.call(cbind, lapply(rses, function(r) {
        m <- match(common, colnames(colData(r)))
        colData(r) <- colData(r)[, m[!is.na(m)]]
        return(r)
    }))
}

rse_merge <- function(rses) {
    res <- rses[[1]]
    for(i in (1 + seq_len(length(rses) - 1))) res <- rse_merge_pair(res, rses[[i]])
    return(res)
}


rse_merge_pair_jx <- function(rse1, rse2) {
    rses <- list(rse1, rse2)
    jxn_raw <- lapply(rses, function(x) { assays(x)$counts })
    all_jxn <- unique(unlist(sapply(jxn_raw, rownames)))
    jxn <- matrix(0, nrow = length(all_jxn), ncol = sum(sapply(jxn_raw, ncol)))

    m1 <- match(all_jxn, rownames(jxn_raw[[1]]))
    jxn[!is.na(m1), seq_len(ncol(jxn_raw[[1]]))] <- jxn_raw[[1]][m1[!is.na(m1)], ]
    m2 <- match(all_jxn, rownames(jxn_raw[[2]]))
    jxn[!is.na(m2), seq_len(ncol(jxn_raw[[2]])) + ncol(jxn_raw[[1]])] <- jxn_raw[[2]][m2[!is.na(m2)], ]

    ## Merge jx range info
    jxn_gr <- rep(rowRanges(rse1)[1], nrow(jxn))
    jxn_gr[!is.na(m1)] <- rowRanges(rse1)[m1[!is.na(m1)]]
    jxn_gr[!is.na(m2)] <- rowRanges(rse2)[m2[!is.na(m2)]]
    stopifnot(sum(jxn_gr == rowRanges(rse1)[1]) == 1)

    ## Now merge the columns
    cols <- sapply(rses, function(x) colnames(colData(x)))
    common <- intersect(cols[[1]], cols[[2]])
    jxn_col <- do.call(rbind, lapply(rses, function(r) {
        m <- match(common, colnames(colData(r)))
        colData(r)[, m[!is.na(m)]]
    }))

    rse_jx <- SummarizedExperiment(assays = list(counts = jxn), rowRanges = jxn_gr, colData = jxn_col)
    # fix junction row names
    rownames(rse_jx) <- paste0(seqnames(rse_jx), ":", start(rse_jx), "-",
        end(rse_jx), "(", strand(rse_jx), ")")
}

rse_merge_jx <- function(rses) {
    res <- rses[[1]]
    for(i in (1 + seq_len(length(rses) - 1))) res <- rse_merge_pair_jx(res, rses[[i]])
    assays(res)$rp10m <- recount::getRPKM(res, 'Length')
    return(res)
}

## Combine across regions
rse_gene <- rse_merge(all[types == 'gene'])
rse_exon <- rse_merge(all[types == 'exon'])
rse_tx <- rse_merge(all[types == 'tx'])
rse_jx <- rse_merge(all[types == 'jx'])

exprs <- list(
    'gene' = assays(rse_gene)$rpkm,
    'exon' = assays(rse_exon)$rpkm,
    'jx' = assays(rse_jx)$rp10m,
    'tx' = assay(rse_tx)
)


## Identify potential cutoffs
seed <- 20180130
seeds <- seed + 0:3
names(seeds) <- names(exprs)

dir.create('pdf', showWarnings = FALSE)

cutoffs <- sapply(names(exprs), function(type) {
    
    message(type)
    pdf(paste0('pdf/suggested_expr_cutoffs_', tolower(type), '.pdf'), width = 12)
    cuts <- expression_cutoff(exprs[[type]], seed = seeds[type])
    message(paste(cuts, collapse = ' '))
    cut <- max(cuts)
    dev.off()
    
    return(cut)

})

cutoffs

means <- lapply(exprs, rowMeans)

## Add the mean expressions, whether it passes the expression cutoff
## and save the data
dir.create('rda', showWarnings = TRUE)
dir.create('rda/unfiltered', showWarnings = TRUE)

rowRanges(rse_gene)$meanExprs <- means[['gene']]
rowRanges(rse_gene)$passExprsCut <- means[['gene']] > cutoffs['gene']
save(rse_gene, file = 'rda/unfiltered/rse_gene_unfiltered.Rdata')
rse_gene <- rse_gene[rowRanges(rse_gene)$passExprsCut]
save(rse_gene, file = 'rda/rse_gene.Rdata')

rowRanges(rse_exon)$meanExprs <- means[['exon']]
rowRanges(rse_exon)$passExprsCut <- means[['exon']] > cutoffs['exon']
save(rse_exon, file = 'rda/unfiltered/rse_exon_unfiltered.Rdata')
rse_exon <- rse_exon[rowRanges(rse_exon)$passExprsCut]
save(rse_exon, file = 'rda/rse_exon.Rdata')

rowRanges(rse_jx)$meanExprs <- means[['jx']]
rowRanges(rse_jx)$passExprsCut <- means[['jx']] > cutoffs['jx']
save(rse_jx, file = 'rda/unfiltered/rse_jx_unfiltered.Rdata')
rse_jx <- rse_jx[rowRanges(rse_jx)$passExprsCut]
save(rse_jx, file = 'rda/rse_jx.Rdata')

rowRanges(rse_tx)$meanExprs <- means[['tx']]
rowRanges(rse_tx)$passExprsCut <- means[['tx']] > cutoffs['tx']
save(rse_tx, file = 'rda/unfiltered/rse_tx_unfiltered.Rdata')
rse_tx <- rse_tx[rowRanges(rse_tx)$passExprsCut]
save(rse_tx, file = 'rda/rse_tx.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

