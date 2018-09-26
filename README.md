# qsva_brain

This repository contains the code for the SCZD case vs control differential expression analysis for the BrainSeq Phase II project (see the [`brainseq_phase2` repo](https://github.com/LieberInstitute/brainseq_phase2)). It was carried out by [Amy Peterson](https://amy-peterson.github.io/) as part of her JHSPH MPH project.

## License

<img src="https://licensebuttons.net/l/by-nc/3.0/88x31.png" alt width="88" height="31" scale="0">
Attribution-NonCommercial: CC BY-NC

This license lets others remix, tweak, and build upon our work non-commercially.

[View License Deed](https://creativecommons.org/licenses/by-nc/4.0) | [View Legal Code](https://creativecommons.org/licenses/by-nc/4.0/legalcode)

## Citation

If you use anything in this repository please cite the [BrainSeq Phase II project](https://github.com/LieberInstitute/brainseq_phase2#citation).

## Scripts

The main script directories follow this order:

1. `expr_data`: scripts for integrating the degradation datasets from multiple brain regions.
2. `means`: calculate the base-pair coverage mean files for the degradation data.
3. `ERs`: identify expressed regions using `derfinder` on the degradation data. Also identify which ERs are strongly associated with degradation time.
4. `brainseq_phase2_qsv`: compute the coverage matrix for the BrainSeq Phase 2 data, identify the quality surrogate variables (qSVs) and then perform the differential expression analysis between schizophrenia cases and non-psychiatric controls.

## Files in this project

* README.md: this file
* qsva_brain.Rproj: RStudio file for organizing the project

### expr_data

* merge_data.R: R script merging the different degradation datasets.
* run_merge_data.sh: bash shell script for running the R script at JHPCE.
* logs: directory with log files
* pdf: directory with image files

### means

* qsva_bws.R: R script for computing the mean base-pair coverage and saving them as BigWig files. Uses the `recount.bwtool` R package.
* qsva_bws.sh: bash shell script for running the previous R script at JHPCE.

### ERs

* make_ERs_stranded.R: R script for identifying the expressed regions, then identifying which ERs are strongly associated with degradation signal and saving the results for later use.
* make_ERs_stranded.sh: bash shell script for running the previous R script at JHPCE.

### brainseq_phase2_qsv

* quantify_top1000.R and quantify_top1000.sh: R and bash script for computing the coverage matrix.
* explore_replicates.R: R script for exploring the coverage matrix results.
* make_qSVs.R: identify the quality surrogate variables following different options. That is, with all 900 samples, with 712 after dropping potential confounding samples, and then for each brain region individually.
* explore_qsvs.R: perform PCA on the BrainSeq gene expression data and explore that data. Also explore the qSVs and their associations with measured covariates.
* casectrl_HIPPO.R and casectrl_DLPFC.R: perform the case-control differential expression analysis for each brain region at the gene expression level.
* casectrl_HIPPO_allFeatures.R and casectrl_DLPFC_allFeatures.R: similar to the above, but for exons, exon-exon junctions and transcript expression levels.
* casectrl_HIPPO_plots.R and casectrl_DLPFC_plots.R: R scripts for making the DEqual plots that assess the performance of the qSVA framework.
* explore_case_control.R: R script for exploring the case-control results, performing gene ontology enrichment analyses, comparing results against BrainSeq Phase I (polyA) and visualizing the top results.
* pdf: directory with image files

### Length of scripts

Number of lines in each script.

```
## R scripts
$ for i in */*.R; do echo $i; wc -l $i; done
ERs/make_ERs_stranded.R
     258 ERs/make_ERs_stranded.R
brainseq_phase2_qsv/casectrl_DLPFC.R
     263 brainseq_phase2_qsv/casectrl_DLPFC.R
brainseq_phase2_qsv/casectrl_DLPFC_allFeatures.R
     164 brainseq_phase2_qsv/casectrl_DLPFC_allFeatures.R
brainseq_phase2_qsv/casectrl_DLPFC_plots.R
      56 brainseq_phase2_qsv/casectrl_DLPFC_plots.R
brainseq_phase2_qsv/casectrl_HIPPO.R
     261 brainseq_phase2_qsv/casectrl_HIPPO.R
brainseq_phase2_qsv/casectrl_HIPPO_allFeatures.R
     165 brainseq_phase2_qsv/casectrl_HIPPO_allFeatures.R
brainseq_phase2_qsv/casectrl_HIPPO_plots.R
      55 brainseq_phase2_qsv/casectrl_HIPPO_plots.R
brainseq_phase2_qsv/explore_case_control.R
     690 brainseq_phase2_qsv/explore_case_control.R
brainseq_phase2_qsv/explore_qsvs.R
    1510 brainseq_phase2_qsv/explore_qsvs.R
brainseq_phase2_qsv/explore_replicates.R
     137 brainseq_phase2_qsv/explore_replicates.R
brainseq_phase2_qsv/make_qSVs.R
     318 brainseq_phase2_qsv/make_qSVs.R
brainseq_phase2_qsv/quantify_top1000.R
      97 brainseq_phase2_qsv/quantify_top1000.R
expr_data/merge_data.R
     264 expr_data/merge_data.R
means/qsva_bws.R
      29 means/qsva_bws.R
      
## bash scripts
$ for i in */*.sh; do echo $i; wc -l $i; done
ERs/make_ERs_stranded.sh
      21 ERs/make_ERs_stranded.sh
brainseq_phase2_qsv/quantify_top1000.sh
      22 brainseq_phase2_qsv/quantify_top1000.sh
expr_data/run_merge_data.sh
      23 expr_data/run_merge_data.sh
means/qsva_bws.sh
      19 means/qsva_bws.sh
```

## LIBD internal

JHPCE location: `/dcl01/ajaffe/data/lab/qsva_brain`
