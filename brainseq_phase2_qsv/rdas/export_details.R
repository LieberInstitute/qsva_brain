library('purrr')
library('knitr')
path <- getwd()
ls_cmd <- 'ls brainseq_phase2_qsvs_age17_noHGold*Rdata dxStats_*_filtered_qSVA_noHGoldQSV_match*.rda /dcl01/lieber/ajaffe/lab/brainseq_phase2/casecontrolint/rda/limma_casecontrol*.Rdata dxStats_*_filtered_qSVA_geneLevel_noHGoldQSV_match*.rda'
ls_files <- system(ls_cmd, intern = TRUE)

ls_files <- map_chr(ls_files, ~ ifelse(grepl('^\\/', .x), .x, file.path(path, .x)))


file_details <- function(x) {
    y <- load(x)
    
    cat('#### `', basename(x), '`\n\n', sep = '')
    cat('* JHPCE path: `', x, '`\n', sep = '')
    cat('* [Script](TODO)\n')
    cat('* Contents:\n')
    cat(print(kable(map_dfr(y, ~ data.frame(object = .x, class = class(get(.x)), description = 'TODO', stringsAsFactors = FALSE)), format = 'markdown')))
    cat('\n* Details:\n\n```\n')
    
    rm(x, y)
    cat(print(ls.str()))
    cat('\n```\n\n')
    
    return(NULL)
}

sink('file_details.md')
walk(ls_files, file_details)
sink()
