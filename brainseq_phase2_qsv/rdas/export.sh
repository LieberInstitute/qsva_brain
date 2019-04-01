#!/bin/bash

date
tar -cvzf SCZDvsControlDEGobjects.tar.gz brainseq_phase2_qsvs_age17_noHGold*Rdata dxStats_*_filtered_qSVA_noHGoldQSV_match*.rda /dcl01/lieber/ajaffe/lab/brainseq_phase2/casecontrolint/rda/limma_casecontrol*.Rdata dxStats_*_filtered_qSVA_geneLevel_noHGoldQSV_match*.rda
date
