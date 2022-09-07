#This data is tabulated in ~/Labbet/Experimental projects/LGI1/CASPR2_CSF_scRNAseq/Results/2021_04_21_Patient_database.csv
csfSce <- readRDS("../../2020/200408_full_scRNAseq_analysis_1/SingleCellExpFiles/csfSce_7_specificity_added.rds")

BSCE <- csfSce[,which(csfSce$proteinCluster == "B_lin")]
dim(BSCE)
#19846   471

table(BSCE$AdaptiveImmuneReceptor, useNA = "always")
#BCR <NA>
#471    0

#So we have 471 complete B cells.

table(paste0(BSCE$donor, "_", BSCE$exp_number))
#0022_1 0051_2 0253_1 1124_1 1166_1 1227_1 1227_2 1284_1
#1      4      3      6    198     78     10    171
