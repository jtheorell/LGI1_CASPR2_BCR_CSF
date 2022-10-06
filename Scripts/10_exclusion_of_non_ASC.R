#Sadly, as can be viewed here and in the output from 9 as well, there is only 
#a very tiny fraction of the specific cells that are non-ASC: 
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_5_BCR.rds")
table(csfSce$Cell_type, csfSce$Specific)

#FALSE Not_tested TRUE
#ASC    31        127  101
#B      13         68    7

#These are distributed among donors 1166 (5 cells) and 1284 (2 cells).
#As the B cells are so clearly skewed towards non-specificity, we will exclude them
#from any further transcriptome analyses here, as we will otherwise get a 
#result downstream just showing that ASCs genes are overrepresented among the specifics
#which we already know now. 

csfSce <- csfSce[,-which(csfSce$Cell_type == "B")]

saveRDS(csfSce, "Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")