library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
library(data.table)

#This time, the data comes from the study: 
#A pathogenic and clonally expanded B celltranscriptome in active multiple sclerosis
#Akshaya Ramesha,b,1, Ryan D. Schuberta,b,1, Ariele L. Greenfielda,b,2, Ravi Dandekara,b,2, 
#Rita Loudermilka,b,Joseph J. Sabatino Jra,b, Matthew T. Koelzerc, Edwina B. Trana,b, K
#anishka Koshala,b, Kicheol Kima,b,Anne-Katrin Pröbstela,b, Debarko Banerjia,b, U
#niversity of California, San Francisco MS-EPIC Team3, 
#Chu-Yueh Guoa,b,Ari J. Greena,b, Riley M. Bovea,b, Joseph L. DeRisid,e, 
#Jeffrey M. Gelfanda,b, Bruce A. C. Creea,b, Scott S. Zamvila,b,f,Sergio E. Baranzinia,b, 
#Stephen L. Hausera,b, and Michael R. Wilson
#PNAS 2022
#The data is from: 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133028

rawMSList <- lapply(list.files("../External/Data/MS_Ramesh/Raw_CSF_files", full.names = TRUE), 
                    fread)
#Here, we get loads of warnings, as the first column contains the row names, but
#the fread function cannot handle that. It is solved here:
rawMSSparseMatList <- lapply(rawMSList, function(x){
    locSMat <- Matrix(as.matrix(x[,2:ncol(x)]))
    row.names(locSMat) <- x$V1
    locSMat
})
fileNameList <- list.files("../External/Data/MS_Ramesh/Raw_CSF_files")

sparseMatListNamed <- lapply(seq_along(rawMSSparseMatList), function(x){
    locMat <- rawMSSparseMatList[[x]]
    colnames(locMat) <- paste0(substr(fileNameList[x],1,10), 
                               "_cell_", seq(1,ncol(locMat)))
    locMat
})

rawMSSparseMat <- do.call("cbind", sparseMatListNamed)


#So now we can create our SingleCellExperiment
msSce <- SingleCellExperiment(assays = list("counts" = rawMSSparseMat))


msSce   <- getBMFeatureAnnos(msSce, filters = "hgnc_symbol", 
                            attributes = c("ensembl_gene_id", "hgnc_symbol",
                                           "chromosome_name", "gene_biotype", 
                                           "start_position", "end_position"), 
                            dataset = "hsapiens_gene_ensembl")

#Sadly, this dataset does not use the ensembl ids as row names. It is however not
#the only dataset with this problem, so we will keep the row names as they are
#and use these names downstream.

dir.create("Data/Comp_to_others/MS_Ramesh/SingleCellExpFiles", recursive = TRUE)
saveRDS(msSce, "../External/Data/MS_Ramesh/0_raw.rds")

mito <- which(rowData(msSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(msSce, subsets=list(Mt=mito))

colData(msSce) <- cbind(colData(msSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

msSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))

#low_lib_size          low_n_features  high_subsets_Mt_percent                 discard 
#0                       0                    1870                    1870 
#This clearly shows that this dataset has already been cleared of most non-useful cells,
#apart from the clearence of the mitochondrially high cells. 

msSce$donor <- substr(colnames(msSce), 1, 10)

dir.create("Diagnostics/Comp_to_others/MS_Ramesh", recursive = TRUE)

gridExtra::grid.arrange(
    plotColData(msSce, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(msSce, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(msSce, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/MS_Ramesh/Exclusion_metrics_per_donor.pdf', height = 5, width = 18)
dev.off()

#1960 cells are excluded here. 
msSce <- msSce[,-which(msSce$discard)]
saveRDS(msSce, "../External/Data/MS_Ramesh/1_lowQual_excluded.rds")

#NORMALIZATION
msSce <- computeSumFactors(msSce) 
msSce <- logNormCounts(msSce)

assayNames(msSce)

#Here, we also add information about the group, as this is relevant downstream
groupCode <- read.csv("../External/Data/MS_Ramesh/Group_code.csv")


msSce$Group <- FALSE
for(i in groupCode$Last_four_digs){
    msSce$Group[which(substr(colnames(msSce), 7,10) == i)] <- 
        groupCode$Group[which(groupCode$Last_four_digs == i)]
}

#CiS Healthy    RRMS 
#5647    5554   36579 
#So all the FALSE is gone. The CiS will be classified as MS in this context. 

#And now over to some singler analysis, to define the cell types. Here, we use the
#reference dataset created with the influensa data that like our dataset has surface
#protein information. 
library(SingleR)

ref <- readRDS("Data/Comp_to_others/Influensa/SingleCellExpFiles/4_cell_types_for_SingleR.rds")
singlerSce <- msSce[which(rowData(msSce)$ensembl_gene_id %in% row.names(ref)),]
refSce <- ref[which(row.names(ref) %in% rowData(msSce)$ensembl_gene_id),]
refSceOrdered <- refSce[order(row.names(refSce)),]
rownames(singlerSce) <- rowData(singlerSce)$ensembl_gene_id
singlerSceOrdered <- singlerSce[order(rownames(singlerSce)),]

identical(row.names(singlerSceOrdered), row.names(refSceOrdered))
#TRUE

#So now over to the fun! This turns into a bit of a marathon though, as it takes
#ages for this to run. 
pred <- SingleR(test=singlerSceOrdered, ref=refSceOrdered, labels=refSceOrdered$cellType)
table(pred$labels)
#This gives us 110 plasmablasts. 
table(pred$labels, singlerSceOrdered$Group)
#CiS Healthy  RRMS
#ASC           5       8    97
#B            86      23   686
#CD4T       1562     789  5661
#CD8Naive   2965    2834 19831
#CD8T        956     895  6392
#Myeloid      26     860  2996
#NK           30      14    78
#pDC           3      12   340
#Precursor    14     119   498
#This shows, as expected, that the number of cells in the healthy setting that fulfill
#the ASC criteria are very few. Therefore, the non-MS ASC will be excluded. 

plotScoreHeatmap(pred)
dev.copy(pdf,'Diagnostics/Comp_to_others/MS_Ramesh/Cell_type_prediction_SingleR.pdf', height = 5, width = 15)
dev.off()
identical(colnames(msSce), row.names(pred))
msSce$SingleR <- pred$labels
#Now, we exclude the rest of the cells
saveRDS(msSce, "../External/Data/MS_Ramesh/2_norm_singler.rds")
msSceASC <- msSce[,which(msSce$SingleR == "ASC" & msSce$Group  == "RRMS")]

#Now, we will head straight for the singler analysis. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- msRameshSce[which(rowData(msRameshSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
saveRDS(pred, "Data/Comp_to_others/MS_Ramesh/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/MS_Ramesh/MonacoSingler.rds")
#Now, what is the overlap here?
table(pred$pruned.labels, msRameshSce$SingleR)
#                               ASC    B CD4T CD8Naive CD8T Myeloid   NK  pDC Precursor
#Central memory CD8 T cells       1    2  191     2494  692       1    1    0         6
#Classical monocytes              0    0    0        0    0     343    0    0         2
#Effector memory CD8 T cells      6    2   91      444 2668       1    2    0         6
#Exhausted B cells                5   70    0        0    0       0    0    0         2
#Follicular helper T cells        0    5  773     3128   10       0    0    0        55
#Intermediate monocytes           3    1    7       21    6    1436    0    1        74
#MAIT cells                       0    0  256      954  392       0    0    0        11
#Myeloid dendritic cells          2    5   16       33    5    2059    0   10       214
#Naive B cells                    0   48    0        0    0       0    0    0         0
#Naive CD4 T cells                0    0   24      388    0       0    0    0         4
#Naive CD8 T cells                0    0    1       98    0       0    0    0         0
#Natural killer cells             3    5   21      247  483       5  115    0        44
#Non classical monocytes          0    0    0        0    0      10    0    0         3
#Non-switched memory B cells      3  311    0        3    0       2    0    0         5
#Non-Vd2 gd T cells               5    0   29      207  386       0    2    0         4
#Plasmablasts                    38    0    0        0    0       0    0    0         1
#Plasmacytoid dendritic cells     1    0    0        0    0       7    0  337         0
#Progenitor cells                 0    1    0        1    0       0    0    0        11
#Switched memory B cells         13  335    3       19    0       0    0    0        10
#T regulatory cells               7    1  756     1707   33       0    0    0        27
#Terminal effector CD4 T cells    1    2  114      195  293       2    0    0         3
#Terminal effector CD8 T cells    0    1   14      125  766       0    0    0         2
#Th1 cells                        5    4 2272     5190  509       0    1    0        47
#Th1/Th17 cells                   2    1 2194     5842  260       0    0    0        36
#Th17 cells                       1    1  442     1397   12       0    0    0        22
#Th2 cells                        1    0  450     1448    7       0    0    0        25
#Vd2 gd T cells                   2    0  346     1660 1720       5    1    0        10

msScePB <- msRameshSce[,which(msRameshSce$SingleR == "ASC" & pred$pruned.labels == "Plasmablasts")]
table(msScePB$Group)
#All these 38 cells are RRMS
saveRDS(msScePB, "Data/Comp_to_others/MS_Ramesh/SingleCellExpFiles/3_PB.rds")


