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

#We have excluded all non-RRMS or CIS (clinically isolated syndrome) patients/controls

rawMSList <- lapply(list.files("../External/Data/MS_Ramesh_CSF/Raw_CSF_files", full.names = TRUE), 
                    fread)
#Here, we get loads of warnings, as the first column contains the row names, but
#the fread function cannot handle that. It is solved here:
rawMSSparseMatList <- lapply(rawMSList, function(x){
    locSMat <- Matrix(as.matrix(x[,2:ncol(x)]))
    row.names(locSMat) <- x$V1
    locSMat
})
fileNameList <- list.files("../External/Data/MS_Ramesh_CSF/Raw_CSF_files")

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

saveRDS(msSce, "../External/Data/MS_Ramesh_CSF/0_raw.rds")

mito <- which(rowData(msSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(msSce, subsets=list(Mt=mito))

colData(msSce) <- cbind(colData(msSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

msSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))

#low_lib_size          low_n_features  high_subsets_Mt_percent                 discard 
#        0                       0                    1736                    1736 
#This clearly shows that this dataset has already been cleared of most non-useful cells,
#apart from the clearence of the mitochondrially high cells. 

msSce$donor <- substr(colnames(msSce), 1, 10)

dir.create("Diagnostics/Comp_to_others/MS_Ramesh_CSF", recursive = TRUE)

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
dev.copy(pdf,'Diagnostics/Comp_to_others/MS_Ramesh_CSF/Exclusion_metrics_per_donor.pdf', height = 5, width = 18)
dev.off()

#1736 cells are excluded here. 
msSce <- msSce[,-which(msSce$discard)]
saveRDS(msSce, "../External/Data/MS_Ramesh_CSF/1_lowQual_excluded.rds")

#NORMALIZATION
msSce <- computeSumFactors(msSce) 
msSce <- logNormCounts(msSce)

assayNames(msSce)

#Here, we make sure that there are no non-MS patients left: 
groupCode <- read.csv("../External/Data/MS_Ramesh_CSF/Group_code.csv")


msSce$Group <- FALSE
for(i in groupCode$Last_four_digs){
    msSce$Group[which(substr(colnames(msSce), 7,10) == i)] <- 
        groupCode$Group[which(groupCode$Last_four_digs == i)]
}

#CiS  RRMS 
#5647 3657
#So all the FALSE is gone. The CiS will be classified as MS in this context. 

#And now over to some singler analysis. We will use the Monaco dataset as a standard
#to define the cell types. 

library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- msSce[which(rowData(msSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
dir.create("Data/Comp_to_others/MS_Ramesh_CSF", recursive = TRUE)
saveRDS(pred, "Data/Comp_to_others/MS_Ramesh_CSF/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/MS_Ramesh/MonacoSingler.rds")
msSce$cellType <- pred$pruned.labels
msSceB <- msSce[,which(msSce$cellType  %in% c("Exhausted B cells",
                                                   "Naive B cells",
                                                   "Non-switched memory B cells",
                                                   "Switched memory B cells",
                                                   "Plasmablasts"))]
msSceB$cellType[grep("B", msSceB$cellType)] <- "B"
msSceB$cellType[grep("Plasmablasts", msSceB$cellType)] <- "ASC"
table(msSceB$cellType)
#ASC   B 
#39 806 

saveRDS(msSceB, "../External/Data/MS_Ramesh_CSF/2_normalised_post_monaco.rds")

#And now, we will further classify these cells using the influensa dataset, to make
#perfectly sure that we have clean groups

ref <- readRDS("Data/Comp_to_others/Influensa/3_cell_types_for_SingleR.rds")
singlerSce <- msSceB[which(rowData(msSceB)$ensembl_gene_id %in% row.names(ref)),]
refSce <- ref[which(row.names(ref) %in% rowData(msSceB)$ensembl_gene_id),]
refSceOrdered <- refSce[order(row.names(refSce)),]
rownames(singlerSce) <- rowData(singlerSce)$ensembl_gene_id
singlerSceOrdered <- singlerSce[order(rownames(singlerSce)),]

identical(row.names(singlerSceOrdered), row.names(refSceOrdered))
#TRUE

pred <- SingleR(test=singlerSceOrdered, ref=refSceOrdered, labels=refSceOrdered$cellType)

table(msSceB$cellType, pred$pruned.labels)
#    ASC   B CD4T Precursor
#ASC  38   0    0         1
#B    15 746   28        15
msSceB$cellTypeKotliarov <- pred$pruned.labels
msSceBRefined <- msSceB[,which(msSceB$cellTypeKotliarov %in% c("ASC", "B"))]
msSceBRefined <- msSceBRefined[,-which(msSceBRefined$cellTypeKotliarov == "ASC" &
                                                     msSceBRefined$cellType != "ASC" )]
table( msSceBRefined$cellType,  msSceBRefined$cellTypeKotliarov)
#    ASC   B
#ASC  38   0
#B     0 746

saveRDS(msSceBRefined, "../External/Data/MS_Ramesh_CSF/3_post_Kotliarov.rds")
