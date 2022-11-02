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

rawMSList <- lapply(list.files("Data/MS_10X/Raw_files", full.names = TRUE), 
                    fread, row.names = 1)
#Here, we get loads of warnings, as the fread function cannot handle row names. Solved below. 
rawMSSparseMat <- do.call("cbind", lapply(rawMSList, function(x){
    locSMat <- Matrix(as.matrix(x[,2:ncol(x)]))
    row.names(locSMat) <- x$V1
    locSMat
}))

#So now we can create our SingleCellExperiment
msSce <- SingleCellExperiment(assays = list("counts" = rawMSSparseMat))


msSce   <- getBMFeatureAnnos(msSce, filters = "hgnc_symbol", 
                            attributes = c("ensembl_gene_id", "hgnc_symbol",
                                           "chromosome_name", "gene_biotype", 
                                           "start_position", "end_position"), 
                            dataset = "hsapiens_gene_ensembl")

#Sadly, this dataset does not use the ensembl ids as row names, so to standardize, 
#we will lose 30% of the transcripts, as they do not have a standardized ensembl id. 
msSceRed <- msSce[-which(is.na(rowData(msSce)$ensembl_gene_id)),]
rownames(msSceRed) <- rowData(msSceRed)$ensembl_gene_id

#And we reorder these, for them to be the same order as the other datasets
msSceRed <- msSceRed[order(row.names(msSceRed)),]

dir.create("Data/MS_10X/SingleCellExpFiles")
saveRDS(msSceRed, "Data/MS_10X/SingleCellExpFiles/0_raw.rds")

msSce <- msSceRed
mito <- which(rowData(msSce)$chromosome_name=="MT")

#For this to work, we need to have a "counts" slot. Therefore, we will create 
#this in a separate sce
stats <- perCellQCMetrics(msSce, subsets=list(Mt=mito))

colData(msSce) <- cbind(colData(msSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

msSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))

#low_lib_size          low_n_features high_subsets_Mt_percent                 discard 
#0                       0                    1960                    1960
#This clearly shows that this dataset has already been cleared of most non-useful cells. 

msSce$donor <- gsub(".+-|_CSF_.+", "", colnames(msSce))

dir.create("Diagnostics/MS_10X", recursive = TRUE)

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
dev.copy(pdf,'Diagnostics/MS_10X/Exclusion_metrics_per_donor.pdf', height = 5, width = 20)
dev.off()

#There are possible batch effects here between the donors, but they will not be
#adressed, as we do not know if the numbers identified corresponds to separate donors.
#48 cells are excluded here. 
msSce <- msSce[,-which(msSce$discard)]
saveRDS(msSce, "Data/MS_10X/SingleCellExpFiles/1_lowQual_excluded.rds")

#NORMALIZATION
msSce <- computeSumFactors(msSce) 
msSce <- logNormCounts(msSce)

assayNames(msSce)

#And now over to some singler analysis, to define the cell types. Again, I would
#not feel confortable doing this if I did not believe the ASC to be so different from
#all other cell types. The original authors also used SingleR. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- msSce[which(rowData(msSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
table(pred$labels)
#THis gives us 105 plasmablasts. This will be one of our final reference dataset.
plotScoreHeatmap(pred)
dev.copy(pdf,'Diagnostics/MS_10X/Cell_type_prediction_fine_SingleR.pdf', height = 5, width = 15)
dev.off()
msSce$SingleR <- pred$labels
#Now, we exclude the rest of the cells
saveRDS(msSce, "Data/MS_10X/SingleCellExpFiles/2_norm_singler.rds")
msPbSce <- msSce[,which(msSce$SingleR == "Plasmablasts")]
saveRDS(msPbSce, "Data/MS_10X/SingleCellExpFiles/3_plasmablasts.rds")


