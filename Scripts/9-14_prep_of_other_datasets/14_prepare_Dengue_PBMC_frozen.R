library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)

#This time, the data comes from the study: 
#CD27hi CD38hi plasmablasts are activated B cells of mixed origin with distinct function
#Angeline Rouers, 1,2,9 Ramapraba Appanna, 1,9 Marion Chevrier,1 Josephine Lum, 1 Mai Chan Lau, 1 Lingqiao Tan, 1
#Thomas Loy, 1,2 Alicia Tay, 1 Raman Sethi, 1 Durgalakshmi Sathiakumar, 1 Kaval Kaur, 1 Julia Bo ̈ hme, 1
#Yee-Sin Leo, 3,4,5,6,7,8 Laurent Renia, 1,2 Shanshan W. Howland,1 Amit Singhal, 1,2 Jinmiao Chen, 1,10
#and Katja Fink 1,10,11,*
#iScience, 2021

#The dataset was downloaded from: 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172176

#This is the SmartSeq2 files of sorted CD27hiCD38hi plasmablasts from patients with primary or secondary
#dengue fever. 
rawMat <- read.delim("../External/Data/Dengue/GSE172176_TPM_E1481DK2_LNA002_006_018_036_PB_890_cells_from_MCLAU.txt", 
                     row.names = 1)

#We have been able to retrieve the metadata directly from the authors: 
metaData <- read.csv("Data/Comp_to_others/Dengue/Metainfo.csv", row.names = 1)

metaDataOrdered <- metaData[match(colnames(rawMat), metaData$SampleID),]

identical(metaDataOrdered$SampleID, colnames(rawMat))
#TRUE

sparseMat <- Matrix(as.matrix(rawMat))

pbSce <- SingleCellExperiment(assays = list("tpm" = sparseMat), 
                              colData = metaDataOrdered)
pbSce$donor <- pbSce$Patient

#For this to work as expected, we need to remove the Ensembl ID version number. 
row.names(pbSce) <- gsub("|\\..+", "", row.names(pbSce))
#And now, we do a conventional analysis on this to see if we find any ASC. 
pbSce  <- getBMFeatureAnnos(pbSce, filters = "ensembl_gene_id", 
                              attributes = c("ensembl_gene_id", "hgnc_symbol",
                                             "chromosome_name", "gene_biotype", 
                                             "start_position", "end_position"), 
                              dataset = "hsapiens_gene_ensembl")

saveRDS(pbSce, "../External/Data/Dengue/0_raw.rds")
#pbSce <- readRDS("Data/Comp_to_others/Dengue/SingleCellExpFiles/0_raw.rds")

mito <- which(rowData(pbSce)$chromosome_name=="MT")

#For this to work, we need to have a "counts" slot. Therefore, we will create 
#this in a separate sce
metricSce <- pbSce
counts(metricSce) <- tpm(metricSce)
stats <- perCellQCMetrics(metricSce, subsets=list(Mt=mito))

colData(pbSce) <- cbind(colData(pbSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

pbSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))
#low_lib_size          low_n_features high_subsets_Mt_percent        discard 
#3                       2                      43                      48
#Here, we investigate if the possible "donors" that are non-verified, seem 
#to separate into batches. 

dir.create("Diagnostics/Comp_to_others/Dengue/SCE", recursive = TRUE)

gridExtra::grid.arrange(
    plotColData(pbSce, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(pbSce, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(pbSce, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/Dengue/SCE/Exclusion_metrics_per_donor.pdf', height = 5, width = 20)
dev.off()

#There are possible batch effects here between the donors, but they will not be
#adressed, as we do not know if the numbers identified corresponds to separate donors.
#48 cells are excluded here. 
pbSce <- pbSce[,-which(pbSce$discard)]
saveRDS(pbSce, "../External/Data/Dengue/1_lowQual_excluded.rds")

#And now, we will run qumi here, as we do not have any count data for this
#experiment, and therefore cannot compare the data to the other datasets otherwise.
pbSceQumiMat <- quminorm(tpm(pbSce), shape = 3, mc.cores = 7)

#And this is now, slightly deciveinly, turned into a counts slot
counts(pbSce) <- pbSceQumiMat

pbSce <- computeSumFactors(pbSce) 
pbSce <- logNormCounts(pbSce)

dir.create("Data/Comp_to_others/Dengue/SingleCellExpFiles")
saveRDS(pbSce, "Data/Comp_to_others/Dengue/SingleCellExpFiles/PB.rds")




