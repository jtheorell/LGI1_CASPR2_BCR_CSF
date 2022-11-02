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
rawMat <- read.delim("Data/PlasmaBlast/GSE172176_TPM_E1481DK2_LNA002_006_018_036_PB_890_cells_from_MCLAU.txt", 
                     row.names = 1)

colNamesRaw <- colnames(rawMat)

#We have not been able to decipher which cells comes from which donor, even if there
#are four entities when taking hte first four signs in the strings: 
table(substr(colNamesRaw, 1,4))
#RHB3 RHB4 RHB5 RHH1 
#296  466  112   16 
#But if this is true, we will still struggle do work with statstistics, as so very
#few cells come from the fourth donor. And we only have three AE patients anyway, so we will
#only compare the cell types as they are.

sparseMat <- Matrix(as.matrix(rawMat))

pbSce <- SingleCellExperiment(assays = list("tpm" = sparseMat))

#For this to work as expected, we need to remove the Ensembl ID version number. 
row.names(pbSce) <- gsub("|\\..+", "", row.names(pbSce))
#And now, we do a conventional analysis on this to see if we find any ASC. 
pbSce  <- getBMFeatureAnnos(pbSce, filters = "ensembl_gene_id", 
                              attributes = c("ensembl_gene_id", "hgnc_symbol",
                                             "chromosome_name", "gene_biotype", 
                                             "start_position", "end_position"), 
                              dataset = "hsapiens_gene_ensembl")

dir.create("Data/PlasmaBlast/SingleCellExpFiles")
saveRDS(pbSce, "Data/PlasmaBlast/SingleCellExpFiles/0_raw.rds")
#pbSce <- readRDS("Data/PlasmaBlast/SingleCellExpFiles/0_raw.rds")

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
pbSce$donor <- substr(colNamesRaw, 1,4)

dir.create("Diagnostics/PlasmaBlast/SCE", recursive = TRUE)

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
dev.copy(pdf,'Diagnostics/PlasmaBlast/SCE/Exclusion_metrics_per_donor.pdf', height = 5, width = 20)
dev.off()

#There are possible batch effects here between the donors, but they will not be
#adressed, as we do not know if the numbers identified corresponds to separate donors.
#48 cells are excluded here. 
pbSce <- pbSce[,-which(pbSce$discard)]
saveRDS(pbSce, "Data/PlasmaBlast/SingleCellExpFiles/1_lowQual_excluded.rds")

#NORMALIZATION
#As this data is actualy tpm data, we create a second slot with that correct name now
tpm(pbSce) <- counts(pbSce)
assay(pbSce, "logcpm") <- log10(tpm(pbSce)+1)

assayNames(pbSce)
#"counts" "tpm" "logcpm"
#Ad with that, we are ready to go on to the comparisons, believe it or not!

#Here, we import the data in question. 
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

#First, we generate the comparable set of values: 
assay(csfSce, "logcpm") <- log10(assay(csfSce, "cpm")+1)
#And now, we use a smaller set of genes that are present in both datasets. Both
#are generated with GRCh38, so they should be reasonably comparable in terms of 
#transcript names. 
#First, it turns out that the transcipt name versions used in the plasmablast dataset
#are more fine-grained. We will see if it is possible to just remove some info. 
rowNamesLessSpecific <- gsub("|\\..+", "", row.names(pbSce))
#Luckliy, they are still all unique:
identical(length(unique(rowNamesLessSpecific)), length(rowNamesLessSpecific))
row.names(pbSce) <- rowNamesLessSpecific
#This means that we can use these to identify overlapping transcripts:
commonNamesCsf <- which(row.names(csfSce) %in% row.names(pbSce))
#All but 1.8% of the transcripts are present in both. 
commonNamesPb <- which(row.names(pbSce) %in% row.names(csfSce))

csfSceCommon <- csfSce[commonNamesCsf,]
pbSceCommon <- pbSce[commonNamesPb,]

identical(row.names(csfSceCommon), row.names(pbSceCommon))
#TRUE

#So no need to change any orders!
csfSceCommon$group <- "encephalitis"
pbSceCommon$group <- "dengue"
groupVec <- c(csfSceCommon$group, pbSceCommon$group)

ascMat <- cbind(assay(csfSceCommon, "logcpm"), assay(pbSceCommon, "logcpm"))

#Now, these experiments are combined. 
ascSce <- SingleCellExperiment(assays=list("logcounts" = ascMat))
ascSce$group <- groupVec
rowData(ascSce) <- rowData(csfSceCommon)

