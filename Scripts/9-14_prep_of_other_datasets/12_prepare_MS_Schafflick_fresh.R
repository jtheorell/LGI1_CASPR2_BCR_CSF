library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
library(data.table)

#This data comes from the study https://www.nature.com/articles/s41467-019-14118-w#Sec10
#Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis

#Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138266
#The CSF cells have been excluded for diversification reasons, as we
#already have a dataset with MS CSF cells in the study. 
matrixNameList <- list.files("../External/Data/MS_Schafflick/GSE138266_RAW", pattern = "matrix", 
                         full.names = TRUE)

matrixList <- lapply(matrixNameList, readMM)

#For the gene names, they should all be identical, so we will check if that is the case for two

identical(fread("../External/Data/MS_Schafflick/GSE138266_RAW/GSM4104122_MS19270_CSF_GRCh38_genes.tsv.gz"),
          fread("../External/Data/MS_Schafflick/GSE138266_RAW/GSM4104124_MS71658_CSF_GRCh38_genes.tsv.gz"))
#TRUE

#So then we combine it all into one big dataframe after giving the cells some 
#colnames and then all get these rownames. 
fileNameList <- list.files("../External/Data/MS_Schafflick/GSE138266_RAW", pattern = "matrix")

matrixListNamed <- lapply(seq_along(matrixList), function(x){
    locMat <- matrixList[[x]]
    colnames(locMat) <- paste0(substr(fileNameList, 1, 10)[x],
                               "_cell_", seq(1,ncol(locMat)))
    locMat
})

matrixComined <- do.call("cbind", matrixListNamed)

rowDataMS <-  fread("../External/Data/MS_Schafflick/GSE138266_RAW/GSM4104122_MS19270_CSF_GRCh38_genes.tsv.gz",
                    header = FALSE)

row.names(matrixComined) <- rowDataMS$V1

msSchafSce<- SingleCellExperiment(assays = list(counts = matrixComined))

rowData(msSchafSce)$original_hgnc <- rowDataMS$V2

msSchafSce  <- getBMFeatureAnnos(msSchafSce, filters = "ensembl_gene_id", 
                             attributes = c("ensembl_gene_id", "hgnc_symbol",
                                            "chromosome_name", "gene_biotype", 
                                            "start_position", "end_position"), 
                             dataset = "hsapiens_gene_ensembl")

dir.create("Data/Comp_to_others/MS_Schafflick/SingleCellExpFiles", recursive = TRUE)
saveRDS(msSchafSce, "../External/Data/MS_Schafflick/0_raw.rds")

############
#EXCLUSIONS
###########
mito <- which(rowData(msSchafSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(msSchafSce, subsets=list(Mt=mito))

#At this stage, we need to remove all completely empty droplets (sum 0), that make up a very considerable
#portion of this dataset (about 82%)
msSchafSceReal <- msSchafSce[,-which(stats$sum == 0)]
stats <- stats[-which(stats$sum == 0),]

colData(msSchafSceReal) <- stats

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

msSchafSceReal$discard <- reasons@listData$discard
colSums(as.matrix(reasons))

#low_lib_size        low_n_features high_subsets_Mt_percent           discard 
#0                   0                     38443                    38443

msSchafSceReal$donor <- substr(colnames(msSchafSceReal), 1, 10)

#Here we make a subset to plot
msSchafSceSub <- msSchafSceReal[,sample(1:ncol(msSchafSceReal), 10000)]

dir.create("Diagnostics/Comp_to_others/MS_Schafflick", recursive = TRUE)

gridExtra::grid.arrange(
    plotColData(msSchafSceSub, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(msSchafSceSub, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(msSchafSceSub, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/MS_Schafflick/Exclusion_metrics_per_donor_automatic.pdf', height = 5, width = 15)
dev.off()

#So this very clearly did not work. Insteat, we are going to introduce manual
#thresholds. 
lowLibSize <- fewFeatures <- highMito <- rep(FALSE, nrow(stats))
lowLibSize[which(stats$sum < 1000)] <- TRUE

fewFeatures[which(stats$detected < 500)] <- TRUE

highMito[which(stats$subsets_Mt_percent > 5)] <- TRUE

reasonsReal <- data.frame(lowLibSize, fewFeatures, highMito)
reasonsReal$discard <- sapply(seq_len(nrow(reasonsReal)), function(x){
    any(reasonsReal[x,1], reasonsReal[x,2], reasonsReal[x,3])
})
colSums(as.matrix(reasonsReal))
#lowLibSize fewFeatures    highMito     discard 
#103536      104197        4076      104818 
msSchafSceReal$discard <- FALSE
msSchafSceReal$discard <- reasonsReal$discard
msSchafSceSub <- msSchafSceReal[,sample(1:ncol(msSchafSceReal), 10000)]

gridExtra::grid.arrange(
    plotColData(msSchafSceSub, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(msSchafSceSub, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(msSchafSceSub, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/MS_Schafflick/Exclusion_metrics_per_donor_real.pdf', height = 5, width = 15)
dev.off()

#104818  excluded here, allmost all from one donor. 
msSchafSce<- msSchafSceReal[,-which(msSchafSceReal$discard)]
saveRDS(msSchafSce, "../External/Data/MS_Schafflick/1_lowQual_excluded.rds")

msSchafSce<- computeSumFactors(msSchafSce)
msSchafSce<- logNormCounts(msSchafSce)

assayNames(msSchafSce)

#Here, a group column is introduced
groupCode <- read.csv("../External/Data/MS_Schafflick/Group_code.csv")

msSchafSce$Group <- FALSE
for(i in groupCode$Four_last){
    msSchafSce$Group[which(substr(colnames(msSchafSce), 7,10) == i)] <- 
        groupCode$Group[which(groupCode$Four_last == i)]
}

table(msSchafSce$Group)
#IIH    MS 
#15014 18409 

#And now over to some singler analysis, to define the cell types. Here, we use the
#reference dataset created with the influensa data that like our dataset has surface
#protein information. 
library(SingleR)

ref <- readRDS("Data/Comp_to_others/Influensa/SingleCellExpFiles/4_cell_types_for_SingleR.rds")
singlerSce <- msSchafSce[which(row.names(msSchafSce) %in% row.names(ref)),]
refSce <- ref[which(row.names(ref) %in% row.names(msSchafSce)),]
refSceOrdered <- refSce[order(row.names(refSce)),]
singlerSceOrdered <- singlerSce[order(rownames(singlerSce)),]

identical(row.names(singlerSceOrdered), row.names(refSceOrdered))
#TRUE

#So now over to the fun! This turns into a bit of a marathon though, as it takes
#ages for this to run. 
pred <- SingleR(test=singlerSceOrdered, ref=refSceOrdered, labels=refSceOrdered$cellType)

table(pred$labels, singlerSceOrdered$Group)
#           IIH   MS
#ASC          2  473
#B           60  442
#CD4T      2974 4702
#CD8Naive  7173 7619
#CD8T      1939 3623
#Myeloid   2424 1056
#NK           1  137
#pDC        114  223
#Precursor  327  134

#Once again reinforcing the fact that healthy people do not have ASC in their CSF. 

plotScoreHeatmap(pred)
dev.copy(pdf,'Diagnostics/Comp_to_others/MS_Schafflick/Cell_type_prediction_SingleR.pdf', height = 5, width = 15)
dev.off()

identical(row.names(pred), colnames(msSchafSce))
#TRUE

msSchafSce$SingleR <- pred$labels
#Now, we exclude the rest of the cells
saveRDS(msSchafSce, "../External/Data/MS_Schafflick/2_norm_singler.rds")

#Now, we will head straight for the singler analysis. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- msSchafSce [which(rowData(msSchafSce )$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
table(pred$pruned.labels)
saveRDS(pred, "Data/Comp_to_others/MS_Schafflick/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/MS_Schafflick/MonacoSingler.rds")
table(pred$pruned.labels, msSchafSce$SingleR)

#                               ASC    B CD4T CD8Naive CD8T Myeloid   NK  pDC Precursor
#Central memory CD8 T cells       1    1  274     1511  295       0    0    0         3
#Classical monocytes              0    1    0        0    0     245    0    0         2
#Effector memory CD8 T cells      4    1   53      225  711       0    0    0         1
#Exhausted B cells               46   41    0        0    0       0    0    0         0
#Follicular helper T cells        0    0  460      800    0       0    0    0         4
#Intermediate monocytes           1   10   10       19    6    1465    0    0       131
#MAIT cells                       6    1  368      684  243       0    0    0         2
#Myeloid dendritic cells         21   17   12       15    3    1735    0    3       215
#Naive B cells                    0   54    0        0    0       0    0    0         0
#Naive CD4 T cells                0    0   37      201    0       0    0    0         0
#Naive CD8 T cells                0    0    5       78    0       0    0    0         0
#Natural killer cells            22    9   31      148  388       5  135    1        17
#Non classical monocytes          0    0    0        0    0       9    0    0         4
#Non-switched memory B cells      1  113    0        2    0       0    0    0         4
#Non-Vd2 gd T cells              16    7   47      219  528       0    2    0         3
#Plasmablasts                   181    0    0        0    0       0    0    0         0
#Plasmacytoid dendritic cells     1    0    0        1    0       8    0  322         0
#Progenitor cells                 0    0    0        2    0       0    0    0         8
#Switched memory B cells        128  232    6        2    4       0    0    0         2
#T regulatory cells               1    2  718      761    9       0    1    0         9
#Terminal effector CD4 T cells    1    0  100      142  172       0    0    0         0
#Terminal effector CD8 T cells    0    1   15      171  408       0    0    0         0
#Th1 cells                        2    0 1899     2895  177       0    0    0        13
#Th1/Th17 cells                   0    0 1747     2758  118       0    0    0         6
#Th17 cells                       0    0  587      906    4       0    0    0         9
#Th2 cells                        2    2  528      956    2       0    0    0         3
#Vd2 gd T cells                  22    5  771     2284 2493       0    0    0         9

#As the predictions are less certain for the CD8 T cells, we will go only for the
#ones with high fidelity
msSchafCD8 <- msSchafSce[,which(grepl("CD8", msSchafSce$SingleR) & pred$pruned.labels == "Central memory CD8 T cells")]
saveRDS(msSchafCD8, "Data/Comp_to_others/MS_Schafflick/SingleCellExpFiles/3_CD8_CM.rds")

#We will also export the plasmablasts here, that are still 181 in number. 
msScePB <- msSchafSce [,which(msSchafSce$SingleR == "ASC" & pred$pruned.labels == "Plasmablasts")]
table(msScePB$Group)
#All 181 cells are from the MS group. 
saveRDS(msScePB, "Data/Comp_to_others/MS_Schafflick/SingleCellExpFiles/3_PB.rds")





