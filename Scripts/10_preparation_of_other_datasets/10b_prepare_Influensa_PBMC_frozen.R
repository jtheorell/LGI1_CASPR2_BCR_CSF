library(scRNAseq)
library(SingleCellExperiment)
library(scran)
library(scater)
library(DropletUtils)
#This dataset is from the study https://www.nature.com/articles/s41591-020-0769-8
#Broad immune activation underlies shared set point signatures for vaccine 
#responsiveness in healthy individuals and disease activity in patients with lupus
#By cotliarov et al. 
influensaSce <- KotliarovPBMCData(mode = c("rna", "adt"), ensembl = TRUE, location = TRUE)

#Here, we are going to use the protein data to define the cells according to our 
#definitions. But first the normal stuff. 
influensaSce  <- getBMFeatureAnnos(influensaSce, filters = "ensembl_gene_id", 
                               attributes = c("ensembl_gene_id", "hgnc_symbol",
                                              "chromosome_name", "gene_biotype", 
                                              "start_position", "end_position",
                                              "transcript_length"), 
                               dataset = "hsapiens_gene_ensembl")

saveRDS(influensaSce, "../External/Data/Influensa/Influensa/0_raw.rds")

mito <- which(rowData(influensaSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(influensaSce, subsets=list(Mt=mito))

colData(influensaSce) <- cbind(colData(influensaSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

colSums(as.matrix(reasons))
#low_lib_size          low_n_features high_subsets_Mt_percent                 discard 
#4392                    4448                    4896                    5891

#So 10% are lost. Reasonable. 
influensaSce$discard <- reasons@listData$discard

dir.create("Diagnostics/Comp_to_others/Influensa", recursive = TRUE)

influensaSce$donor <- gsub(".+_.....|", "", colnames(influensaSce))

gridExtra::grid.arrange(
    plotColData(influensaSce, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(influensaSce, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(influensaSce, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/Influensa/Exclusion_metrics_per_donor.pdf', height = 5, width = 20)
dev.off()

influensaSceRed <- influensaSce[,-which(influensaSce$discard)]

#NORMALIZATION
influensaSce <- computeSumFactors(influensaSce) 
influensaSce <- logNormCounts(influensaSce)

#Now, we exclude the rest of the cells
saveRDS(influensaSce, "../External/Data/Influensa/1_logNorm.rds")

############
#Creation of dataset for secondary definition of cell types
############
#First, we remove all cells that lack expression of any 
controls <- grep("IgG", rownames(altExp(influensaSce))) 
qc.stats <- cleanTagCounts(altExp(influensaSce), controls=controls)
summary(qc.stats$zero.ambient) # libraries removed with no ambient contamination
#Mode   FALSE    TRUE 
#logical   58625      29 
summary(qc.stats$high.controls)
#Mode   FALSE    TRUE 
#logical   58439     215

hist(log10(qc.stats$sum.controls + 1), col='grey', breaks=50,
     main="", xlab="Log-total count for controls per cell")

thresholds <- attr(qc.stats$high.controls, "thresholds")
abline(v=log10(thresholds["higher"]+1), col="red", lty=2)
dev.copy(pdf,'Diagnostics/Comp_to_others/Influensa/ADT_control_exclusion_threshold.pdf', height = 5, width = 6)
dev.off()

discard <- qc.stats$zero.ambient
discard[which(qc.stats$high.controls)] <- TRUE
influensaSce <- influensaSce[,-which(discard)]

#Now, here comes the normalisation in a more standardized version than what I have used before. 
baseline <- ambientProfileBimodal(altExp(influensaSce))

sf.amb <- medianSizeFactors(altExp(influensaSce), reference=baseline)
summary(sf.amb)
length(which(sf.amb == 0))
#0
sizeFactors(altExp(influensaSce)) <- sf.amb

#This one below still generates a lot of zero values, probablye due to a lack 
#of control data. So we will stick to the one above. 
#sf.control <- librarySizeFactors(altExp(influensaSce), subset_row=controls) 
#
#sizeFactors(altExp(influensaSce)) <- sf.control

altExp(influensaSce) <- logNormCounts(altExp(influensaSce))

# Checking that we have normalized values:
assayNames(altExp(influensaSce))

saveRDS(influensaSce, "../External/Data/Influensa/2_ADT_logNorm.rds")

#Now, we cluster the data. 
set.seed(101)
altExp(influensaSce) <- runPCA(altExp(influensaSce), ncomponents=25)
altExp(influensaSce) <- runUMAP(altExp(influensaSce), dimred = "PCA")

reducedDimNames(altExp(influensaSce))

#Now, we need to export all markers to orient ourselves and identify the plamablasts. 
dir.create("Diagnostics/Comp_to_others/Influensa/Marker_UMAPs")
lapply(row.names(altExp(influensaSce)), function(x){
    plotUMAP(altExp(influensaSce), colour_by=x)
    ggsave(paste0("Diagnostics/Comp_to_others/Influensa/Marker_UMAPs/UMAP_of_", x, ".pdf"), 
           height = 5, width = 6)
})

#What we will do first is to create a vector of major cell types, including of course
#the ASC, and then save this to be used to separate the cells in the pure transcriptomic
#datasets.


#So this clearly separates out the B-lineage cells from the rest. We start by
#removing all the others. 
inflUmap <- reducedDim(altExp(influensaSce), "UMAP")
cellTypes <- rep("CD8T", ncol(influensaSce))
cellTypes[which(inflUmap[,1] < -5 & inflUmap[,2] > 0)] <- "B"
cellTypes[which(inflUmap[,1] > 5 & inflUmap[,2] < -10)] <- "pDC"
cellTypes[which(inflUmap[,1] < -5 & inflUmap[,2] < -7.5)] <- "Precursor"
cellTypes[which(inflUmap[,1] > -5 & inflUmap[,1] < 5 & inflUmap[,2] < -7.5)] <- "Myeloid"
cellTypes[which(inflUmap[,1] < -5 & inflUmap[,2] > -7.5 & inflUmap[,2] < 0)] <- "NK"
cellTypes[which(inflUmap[,1] > 2.5 & inflUmap[,2] > 0)] <- "CD4T"

#And now, we zoom in on the B-cell compartment to identify the ASC. 

influensaSceB <- influensaSce[,which(cellTypes == "B")]

#Now, we will check if there are any typical ASC here. 
plot(logcounts(altExp(influensaSceB))["CD20_PROT",], 
     logcounts(altExp(influensaSceB))["CD38_PROT",])

#There are indeed such cells. These are identified here as ASC.
cellTypes[which(cellTypes == "B" & ((logcounts(altExp(influensaSce))["CD20_PROT",] < 3 &
                                         logcounts(altExp(influensaSce))["CD38_PROT",] > 6) |
                                        (logcounts(altExp(influensaSce))["CD20_PROT",] < 4 &
                                             logcounts(altExp(influensaSce))["CD38_PROT",] > 7)))] <- "ASC"

#THis is added to the dataset
influensaSce$cellType <- cellTypes

#And now we plot this in two ways. First the specific plot for these cells. 
influensaSceB <- influensaSce[,which(influensaSce$cellType %in% c("B", "ASC"))]

cellTypeColor <- rep("black", ncol(influensaSceB))
cellTypeColor[which(influensaSceB$cellType == "ASC")] <- "red"
    
pdf("Diagnostics/Comp_to_others/Influensa/CD20_vs_CD38_on_B-cells.pdf")
plot(logcounts(altExp(influensaSceB))["CD20_PROT",], 
     logcounts(altExp(influensaSceB))["CD38_PROT",], col = cellTypeColor)
dev.off()

#And now the whole shebang on the umap
altExp(influensaSce)$cellType <- influensaSce$cellType
plotUMAP(altExp(influensaSce), colour_by="cellType")
ggsave(paste0("Diagnostics/Comp_to_others/Influensa/Cell_types_on_UMAP.pdf"), 
       height = 5, width = 6)

#This works nicely. So now we save these two datasets. Before, though, we will
#scale down all the non-interesting cell types and save that data for classificaiton
#of the other datasets. 
infRedvec <- unlist(sapply(unique(influensaSce$cellType), function(x){
    locVec <- which(influensaSce$cellType == x)
    if(length(locVec) > 500){
        locVec[sample(1:length(locVec), 500)]
    } else {
        locVec
    }
}))

infRed <- influensaSce[,infRedvec]

#Now, these are saved. 
saveRDS(infRed, "Data/Comp_to_others/Influensa/3_cell_types_for_SingleR.rds")

#Now, we will head straight for the singler analysis. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- influensaSce[which(rowData(influensaSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
dir.create("Data/Comp_to_others/Influensa", recursive = TRUE)

saveRDS(pred, "Data/Comp_to_others/Influensa/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/Influensa/MonacoSingler.rds")

unique(pred$pruned.labels)
#[1] "Natural killer cells"          "Naive CD8 T cells"             "Naive CD4 T cells"            
#[4] "Classical monocytes"           "Terminal effector CD8 T cells" "Th17 cells"                   
#[7] "Naive B cells"                 "Th1/Th17 cells"                "Th1 cells"                    
#[10] "Vd2 gd T cells"                "MAIT cells"                    "Intermediate monocytes"       
#[13] "Effector memory CD8 T cells"   "Central memory CD8 T cells"    "Terminal effector CD4 T cells"
#[16] NA                              "T regulatory cells"            "Th2 cells"                    
#[19] "Exhausted B cells"             "Non-switched memory B cells"   "Non classical monocytes"      
#[22] "Myeloid dendritic cells"       "Follicular helper T cells"     "Non-Vd2 gd T cells"           
#[25] "Plasmacytoid dendritic cells"  "Low-density basophils"         "Switched memory B cells"      
#[28] "Progenitor cells"              "Plasmablasts"                  "Low-density neutrophils"   

influensaSce$cellTypeMonaco <- pred$pruned.labels
influensaSceB <- influensaSce[,which(influensaSce$cellType %in% c("ASC", "B"))]
table(influensaSceB$cellTypeMonaco, influensaSceB$cellType)
#                               ASC    B
#Central memory CD8 T cells       0    3
#Classical monocytes              0    5
#Effector memory CD8 T cells      0    1
#Exhausted B cells                0  865
#Myeloid dendritic cells          0    1
#Naive B cells                    0 3688
#Naive CD4 T cells                0    4
#Natural killer cells             0    2
#Non classical monocytes          0    1
#Non-switched memory B cells      0  993
#Non-Vd2 gd T cells               0    3
#Plasmablasts                    47    5
#Progenitor cells                 0    2
#Switched memory B cells          0  516
#T regulatory cells               1    3
#Terminal effector CD8 T cells    0    1
#Th1 cells                        0    4
#Th1/Th17 cells                   0    2
#Th2 cells                        0    3
#Vd2 gd T cells                   0   16

#These two methods are very strongly overlapping when it comes to the ASC, and 
#reasonably well so also for the B-cells. 

influensaSceBMon <- influensaSce[,which(influensaSce$cellTypeMonaco %in% c("Exhausted B cells",
                                                                "Naive B cells",
                                                                "Non-switched memory B cells",
                                                                "Switched memory B cells",
                                                                "Plasmablasts"))]
table(influensaSceBMon$cellTypeMonaco, influensaSceBMon$cellType)
#                             ASC    B CD4T CD8T Myeloid   NK  pDC Precursor
#Exhausted B cells              0  865    7    5       4    1    0         0
#Naive B cells                  0 3688    0    0       3    0    0         1
#Non-switched memory B cells    0  993    1    0       0    0    0         0
#Plasmablasts                  47    5    5    4       3    0    1         0
#Switched memory B cells        0  516   11    9       0    3    0         0
#Here, it is however clear that Monaco only will contaminate the plasmablast
#compartment quite significantly. Therefore, a double method is necessary. 
influensaSceBRefined <- influensaSceB[,which(influensaSceB$cellTypeMonaco %in% 
                                                         c("Exhausted B cells",
                                                           "Naive B cells",
                                                           "Non-switched memory B cells",
                                                           "Switched memory B cells",
                                                           "Plasmablasts"))]

influensaSceBRefined <- influensaSceBRefined[,-which(influensaSceBRefined$cellTypeMonaco == "Plasmablasts" &
                                                         influensaSceBRefined$cellType != "ASC")]

saveRDS(influensaSceBRefined, "../External/Data/Influensa/4_post_monaco.rds")



