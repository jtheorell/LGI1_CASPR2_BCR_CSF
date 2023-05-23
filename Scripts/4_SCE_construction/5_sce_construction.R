library(ensembldb)
library(tximport)
library(SingleCellExperiment)
library(scran)
library(scater)
library(uwot)
library(DepecheR)
library(bluster)
library(ggplot2)
library(scry)

#I use the file Homo_sapiens.GRCh38.98.gtf.gz and combine it with the ERCC information
#In the following way.
#If not run before: 
#DB <- ensDbFromGtf(gtf="~/Labbet/Saved_transcriptome_files/HomoS_GRCh38.98_cdna_ncrna_ERCC/Homo_sapiens.GRCh38.98.gtf.gz")

EDB <- EnsDb("~/Labbet/Saved_transcriptome_files/HomoS_GRCh38.98_cdna_ncrna_ERCC/Homo_sapiens.GRCh38.98.sqlite")
tx2geneId <- transcripts(EDB, columns=c("tx_id", "gene_id"), return.type="DataFrame")
#Now add the ERCC spike-ins
ERCCgtf <- read.table('~/Labbet/Saved_transcriptome_files/HomoS_GRCh38.98_cdna_ncrna_ERCC/ERCC92.gtf.gz', 
                      header = FALSE, sep = '\t')
ERCCNamed <- data.frame("tx_id"=as.character(ERCCgtf[,1]), 
                        "gene_id"=as.character(ERCCgtf[,1]), 
                        stringsAsFactors = FALSE)
tx2geneComplete <- rbind(tx2geneId, ERCCNamed)

#And now we import this.
dirs <- list.files("Salmon_results")
fullDirs <- file.path("Salmon_results", dirs, "quant.sf") 

txiFileGenes <- tximport(files=fullDirs, type = "salmon",  
                         tx2gene = tx2geneComplete,
                         ignoreTxVersion=TRUE)

dirs <- list.files("Salmon_results")
fullDirs <- file.path("Salmon_results", dirs, "quant.sf") 

txiFileGenes <- tximport(files=fullDirs, type = "salmon",  
                         tx2gene = tx2geneComplete,
                         ignoreTxVersion=TRUE)

csfSce <- SingleCellExperiment(assays = 
                                   list(counts = txiFileGenes$counts, 
                                        cpm = txiFileGenes$abundance))

colnames(csfSce) <- dirs

#Add some meta information embedded in the cell names
csfSce$donor <- gsub("|_._._...", "\\1", colnames(csfSce))
csfSce$exp_number <- gsub("...._|_._...", "\\1", colnames(csfSce))
csfSce$plate <- gsub("...._._|_...", "\\1", colnames(csfSce))
csfSce$plate_position <- gsub("...._._._|", "\\1", colnames(csfSce))

#Here, rowData is added, to get information not only about the gene symbol, but
#also the chromosome position, gene biotype, etc. 
csfSce <- getBMFeatureAnnos(csfSce, filters = "ensembl_gene_id", 
                            attributes = c("ensembl_gene_id", "hgnc_symbol",
                                           "chromosome_name", "gene_biotype", 
                                           "start_position", "end_position"), 
                            dataset = "hsapiens_gene_ensembl")

#Now, the ERCC are defined as such
is.spike <- grepl("ERCC", rownames(csfSce))
csfSce <- splitAltExps(csfSce, ifelse(is.spike,"ERCC", "gene"), ref = "gene")

#csfSce <- readRDS("../External/Data/csfSce_1b_preQC_1166_1227_1284_exp_1.rds")

#############
#Inclusion of stats
#############
#Identification of mitochondrial genes
mito <- which(rowData(csfSce)$chromosome_name=="MT")

#Now the preCellQCMetrcis are calculated. 
stats <- perCellQCMetrics(csfSce, subsets=list(Mt=mito))
colData(csfSce) <- cbind(colData(csfSce), stats)

#############
#Inclusion of BCR data
#############

BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")
BCR_H <- BCR_all[which(BCR_all$LOCUS == "H"),]

#Now, we will create dummy rows for all that are lacking a BCR. 
dummyRows <- colnames(csfSce)[-which(colnames(csfSce) %in% BCR_H$CELL)]
dummyCols <- colnames(BCR_H)

dummyMat <- as.data.frame(matrix(NA, nrow = length(dummyRows),
                   ncol = length(dummyCols), 
                   dimnames = list(dummyRows, dummyCols)))
dummyMat$CELL <- dummyRows
fullBCRMat <- rbind(BCR_H, dummyMat)

fullBCRMatOrd <- fullBCRMat[match(colnames(csfSce), fullBCRMat$CELL),]

identical(fullBCRMatOrd$CELL, colnames(csfSce))
#TRUE
#So now we can add this to the csfSce
colData(csfSce) <- cbind(colData(csfSce), fullBCRMatOrd[,c("V_CALL","JUNCTION",
                                                                 "JUNCTION_LENGTH",
                                                                 "ISOTYPE","HamClone",
                                                                 "light_type","Clonal",
                                                                 "All_mutations",
                                                                 "Non_silent_mutations",
                                                                 "Specific", "LRR", "EPTP",
                                                                 "Specific_UCA",
                                                                 "SubClone")])
#And for the sake of it: we will change all the NAs in the Specific column to "not tested"
csfSce$Specific[which(is.na(csfSce$Specific))] <- "Not_tested"

#############
#Inclusion of cytometry data
#############
#To be able to do this, we need to exclude all A1 wells. 
csfSce <- csfSce[,-grep("A01", colnames(csfSce))]
#Here, the flow data is integrated into the SCE
flowData <- read.csv("Data/Cytometry/flowDataPlusIndex.csv")

flowDataHere <- flowData[which(flowData$Cell %in% colnames(csfSce)),]
#Cells that have been excluded already at the flow step should of course
#not be included here.
table(csfSce$Specific[-which(colnames(csfSce) %in% flowDataHere$Cell)])
#Not_tested 
#         30
csfSce <- csfSce[,which(colnames(csfSce) %in% flowDataHere$Cell)]

flowDataHereOrdered <- flowDataHere[match(colnames(csfSce),flowDataHere$Cell),]

flowDataDf <- flowDataHereOrdered[,c("FSC.A", "SSC.A", "CD3", "CD4", "CD8",
                                     "CD16", "CD19", "CD20", "CD27", "CD38",
                                     "CD56", "CD138", "IgD")]

flowDataMat <- t(flowDataDf)

colnames(flowDataMat) <- flowDataHereOrdered$Cell

identical(colnames(flowDataMat), colnames(csfSce))
#TRUE. And with that we are ready to integrate this. 

flowDataSce <- 
    SingleCellExperiment(assays = list("normcounts" = flowDataMat))
altExp(csfSce, "flowData") <- flowDataSce

###############
#GLM-PCA.
###############
#Now, we generate a GLM-PCA. 
csfSceTownes <- devianceFeatureSelection(csfSce, assay="counts", sorted=TRUE)

dir.create("Diagnostics/GLMPCA")

pdf("Diagnostics/GLMPCA/Binomial_deviance_selection.pdf")
plot(rowData(csfSceTownes)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")
dev.off()

#We select these 2000 genes. 
selectSCE <- csfSceTownes[1:2000,]

Sys.time()
set.seed(101)
glmpcaSCE <- GLMPCA(selectSCE, L=30, assay="counts")
Sys.time()

#Four minutes. 
saveRDS(glmpcaSCE, "../External/Data/GLMPCA_model.rds")
#glmpcaSCE <- readRDS("../External/Data/GLMPCA_model.rds")

#And this is added back: 
reducedDim(csfSce, "GLMPCA") <- reducedDim(glmpcaSCE)

#############
#Cell type definitions
#############

#We will here use the nice euclidean distance method for selection of the gates. 
#The euclidean space that we will use is the transcriptomic one, and we will
#use the protein to guide us. 

#One of the first things we need to do here, is to make sure that we do not
#lose any ASC because CD138 leaks into CD3. 
flowDat <- as.data.frame(t(normcounts(altExp(csfSce, "flowData"))))
flowDat$Specific <- csfSce$Specific
flowDat$sum <- csfSce$sum
dir.create("Diagnostics/Cell_type_definitions")
ggplot(flowDat, aes(x = CD3, y = CD138, color = Specific)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Cell_type_definitions/CD3_vs_CD138_with_specificity.pdf")

#And here, the definitions start, as the flow plot gives great guidance. 
cellType <- rep("other", nrow(flowDat))
cellType[which(flowDat$CD3 < 2.5 & flowDat$CD138 > 2.5)] <- "ASC"
cellType[which(flowDat$CD3 > 2.5 & flowDat$CD138 > 2.5)] <- "Doublet"
cellType[which(flowDat$CD3 > 1.25 & flowDat$CD138 < 2.5)] <- "T.cell"

ggplot(flowDat[which(cellType == "other"),], aes(x = CD3, y = CD19, color = Specific)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Cell_type_definitions/CD3_vs_CD19_on_non-ASC_non_T_with_specificity.pdf")

cellType[which(cellType == "other" & flowDat$CD19 > 1.5)] <- "B"

ggplot(flowDat[which(cellType == "B"),], aes(x = CD38, y = CD20, color = Specific)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Cell_type_definitions/CD38_vs_CD20_on_B_with_specificity.pdf")

cellType[which(cellType == "B" & flowDat$CD20 < 1 & flowDat$CD38 > 2)] <- "ASC"

ggplot(flowDat[which(cellType == "B"),], aes(x = CD138, y = CD20, color = Specific)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Cell_type_definitions/CD138_vs_CD20_on_B_with_specificity.pdf")

cellType[which(cellType == "B" & flowDat$CD20 < 2 & flowDat$CD138 > 0.5)] <- "ASC"

csfSce$cellType <- cellType

#Complete_plus_cell_type
saveRDS(csfSce, "../External/Data/csfSce_1c_flow_included_preQC.rds")
#csfSce <- readRDS("../External/Data/csfSce_1c_flow_included_preQC.rds")


#############
#Cell exclusions based on qc params
#############

#There are considerable differences in the transcriptomes
#of the B-cells and the ASC, to no-ones surprise. This means that there are
#reasons to put the thresholds very high, to avoid losing the B-cells alltogheter. 
#This is the most clear for ERCC, where especially 1166 has up to 20% ERCC in 
#the B-cell compartment, but also has a higher level for the ASC. Another
#reason that this is acceptable is that the cells in here of course have been
#heavily pre-selected. 
lowCount <- fewFeatures <- highMito <- highERCC <- rep(FALSE, ncol(csfSce))

lowCount[which(csfSce$sum < 100000)] <- TRUE
fewFeatures[which(csfSce$detected < 1000)] <- TRUE
highMito[which(csfSce$subsets_Mt_percent > 10)] <- TRUE
#For the ERCC, there is a considerably larger population of cells with high percentages 
#in the JR1166 donor. 
highERCC[which(csfSce$altexps_ERCC_percent > 30)] <- TRUE

csfSce$lowQual <- FALSE
csfSce$lowQual[which(lowCount | fewFeatures | highMito | highERCC)] <- TRUE

#This would lead to the exclusion of 15 cells
dir.create("Diagnostics/SCE")

gridExtra::grid.arrange(
    plotColData(csfSce, x="donor", 
                y="sum", colour_by = "lowQual") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(csfSce, x="donor", 
                y="detected", colour_by = "lowQual") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(csfSce, x="donor", 
                y="subsets_Mt_percent", colour_by = "lowQual") + 
        ggtitle("Mito percent"),
    plotColData(csfSce, x="donor", 
                y="altexps_ERCC_percent", colour_by = "lowQual") + 
        ggtitle("ERCC percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/SCE/Exclusion_metrics_per_donor.pdf', height = 5, width = 10)
dev.off()

#We also need to investigate this from a specificity point of view, as each cell counts
gridExtra::grid.arrange(
    plotColData(csfSce, x="Specific", 
                y="sum", colour_by = "lowQual") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(csfSce, x="Specific", 
                y="detected", colour_by = "lowQual") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(csfSce, x="Specific", 
                y="subsets_Mt_percent", colour_by = "lowQual") + 
        ggtitle("Mito percent"),
    plotColData(csfSce, x="Specific", 
                y="altexps_ERCC_percent", colour_by = "lowQual") + 
        ggtitle("ERCC percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/SCE/Exclusion_metrics_per_specificity.pdf', height = 5, width = 10)
dev.off()

#Looks entirely reasonable
#How do the excluded distribute over the specificity groups?
table(csfSce$Specific[which(csfSce$lowQual)])
#FALSE Not_tested       TRUE 
#.   1         152          4 

csfSce <- csfSce[,!csfSce$lowQual]

#Here, 157 cells, or 9%, are lost. 
#We are left with 1674 cells.
#We now exclude all non-B-lineage cells
csfSce <- csfSce[,which(csfSce$cellType %in% c("ASC", "B"))]

#Here, we also export a file with the B-lineage cell names, to be used for RNA velocity
write.csv(colnames(csfSce), 
          "../External/Data/All_full_cells.csv", row.names = FALSE)

csfSce <- computeSumFactors(csfSce) 
csfSce <- logNormCounts(csfSce)
#As we are using the monaco data downstream for the other datasets, it is good to show 
#the overlap here: 

library(SingleR)
ref <- MonacoImmuneData()
singlerSce <- csfSce[which(rowData(csfSce)$hgnc_symbol %in% row.names(ref)),]
rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
table(pred$pruned.labels, csfSce$cellType)


#                            ASC   B
#Exhausted B cells            16 154
#Naive B cells                 0   7
#Non-switched memory B cells   0  16
#Plasmablasts                251   0
#Switched memory B cells       1  17

#So the 16 cells that are classified as exhausted and at the same time as ASC 
#need to be further scrutinised
flowDat <- as.data.frame(t(normcounts(altExp(csfSce, "flowData"))))
flowDat$Monaco <- pred$pruned.labels
flowDat$Mismatch <- FALSE
flowDat$Mismatch[which(pred$pruned.labels == "Exhausted B cells" &
                           csfSce$cellType == "ASC")] <- TRUE
ggplot(flowDat, aes(x = CD138, y = CD20, color = Mismatch)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Cell_type_definitions/CD138_vs_CD20_on_B_with_Monaco_cellType_mismatch.pdf")

#This shows that all the cells that are registered as exhausted but that are classified
#as ASC express high levels of CD138. Three of them do however also express CD20, 
#but this is considered of secondary importance. 

#Here, the file is saved
dir.create("SingleCellExpFiles")

#Here, we exclude all cells that do not have a complete BCR
#Three cells are not represented among the transcriptome-containing cells, 
#but these are likely to disappear downstream in the filtering steps to come anyway.
saveRDS(csfSce, file="Data/SingleCellExpFiles/csfSce_1_plus_flow_and_GLMPCA.rds")

#Now, as we have more information, including the cell type, we can identify and exclude cells
#that either have a TCR but that are B-cells phenotypically, or the opposite. 

tcrContainers <- read.csv("Data/BCR_auxiliaries/tcr_bcr_doublets.csv")[,2]

length(which(colnames(csfSce) %in% tcrContainers))
#1
#What group does this one belong to?
csfSce$Specific[which(colnames(csfSce) %in% tcrContainers)]
#None

#These are of course excluded. There are none, currently.
#csfSce <- csfSce[,-which(colnames(csfSce) %in% tcrContainers)]


#That still leaves us with an impressive 754 cells. 
#"counts"    "cpm"       "logcounts"

#And now, before we go on to model the genes of interest, we remove genes
#that we know are not going to be of interest downstream. 
#Now, we exclude everything non-protein coding
sum(colSums(counts(csfSce)))
#553855763
csfSce <- csfSce[which(rowData(csfSce)$gene_biotype == "protein_coding"),]
sum(colSums(counts(csfSce)))
#405869883

#Here, we also remove all the mitochondrial transcripts
csfSce <- csfSce[-which(rowData(csfSce)$chromosome_name == "MT"),]
sum(colSums(counts(csfSce)))
#390400602

#When we model the gene variation, we have to include the cell types
donorCell <- paste0(csfSce$donor, "_", csfSce$cellType)
table(donorCell)
#1166_ASC   1166_B 1227_ASC   1227_B 1284_ASC   1284_B 
#.    181       26       39       42       51      128 

#Now, how many genes and reads per cell do we have on median?

#And with that- a ready model. 
saveRDS(csfSce, "Data/SingleCellExpFiles/csfSce_2_norm.rds")


#And we also create a file for only asc
saveRDS(csfSce[,which(csfSce$cellType == "ASC")], "Data/SingleCellExpFiles/csfSce_3_ASC_only.rds")




