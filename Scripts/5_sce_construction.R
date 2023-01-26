library(ensembldb)
library(tximport)
library(SingleCellExperiment)
library(scran)
library(scater)
library(uwot)
library(DepecheR)
library(bluster)

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
#Cell exclusions based on qc params
#############

#Identification of mitochondrial genes
mito <- which(rowData(csfSce)$chromosome_name=="MT")

#Now the preCellQCMetrcis are calculated. 
stats <- perCellQCMetrics(csfSce, subsets=list(Mt=mito))
colData(csfSce) <- cbind(colData(csfSce), stats)
csfSce$plate <- factor(csfSce$plate)

#As shown in 5a, there are considerable differences in the transcriptomes
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

csfSce <- csfSce[,!csfSce$lowQual]

#Here, 188 cells, or 10%, are lost. 
#We are left with 1684 cells.
#We of course also exclude the A1, as they are positive controls with 10 cells
#each
csfSce <- csfSce[,-grep("A01", colnames(csfSce))]
#This removes another 10. 

#Now, we generate a GLM-PCA. 
library(scry)
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

#And this is added back: 
reducedDim(csfSce, "GLMPCA") <- reducedDim(glmpcaSCE)

#Now, we are going to integrate this data with the flow cytometry data, 
#to generate both umaps and clusters. 

#Here, the flow data is integrated into the SCE
cellTypeData <- read.csv("Data/Cytometry/flowDataPlusIndex.csv")

flowDataHere <- cellTypeData[which(cellTypeData$Cell %in% colnames(csfSce)),]
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

#Now; is the flow data and the GLM-PCA on the same scale?
summary(as.vector(as.matrix(reducedDim(csfSce, "GLMPCA"))))
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-9.604e-04 -1.030e-08 -5.000e-10  0.000e+00  1.740e-08  3.008e-04 

#So that will not be the case. We will take this whole thing times 10000 just to 
#make the scale more reasonable. 
glmpca10000 <- reducedDim(csfSce, "GLMPCA")*10000

#And now, we check the flow data. 
summary(as.vector(normcounts(altExp(csfSce, "flowData"))[3:nrow(altExp(csfSce, "flowData"))]))
##Here, we had to start by removing the FSC and SSC, that clearly impact the results
#So that by itself shows that we need to expand here. So what we will do is to use
#The dScale function in DepecheR. 
scaledFlowData <- dScale(flowDataDf)*2

#These are now combined and used to generate a new umap and clusterings. 
set.seed(200)
combinedUmap <- uwot::umap(cbind(scaledFlowData, glmpca10000))

#And we display the flow markers on this
dir.create("Diagnostics/Protein_plus_transcriptome")
dColorPlot(scaledFlowData, xYData = combinedUmap, 
           plotDir = "Diagnostics/Protein_plus_transcriptome/Prot_on_UMAP")

#Now, this is integrated into the sce, and then we will run some clustering.
reducedDim(csfSce, "flowAndGLMPCA") <- cbind(scaledFlowData, glmpca10000)

reducedDim(csfSce, "combinedUmap") <- combinedUmap

#Now, we cluster the cells. 
set.seed(675)
cellClusters <- clusterCells(csfSce, use.dimred = "flowAndGLMPCA", 
                             BLUSPARAM = NNGraphParam(cluster.fun="louvain" ))
dColorPlot(cellClusters, xYData = combinedUmap,
           plotDir = "Diagnostics/Protein_plus_transcriptome",
           plotName = "Louvain_clusters",
           colorScale = "dark_rainbow")

#These clusters are interpretable. We will not need to go into the CD4 clusters, 
#So the definitions are here: 

clustersNamed <- sapply(as.character(cellClusters), switch, "1" = "CD4T", "2" = "CD8T",
                        "3" = "CD4T", "4" = "B", "5" = "ASC", 
                        "6" = "NK")

dColorPlot(clustersNamed, xYData = combinedUmap,
           plotDir = "Diagnostics/Protein_plus_transcriptome",
           plotName = "Louvain_clusters_named",
           colorScale = "dark_rainbow")

csfSce$clustersLouvain <- clustersNamed
#saveRDS(csfSce, "../External/Data/csfSce_1c_GLMPCA_plus_cell_types_exp_1.rds")
#Here, we exclude all non-ASC, B or CD8 T
csfSce <- csfSce[,grep("8|B|A", csfSce$clustersLouvain)]
table(csfSce$clustersLouvain)
#ASC       B    CD8T 
#297      235   230

#Here, the file is saved
dir.create("SingleCellExpFiles")

#Here, we exclude all cells that do not have a complete BCR
#Three cells are not represented among the transcriptome-containing cells, 
#but these are likely to disappear downstream in the filtering steps to come anyway.
saveRDS(csfSce, file="Data/SingleCellExpFiles/csfSce_1_plus_flow_and_GLMPCA.rds")

#Now, as we have more information, including the cell type, we can identify and exclude cells
#that either have a TCR but that are B-cells phenotypically, or the opposite. 

tcrContainers <- read.csv("Data/BCR_auxiliaries/tcr_bcr_doublets.csv")[,2]

length(which(colnames(csfSce)[which(csfSce$clustersLouvain != "CD8T")] %in% tcrContainers))
#6

#These are of course excluded. 
csfSce <- csfSce[,-which(colnames(csfSce) %in% tcrContainers & 
                                csfSce$clustersLouvain != "CD8T")]

#Now, the same for the BCRs
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

CD8T_BCR <- which(colnames(csfSce) %in% BCR_all$CELL & csfSce$clustersLouvain == "CD8T")
length(CD8T_BCR)
#2
csfSce <- csfSce[,-CD8T_BCR]

#That still leaves us with an impressive 754 cells. 

#Now, we normalise the data. 
csfSce <- computeSpikeFactors(csfSce, spikes = "ERCC") 
csfSce <- logNormCounts(csfSce)

assayNames(csfSce)
#"counts"    "cpm"       "logcounts"

#And now, before we go on to model the genes of interest, we remove genes
#that we know are not going to be of interest downstream. 
#Now, we exclude everything non-protein coding
sum(colSums(counts(csfSce)))
#838624479
csfSce <- csfSce[which(rowData(csfSce)$gene_biotype == "protein_coding"),]
sum(colSums(counts(csfSce)))
#669217972

#Here, we also remove all the mitochondrial transcripts
csfSce <- csfSce[-which(rowData(csfSce)$chromosome_name == "MT"),]
sum(colSums(counts(csfSce)))
#643639381

#When we model the gene variation, we have to include the cell types
donorCell <- paste0(csfSce$donor, "_", csfSce$clustersLouvain)
table(donorCell)
#1166_ASC    1166_B 1166_CD8T  1227_ASC    1227_B 1227_CD8T  1284_ASC    1284_B 1284_CD8T 
#190        36        46        49        48        55        54       149       127 

dec.Sce <- modelGeneVarWithSpikes(csfSce, "ERCC", block = donorCell)
chosen.hvgs <- getTopHVGs(dec.Sce)

#And with that- a ready model. 
saveRDS(csfSce, "Data/SingleCellExpFiles/csfSce_2_norm.rds")

