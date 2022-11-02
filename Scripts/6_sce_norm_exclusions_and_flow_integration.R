library(SingleCellExperiment)
library(scran)
library(scater)

csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_1_preQC.rds")
dim(csfSce)
#58565   381
#And this is good, as it is the same number as the total number of licensed BCRs

#Before we go into any other exclusions, we use the information from the cell
#surface and the separate TraCeR analysis, to exclude cells with a BCR but a non-
#BCR phenotype, as well as cells containing a TCR and a BCR. 

#Exclusion of cells with a TCR. TraCeR has been used to identify these: 
tcrContainers <- read.csv("Data/BCR_auxiliaries/tcr_bcr_doublets.csv")[,2]

csfSce <- csfSce[,-which(colnames(csfSce) %in% tcrContainers)]
#In total, this makes us lose 10 cells. 

#Now, we import information about the cell type, and exclude the discordant cells
#i.e., cells that carry a BCR but that does not have a B-cell phenotype. 
cellTypeData<- read.csv("Data/Cytometry/flowDataPlusIndexAndcellType.csv", row.names = 1)

csfSce$Cell_type <- NA
for(i in colnames(csfSce)){
    csfSce$Cell_type[which(colnames(csfSce) == i)] <- 
        cellTypeData$Cell_type[which(cellTypeData$Cell == i)]
}

table(csfSce$Cell_type)
#ASC       B    CD4T    CD8T doublet      NK 
#267      95       4       2       2       1

#So these 9 non-B are excluded here. 
csfSce <- csfSce[,-which(csfSce$Cell_type %in% c("CD4T", "CD8T", "doublet", "NK"))]

#Now, we get over into the normal exclusion steps. 

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

#So after losing these 15 cells, we are down to 347. 

#Normalization. There are 100s of methods, but according to this paper, 
#https://www.frontiersin.org/articles/10.3389/fgene.2020.00041/full 
#the scran method is not too bad. Also, we are only going to make
#comparisons within individual donors between the antibody-positive and 
#-negative groups, there is no reason to believe that we would miss anything
#greatly important by choosing a slightly worse method, or for example not adjusting
#for amplification biases, etc, that could provide a much bigger problem if the
#main focus was on separation cellular subsets. 
#We save this in case something goes wrong: 
preQCSce <- csfSce
#We will take advantage of the spike-ins here. However, it is clear that 
#there are different amounts in the different batches, and that can sadly
#not any longer be taken into account in this algorithm. But this might be 
#a minor problem anyway, as all comparisons will be conducted internally to each
#donor. 
csfSce <- computeSpikeFactors(csfSce, spikes = "ERCC") 
csfSce <- logNormCounts(csfSce)

summary(sizeFactors(csfSce))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03379 0.33310 0.58475 1.00000 1.07073 6.06699 

assayNames(csfSce)
#"counts"    "cpm"       "logcounts"

#And now, before we go on to model the genes of interest, we remove genes
#that we know are not going to be of interest downstream. 
#Now, we exclude everything non-protein coding
sum(colSums(counts(csfSce)))
#398634328
csfSce <- csfSce[which(rowData(csfSce)$gene_biotype == "protein_coding"),]
sum(colSums(counts(csfSce)))
#282369946
#Not very unexpectedly, this removes a considerable portion (30%) of the data. 

#Here, we also remove all the mitochondrial transcripts
csfSce <- csfSce[-which(rowData(csfSce)$chromosome_name == "MT"),]
sum(colSums(counts(csfSce)))
#271879831

#When we model the gene variation, we have to include the ASC/B cell definitions

donorCell <- paste0(csfSce$donor, "_", csfSce$Cell_type)
table(donorCell)
#1166_ASC   1166_B 1227_ASC   1227_B 1284_ASC   1284_B 
#175       29       41       16       43       43 

dec.Sce <- modelGeneVarWithSpikes(csfSce, "ERCC", block = donorCell)
chosen.hvgs <- getTopHVGs(dec.Sce)

#Here, the flow data is integrated into the SCE
flowDataHere <- cellTypeData[which(cellTypeData$Cell %in% colnames(csfSce)),]
flowDataHereOrdered <- flowDataHere[match(colnames(csfSce),flowDataHere$Cell),]

flowDataMat <- t(flowDataHereOrdered[,1:19])
colnames(flowDataMat) <- flowDataHereOrdered$Cell

identical(colnames(flowDataMat), colnames(csfSce))
#TRUE. And with that we are ready to integrate this. 

flowDataSce <- 
    SingleCellExperiment(assays = list("normcounts" = flowDataMat))
altExp(csfSce, "flowData") <- flowDataSce

#And here it feels reasonable to start a new chapter. 
saveRDS(csfSce, "Data/SingleCellExpFiles/csfSce_2_norm_and_flow_integration.rds")
