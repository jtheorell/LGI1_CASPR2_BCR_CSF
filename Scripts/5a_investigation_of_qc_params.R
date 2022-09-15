#Here, we check if it is necessary to exclude the transcriptomes with
#a high ERCC content, or if they are representative in all other respects. 
library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)
library(uwot)

#Here the data from step 1 is added.
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_1_preQC.rds")
flowData <- read.csv("Data/Flow_data/flowDataPlusIndexAndcellType.csv")

#Now, we add a column of cell type information
csfSce$cell_type <- "ASC"
csfSce$cell_type[which(colnames(csfSce) %in% 
                           flowData$Cell[which(flowData$Cell_type == "B")])] <- "B"
#Identification of mitochondrial genes
mito <- which(rowData(csfSce)$chromosome_name=="MT")

#Now the preCellQCMetrcis are calculated. 
stats <- perCellQCMetrics(csfSce, subsets=list(Mt=mito))
colData(csfSce) <- cbind(colData(csfSce), stats)
csfSce$plate <- factor(csfSce$plate)

#Now, we are making a very crude analysis here. We select the 1000 most variable
#genes. 
csfSce <- computeSumFactors(csfSce, cluster=csfSce$donor) 
csfSce <- logNormCounts(csfSce)


dec.Sce <- modelGeneVarWithSpikes(csfSce, "ERCC", block = csfSce$donor)
chosen.hvgs <- getTopHVGs(dec.Sce, n = 1000)

csfSce$umap1d <- uwot::umap(t(logcounts(csfSce)[chosen.hvgs,]),
                          n_components = 1)

ggData <- as.data.frame(colData(csfSce))
ggplot(ggData, aes(x=altexps_ERCC_percent, y=umap1d, shape=cell_type, color=donor)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/SCE/ERCC_vs_1DUmap.pdf")
#This shows that theB cells generally have considerably higher ERCC counts, 
#which is expected, as their transcriptomes are smaller, and that  
#donor 1166 has higher counts both for the ASC and the B cells, so the threshold
#needs to be high there. All in all, this plot would speak for a threshold
#of 30%. 

ggplot(ggData, aes(x=subsets_Mt_percent, y=umap1d, shape=cell_type, color=donor)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/SCE/Mito_vs_1DUmap.pdf")
#Once again, there is more spread in the B cell pool, but it seems like
#all but one non-extreme outlier is below 10%, so this can be kept as a threshold

ggplot(ggData, aes(x=sum, y=umap1d, shape=cell_type, color=donor)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/SCE/Total_count_vs_1DUmap.pdf")

#Even if this graph only shows a normal distribution at the low tail, it is still
#so that it is hard to use cells with counts below 100000 counts, so these are
#excluded. 

ggplot(ggData, aes(x=detected, y=umap1d, shape=cell_type, color=donor)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/SCE/Number_of_detected_features_vs_1DUmap.pdf")

#Here, it seems like 1000 genes is a sensible cutoff. 




