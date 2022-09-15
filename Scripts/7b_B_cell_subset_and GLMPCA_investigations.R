library(SingleCellExperiment)
library(uwot)
library(ggplot2)
library(BiocParallel)
#Here, we investigate the 
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_3_GLMPCA.rds")


#Now, we will investigate how many classical B cell subsets that are still 
#present among these IgG+ cells. 

#What we will do is that we will identify
#prototypical cells from the five cell types: ASC, naive, memory, double negative
#and non-switched memory, and then we will use the transcriptome information to
#identify which of these other groups that they best fit into.
#We start by extracting the dataset for ease of use: 

flowDataB <- as.data.frame(t(normcounts(altExp(csfSce, "flowData"))))[which(csfSce$Cell_type == "B"),]
plot(flowDataB$IgD, flowDataB$CD27)


#When looking at this data above, it is clear that there are no truly IgD positive
#cells in this dataset any longer (the positive values go up to around 6 in the full
#dataset), which is not strange, considering the exclusion
#of all IgM and IgD positive csells on the BCR level. 
#What we do have is a difference along the CD27 axis. But first, it might be
#prudent to investigate if this at all correlates with any transcriptomic differences
#To check this, we will make an umap based on the GLMPCA. 
glmpcaUmap <- umap(reducedDim(csfSce))

#And now, we plot CD27 on top of this.
ggplotData <- data.frame(glmpcaUmap, "Cell_type" = csfSce$Cell_type,
                         "CD27" = as.data.frame(t(normcounts(altExp(csfSce, "flowData"))))$CD27,
                         "CD38" = as.data.frame(t(normcounts(altExp(csfSce, "flowData"))))$CD38,
                         "Donor" = csfSce$donor)

dir.create("Diagnostics/Integrated")
ggplot(ggplotData, aes(x=X1, y=X2, shape=Cell_type, size=CD27, color = Donor)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Integrated/GLMPCAUMAP_CD27_Donor_Cell_type.pdf")
ggplot(ggplotData, aes(x=X1, y=X2, shape=Cell_type, size=CD38, color = CD27)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Integrated/GLMPCAUMAP_CD27_CD38.pdf")
ggplot(ggplotData, aes(x=X1, y=X2, shape=Cell_type, size=CD27, color = CD38)) +
    geom_point() + theme_bw()
ggsave("Diagnostics/Integrated/GLMPCAUMAP_CD38_CD27.pdf")

#It turns out that CD27 not seems to correlate at all to the transcriptomic signal. 
#But CD38 instead does correlate very strongly. This speaks in favour of not making
#any further separation of the cells at this stage. Here, the correlations are 
#formalised. 

corrListProtX1 <- apply(normcounts(altExp(csfSce, "flowData")), 
                        1, function(x) cor(glmpcaUmap[,1], x))
corrListProtX2 <- apply(normcounts(altExp(csfSce, "flowData")), 
                        1, function(x) cor(glmpcaUmap[,2], x))

#So, on the protein side, as expected, it is CD38 and CD138 that are most correlated, trailed by a nagative correlation to CD20. 
#This is in other words the ASC/B axis. 

#Another factor that is highly interesting however is correlating the 
#PCA signal to the transcriptome to identify what drives the separation of the
#ASC into two clusters. In other words; what explains X2?
corrVecSCEX2 <- unlist(bplapply(row.names(csfSce), function(x){
    cor(glmpcaUmap[,2], logcounts(csfSce)[x,])
}))

names(corrVecSCEX2) <- rowData(csfSce)$hgnc_symbol
corrVecSCEX2Ordered <- corrVecSCEX2[order(abs(corrVecSCEX2), decreasing = TRUE)]

#So deeming from this, it seems like this separate population of ASC are 
#driven by a large set of upregulated genes. We might look deeper into this later. 

dir.create("Data/Integrated")
corrDfX2 <- data.frame("names" = names(corrVecSCEX2Ordered), 
                       "Pearson_correlation" = corrVecSCEX2Ordered)
write.csv(corrDfX2, "Data/Integrated/GLMPCAUMAP_V2_vs_all_transcripts.csv",
          row.names = FALSE)

test <- read.csv("Data/Integrated/GLMPCAUMAP_V2_vs_all_transcripts.csv")
