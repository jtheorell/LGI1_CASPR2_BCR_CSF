library(scran)
library(scater)
#We will also make a specific investigation of a gene set constructed by Professor Irani
#that relates to currently investigated treatment targets. 
#HGNC-symbols for 
#The CD20 gene is called MS4A1
#CD138 is called syndecan 1 (SDC1)
#BAFF is called TNFSF13B
#BAFF-R is called TNFRSF13C
#APRIL is called ANP32B
#CD267/TACI is called TNFRSF13B
#GM-CSF is called CSF2 and its receptors called thereafter
#I did not include SIPR, as it is unclear to me which gene that is meant
#by this. 
library(ggplot2)
resList <- readRDS("Results/Comp_to_others/All_DE_integrated_downsamp.rds")
#We start by removing the specific-vs-non, as that comparison is completely negative
resList <- resList[-which(names(resList) == "Median")]
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")
targetsAE <- c("CD27", "CD38", "CXCR4", "CXCR5", "MS4A1", "CD19", "SDC1", "CD22", "TNFSF13B",
               "TNFRSF13C", "ANP32B", "TNFRSF13B", "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IFNG",
               "IFNGR1","IFNGR2", "TNF", "TNFRSF1A", "TNFRSF1B", "IL10", "IL10RA","IL10RB", 
               "CSF2", "CSF2RA", "CSF2RB", "IL6", "IL6R", "BTK", "ITGA4", "ITGB1", "CD40")

#Now we are going to investigate which positions these arrive in for the three datasets
#But before that: how many are present among the common genes?
length(which(row.names(resList[[1]]) %in% targetsAE))/length(targetsAE)
#0.6470588
#Excluded genes are: 
targetsAE[-which(targetsAE %in% row.names(resList[[1]]))]
# "CXCR5"    "TNFSF13B" "IL12A"    "IL12B"    "IL12RB2"  "IFNG"     "TNF"      "TNFRSF1A" "IL10"    
# "CSF2"     "CSF2RA"   "IL6"  

#Now, we create a dataframe with all the positions and related data. 
groupList <- list("")
targetList <- lapply(1:length(resList), function(x){
    locPos <- which(row.names(resList[[x]]) %in% targetsAE)
    locDf <- data.frame(row.names(resList[[x]])[locPos],
                        locPos,
                        resList[[x]]$logFC[locPos])
    colnames(locDf) <- c("Name", paste0(names(resList)[[x]], "_position"), 
                         paste0(names(resList)[[x]], "_logFC"))
    locDf <- locDf[order(locDf$Name),]
})

all(sapply(targetList, function(x){
    identical(targetList[[1]]$Name, x$Name)
}))
#TRUE
#So, now we can combine them, and keep the names as rownames in the new object. 
targetDf <- data.frame(do.call("cbind", lapply(targetList, "[", 2:3)),
                            row.names = targetList[[1]]$Name)

targetDf$medianPosition <- rowMedians(as.matrix(targetDf[,grep("position", colnames(targetDf))]))
targetDf$medianFC <- rowMedians(as.matrix(targetDf[,grep("logFC", colnames(targetDf))]))
targetDfOrdered <- targetDf[order(targetDf$medianF, decreasing = TRUE),]

#Now, we are going to establish "confidence intervals", or more correctly, determine
#the average AUC values that correspond to 95 and 99% of the data. This will be 
#calculated by randomly selecting 32 of the top expressed genes 1000 times, and after 
#this identifying
#which average AOC values that are on the 95th and 99th percentile of the total data. 

csfSceRed <- csfSce[which(rowData(csfSce)$hgnc_symbol %in% row.names(resList[[1]])),]
row.names(csfSceRed) <- rowData(csfSceRed)$hgnc_symbol
dec.Sce <- modelGeneVarWithSpikes(csfSceRed, "ERCC")
chosen.hvgs <- getTopHVGs(dec.Sce)
#This selects 2069 genes. 

fcVec <- unlist(lapply(1:1000, function(x){
    locSamp <- sample(row.names(resList[[1]])[which(row.names(resList[[1]]) %in% chosen.hvgs)], 34)
    locMat <- do.call("cbind", lapply(1:length(resList), function(y){
        locPos <- which(row.names(resList[[y]]) %in% locSamp)
        locDf <- data.frame(row.names(resList[[y]])[locPos],
                            resList[[y]]$logFC[locPos])
        colnames(locDf) <- c("Name", "logFC")
        locDf <- locDf[order(locDf$Name),]
        locDf$logFC
    }))
    
    rowMeans(locMat)
}))

quantVals <- quantile(fcVec, c(0.01, 0.99))
#    1%       99% 
#    -2.385547  4.228274 

#Now, we plot this. 


# Color by groups
ggData <- data.frame("Median_logFC" = c(targetDfOrdered$medianFC, 
                                        fcVec),
                     "Group" = c(rep("Target_genes", nrow(targetDfOrdered)),
                                 rep("All_genes", length(fcVec))))
ggplot(ggData, aes(x=Median_logFC, color=Group, fill=Group)) + 
    geom_histogram(aes(y=after_stat(density)), alpha=0.5, 
                   position="identity")+ geom_vline(xintercept = quantVals[1]) +
    geom_vline(xintercept = quantVals[2]) +
    geom_density(alpha=.2) + theme_bw()
ggsave("Results/Comp_to_others/Drug_target_genes_vs_all_downsamp.pdf", width = 7, height = 5)

#Are there any that fall outside of the margins?
row.names(targetDfOrdered)[which(targetDfOrdered$medianFC > quantVals[2])]
#None
row.names(targetDfOrdered)[which(targetDfOrdered$medianFC < quantVals[1])]
#"ITGB1"

#Now, the non-MS. 
targetDfNoMS <- targetDf[,-grep("MS|Specific|median", colnames(targetDf))]
targetDfNoMS$meanPosition <- rowMeans(targetDfNoMS[,c(1,3)])
targetDfNoMS$meanAUC <- rowMeans(targetDfNoMS[,c(2,4)])

targetDfNoMS$medianPosition <- rowMedians(as.matrix(targetDfNoMS[,grep("position", colnames(targetDfNoMS))]))
targetDfNoMS$medianFC <- rowMedians(as.matrix(targetDfNoMS[,grep("logFC", colnames(targetDfNoMS))]))

#We will now repeat the above but without the MS. 
noMsResList <- resList[-grep("MS", names(resList))]
set.seed(101)
cohenVecNoMS <- unlist(lapply(1:1000, function(x){
    locSamp <- sample(row.names(noMsResList[[1]])[which(row.names(noMsResList[[1]]) %in% chosen.hvgs)], 34)
    locMat <- do.call("cbind", lapply(1:length(noMsResList), function(y){
        locPos <- which(row.names(noMsResList[[y]]) %in% locSamp)
        locDf <- data.frame(row.names(noMsResList[[y]])[locPos],
                            noMsResList[[y]]$logFC[locPos])
        colnames(locDf) <- c("Name", "logFC")
        locDf <- locDf[order(locDf$Name),]
        locDf$logFC
    }))
    
    rowMeans(locMat)
}))

quantValsNoMs <- quantile(cohenVecNoMS, c(0.01, 0.99))

#    1%       99% 
#    -2.669788  4.093346 

# Color by groups
ggData <- data.frame("Median_logFC" = c(targetDfNoMS$medianFC, 
                                        cohenVecNoMS),
                     "Group" = c(rep("Target_genes", nrow(targetDfNoMS)),
                                 rep("All_genes", length(cohenVecNoMS))))
ggplot(ggData, aes(x=Median_logFC, color=Group, fill=Group)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, 
                   position="identity") + geom_vline(xintercept = quantVals[1]) +
    geom_vline(xintercept = quantVals[2]) +
    geom_density(alpha=.2) + theme_bw()
ggsave("Results/Comp_to_others/Drug_target_genes_vs_no_MS_data_downsamp.pdf", width = 7, height = 5)

#Which are the tlow ones?
row.names(targetDfNoMS)[which(targetDfNoMS$medianFC < quantValsNoMs[1])]
#"ITGB1"
#ANd hte high one?
row.names(targetDfNoMS)[which(targetDfNoMS$medianFC > quantValsNoMs[2])]
#"SDC1", i.e. CD1138


#NOw, what about the MS only?
targetDfMS <- targetDf[,grep("MS", colnames(targetDf))]
targetDfMS$meanPosition <- rowMeans(targetDfMS[,c(1,3)])
targetDfMS$meanAUC <- rowMeans(targetDfMS[,c(2,4)])

targetDfMS$medianPosition <- rowMedians(as.matrix(targetDfMS[,grep("position", colnames(targetDfMS))]))
targetDfMS$medianCohen <- rowMedians(as.matrix(targetDfMS[,grep("logFC", colnames(targetDfMS))]))

msResList <- resList[grep("MS", names(resList))]
set.seed(101)
cohenVecMS <- unlist(lapply(1:1000, function(x){
    locSamp <- sample(row.names(msResList[[1]])[which(row.names(msResList[[1]]) %in% chosen.hvgs)], 34)
    locMat <- do.call("cbind", lapply(1:length(msResList), function(y){
        locPos <- which(row.names(msResList[[y]]) %in% locSamp)
        locDf <- data.frame(row.names(msResList[[y]])[locPos],
                            msResList[[y]]$logFC[locPos])
        colnames(locDf) <- c("Name", "logFC")
        locDf <- locDf[order(locDf$Name),]
        locDf$logFC
    }))
    
    rowMeans(locMat)
}))

quantValsMs <- quantile(cohenVecMS, c(0.01, 0.99))

#1%      99% 
#-2.074742  5.290023

# Color by groups
ggData <- data.frame("Median_Cohen" = c(targetDfMS$medianCohen, 
                                        cohenVecMS),
                     "Group" = c(rep("Target_genes", nrow(targetDfMS)),
                                 rep("All_genes", length(cohenVecMS))))
ggplot(ggData, aes(x=Median_Cohen, color=Group, fill=Group)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, 
                   position="identity") + geom_vline(xintercept = quantValsMs[1]) +
    geom_vline(xintercept = quantValsMs[2]) + 
    geom_density(alpha=.2) + theme_bw()
ggsave("Results/Comp_to_others/Drug_target_genes_vs_MS_data_downsamp.pdf")

#Which are the low ones? 
row.names(targetDfMS)[which(targetDfMS$medianCohen < quantValsMs[1])]
#"MS4A1", i.e. CD20
row.names(targetDfMS)[which(targetDfMS$medianCohen > quantValsMs[2])]
#None
