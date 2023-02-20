#In this analysis, we will complement our specific-vs-not analysis with also
#comparing our specific cells to other ASCs from other datasets. 
library(SingleCellExperiment)
library(scran)
library(scater)
library(BiocParallel)
library(quminorm)
library(edgeR)

sceList <- readRDS("Data/Comp_to_others/All_sce_common_genes_pre_all.rds")

#We start by removing the non-specific ASC from our datset and the B-cells from all datasets
sceListAsc <- lapply(sceList, function(x){
    if(length(which("Specific" %in% colnames(colData(x)) > 0))){
        x <- x[,which(x$Specific == "TRUE")]
    }
    x[,which(x$cellType == "ASC")]
})

dir.create("Results/Specific_vs_others")

#Now, before starting any downsampling, we need to save a picture showing why 
#it is needed. 
plotSce <- SingleCellExperiment(
    assays = list("logcounts" = do.call("cbind",lapply(names(sceListAsc), function(x){
        locSce <- sceListAsc[[which(names(sceListAsc) == x)]]
        locSce <- locSce[which(row.names(locSce) %in% c("XBP1", "PRDM1", "IRF4",
                                                        "B2M", "ACTB", "PTPRC")),]
        locSce$Dataset <- x
        logcounts(locSce)
    }))))
plotSce$group <- unlist(sapply(names(sceListAsc), function(x){
    rep(x, ncol(sceListAsc[[x]]))
}))
plotSce$cellType <- "ASC"

plotExpression(plotSce, features=row.names(plotSce), 
               x="group",  colour_by="cellType", ncol = 3)
dev.copy(pdf,'Results/Specific_vs_others/High_genes_pre_qumi_and_downsamp.pdf', height = 8, width = 8)
dev.off()

plotExpression(plotSce, features="XBP1", 
               x="group",  colour_by="cellType", ncol = 3)
dev.copy(pdf,'Results/Specific_vs_others/XBP1_pre_qumi_and_downsamp.pdf', height = 8, width = 8)
dev.off()

#Now, we will generate quasi-UMIs for the aeSce, using the quminorm package. To do this
#in a correct way, we will need the full data, including the non-protein coding 
#transcripts, to get the correct distribution. 
rawAe <- readRDS("Data/SingleCellExpFiles/csfSce_1_plus_flow_and_GLMPCA.rds")
aeSceTpm <- calculateTPM(rawAe)
#Here, we have tested different values, to reach an optimal similarity in the expression
#of B2M, that in the original setting is too harshly affected. 
aeSceQumiMat <- quminorm(aeSceTpm, shape = 3, mc.cores = 7)

#Now, we remove all rows and columns not included in the asc dataset here used
aeSceQumiMatRed <- aeSceQumiMat[which(row.names(aeSceQumiMat) %in% rowData(sceListAsc$AE_CSF)$ensembl_gene_id),
                                which(colnames(aeSceQumiMat) %in% colnames(sceListAsc$AE_CSF))]
aeSceQumiMatOrd <- aeSceQumiMatRed[match(rowData(sceListAsc$AE_CSF)$ensembl_gene_id, row.names(aeSceQumiMatRed)),]

all(identical(rownames(aeSceQumiMatOrd), rowData(sceListAsc$AE_CSF)$ensembl_gene_id),
    identical(colnames(aeSceQumiMatOrd), colnames(sceListAsc$AE_CSF)))
#TRUE
rownames(aeSceQumiMatOrd) <- rownames(sceListAsc$AE_CSF)

#And here the qumi sum
median(colSums(aeSceQumiMatOrd))
#25431

aeSceQumi <- sceListAsc$AE_CSF
assays(aeSceQumi) <- list("counts" = aeSceQumiMatOrd)

#And for systematics sake, we do the same for the smartseq3 dataset
rawHealthy <- readRDS("../External/Data/Healthy/3_post_Kotliarov.rds")

healthySceTpm <- calculateTPM(rawHealthy)
healthySceQumiMat <- quminorm(healthySceTpm, shape = 3, mc.cores = 7)
healthySceQumiRed <- healthySceQumiMat[which(row.names(healthySceQumiMat) %in% rowData(sceListAsc$Healthy_PBMC)$ensembl_gene_id),
                                which(colnames(healthySceQumiMat) %in% colnames(sceListAsc$Healthy_PBMC))]

healthySceQumiOrd <- healthySceQumiRed[match(rowData(sceListAsc$Healthy_PBMC)$ensembl_gene_id, row.names(healthySceQumiRed)),]

all(identical(rownames(healthySceQumiOrd), rowData(sceListAsc$Healthy_PBMC)$ensembl_gene_id),
    identical(colnames(healthySceQumiOrd), colnames(sceListAsc$Healthy_PBMC)))
#TRUE
rownames(healthySceQumiOrd) <- rownames(sceListAsc$Healthy_PBMC)

median(colSums(healthySceQumiOrd))
#7967
healthySceQumi <- sceListAsc$Healthy_PBMC
assays(healthySceQumi) <- list("counts" = healthySceQumiOrd)

#And these are integrated
sceListAsc$AE_CSF <- aeSceQumi
sceListAsc$Healthy_PBMC <- healthySceQumi

#Here, we plot an intermediate version
sceListAsc$AE_CSF <- computeSumFactors(sceListAsc$AE_CSF)
sceListAsc$AE_CSF <- logNormCounts(sceListAsc$AE_CSF)
sceListAsc$Healthy_PBMC <- computeSumFactors(sceListAsc$Healthy_PBMC)
sceListAsc$Healthy_PBMC <- logNormCounts(sceListAsc$Healthy_PBMC)

#Now, we are going to downsample all datasets
saveRDS(sceListAsc, "Data/Comp_to_others/All_sce_common_genes_qumi.rds")

plotSce <- SingleCellExperiment(
    assays = list("logcounts" = do.call("cbind",lapply(names(sceListAsc), function(x){
        locSce <- sceListAsc[[which(names(sceListAsc) == x)]]
        locSce <- locSce[which(row.names(locSce) %in% c("XBP1", "PRDM1", "IRF4",
                                                        "B2M", "ACTB", "PTPRC"))]
        locSce$Dataset <- x
        logcounts(locSce)
    }))))
plotSce$group <- unlist(sapply(names(sceListAsc), function(x){
    rep(x, ncol(sceListAsc[[x]]))
}))
plotSce$cellType <- "ASC"
plotExpression(plotSce, features=row.names(plotSce),
               x="group", colour_by="cellType", ncol = 3)
dev.copy(pdf,'Results/Specific_vs_others/High_genes_post_qumi_pre_downsamp.pdf', height = 8, width = 8)
dev.off()

plotExpression(plotSce, features="XBP1",
               x="group", colour_by="cellType", ncol = 3)
dev.copy(pdf,'Results/Specific_vs_others/XBP1_post_qumi_pre_downsamp.pdf', height = 8, width = 8)
dev.off()

medianSums <- lapply(sceListAsc, function(x) median(colSums(as.matrix(counts(x)))))

minMedSum <- round(min(unlist(medianSums)))
#4438
madVal <- mad(colSums(as.matrix(counts(sceListAsc[[which.min(unlist(medianSums))]]))))
#1822.115

minVal <- min(colSums(as.matrix(counts(sceListAsc[[which.min(unlist(medianSums))]]))))
#1840
#So this is the amount we will downsample to. There are of course multiple ways of 
#doing this, but what I will do, as the matrices are not overwhelmingly large, is
#to expand each cell to a 0-1 vector, and then sample this vector. However, adding
#another layer of complexity, we want the median within each dataset to have the
#same number of counts,but that is not the same as saying that all cells should
#have this count. For this reason, we will first identify the median cell and then scale
#the values based on the MAD around this. 
#We are going to run this 100 times
#Before doing this, we will import the genes that are expressed differentially
#in our ASC dataset, as this will be the basis for downstream analyses. 
ascExprs <- row.names(readRDS("Results/Specific_vs_not/EdgeR_result.rds"))

sceListRed <- lapply(sceListAsc, function(x) x[which(row.names(x) %in% ascExprs),])

#Are they all identical now, row-wise? 
all(sapply(sceListRed, function(x) identical(rownames(sceListRed$AE_CSF), rownames(x))))
#TRUE

#For this to work, we also need to name the cells of the covid dataset
colnames(sceListRed$Covid_PBMC) <- paste0("Cell_", seq_len(ncol(sceListRed$Covid_PBMC)))

dir.create("Results/Specific_vs_others/Bootstrap_DE")
for(i in 1:100){
    print(paste0("Now, iteration ", i, " is starting"))
    allCountListDown <- lapply(sceListRed, function(x){
        locColSums <- colSums(as.matrix(counts(x)))
        locMedian <- median(locColSums)
        locMad <- mad(locColSums)
        locColSumsCenter <- locColSums-locMedian
        locColSumsScale <- locColSumsCenter/locMad
        locColSumsRefit <- round((locColSumsScale*madVal)+minMedSum)
        #And now, we let all values that are below the lowest value in the sample
        #with the lowest median be set to that value.
        locColSumsRefit[which(locColSumsRefit < minVal)] <- minVal
        #Now, we have the single-cell resolution values to rescale to. Now, we
        #do the actual sampling of the data. There are a few cells that have fewer
        #than the expected number of transcripts. For this reason, and to keep
        #this statistically sound, we will allow for resampling.
        set.seed(i)
        xSampled <- do.call("cbind",bplapply(names(locColSumsRefit), function(y){
            print(y)
            locDat <- counts(x)[,y]
            posValVec <- locDat[which(locDat > 0)]
            posValVecLong <- unlist(sapply(names(posValVec), function(z){
                locVal <- unname(posValVec[which(names(posValVec) == z)])
                if(locVal > 1){
                    innerVal <- rep(1, locVal)
                    names(innerVal) <- paste0("_", seq(1,locVal))
                    innerVal
                } else {
                    innerVal <- locVal
                    innerVal
                }
            }))
            posValVecSampled <- posValVecLong[sample(seq_along(posValVecLong),
                                                     locColSumsRefit[y],replace = TRUE)]
            #Here, we de-uniqueify the names
            names(posValVecSampled) <- gsub("\\._.+", "", names(posValVecSampled))
            posValVecShort <- sapply(unique(names(posValVecSampled)), function(z){
                sum(posValVecSampled[which(names(posValVecSampled) == z)])
            })
            #ANd now, we add this back together
            nullGenes <- names(locDat)[-which(names(locDat) %in% names(posValVecShort))]
            negValTotal <- rep(0, length(nullGenes))
            names(negValTotal) <- nullGenes
            valVecResurrected <- c(posValVecShort, negValTotal)
            valVecOrdered <- valVecResurrected[match(names(locDat), names(valVecResurrected))]
        }))
        colnames(xSampled) <- names(locColSumsRefit)
        xSampled
    })
    #Now, that the data is all downsampled and more comparable, we will add a new lognormcount calculation
    #here
    sceListAscDown <- lapply(names(sceListRed), function(x){
        locSce <- sceListRed[[which(names(sceListRed) == x)]]
        assays(locSce) <- list("counts" = allCountListDown[[which(names(allCountListDown) == x)]])
        locSce <- computeSumFactors(locSce)
        locSce <- logNormCounts(locSce)
        locSce$group <- x
        locSce
    })
    names(sceListAscDown) <- names(sceListRed)
    
    ###################
    #DATASET CONSTRUCTION
    ###################
    datList <- sceListAscDown[-grep("AE", names(sceListAscDown))]
    nameList <- names(datList)
    aeSceSpec <- sceListAscDown$AE_CSF
    sceList <- lapply(seq_along(datList), function(x){
        locDat <- datList[[x]]
        locCounts <- cbind(counts(aeSceSpec), counts(locDat))
        locSce <- SingleCellExperiment(assays = list("counts" = locCounts),
                                       rowData = rowData(sceListAscDown$Influensa))
        locSce$group <- c(aeSceSpec$group, locDat$group)
        locSce$donor <- c(aeSceSpec$donor, locDat$donor)
        locSce
    })
    names(sceList) <- nameList
    #############
    #ANALYSIS
    #############
    #NOw, we do the pseudo-bulk and DE analysis
    deResList <- lapply(sceList, function(x){
        current <- aggregateAcrossCells(x,
                                        id=colData(x)[,c("group", "donor")])
        current
        # Creating up a DGEList object for use in edgeR:
        y <- DGEList(counts(current), samples=colData(current))
        y
        #Now, we controversially only exclude donors with less than 4 cells
        discarded <- current$ncells < 4
        y <- y[,!discarded]
        #summary(discarded)
        #THis is for iteration 1, dataset 1
        #   Mode   FALSE    TRUE
        # logical      16      16
        y <- calcNormFactors(y)
        design <- model.matrix(~factor(current$group[!discarded]), y$donor)
        #design
        y <- estimateDisp(y, design)
        #summary(y$trended.dispersion)
        fit <- glmQLFit(y, design, robust=TRUE)
        #summary(fit$var.prior)
        res <- glmQLFTest(fit, coef=ncol(design))
        #summary(decideTests(res))
        #And here this result is saved.
    })
    saveRDS(deResList, paste0("Results/Specific_vs_others/Bootstrap_DE/All_DE_results_", i, ".rds"))
}

#So now, we have an as internally comparable dataset as possible.
#This last bootstrap is saved for later use. Here, we both save a full version, and a version where
#the AE data has been split and the ones without known specificity are excluded. 

saveRDS(sceListAscDown,
        "Data/Comp_to_others/Downsampled_specificity_included_boot_100.rds")

plotSce <- SingleCellExperiment(
    assays = list("logcounts" = do.call("cbind",lapply(names(sceListAscDown), function(x){
        locSce <- sceListAscDown[[which(names(sceListAscDown) == x)]]
        locSce <- locSce[which(row.names(locSce) %in% c("XBP1", "PRDM1", "IRF4",
                                                        "B2M", "ACTB", "PTPRC"))]
        locSce$Dataset <- x
        logcounts(locSce)
    }))))
plotSce$group <- unlist(lapply(sceListAscDown, function(x) x$group))
plotSce$cellType <- "ASC"
plotExpression(plotSce, features=row.names(plotSce),
               x="group", colour_by="cellType", ncol = 3)
dev.copy(pdf,'Results/Specific_vs_others/High_genes_post_qumi_and_downsamp.pdf', height = 8, width = 8)
dev.off()

plotExpression(plotSce, features="XBP1",
               x="group", colour_by="cellType", ncol = 3)
dev.copy(pdf,'Results/Specific_vs_others/XBP1_post_qumi_and_downsamp.pdf', height = 8, width = 8)
dev.off()


