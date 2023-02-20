library(SingleCellExperiment)
library(scuttle)
library(scater)
library(scran)
library(quminorm)
library(BiocParallel)

sceList <- readRDS("Data/Comp_to_others/All_sce_common_genes_pre_all.rds")

#Here, we get the right rows
aeSceRows <- sceList$AE_CSF

#Here, we get the columns
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")
aeZoomed <- aeSce[,which(aeSce$clustersLouvain != "CD8T" | aeSce$Specific != "Not_tested")]
rawAe <- readRDS("Data/SingleCellExpFiles/csfSce_1_plus_flow_and_GLMPCA.rds")
rawAeAsc <- rawAe[,which(colnames(rawAe)%in% colnames(aeZoomed))]
aeSceTpm <- calculateTPM(rawAeAsc)
#Here, we have tested different values, to reach an optimal similarity in the expression
#of B2M, that in the original setting is too harshly affected. 
aeSceQumiMat <- quminorm(aeSceTpm, shape = 3, mc.cores = 7)

#Now, we remove all rows and columns not included in the asc dataset here used
aeSceQumiMatRed <- aeSceQumiMat[which(row.names(aeSceQumiMat) %in% rowData(aeSceRows)$ensembl_gene_id),]
aeSceQumiMatOrd <- aeSceQumiMatRed[match(rowData(aeSceRows)$ensembl_gene_id, row.names(aeSceQumiMatRed)),]

all(identical(rownames(aeSceQumiMatOrd), rowData(aeSceRows)$ensembl_gene_id),
    identical(colnames(aeSceQumiMatOrd), colnames(aeZoomed)))
#TRUE
rownames(aeSceQumiMatOrd) <- rownames(aeSceRows)

aeSceQumi <- SingleCellExperiment(assays = list("counts" = aeSceQumiMatOrd),
                                  rowData = rowData(aeSceRows), colData = colData(aeZoomed))
aeSceQumi$cellType <- aeSceQumi$clustersLouvain

#Here, it is included
sceList$AE_CSF <- aeSceQumi

#The same is now done for the healthy cells. 
rawHealthy <- readRDS("../External/Data/Healthy/3_post_Kotliarov.rds")

healthySceTpm <- calculateTPM(rawHealthy)
healthySceQumiMat <- quminorm(healthySceTpm, shape = 3, mc.cores = 7)
healthySceQumiRed <- healthySceQumiMat[which(row.names(healthySceQumiMat) %in% rowData(sceList$Healthy_PBMC)$ensembl_gene_id),
                                       which(colnames(healthySceQumiMat) %in% colnames(sceList$Healthy_PBMC))]

healthySceQumiOrd <- healthySceQumiRed[match(rowData(sceList$Healthy_PBMC)$ensembl_gene_id, row.names(healthySceQumiRed)),]

all(identical(rownames(healthySceQumiOrd), rowData(sceList$Healthy_PBMC)$ensembl_gene_id),
    identical(colnames(healthySceQumiOrd), colnames(sceList$Healthy_PBMC)))
#TRUE
rownames(healthySceQumiOrd) <- rownames(sceList$Healthy_PBMC)

median(colSums(healthySceQumiOrd))
#14235
healthySceQumi <- sceList$Healthy_PBMC
assays(healthySceQumi) <- list("counts" = healthySceQumiOrd)

#And these are integrated
sceList$Healthy_PBMC <- healthySceQumi

#And this data is saved. 
saveRDS(sceList, "Data/Comp_to_others/All_sce_common_genes_post_qumi_including_B-cells.rds")

medianSums <- lapply(sceList, function(x) median(colSums(as.matrix(counts(x)))))

minMedSum <- round(min(unlist(medianSums)))
#1449
madVal <- mad(colSums(as.matrix(counts(sceList[[which.min(unlist(medianSums))]]))))
#418.0932

minVal <- min(colSums(as.matrix(counts(sceList[[which.min(unlist(medianSums))]]))))
#63

#Here, we are now going to exclude a number of genes, to remove the less interesting
#part of the data. It should of course be done before hte median calculations, but
#it is not in 15a, so it will not be here eihter. 
ascExprs <- row.names(readRDS("Results/Specific_vs_not/EdgeR_result.rds"))

sceListRed <- lapply(sceList, function(x) x[which(row.names(x) %in% ascExprs),])

#Are they all identical now, row-wise? 
all(sapply(sceListRed, function(x) identical(rownames(sceListRed$AE_CSF), rownames(x))))
#TRUE

#Now, we will loop over all the dataset, essentially reproducing the analysis
#from the specific vs others downsampling, but including the B-cells. 
#We need colnames for the covid dataset again. 
colnames(sceListRed$Covid_PBMC) <- paste0("Cell_", seq_len(ncol(sceListRed$Covid_PBMC)))
set.seed(101)
sceListDown <- lapply(sceListRed, function(x){
    locColSums <- colSums(as.matrix(counts(x)))
    locMedian <- median(locColSums)
    locMad <- mad(locColSums)
    locColSumsCenter <- locColSums-locMedian
    locColSumsScale <- locColSumsCenter/locMad
    locColSumsRefit <- round((locColSumsScale*madVal)+minMedSum)
    #And now, we let all values that are below the lowest value in the sample
    #with the lowest median be set to that value.
    locColSumsRefit[which(locColSumsRefit < minVal)] <- minVal
    locCountDat <- do.call("cbind",bplapply(names(locColSumsRefit), function(y){
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
    colnames(locCountDat) <- names(locColSumsRefit)
    
    locSceDown <- SingleCellExperiment(assays = list("counts" = locCountDat), 
                                      colData = colData(x))
    
    locSceDown <- computeSumFactors(locSceDown)
    locSceDown <- logNormCounts(locSceDown)
    locSceDown
})

#Now, we will separate the non-specific in our data. 
sceListDown$AE_pos <- sceListDown$AE_CSF[,which(sceListDown$AE_CSF$Specific == "TRUE")]
sceListDown$AE_neg <- sceListDown$AE_CSF[,which(sceListDown$AE_CSF$Specific == "FALSE")]
sceListDown$AE_CSF <- NULL

saveRDS(sceListDown,"Data/Comp_to_others/Downsampled_B_and_LGI1-C2_neg_included.rds")


