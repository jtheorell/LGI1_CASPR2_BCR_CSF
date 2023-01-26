#Here comes now a combined analysis. Previous analyses have shown that, at least
#when using Batchelor, the batch correction completely removes all differences 
#between the ASC in the CSF and in the periphery, possibly as they are so comparatively
#homogenous that any nearest neighbor analysis will anchor the whole dataset to itself.
#Now the strategy here is instead to use two separate peripheral dataset: one from
#healthy individuals and one from Dengue infected ditto, and not do any batch
#correction but only consider the hits that are consistent between the two 
#datasets. This of course builds on the assumption that there are as large batch
#differences between the dengue and healthy datasets as it is between our dataset
#and the two others, but that will have to do.
library(SingleCellExperiment)
library(scran)
library(scater)
library(BiocParallel)
library(quminorm)
library(edgeR)

#In these analyses, we will stop calling our data "csfSce" as there are two such
#datasets now. Instead aeSce for "Autoimmune encephalitis"
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

#As can be seen here, we do not have any ASCs from 1227 with any other specificty than LGI1
#     FALSE Not_tested TRUE
#1166    14         92   75
#1227     0         30   13
#1284     4         10   31
#This means that all comparisons between the specific and the non-specific fall for the
#internal analysis, and it motivates our jump into the unknown and all these external
#datasets
covidSce <- readRDS("Data/Comp_to_others/Covid/SingleCellExpFiles/2_plasmablasts.rds")
healthySce <- readRDS("Data/Comp_to_others/Healthy/SingleCellExpFiles/3_PB.rds")
influensaSce <- readRDS("Data/Comp_to_others/Influensa/SingleCellExpFiles/3_PB.rds")
msRameshSce <- readRDS("Data/Comp_to_others/MS_Ramesh/SingleCellExpFiles/3_PB.rds")
msSchafSce <- readRDS("Data/Comp_to_others/MS_Schafflick/SingleCellExpFiles/3_PB.rds")
dengueSce <- readRDS("Data/Comp_to_others/Dengue/SingleCellExpFiles/PB.rds")

##################
#SHARED GENES
##################
#We will now have to do something in a bit of a convoluted way. First, we will aggregate
#the cells into the donors, do the exclusions and after that we will identify the
#commonly expressed genes. After this, we will go back to the single-cell resolution, 
#do the subsampling and then the actual DE analysis. 

rawSceList <- list(aeSce, covidSce, healthySce, influensaSce, msRameshSce, msSchafSce, dengueSce)

#First, we filter out genes that are not of common interest. These are the ones that
#are lowly expressed. However, as we have very different counts in the datasets, 
#we will start by defining a "standard" count that we expect from an interesting
#gene in any of the datasets. 
#The smallest dataset has the following number of counts: 
#min(unlist(lapply(rawSceList, function(x){sum(as.vector(counts(x)))})))
#335598
#If we were to pick any gene that has a total number of 3 counts in this dataset, 
#we of course need to scale this for the other datasets, to get a similar sensitivity
#THe highest count is namely: 
#max(unlist(lapply(rawSceList, function(x){sum(as.vector(counts(x)))})))
#197525048, or 588 times as large. 
#This means that if we for the smallest datasets require a count of 3, we will require a count
#of 588 for the largest datasets. 
#We will apply this rule, but include all genes that are identified, as we are 
#interested in diversity. 
scaleVal <- 3/335598

allFilteredNames <- lapply(rawSceList, function(x){
    sumMat <- do.call("cbind", lapply(unique(x$donor), function(y) 
        rowSums(counts(x[,which(x$donor == y)]))))
    #Here, we select all genes with expression in at least three donors
    #and with a sum for at least one donor that exceeds the experiment-specific
    #threshold
    keep <- which(rowMax(sumMat) > sum(as.vector(counts(x)))*scaleVal &
                      rowSums(sumMat != 0) > 2)
    rowData(x)$hgnc_symbol[keep]
})

allExpressedHgnc <- unique(unlist(allFilteredNames))
length(allExpressedHgnc)

#Interestingly, the number of shared transcripts is considerably larger with the 
#hgnc symbols. It is however also so, that the hgnc symbols are not necessarily
#unique, but can map to multiple ensembl ids. If we only look at the 
#ensembl ids, however, we lose genes like TNF and IFNGR2, which we would like to include for 
#downstream analyses. We will now check how many genes we would lose if we exlude
#all tuplicated hgnc symbols. 
duplicatedHgnc <- unique(unlist(lapply(rawSceList, function(x) 
    rowData(x)$hgnc_symbol[which(duplicated(rowData(x)$hgnc_symbol))])))

length(duplicatedHgnc)
#20
ducplicatedSharedHgnc <- duplicatedHgnc[which(duplicatedHgnc %in% allExpressedHgnc)]
ducplicatedSharedHgnc
# ""             "PINX1"        "TBCE"         "ABCF2"        "POLR2J3"      "GGT1"         NA            
# "TNFRSF10A-DT" "PRPF31" 
#These are removed from the list
allExpressedHgnc <- allExpressedHgnc[-which(allExpressedHgnc %in% ducplicatedSharedHgnc)]

#Now, which of these are common to all datasets?
allCommonHgnc <- Reduce(intersect, lapply(rawSceList, function(x) rowData(x)$hgnc_symbol))
allCommonExprHgnc <- allExpressedHgnc[which(allExpressedHgnc %in% allCommonHgnc)]

write.csv(allCommonExprHgnc , "Data/Comp_to_others/commonExpressedHGNC.csv", row.names = FALSE)

#Now, we reduce all seven datasets to these genes, change the row names
#and order them. 
#AE
aeSceTemp <- aeSce[which(rowData(aeSce)$hgnc_symbol %in% allCommonExprHgnc),]
row.names(aeSceTemp) <- rowData(aeSceTemp)$hgnc_symbol
aeSceTemp <- aeSceTemp[order(row.names(aeSceTemp)),]

#Covid
covidSceTemp <- covidSce[which(rowData(covidSce)$hgnc_symbol %in% allCommonExprHgnc),]
row.names(covidSceTemp) <- rowData(covidSceTemp)$hgnc_symbol
covidSceTemp <- covidSceTemp[order(row.names(covidSceTemp)),]

#Dengue
dengueSceTemp <- dengueSce[which(rowData(dengueSce)$hgnc_symbol %in% allCommonExprHgnc),]
row.names(dengueSceTemp) <- rowData(dengueSceTemp)$hgnc_symbol
dengueSceTemp <- dengueSceTemp[order(row.names(dengueSceTemp)),]

#Healthy
healthySceTemp <- healthySce[which(rowData(healthySce)$hgnc_symbol %in% allCommonExprHgnc),]
row.names(healthySceTemp) <- rowData(healthySceTemp)$hgnc_symbol
healthySceTemp <- healthySceTemp[order(row.names(healthySceTemp)),]

#Influensa
influensaSceTemp <- influensaSce[which(rowData(influensaSce)$hgnc_symbol %in% allCommonExprHgnc),]
row.names(influensaSceTemp) <- rowData(influensaSceTemp)$hgnc_symbol
influensaSceTemp <- influensaSceTemp[order(row.names(influensaSceTemp)),]

#MS-Ramesh
msRameshSceTemp <- msRameshSce[which(rowData(msRameshSce)$hgnc_symbol %in% allCommonExprHgnc),]
row.names(msRameshSceTemp) <- rowData(msRameshSceTemp)$hgnc_symbol
msRameshSceTemp <- msRameshSceTemp[order(row.names(msRameshSceTemp)),]

#MS-Schafflick
msSchafSceTemp <- msSchafSce[which(rowData(msSchafSce)$hgnc_symbol %in% allCommonExprHgnc),]
row.names(msSchafSceTemp) <- rowData(msSchafSceTemp)$hgnc_symbol
msSchafSceTemp <- msSchafSceTemp[order(row.names(msSchafSceTemp)),]

#Are they now all identical?
all(sapply(list(row.names(covidSceTemp),
                row.names(dengueSceTemp),
                row.names(healthySceTemp),
                row.names(influensaSceTemp),
                row.names(msRameshSceTemp),
                row.names(msSchafSceTemp)), 
           FUN = identical, row.names(aeSceTemp)))
#TRUE
#So now, not only the genes are the same, but also the order of them. 
aeSce <- aeSceTemp
covidSce <- covidSceTemp
dengueSce <- dengueSceTemp
healthySce <- healthySceTemp
influensaSce <- influensaSceTemp
msRameshSce <- msRameshSceTemp
msSchafSce <- msSchafSceTemp 

allDatListRaw <- list(aeSce[,which(aeSce$Specific == "TRUE")], 
                      aeSce[,which(aeSce$Specific == "FALSE")],
                      covidSce, dengueSce, healthySce, influensaSce, msRameshSce, msSchafSce)
names(allDatListRaw) <- c("AE_specific", "AE_not_specific", "Covid", "Dengue", "Healthy", 
                       "Influensa", "MS_Ramesh", "MS_Schafflick")

saveRDS(allDatListRaw, "Data/Comp_to_others/All_sce_common_genes_pre_all.rds")

#Now, before starting any downsampling, we need to save a picture showing why 
#it is needed. 
plotSce <- SingleCellExperiment(
    assays = list("logcounts" = do.call("cbind",lapply(names(allDatListRaw), function(x){
        locSce <- allDatListRaw[[which(names(allDatListRaw) == x)]]
        locSce <- locSce[which(row.names(locSce) %in% c("XBP1", "PRDM1", "IRF4",
                                                        "B2M", "ACTB", "PTPRC")),]
        locSce$Dataset <- x
        logcounts(locSce)
    }))))
plotSce$group <- unlist(sapply(names(allDatListRaw), function(x){
    rep(x, ncol(allDatListRaw[[x]]))
}))

plotExpression(plotSce, features=row.names(plotSce), 
               x="group", colour_by="group", ncol = 3)
dev.copy(pdf,'Results/Comp_to_others/High_genes_pre_qumi_and_downsamp.pdf', height = 8, width = 8)
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
aeSceQumiMatRed <- aeSceQumiMat[which(row.names(aeSceQumiMat) %in% rowData(aeSce)$ensembl_gene_id),
                                which(colnames(aeSceQumiMat) %in% colnames(aeSce))]
aeSceQumiMatOrd <- aeSceQumiMatRed[match(rowData(aeSce)$ensembl_gene_id, row.names(aeSceQumiMatRed)),]

all(identical(rownames(aeSceQumiMatOrd), rowData(aeSce)$ensembl_gene_id),
    identical(colnames(aeSceQumiMatOrd), colnames(aeSce)))
#TRUE
rownames(aeSceQumiMatOrd) <- rownames(aeSce)

#And here the qumi sum
median(colSums(aeSceQumiMatOrd))
#25060

aeSceQumi <- aeSce
assays(aeSceQumi) <- list("counts" = aeSceQumiMatOrd)

#And for systematics sake, we do the same for the smartseq3 dataset
rawHealthy <- readRDS("Data/Comp_to_others/Healthy/SingleCellExpFiles/3_PB.rds")

healthySceTpm <- calculateTPM(rawHealthy)
healthySceQumiMat <- quminorm(healthySceTpm, shape = 3, mc.cores = 7)
healthySceQumiRed <- healthySceQumiMat[which(row.names(healthySceQumiMat) %in% rowData(healthySce)$ensembl_gene_id),
                                which(colnames(healthySceQumiMat) %in% colnames(healthySce))]

healthySceQumiOrd <- healthySceQumiRed[match(rowData(healthySce)$ensembl_gene_id, row.names(healthySceQumiRed)),]

all(identical(rownames(healthySceQumiOrd), rowData(healthySce)$ensembl_gene_id),
    identical(colnames(healthySceQumiOrd), colnames(healthySce)))
#TRUE
rownames(healthySceQumiOrd) <- rownames(healthySce)

median(colSums(healthySceQumiOrd))
#7967
healthySceQumi <- healthySce
assays(healthySceQumi) <- list("counts" = healthySceQumiOrd)

#Now, we are going to downsample all datasets
allDatList <- list(aeSceQumi, covidSce, dengueSce, healthySceQumi, influensaSce, msRameshSce, msSchafSce)
names(allDatList) <- c("AE", "Covid", "Dengue", "Healthy", 
                       "Influensa", "MS_Ramesh", "MS_Schafflick")
saveRDS(allDatList, "Data/Comp_to_others/All_sce_common_genes_qumi.rds")

#Here, we plot an intermediate version
allDatListQumi <- allDatList
allDatListQumi$AE <- computeSumFactors(allDatListQumi$AE)
allDatListQumi$AE <- logNormCounts(allDatListQumi$AE)
allDatListQumi$Healthy <- computeSumFactors(allDatListQumi$Healthy)
allDatListQumi$Healthy <- logNormCounts(allDatListQumi$Healthy)

allDatListQumi$AE_specific <- allDatListQumi$AE[,which(allDatListQumi$AE$Specific == "TRUE")]
allDatListQumi$AE_not_specific <- allDatListQumi$AE[,which(allDatListQumi$AE$Specific == "FALSE")]

allDatListQumi <- allDatListQumi[-1]

plotSce <- SingleCellExperiment(
    assays = list("logcounts" = do.call("cbind",lapply(names(allDatListQumi), function(x){
        locSce <- allDatListQumi[[which(names(allDatListQumi) == x)]]
        locSce <- locSce[which(row.names(locSce) %in% c("XBP1", "PRDM1", "IRF4",
                                                        "B2M", "ACTB", "PTPRC"))]
        locSce$Dataset <- x
        logcounts(locSce)
    }))))
plotSce$group <- unlist(sapply(names(allDatListQumi), function(x){
    rep(x, ncol(allDatListQumi[[x]]))
}))
plotExpression(plotSce, features=row.names(plotSce),
               x="group", colour_by="group", ncol = 3)
dev.copy(pdf,'Results/Comp_to_others/High_genes_post_qumi_pre_downsamp.pdf', height = 8, width = 8)
dev.off()


medianSums <- lapply(allDatList, function(x) median(colSums(as.matrix(counts(x)))))

minMedSum <- round(min(unlist(medianSums)))
#4334
madVal <- mad(colSums(as.matrix(counts(allDatList[[which.min(unlist(medianSums))]]))))
#1794.687

minVal <- min(colSums(as.matrix(counts(allDatList[[which.min(unlist(medianSums))]]))))
#1811
#So this is the amount we will downsample to. There are of course multiple ways of 
#doing this, but what I will do, as the matrices are not overwhelmingly large, is
#to expand each cell to a 0-1 vector, and then sample this vector. However, adding
#another layer of complexity, we want the median within each dataset to have the
#same number of counts,but that is not the same as saying that all cells should
#have this count. For this reason, we will first identify the median cell and then scale
#the values based on the MAD around this. 
#We are going to run this 100 times
dir.create("Results/Comp_to_others/Bootstrap")
dir.create("Results/Comp_to_others/Bootstrap_DE")
dir.create("Data/Comp_to_others/SCEs")
for(i in 1:100){
    print(paste0("Now, iteration ", i, " is starting"))
    allCountListDown <- lapply(allDatList, function(x){
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
    allDatListDown <- lapply(names(allDatList), function(x){
        locSce <- allDatList[[which(names(allDatList) == x)]]
        assays(locSce) <- list("counts" = allCountListDown[[which(names(allCountListDown) == x)]])
        locSce <- computeSumFactors(locSce)
        locSce <- logNormCounts(locSce)
        locSce
    })
    names(allDatListDown) <- names(allDatList)
    #Now, we split out the specific and the not specific cells
    allDatListDivided <- c("AE_specific" = allDatListDown[[1]][,which(allDatListDown[[1]]$Specific == "TRUE")],
                           "AE_not_specific" = allDatListDown[[1]][,which(allDatListDown[[1]]$Specific == "FALSE")],
                           allDatListDown[2:length(allDatListDown)])
    allDatListDivided <- lapply(names(allDatListDivided), function(x){
        locSce <- allDatListDivided[[which(names(allDatListDivided) == x)]]
        locSce$group <- x
        locSce
    })
    names(allDatListDivided) <- c("AE_specific", "AE_not_specific", names(allDatListDown[2:length(allDatListDown)]))
    ###################
    #DATASET CONSTRUCTION
    ###################
    datList <- allDatListDivided[-grep("AE", names(allDatListDivided))]
    nameList <- names(datList)
    aeSceSpec <- allDatListDivided$AE_specific
    sceList <- lapply(seq_along(datList), function(x){
        locDat <- datList[[x]]
        locCounts <- cbind(counts(aeSceSpec), counts(locDat))
        locSce <- SingleCellExperiment(assays = list("counts" = locCounts),
                                       rowData = rowData(allDatListDivided$Influensa))
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
    saveRDS(deResList, paste0("Results/Comp_to_others/Bootstrap_DE/All_DE_results_", i, ".rds"))
}

#So now, we have an as internally comparable dataset as possible.
#This last bootstrap is saved for later use. Here, we both save a full version, and a version where
#the AE data has been split and the ones without known specificity are excluded. 

saveRDS(allDatListDivided,
        paste0("Data/Comp_to_others//Downsampled_specificity_included_boot_100.rds"))

plotSce <- SingleCellExperiment(
    assays = list("logcounts" = do.call("cbind",lapply(names(allDatListDivided), function(x){
        locSce <- allDatListDivided[[which(names(allDatListDivided) == x)]]
        locSce <- locSce[which(row.names(locSce) %in% c("XBP1", "PRDM1", "IRF4",
                                                        "B2M", "ACTB", "PTPRC"))]
        locSce$Dataset <- x
        logcounts(locSce)
    }))))
plotSce$group <- unlist(lapply(allDatListDivided, function(x) x$group))
plotExpression(plotSce, features=row.names(plotSce),
               x="group", colour_by="group", ncol = 3)
dev.copy(pdf,'Results/Comp_to_others/High_genes_post_qumi_and_downsamp.pdf', height = 8, width = 8)
dev.off()

#It is annoying that the ACTB looks like it has been choped off for our data, but
#that it a cosequence of the fact that ACTB seems to be a dominant feature in our
#dataset and thus it is often #penalised" in the qumi procedure. 




