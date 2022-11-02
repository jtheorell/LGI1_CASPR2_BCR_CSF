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

healthySce <- readRDS("Data/SmartSeq3_PBMC/SingleCellExpFiles/3_plasmablasts.rds")
dengueSce <- readRDS("Data/PlasmaBlast/SingleCellExpFiles/1_lowQual_excluded.rds")
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

#Now, the comparisons will be slightly different in nature, as the dengue dataset
#only has TPM. We will therefore have to calculate tpm for our dataset first
#(we only have CPM, which is different), and this needs to be done for the full
#dataset, before exclusions of non-protein coding transcripts. After this, we can
#log-transform. We will also have to call the data "counts" and "logcounts" for
#the comparative algorithm to work. 
#First, we import the early data
csfSceFull <- readRDS("Data/SingleCellExpFiles/csfSce_1_preQC.rds")
tpmAssay <- calculateTPM(csfSceFull)

tpmAssayReduced <- tpmAssay[which(rownames(tpmAssay) %in% rownames(csfSce)),
                            which(colnames(tpmAssay) %in% colnames(csfSce))]
identical(rownames(tpmAssayReduced), rownames(csfSce))
#TRUE
identical(colnames(tpmAssayReduced), colnames(csfSce))
#TRUE

#And then we can integrate this
tpm(csfSce) <- tpmAssayReduced

#Now, we will identify the transcripts that are shared among the three datasets.

dengueCsfGenes <- row.names(csfSce)[which(row.names(csfSce) %in% row.names(dengueSce))]
#We lose 352, or about 1.8% of the transcripts here. 
dengueHealthyCsfGenes <- dengueCsfGenes[which(dengueCsfGenes %in% row.names(healthySce))]
#Here another 76 transcripts, or 0.4% are lost. 

#Now, we reduce all three datasets to these transcripts. 
healthySce <- healthySce[which(row.names(healthySce) %in% dengueHealthyCsfGenes),]
dengueSce <- dengueSce[which(row.names(dengueSce) %in% dengueHealthyCsfGenes),]
csfSce <- csfSce[which(row.names(csfSce) %in% dengueHealthyCsfGenes),]

all(identical(row.names(healthySce), row.names(dengueSce)), 
    identical(row.names(csfSce), row.names(dengueSce)))
#TRUE
#So now, not only the genes are the same, but also the order of them. 
#Ready to take on the analyses: 

##################
#DENGUE VS CSF
##################
#Here, we are going to make a new assay; logTpm, for the dengue and csf datasets
assay(csfSce, "logTpm") <- log10(tpm(csfSce)+1)
assay(dengueSce, "logTpm") <- log10(tpm(dengueSce)+1)

#Now, we put the dataset together. We will have to name the tpm counts, and 
#the logTpm logcounts, for the algorithm downstream to work. 
csfDengueCounts <- cbind(tpm(csfSce), tpm(dengueSce))
csfDengueLog <- cbind(assay(csfSce, "logTpm"), assay(dengueSce, "logTpm"))
csfDengueSce <- SingleCellExperiment(assays = list("counts" = csfDengueCounts,
                                                   "logcounts" = csfDengueLog),
                                     rowData = rowData(csfSce))
csfDengueSce$group <- c(rep("AE", ncol(csfSce)), rep("dengue", ncol(dengueSce)))

#Now over to the analysis
marker.info <- scoreMarkers(csfDengueSce, csfDengueSce$group)
marker.info

#This can be done with any number of cell types, and as we in this case
#have only two, we can go for the specific one, and take the analysis from there
chosen <- marker.info[["AE"]]

##Here, we change the axis to more informative names. 
identical(row.names(chosen), row.names(csfDengueSce))
#TRUE

row.names(chosen) <- paste0(seq(1,nrow(csfDengueSce)),"_",rowData(csfDengueSce)$hgnc_symbol)
csfDengueSceNewNames <- csfDengueSce
row.names(csfDengueSceNewNames) <- row.names(chosen)

#Now, as a value in 0.5 corresponds to no signal, we make it possible to order 
#according to the largest divergence from this in both directions. As the
#values here generated are made considering multiple clusters at the same time, 
#and we only have two, median, min, max, mean etc are all giving the same value.
normalizedAUC <- chosen$median.AUC-0.5

summary(normalizedAUC)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.201693 -0.006547  0.000000  0.004721  0.017886  0.195465


orderedDengue <- chosen[order(abs(normalizedAUC), decreasing=TRUE),]
dir.create("Results/Comp_to_others/", recursive = TRUE)
plotExpression(csfDengueSceNewNames, features=rownames(orderedDengue)[1:20], 
               x="group", colour_by="group")
dev.copy(pdf,'Results/Comp_to_others/Top_20_features_by_AUC_dengue_vs_CSF.pdf', height = 20, width = 6)
dev.off()

#And now over to the other analysis: 
##################
#HEALTHY VS CSF
##################

#Here, things are easier, as we can just combine the two datasets as they are. 
csfHealthCounts <- cbind(counts(csfSce), counts(healthySce))
csfHealthLogcounts <- cbind(logcounts(csfSce), logcounts(healthySce))
csfHealthSce <- SingleCellExperiment(assays = list("counts" = csfHealthCounts,
                                                   "logcounts" = csfHealthLogcounts),
                                     rowData = rowData(csfSce))
csfHealthSce$group <- c(rep("AE", ncol(csfSce)), rep("healthy", ncol(healthySce)))

marker.info <- scoreMarkers(csfHealthSce, csfHealthSce$group)
marker.info

chosen <- marker.info[["AE"]]

identical(row.names(chosen), row.names(csfHealthSce))
#TRUE

row.names(chosen) <- paste0(seq(1,nrow(csfHealthSce)),"_",rowData(csfHealthSce)$hgnc_symbol)
csfHealthSceNewNames <- csfHealthSce
row.names(csfHealthSceNewNames) <- row.names(chosen)

#Now, as a value in 0.5 corresponds to no signal, we make it possible to order 
#according to the largest divergence from this in both directions. As the
#values here generated are made considering multiple clusters at the same time, 
#and we only have two, median, min, max, mean etc are all giving the same value.
normalizedAUC <- chosen$median.AUC-0.5

summary(normalizedAUC)
orderedHealthy <- chosen[order(abs(normalizedAUC), decreasing=TRUE),]

plotExpression(csfHealthSceNewNames, features=rownames(orderedHealthy)[1:20], 
               x="group", colour_by="group")
dev.copy(pdf,'Results/Comp_to_others/Top_20_features_by_AUC_healthy_vs_CSF.pdf', height = 20, width = 6)
dev.off()

#So now, we will find the intersection between these two datasets. What we will 
#do is to identify which of the 100 top genes that are present among the 100 top 
#genes of the other. 

topHealthy <- rownames(orderedHealthy)[1:100]
topDengue <- rownames(orderedDengue)[1:100]

dengueHealthyIntersection <- topHealthy[which(topHealthy %in% topDengue)]
denHelGeneNames <- gsub(".+_|", "", dengueHealthyIntersection)
#This includes 47 genes. 

##################m#
#MS
####################

#Now on to the MS dataset, which is not smartSeq and 
#therefore greater risk of generating spurious hits. 
msSce <- readRDS("Data/MS_10X/SingleCellExpFiles/3_plasmablasts.rds")
#We start by checking if the 47 genes are present here. 

length(which(rowData(msSce)$hgnc_symbol %in% denHelGeneNames))

#All are present. So then we can just run with the genes here that are common
#to the other datasets. 

msAndCsfCommon <- row.names(msSce)[which(row.names(msSce) %in% dengueHealthyCsfGenes)]
#Here, we go from a total of 19418 common genes down to 17489, which is a loss of
#1929 genes, or 10%. Still decent. 

csfSceMs <- csfSce[which(row.names(csfSce) %in% msAndCsfCommon),]
msSceCsf <- msSce[which(row.names(msSce) %in% msAndCsfCommon),]

identical(row.names(csfSceMs), row.names(msSceCsf))
#TRUE


csfMsCounts <- cbind(counts(csfSceMs), counts(msSceCsf))
csfMsLogcounts <- cbind(logcounts(csfSceMs), logcounts(msSceCsf))
csfMsSce <- SingleCellExperiment(assays = list("counts" = csfMsCounts,
                                                   "logcounts" = csfMsLogcounts),
                                     rowData = rowData(csfSceMs))
csfMsSce$group <- c(rep("AE", ncol(csfSceMs)), rep("MS", ncol(msSceCsf)))

marker.info <- scoreMarkers(csfMsSce, csfMsSce$group)
marker.info

chosen <- marker.info[["AE"]]

identical(row.names(chosen), row.names(csfMsSce))
#TRUE

row.names(chosen) <- paste0(seq(1,nrow(csfMsSce)),"_",rowData(csfMsSce)$hgnc_symbol)
csfMsSceNewNames <- csfMsSce
row.names(csfMsSceNewNames) <- row.names(chosen)

#Now, as a value in 0.5 corresponds to no signal, we make it possible to order 
#according to the largest divergence from this in both directions. As the
#values here generated are made considering multiple clusters at the same time, 
#and we only have two, median, min, max, mean etc are all giving the same value.
normalizedAUC <- chosen$median.AUC-0.5

summary(normalizedAUC)
orderedMS <- chosen[order(abs(normalizedAUC), decreasing=TRUE),]

plotExpression(csfMsSceNewNames, features=rownames(orderedMs)[1:20], 
               x="group", colour_by="group")
dev.copy(pdf,'Results/Comp_to_others/Top_20_features_by_AUC_MS_vs_CSF.pdf', height = 20, width = 6)
dev.off()

#And now to the fun: how many are still kept in this analysis?
#We need to start by chaging the names back. 
topMs <- row.names(orderedMs)[1:100]
msGeneNames <- gsub(".+_|", "", topMs)
msDengueHealthyIntersection <- msGeneNames[which(msGeneNames %in% denHelGeneNames)]


#So that is a wonderful 31 markers! 
#Now, before going on to trying to see patterns here, we will rank them according
#to their order of influence for all the three datasets. 
#We need to rename the full lists of gene names first.
allMSGeneNames <- gsub(".+_|", "", row.names(orderedMS))
allHealthyGeneNames <- gsub(".+_|", "", row.names(orderedHealthy))
allDengueGeneNames <- gsub(".+_|", "", row.names(orderedDengue))

topGenesList <- list("MS"= data.frame(orderedMS$median.AUC[which(allMSGeneNames %in% msDengueHealthyIntersection)],
                            row.names = allMSGeneNames[which(allMSGeneNames %in% msDengueHealthyIntersection)]),
                     "Healthy" = data.frame(orderedHealthy$median.AUC[which(allHealthyGeneNames %in% msDengueHealthyIntersection)],
                            row.names = allHealthyGeneNames[which(allHealthyGeneNames %in% msDengueHealthyIntersection)]),
                     "Dengue" = data.frame(orderedDengue$median.AUC[which(allDengueGeneNames %in% msDengueHealthyIntersection)],
                            row.names = allDengueGeneNames[which(allDengueGeneNames %in% msDengueHealthyIntersection)]))

#And now, we make the order the same for all three. 
topGenesListOrdered <- lapply(topGenesList, function(x){
    locX <- x[order(row.names(x)),]
    names(locX) <- row.names(x)[order(row.names(x))]
    locX
})

all(identical(names(topGenesListOrdered[[1]]), names(topGenesListOrdered[[2]])),
    identical(names(topGenesListOrdered[[3]]), names(topGenesListOrdered[[2]])))
#TRUE
#So then they can be recombined. 
topGenesDf <- as.data.frame(do.call("cbind", topGenesListOrdered))
row.names(topGenesDf) <- names(topGenesListOrdered[[1]])
topGenesDf$average <- rowMeans(topGenesDf)
topGenesDfOrdered <- topGenesDf[order(topGenesDf$average, decreasing = TRUE),]
write.csv(topGenesDfOrdered, "Results/Comp_to_others/TopGenes.csv")
