library(flowCore)
library(flowSpecs)
library(ggplot2)
library(DepecheR)
library(BiocParallel)
library(uwot)
library(FNN)
library(robustbase)

fsList <- readRDS(fsList, "Data/Cytometry/raw_flow_cytometry_data.rds")

#Now, we test what a reasonable transformation cofactor is
transCoFacs <- rep(1000, times = length(5:16))
names(transCoFacs) <- colnames(fsListConsistent[[1]])[5:16]
transCoFacs[c("CD3")] <- 2000
transCoFacs[c("CD27")] <- 4000
testFile <- fsListConsistent[[2]][[3]]
testFileTrans <- arcTrans(testFile, transNames = colnames(testFile)[5:16], 
                          transCoFacs = transCoFacs)
hist(exprs(testFileTrans)[,"CD27"], breaks = 200)

#And now, we normalize the peaks. 
normOrNot <- rep(TRUE, length(colnames(fsListConsistent[[1]])))
names(normOrNot) <- colnames(fsListConsistent[[1]])
normOrNot[c("SSC-A", "SSC-H",  "TIME", "CD3")] <- FALSE
fsListNorm <- lapply(fsListConsistent, function(x){
    peakNorm(fs = x, ctrlPos = 4, standFF = fsListConsistent[[2]][[4]], 
             transNames = colnames(x)[5:16], transCoFacs = transCoFacs, 
             normOrNot = normOrNot)
})

fsListTrans <- lapply(fsListNorm, arcTrans, 
                      transNames = colnames(fsListNorm[[1]])[5:16], 
                      transCoFacs = transCoFacs)

#And here, we view the standards from all batches. 
standFs <- flowSet(lapply(fsListConsistent, "[[", 4))
standFsTrans <- arcTrans(standFs, transNames = colnames(standFs)[5:16], 
                         transCoFacs = transCoFacs)

dir.create("Diagnostics/Norm", recursive = TRUE)

subDf <- flowSet2LongDf(standFsTrans)

for(i in seq_len(ncol(subDf)-3)){
    locName <- colnames(subDf)[i]
    locDf <- subDf[,c(i, which(colnames(subDf) == "acqDate"))]
    colnames(locDf)[1] <- "var"
    ggplot(locDf, aes(x=var, fill=acqDate)) +
        geom_density(alpha = 0.8) + xlab(locName) +ylab("")
    ggsave(paste0("Diagnostics/Norm/", locName, "_preNorm.pdf"))
}

standFs <- flowSet(lapply(fsListTrans, "[[", 4))

subDf <- flowSet2LongDf(standFs)

for(i in seq_len(ncol(subDf)-3)){
    locName <- colnames(subDf)[i]
    locDf <- subDf[,c(i, which(colnames(subDf) == "acqDate"))]
    colnames(locDf)[1] <- "var"
    ggplot(locDf, aes(x=var, fill=acqDate)) +
        geom_density(alpha = 0.8) + xlab(locName) +ylab("")
    ggsave(paste0("Diagnostics/Norm/", locName, "_afterNorm.pdf"))
}

##########################

#Here, we are now creating two new parameters, namely CD4 and IgD, from the 
#combined data currently at hand. We will simply define any cell positive for 
#the double parameter and CD19 as IgD, and any that is negative for CD19 as CD4.
#To get as physiological as possible, we will gate based on the PBMC control
#from the same individual in all cases. 

CD19Peaks <- lapply(fsListTrans, function(x) flowSpecs:::peakIdenti(
    exprs(x[[3]][,"CD19"]), volThresh = 0.02, returnStats = TRUE))

#In all cases, two peaks are found with these settings. We use the upper 
#extension of the lower peak to define what is positive and negative for CD19 and
#this which signal that comes from CD4 or IgD.
#In this work, we also need to deal with the fact that we get artificially 
#low noise in these two parameters. This is especially acute for IgD, where only
#a minor fraction of all events will get a value, and therefore the majority of
#the data will have the value 0. For this reason, we are randomly copying values
#from the negative peak in the PBMC file from the same individual. Hence this
#second set of peaks:
CD4_IgDPeaks <- lapply(fsListTrans, function(x) flowSpecs:::peakIdenti(
    exprs(x[[3]][,"CD4_IgD"]), returnStats = TRUE))

all(unlist(lapply(lapply(CD4_IgDPeaks, "[[", 1), length)) == 2)
#TRUE
#This shows that two peaks were found in all cases, not unexpectedly.

fsListTransAppended <- lapply(seq_along(fsListTrans), function(x){
    locRefFF <- exprs(fsListTrans[[x]][[3]])
    locCD4_IgDNeg <- locRefFF[which(
        locRefFF[,"CD4_IgD"] <= CD4_IgDPeaks[[x]]$Width[[1]][2]),"CD4_IgD"]
    locFs <- fsApply(fsListTrans[[x]], function(y){
        locExprs <- exprs(y)
        #Here, we first create a random negative peak modeled on the negative
        #peak in the PBMC CD4_IgD vector from the particular individual. 
        #To be very certain not to make the noise mx with the data, we create it
        #with a standard deviation being a fifth of the distance between the 
        #negative peak and upper border.
        
        locCD4 <- rnorm(nrow(y), mean = CD4_IgDPeaks[[x]]$PeakPos[1],
                        sd = (CD4_IgDPeaks[[x]]$Width[[1]][2]-
                            CD4_IgDPeaks[[x]]$PeakPos[1])/5)
        locCD4[which(locExprs[,"CD19"] <= CD19Peaks[[x]]$Width[[1]][2])] <- 
            locExprs[which(locExprs[,"CD19"] <= CD19Peaks[[x]]$Width[[1]][2]),"CD4_IgD"]
        locIgD <- rnorm(nrow(y), mean = CD4_IgDPeaks[[x]]$PeakPos[1],
                        sd = (CD4_IgDPeaks[[x]]$Width[[1]][2]-
                                  CD4_IgDPeaks[[x]]$PeakPos[1])/5)
        locIgD[which(locExprs[,"CD19"] > CD19Peaks[[x]]$Width[[1]][2])] <- 
            locExprs[which(locExprs[,"CD19"] > CD19Peaks[[x]]$Width[[1]][2]),"CD4_IgD"]
        newCols <- data.frame("CD4" = locCD4, "IgD" = locIgD)
        appendedFF <- flowSpecs:::appendFFCols(y, as.matrix(newCols))
    })
})

##########################
#Now we will remove the non-sorted files from the ungated CSF files, before
#starting to gate. 

#And now we are ready for the next step: to identify the cells we sorted. 
#Import the index files
indexList <- readRDS("Data/Cytometry/indexList.rds")
#Now, we check that the number of rows in the index files correspond to the 
#number of rows in the CSF flow files. It has previously been checked that the
#index files, which are made up of multiple files from the acquisition, and the
#fcs files, which are concatenated for the same reason, have the same order.
identical(unname(unlist(lapply(indexList, nrow))), 
          unlist(unname(lapply(fsListTransAppended, function(x) fsApply(x, nrow)[[1]]))))
#TRUE


csfAllList <- lapply(fsListTransAppended, "[[", 1)
names(csfAllList) <- lapply(lapply(fsListTransAppended, sampleNames), "[[", 1)

csfAllMatList <- lapply(csfAllList, exprs)

#Here, we do not only exclude all non-index events and bake the data together, 
#we also create a name vector which corresponds to the names in the other datasets
csfIndexMatList <- lapply(seq_along(csfAllMatList), function(x){
    locNonRows <- which(indexList[[x]]$Index == "")
    locResIndexed <- csfAllMatList[[x]][-locNonRows,]
    locIndexIndexed <- indexList[[x]][-locNonRows,]
    
    #Now, we create the names, which is not totally trivial, since the indexes need
    #to be converted from a "A1" to a "A01" format.
    locIndexes <- locIndexIndexed$Index
    lengthOfTwo <- which(nchar(locIndexes) == 2)
    locIndexes[lengthOfTwo] <- paste0(substr(locIndexes[lengthOfTwo], 1, 1), 0, 
                                      substr(locIndexes[lengthOfTwo], 2, 2))
    locNames <- paste0(substr(names(csfAllList)[[x]],3,9),locIndexIndexed$Plate,"_",locIndexes)
    locRes <- cbind(locResIndexed, locIndexIndexed, "Cell" = locNames)
})
names(csfIndexMatList) <- names(csfAllList)

##################
#Now, we are going to cluster the data, based on the gated PBMC and CSF files. 

#These are now, for convenience, reconverted to fcs files, and then immediately
#converted again to the long format used by depeche
csfIndexFS <- flowSet(lapply(csfIndexMatList, function(x) 
    flowFrame(as.matrix(x[,1:19]))))
csfIndexDf <- flowSet2LongDf(csfIndexFS, idInfo = list(
    "Tissue" = "JR...._1_|_.{3,5}\\.fcs", "Donor" = "|_1_.+"))
csfIndexDf$Type <- "Index"


lapply(fsListTransAppended, fsApply, dim)

#This shows that the lowest number of cells in any of the gated files is 5000
#So now we make a new, large dataframe based on 5000 CSF and 5000 PBMC events
#from each donor.
set.seed(101)
clusteringFile <- do.call("rbind", lapply(fsListTransAppended, function(x){
    locFS5000 <- fsApply(x[2:3], function(y) 
        y[sample(1:nrow(y),5000),])
    flowSet2LongDf(locFS5000, idInfo = list("Tissue" = "JR...._1_|_.{3,5}\\.fcs",
                                        "Donor" = "|_1_.+"))
}))
clusteringFile$Type <- "Gated"

combinedFile <-rbind(clusteringFile, csfIndexDf)

#So now, we try to cluster this using Depeche. To do so, we need to
#rescale the FSC and SSC parameters. They get the max and min values
#based on the median of the rest of the data. 
depInData <- combinedFile[,c("FSC.A","SSC.A","CD20","CD56","CD38",
                               "CD138","CD19","CD16","CD8","CD3","CD27",
                               "CD4","IgD")]
medMax <- median(apply(depInData[,3:ncol(depInData)], 2, max))
medMin <- median(apply(depInData[,3:ncol(depInData)], 2, min))

scaledSC <- dScale(depInData[,1:2], control = c(medMin, medMax))
depInData[,1:2] <- scaledSC

#As CD19 plays a central part in separating the cells of interest from the rest
#here, we will increase its weight twofold. Even without this, it works reasonably
#well, but 
depInData[,"CD19"] <- depInData[,"CD19"]*2

#In this setup, there is considerable leakage between the CD3 and CD138 channels. 
#For this reason, we will create a gate that separates them.
plot(depInData$CD3, depInData$CD138)
abline(h = 2, v= 2)
abline(h = 3, v= 3)
abline(h = 3.8, v= 4)
abline(h = 4.5, v= 4.8)
abline(h = 6)
CD3orCD138 <- rep("CD3", nrow(depInData))
CD3 <- depInData$CD3
CD138 <- depInData$CD138

CD3orCD138[which((CD138 > 2 & CD3 < 2) |
                     (CD138 > 3 & CD3 < 3) |
                     (CD138 > 3.8 & CD3 < 4) |
                     (CD138 > 4.5 & CD3 < 4.8) |
                     (CD138 > 6))] <- "CD138"
CD3orCD138[which(CD138 < 2 & CD3 < 1.3)] <- "Other"
dir.create("Diagnostics/Cytometry")
dColorPlot(CD3orCD138, xYData = data.frame(depInData$CD3, depInData$CD138),
           plotName = "Diagnostics/Cytometry/CD3_vs_CD138_with_positivity")

#This worked very well. So now, we will define all CD3+ as negative
#for CD138 and vise versa. 
depInData$CD3[which(CD3orCD138 == "CD138")] <- min(depInData$CD3)
depInData$CD138[which(CD3orCD138 == "CD3")] <- min(depInData$CD138)

depModel <- depeche(depInData,
                    plotDir = "Diagnostics/Depeche")
saveRDS(depModel, "Data/Cytometry/DepModel.rds")

flowUmap <- umap(depInData)
saveRDS(flowUmap, "Data/Cytometry/flowUmap.rds")

dColorPlot(depModel[[1]], xYData = flowUmap, plotDir = "Diagnostics/Cytometry")
dColorPlot(depInData, xYData = flowUmap, plotDir = "Diagnostics/Cytometry")
#Here, we rename the clusters: 
clustersNamed <- sapply(depModel[[1]], switch, "1" = "CD4T", "2" = "NK",
                        "3" = "CD8T", "4" = "B", "5" = "CD4T", "6" = "ASC")
dColorPlot(clustersNamed, xYData = flowUmap, 
           plotName = "Diagnostics/Cytometry/Clusters")
dDensityPlot(flowUmap, 
             idsVector = depModel[[1]], 
             plotDir = "Diagnostics/Cytometry/ClusterDensities")
dDensityPlot(flowUmap, 
             idsVector = paste0(clusteringFile$Donor, "_", clusteringFile$Tissue), 
             plotDir = "Diagnostics/Cytometry/SampleDensities")

##################
#And this now means that we can export the cluster information to the index 
#files. 

csfIndexDf <- do.call("rbind", csfIndexMatList)

#Before it is included, we just check that all is right: 
identical(substr(csfIndexDf$Cell, 1, 4), 
          substr(combinedFile$Donor[which(combinedFile$Type == "Index")],3,6))
#TRUE, meaning that at least the donors are in the same order, and as the 
#number of cells also are identical, there is no reason to believe that #
#there should be a change of order internal to each donor, as no such thing
#has taken place. 

csfIndexDf$Cell_type <- clustersNamed[which(combinedFile$Type == "Index")]

table(csfIndexDf$Cell_type)
#ASC    B CD4T CD8T   NK 
#316  277 1069  264   89

#Which looks simply beautiful. 
write.csv(csfIndexDf, "Data/Cytometry/flowDataPlusIndexAndcellType.csv")
