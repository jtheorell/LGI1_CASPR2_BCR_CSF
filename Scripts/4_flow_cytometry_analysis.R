library(flowCore)
library(flowSpecs)

fsList <- readRDS("Data/Cytometry/raw_flow_cytometry_data.rds")

#Now, we test what a reasonable transformation cofactor is
transCoFacs <- rep(1000, times = length(5:16))
names(transCoFacs) <- colnames(fsList[[1]])[5:16]
transCoFacs[c("CD3")] <- 2000
transCoFacs[c("CD27")] <- 4000
testFile <- fsList[[2]][[3]]
testFileTrans <- arcTrans(testFile, transNames = colnames(testFile)[5:16], 
                          transCoFacs = transCoFacs)
hist(exprs(testFileTrans)[,"CD27"], breaks = 200)

#And now, we normalize the peaks. 
normOrNot <- rep(TRUE, length(colnames(fsList[[1]])))
names(normOrNot) <- colnames(fsList[[1]])
normOrNot[c("SSC-A", "SSC-H",  "TIME", "CD3")] <- FALSE
fsListNorm <- lapply(fsList, function(x){
    peakNorm(fs = x, ctrlPos = 4, standFF = fsList[[2]][[4]], 
             transNames = colnames(x)[5:16], transCoFacs = transCoFacs, 
             normOrNot = normOrNot)
})

fsListTrans <- lapply(fsListNorm, arcTrans, 
                      transNames = colnames(fsListNorm[[1]])[5:16], 
                      transCoFacs = transCoFacs)

#And here, we view the standards from all batches. 
standFs <- flowSet(lapply(fsList, "[[", 4))
standFsTrans <- arcTrans(standFs, transNames = colnames(standFs)[5:16], 
                         transCoFacs = transCoFacs)

dir.create("Diagnostics/Cytometry/Norm", recursive = TRUE)

subDf <- flowSet2LongDf(standFsTrans)

for(i in seq_len(ncol(subDf)-3)){
    locName <- colnames(subDf)[i]
    locDf <- subDf[,c(i, which(colnames(subDf) == "acqDate"))]
    colnames(locDf)[1] <- "var"
    ggplot(locDf, aes(x=var, fill=acqDate)) +
        geom_density(alpha = 0.8) + xlab(locName) +ylab("")
    ggsave(paste0("Diagnostics/Cytometry/Norm/", locName, "_preNorm.pdf"))
}

standFs <- flowSet(lapply(fsListTrans, "[[", 4))

subDf <- flowSet2LongDf(standFs)

for(i in seq_len(ncol(subDf)-3)){
    locName <- colnames(subDf)[i]
    locDf <- subDf[,c(i, which(colnames(subDf) == "acqDate"))]
    colnames(locDf)[1] <- "var"
    ggplot(locDf, aes(x=var, fill=acqDate)) +
        geom_density(alpha = 0.8) + xlab(locName) +ylab("")
    ggsave(paste0("Diagnostics/Cytometry/Norm/", locName, "_afterNorm.pdf"))
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

csfIndexDf <- do.call("rbind", csfIndexMatList)

write.csv(csfIndexDf, "Data/Cytometry/flowDataPlusIndex.csv")

