library(DAseq)
library(DepecheR)
library(ggplot2)

csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_5_BCR.rds")

#We are now going to make an attempt to run a neighbor smoothing thing based on
#specific and non-specific cells. 
#As we have more positive cells than negative, and we want to treat both groups as
#similar as possible, we are going to run 100 bootstraps with 44 cells from the 
#negative (which is the number of negative cells) and 44 positive cells, but with 
#re-sampling, i.e. an individual cell can occur more than once. This is a more correct
#way of treating th data, than running the sampling without re-sampling, and thus
#making the negative subset always containing the same set of cells. 

specTab <- table(colData(csfSce)$donor, colData(csfSce)$Specific)
specTab
#FALSE Not_tested TRUE
#1166    20        110   74
#1227     4         42   11
#1284    20         43   23

#This shows that almost the whole imbalance comes from donor 1. This means
#that we will make the bootstraps on a per-donor basis. 

posRows <- seq_len(ncol(csfSce))[which(colData(csfSce)$Specific == "TRUE")]
negRows <- seq_len(ncol(csfSce))[which(colData(csfSce)$Specific == "FALSE")]

posRowsSplit <- split(posRows, f = colData(csfSce)$donor[which(colData(csfSce)$Specific == "TRUE")])
negRowsSplit <- split(negRows, f = colData(csfSce)$donor[which(colData(csfSce)$Specific == "FALSE")])

negLength <- specTab[,1]

#And here comes the definition of all the subsamplings
rowList <- list()
for(i in 1:100){
    set.seed(i+100)
    rowList[[i]] <- unlist(lapply(names(posRowsSplit), function(x){
        locVec <- c(sample(posRowsSplit[[x]], negLength[x], replace = TRUE),
                    sample(negRowsSplit[[x]], negLength[x], replace = TRUE))
        
    }))
}

#And now we check that all values are represented here. 
neighRows <- seq_len(ncol(csfSce))[-which(colData(csfSce)$Specific == "Not_tested")]
identical(length(unique(unlist(rowList))), length(neighRows))
#TRUE

#The distribution of the values
hist(table(unlist(rowList)), breaks = 50)
#Not unexpectedly, we see two groups, where the top one is likely to be dominated by
#negative values as well as non-JR1166 values, as these are the ones that we have
#the most of and are thus most likely to use only once. 

#We now create a dummy variable with 0 for false and 1 for true. 
specDummy <- rep(NA, ncol(csfSce))
specDummy[which(colData(csfSce)$Specific == "FALSE")] <- 0
specDummy[which(colData(csfSce)$Specific == "TRUE")] <- 1

#Now to the fun: all cells are of course included here. 
specSmoothMat <- do.call("rbind", lapply(seq_along(rowList), function(x){
    locNeigh <- rowList[[x]]
    neighSmooth(specDummy, reducedDim(csfSce, "GLMPCA"), locNeigh, 
                kNeighK = 9, kMeansK = 1)
}))

specSmooth <- colMeans(specSmoothMat)
dir.create("Diagnostics/Neighbors")
pdf("Diagnostics/Neighbors/SpecSmooth_histogram.pdf")
hist(specSmooth, breaks = 50)
dev.off()

ggDat <- data.frame("Specific" = specSmooth, "Cell_type" = csfSce$Cell_type,
                        reducedDim(csfSce, "GLMPCA_UMAP"))
p <- ggplot(ggDat, aes(x = X1, y = X2, color = Specific, shape = Cell_type, 
                       size = 3)) + 
    geom_point() + scale_color_viridis(option = "viridis") + theme_bw() +
    coord_equal()
p
ggsave("Diagnostics/Neighbors/UMAP_smooth_specificity.pdf")

#Interestingly, there are no transcriptomic traits that are specific for the
#specific cells, but there are a considerable fraction of cells only being surrounded
#by non-specific cells. This will be further investigated. First: are these
#a donor-specific anomaly, or are they true for all three? We run a similar
#repeated version here, to avoid just going for one k-nearest neighbor definition.
table(csfSce$donor)
#1166 1227 1284 
#204   57   86
#So we will now generate a similar set of balanced, re-sampled neighbor datasets
#with 57 cells from each donor. 
donRows <- split(seq_len(ncol(csfSce)), f = csfSce$donor)

neighDonRows <- lapply(1:100, function(x){
    unlist(lapply(donRows, sample, size = 57, replace = FALSE))
})

neighDons <- do.call("rbind", lapply(neighDonRows, function(x) {
    nUniqueNeighDons(csfSce$donor, reducedDim(csfSce, "GLMPCA"), neighRows = x,
                     kNeighK = 9,kMeansK = 1)
}))

averageNeighDons <- colMeans(neighDons)

ggDat <- cbind(ggDat, averageNeighDons)
p <- ggplot(ggDat, aes(x = X1, y = X2, color = averageNeighDons, shape = Cell_type, 
                       size = 3)) + 
    geom_point() + scale_color_viridis(option = "viridis") + theme_bw() +
    coord_equal()
p
ggsave("Diagnostics/Neighbors/UMAP_smooth_nDons.pdf")


#So now, cells surrounded by less than all donors are excluded. 
specSmoothRepresentative <- specSmooth[which(averageNeighDons == 3)]

pdf("Diagnostics/Neighbors/SpecSmooth_histogram_3_donors_around.pdf")
hist(specSmoothRepresentative, breaks = 50)
dev.off()

#Interestingly, this plot looks very much like the previous one




