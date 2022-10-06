#Here, the idea is to start by creating a trajectory, from the most
#to the least likely specific cell, and then use this newly created dimension to
#run a partial least squares analysis to identify targets. 
library(DepecheR)
library(mixOmics)
library(ggplot2)

csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

#What we will now do, is to identify the number of specific and non-specific
#cells, and then run bootstrapped neighrbor smoothings for the full dataset
#with replacement, to even out the difference in likelihood to have neighbors from
#the specific and non-specific cells

table(csfSce$Specific)

posRows <- seq_len(ncol(csfSce))[which(colData(csfSce)$Specific == "TRUE")]
negRows <- seq_len(ncol(csfSce))[which(colData(csfSce)$Specific == "FALSE")]

neighRows <- lapply(1:100, function(x){
    locPos <- sample(posRows, length(negRows), replace = TRUE)
    locNeg <- sample(negRows, length(negRows), replace = TRUE)
    c(locPos, locNeg)
})

identical(length(c(posRows, negRows)), length(unique(unlist(neighRows))))
#TRUE

#The distribution of the values
hist(table(unlist(neighRows)), breaks = 50)
#Not unexpectedly, we see two groups, where the top one is likely to be the
#negative values, and the lower are the positive values. The relation 3:1
#fits with the number of cells. 

#We now create a dummy variable with 0 for false and 1 for true. 
specDummy <- rep(NA, ncol(csfSce))
specDummy[which(colData(csfSce)$Specific == "FALSE")] <- 0
specDummy[which(colData(csfSce)$Specific == "TRUE")] <- 1

#Now to the fun: all cells are of course included here. 
specSmoothMat <- do.call("rbind", lapply(seq_along(neighRows), function(x){
    locNeigh <- neighRows[[x]]
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

#Showing some uneven distribution over the UMAP. 

#We also check how well it works, i.e. how many of the negatives that
#are towards the negative side, etc
hist(specSmooth[which(colData(csfSce)$Specific == "FALSE")])
hist(specSmooth[which(colData(csfSce)$Specific == "TRUE")])

#The negatives are almost equally distributed here, whereas the positives are
#localized to an upper area. 
#We start by excluding all variables without any positive points. 
plsData <- t(logcounts(csfSce))

plsData <- plsData[,-which(colSums(plsData) < 5)]
plsResults <- pls(plsData, specSmooth, ncomp = 1)

xLoadings <- plsResults$loadings$X[,1]
xLoadingsOrdered <- xLoadings[order(abs(xLoadings), decreasing = TRUE)]


plsTuning <- tune.spls(plsData, specSmooth, test.keepX = 2^(1:10),ncomp = 1, folds = 10, measure = "MSE")
plsResults <- spls(plsData, specSmooth, keepX = plsTuning$choice.keepX, 
                   ncomp = 1)

#Here is the sPLS vector. 
splsVec <- plsResults$variates$X

dir.create("Results/Partial_least_squares")
pdf("Results/Partial_least_squares/PLS_vs_smooth_specificity.pdf")
plot(splsVec, specSmooth)
dev.off()
cor(splsVec, specSmooth)
#0.7348311
#This is actually the best correlation we can get. I have tried to include the
#full 14000 genes, or the top 100, but this model supersedes them. 

#And here are the loadings: 
splsGenes <- plsResults$loadings$X
splsGenes <- as.data.frame(splsGenes[-which(splsGenes == 0),])
colnames(splsGenes)[1] <- "spls_loadings"
splsGenes$names <- sapply(row.names(splsGenes), 
                          function(x) 
                              rowData(csfSce)$hgnc_symbol[which(row.names(csfSce) == x)])

splsGenes <- splsGenes[order(splsGenes$spls_loadings),]
splsGenes
#                spls_loadings names
#ENSG00000180879   -0.70963260  SSR4
#ENSG00000164687    0.04798408 FABP5
#ENSG00000128708    0.12855863  HAT1
#ENSG00000171848    0.20415405  RRM2
#ENSG00000151725    0.29023220 CENPU
#ENSG00000176890    0.29599138  TYMS
#ENSG00000276345    0.34624430      
#ENSG00000166451    0.37971326 CENPN

#As number 7 lacks a name, we will substitute the entreze code there. 


splsSce <- csfSce[which(row.names(csfSce) %in% row.names(splsGenes)),]
splsGenesOrdered <- splsGenes[match(row.names(splsSce), row.names(splsGenes)),]

identical(row.names(splsGenesOrdered), 
          row.names(splsSce))

row.names(splsSce) <- splsGenesOrdered$names
plotExpression(splsSce, features=row.names(splsSce), 
               x="Specific", colour_by="Specific")
dev.copy(pdf,'Results/Partial_least_squares/sPLS_selected_genes.pdf', height = 8, width = 6)
dev.off()
