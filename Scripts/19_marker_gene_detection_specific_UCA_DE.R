#Here, we are doing a basic marker gene detection analysis of the cells with
#known specificities. We will not take the donor information into consideration
#as the donors and the cells are so few. Instead, we will do the downsteam
#analyses based on the idea that the specific and non-specific cells represent
#discernably different "cell types"

#We are, as always, following the OSCA project. This time: 
#http://bioconductor.org/books/3.15/OSCA.basic/marker-detection.html
library(scran)
library(scater)
library(BiocParallel)
library(edgeR)
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

#To increase interpretability, we remove all genes that lack a hgnc symbol
csfSce <- csfSce[-which(is.na(rowData(csfSce)$hgnc_symbol) | 
                            rowData(csfSce)$hgnc_symbol == ""),]

#A few entries have two rows, i.e. two gene products associated to the same
#hgnc symbol: 
dupDat <- rowData(csfSce)$hgnc_symbol[which(duplicated(rowData(csfSce)$hgnc_symbol))]

#For these, we only include the one with the highest count. 
for(i in dupDat){
    locDat <- csfSce[which(rowData(csfSce)$hgnc_symbol == i),]
    locSums <- rowSums(counts(locDat))
    csfSce <- csfSce[-which(rowData(csfSce)$hgnc_symbol == i)[-which.max(locSums)],]
}

#Now, we have unique hgnc symbols for all, so we now change the row names
row.names(csfSce) <- rowData(csfSce)$hgnc_symbol
#We now remove all cells without known specificity.
specSce <- csfSce[,-which(csfSce$Specific == "Not_tested")]

table(specSce$Specific_UCA, specSce$donor)
#      1166 1227 1284
#FALSE   70    5   31
#TRUE    19    8    4

#This shows that we have at least four cells in each group, which is the level
#that is defined as the lowest for hte other experiment. 

#Now, we set this up as a DE experiment using EdgeR. 

#NOw, we do the pseudo-bulk and DE analysis
nameVec <- paste0(specSce$donor, specSce$Specific_UCA)
sumMat <- do.call("cbind", lapply(unique(nameVec), function(x) 
    rowSums(counts(specSce[,which(nameVec == x)]))))


# Creating up a DGEList object for use in edgeR:
y <- DGEList(sumMat, samples=data.frame("donor" = substr(unique(nameVec), 1,4),
                                        "Specific_UCA" = substr(unique(nameVec), 5,9)))
y    
    
#Now, we controversially only exclude donors with less than 4 cells

keep <- filterByExpr(y, group=y$samples$Specific_UCA)
y <- y[keep,]
#summary(keep)
#Mode   FALSE    TRUE 
#logical   10130    9097 

y <- calcNormFactors(y)
design <- model.matrix(~factor(y$samples$Specific_UCA), y$donor)
#design
y <- estimateDisp(y, design)
#summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
#summary(fit$var.prior)
res <- glmQLFTest(fit, coef=ncol(design))
#summary(decideTests(res))
#factor(y$samples$Specific_UCA)TRUE
#Down                                    0
#NotSig                               9097
#Up                                      0

#So none are significant. This should suffice. 
