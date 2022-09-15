#Now, we create the elastic-net regularized principal component analysis data needed 
#for the downstream analyses. 
library(scry)

#We start bby selecting genes and making a townes pca. 
#We follow: 
#https://bioconductor.org/packages/release/bioc/vignettes/scry/inst/doc/scry.html
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm_and_flow_integration.rds")

csfSceTownes <- devianceFeatureSelection(csfSce, assay="counts", sorted=TRUE)

dir.create("Diagnostics/GLMPCA")

pdf("Diagnostics/GLMPCA/Binomial_deviance_selection.pdf")
plot(rowData(csfSceTownes)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")
dev.off()

#We select these 2000 genes. 
selectSCE <- csfSceTownes[1:2000,]

Sys.time()
set.seed(101)
glmpcaSCE <- GLMPCA(selectSCE, L=100, assay="counts")
Sys.time()

#And this is added back: 
reducedDim(csfSce, "GLMPCA") <- reducedDim(glmpcaSCE)

saveRDS(csfSce, "Data/SingleCellExpFiles/csfSce_3_GLMPCA.rds")

