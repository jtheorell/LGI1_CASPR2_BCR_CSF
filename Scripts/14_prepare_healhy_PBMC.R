library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
#Here, we have downloaded the read count data and metadata from: 
#https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11452
#This dataset was generated for the study: 
#https://www.nature.com/articles/s41587-022-01311-4

colDataMat <- read.delim("Data/SmartSeq3_PBMC/Raw/PBMCs.allruns.barcode_annotation.txt")
#colDataMat <- 
#read.delim(url("https://www.ebi.ac.uk/biostudies/files/E-MTAB-11452/PBMCs.allruns.barcode_annotation.txt"))
#This one below takes a lot of time. The file is 9 GB, as it is stored as a conventional matrix, and not 
#a sparse one. 
rawData <- read.delim("Data/SmartSeq3_PBMC/Raw/PBMCs.allruns.readcounts_intronexon.txt")
#rawData <- 
#read.delim(url(https://www.ebi.ac.uk/biostudies/files/E-MTAB-11452/PBMCs.allruns.readcounts_intronexon.txt"))
rawDataSparse <- Matrix(as.matrix(rawData), sparse = TRUE)

pbmcSce <- SingleCellExperiment(assays = list("counts" = rawDataSparse), colData = colDataMat)


#And now, we do a conventional analysis on this to see if we find any ASC. 
pbmcSce  <- getBMFeatureAnnos(pbmcSce, filters = "ensembl_gene_id", 
                              attributes = c("ensembl_gene_id", "hgnc_symbol",
                                             "chromosome_name", "gene_biotype", 
                                             "start_position", "end_position"), 
                              dataset = "hsapiens_gene_ensembl")

dir.create("Data/SmartSeq3_PBMC/SingleCellExpFiles")
saveRDS(pbmcSce, "Data/SmartSeq3_PBMC/SingleCellExpFiles/0_raw.rds")

#Now, we are going to make a very light exclusion of cells, and then go on to categorize
#them using singler, that should work well for plasmablasts. After that, we can decide
#If we want to make more through exclusions, etc. 

mito <- which(rowData(pbmcSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(pbmcSce, subsets=list(Mt=mito))
colData(pbmcSce) <- cbind(colData(pbmcSce), stats)
pbmcSce$donor <- gsub("|_.+", "", colnames(pbmcSce))

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

pbmcSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))
#low_lib_size          low_n_features high_subsets_Mt_percent                 discard 
#134                    4697                    1221                    5657

#So this means loosing about 10% of the cells, which is considerably less than the 40% lost
#in the original QC: 
table(colData(pbmcSce)$QC_status)
#QCfail QCpass 
#9532  15887

pbmcSce <- pbmcSce[,-which(pbmcSce$discard)]

saveRDS(pbmcSce, "Data/SmartSeq3_PBMC/SingleCellExpFiles/1_lowQual_excluded.rds")
#pbmcSce <- readRDS("Data/SmartSeq3_PBMC/SingleCellExpFiles/1_lowQual_excluded.rds")
#Now, after this very careful exclusion, we go on to transform and identify cell types. 


#NORMALIZATION
pbmcSce <- computeSumFactors(pbmcSce) 
pbmcSce <- logNormCounts(pbmcSce)

assayNames(pbmcSce)

#And now over to some singler analysis, to define the cell types. Again, I would
#not feel confortable doing this if I did not believe the ASC to be so different from
#all other cell types. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- pbmcSce[which(rowData(pbmcSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
table(pred$labels)
#THis gives us 38 plasmablasts. This will be one of our reference datasets.
plotScoreHeatmap(pred)
dev.copy(pdf,'Diagnostics/SmartSeq3_PBMC/SCE/Cell_type_prediction_fine_SingleR.pdf', height = 5, width = 15)
dev.off()
pbmcSce$SingleR <- pred$labels
#Now, we exclude the rest of the cells
saveRDS(pbmcSce, "Data/SmartSeq3_PBMC/SingleCellExpFiles/2_norm_singler.rds")
ss3PbSce <- pbmcSce[,which(pbmcSce$SingleR == "Plasmablasts")]
saveRDS(ss3PbSce, "Data/SmartSeq3_PBMC/SingleCellExpFiles/3_plasmablasts.rds")

#Now, we will go on to a combined analuysis, where our data is compared both
#to this dataset and to the dengue plasmacell dataset, to identify more likely true
#hits. 


