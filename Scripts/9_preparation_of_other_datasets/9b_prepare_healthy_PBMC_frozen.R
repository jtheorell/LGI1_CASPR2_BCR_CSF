library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
library(SingleR)
library(celldex)
#Here, we have downloaded the read count data and metadata from: 
#https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11452
#This dataset was generated for the study: 
#https://www.nature.com/articles/s41587-022-01311-4

colDataMat <- read.delim("../External/Data/Healthy/PBMCs.allruns.barcode_annotation.txt")
#colDataMat <- 
#read.delim(url("https://www.ebi.ac.uk/biostudies/files/E-MTAB-11452/PBMCs.allruns.barcode_annotation.txt"))
#This one below takes a lot of time. The file is 9 GB, as it is stored as a conventional matrix, and not 
#a sparse one. 
rawData <- read.delim("../External/Data/Healthy/PBMCs.allruns.readcounts_intronexon.txt")
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

saveRDS(pbmcSce, "../External/Data/Healthy/0_raw.rds")

#Now, we are going to make a very light exclusion of cells, and then go on to categorize
#them using singler, that should work well for plasmablasts. After that, we can decide
#If we want to make more through exclusions, etc. 

mito <- which(rowData(pbmcSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(pbmcSce, subsets=list(Mt=mito))
colData(pbmcSce) <- cbind(colData(pbmcSce), stats)

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
#17411  26260 

pbmcSce <- pbmcSce[,-which(pbmcSce$discard)]

saveRDS(pbmcSce, "../External/Data/Healthy/1_lowQual_excluded.rds")
#pbmcSce <- readRDS("../External/Data/Healthy/1_lowQual_excluded.rds")
#Now, after this very careful exclusion, we go on to transform and identify cell types. 


#NORMALIZATION
pbmcSce <- computeSumFactors(pbmcSce) 
pbmcSce <- logNormCounts(pbmcSce)

assayNames(pbmcSce)

#And now over to some singler analysis, to define the cell types. We will here use
#the monaco standard dataset to define B-cells and ASC. 

#Now, we will head straight for the singler analysis. 
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- pbmcSce[which(rowData(pbmcSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
table(pred$pruned.labels)
dir.create("Data/Comp_to_others/Healthy")
saveRDS(pred, "Data/Comp_to_others/Healthy/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/Healthy/MonacoSingler.rds")
table(pred$pruned.labels)
#Central memory CD8 T cells           Classical monocytes   Effector memory CD8 T cells 
#                       828                           341                           659 
#Exhausted B cells     Follicular helper T cells        Intermediate monocytes 
#             1185                           151                          3631 
#Low-density basophils       Low-density neutrophils                    MAIT cells 
#                6                            75                          1010 
#Myeloid dendritic cells                 Naive B cells             Naive CD4 T cells 
#              952                          2161                          1909 
#Naive CD8 T cells          Natural killer cells       Non classical monocytes 
#             1161                          8783                           515 
#Non-switched memory B cells            Non-Vd2 gd T cells                  Plasmablasts 
#                        555                          2173                            35 
#Plasmacytoid dendritic cells              Progenitor cells       Switched memory B cells 
#                         301                           241                           119 
#T regulatory cells Terminal effector CD4 T cells Terminal effector CD8 T cells 
#             1367                           473                          1232 
#Th1 cells                Th1/Th17 cells                    Th17 cells 
#     2316                          1037                           304 
#Th2 cells                Vd2 gd T cells 
#.     536                          3699 

pbmcSce$cellType <- pred$pruned.labels
healthySceB <- pbmcSce[,which(pbmcSce$cellType %in% c("Exhausted B cells",
                                                        "Naive B cells",
                                                        "Non-switched memory B cells",
                                                        "Switched memory B cells",
                                                        "Plasmablasts"))]
healthySceB$cellType[grep("B", healthySceB$cellType)] <- "B"
healthySceB$cellType[grep("Plasmablasts", healthySceB$cellType)] <- "ASC"
table(healthySceB$cellType)
#ASC    B 
#35 4020 

saveRDS(healthySceB, "../External/Data/Healthy/2_normalised_post_monaco.rds")


