library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
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

dir.create("Data/Comp_to_others/Healthy/SingleCellExpFiles")
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

#And now over to some singler analysis, to define the cell types. Here, we use the
#reference dataset created with the influensa data that like our dataset has surface
#protein information. 
library(SingleR)

ref <- readRDS("Data/Comp_to_others/Influensa/SingleCellExpFiles/4_cell_types_for_SingleR.rds")
singlerSce <- pbmcSce[which(rownames(pbmcSce) %in% row.names(ref)),]
refSce <- ref[which(row.names(ref) %in% rownames(pbmcSce)),]
refSceOrdered <- refSce[order(row.names(refSce)),]
singlerSceOrdered <- singlerSce[order(rownames(singlerSce)),]

identical(row.names(singlerSceOrdered), row.names(refSceOrdered))
#TRUE

#So now over to the fun! This turns into a bit of a marathon though, as it takes
#ages for this to run. 
pred <- SingleR(test=singlerSceOrdered, ref=refSceOrdered, labels=refSceOrdered$cellType)
table(pred$labels)
#ASC         B      CD4T      CD8T   Myeloid        NK       pDC Precursor 
#101      6536     11053     11074      3391      5585       157       117
dir.create("Diagnostics/Comp_to_others/Healthy")

plotScoreHeatmap(pred)
dev.copy(pdf,'Diagnostics/Comp_to_others/Healthy/SCE/Cell_type_prediction_SingleR.pdf', height = 5, width = 15)
dev.off()

pbmcSce$SingleR <- pred$labels

#Now, we exclude the rest of the cells
saveRDS(pbmcSce, "../External/Data/Healthy/2_norm_singler.rds")

#Now, we will head straight for the singler analysis. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- pbmcSce[which(rowData(pbmcSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
table(pred$pruned.labels)
saveRDS(pred, "Data/Comp_to_others/Healthy/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/Healthy/MonacoSingler.rds")
table(pred$pruned.labels, pbmcSce$SingleR)
#                               ASC    B CD4T CD8Naive CD8T Myeloid   NK  pDC Precursor
#Central memory CD8 T cells       2    3   25      755   42       0    0    0         1
#Classical monocytes              0   54   15       27   17     228    0    0         0
#Effector memory CD8 T cells      5   14   33      117  486       0    3    0         1
#Exhausted B cells               30 1084    7       41   16       7    0    0         0
#Follicular helper T cells        0    1   40      109    1       0    0    0         0
#Intermediate monocytes           0  583  125      385  436    2082   18    0         2
#Low-density basophils            0    1    0        0    1       3    0    0         1
#Low-density neutrophils          0    7   15        7   31      14    1    0         0
#MAIT cells                       0   39   98      361  510       0    1    0         1
#Myeloid dendritic cells          0  205   28      121   35     547    1    0        15
#Naive B cells                    0 1997   12       82   57      13    0    0         0
#Naive CD4 T cells                0    0   13     1895    0       0    0    0         1
#Naive CD8 T cells                0    0    0     1161    0       0    0    0         0
#Natural killer cells             2  672  129      838 3273     103 3744    0        22
#Non classical monocytes          0   16    1        3    5     490    0    0         0
#Non-switched memory B cells      0  555    0        0    0       0    0    0         0
#Non-Vd2 gd T cells               0  197   68      606 1152       3  144    0         3
#Plasmablasts                    35    0    0        0    0       0    0    0         0
#Plasmacytoid dendritic cells     0   69    2       56   15       1    0  157         1
#Progenitor cells                 0   79    6       51   35      16    1    0        53
#Switched memory B cells          4   98    2       10    4       1    0    0         0
#T regulatory cells              14   12  262     1067   10       0    1    0         1
#Terminal effector CD4 T cells    0   49   23      179  221       1    0    0         0
#Terminal effector CD8 T cells    2   32   10      161 1023       0    3    0         1
#Th1 cells                        5   16  356     1888   49       0    0    0         2
#Th1/Th17 cells                   0    1  424      576   36       0    0    0         0
#Th17 cells                       0    2  164      134    4       0    0    0         0
#Th2 cells                        3    3  176      352    1       0    0    0         1
#Vd2 gd T cells                   6  180  247     1371 1875       4   10    0         6

#We will also export the plasmablasts here. 
healthyScePB <- pbmcSce[,which(pbmcSce$SingleR == "ASC" & pred$pruned.labels == "Plasmablasts")]
saveRDS(healthyScePB, "Data/Comp_to_others/Healthy/SingleCellExpFiles/3_PB.rds")


