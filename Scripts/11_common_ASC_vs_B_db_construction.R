#Downstream, we will focus on especially our cells and how they compare to the other
#datasets. To do this, we need to make sure that our data is comparable to the others.
#Therefore, we are here making sure that we are looking at the same genes. 

aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
aeSce <- aeSce[,-which(aeSce$cellType == "ASC" & aeSce$Specific != "TRUE")]
aeSce <- aeSce[,-which(aeSce$cellType == "B" & aeSce$Specific == "TRUE")]

#Before going on, we will change the row names to HGNC symbols and identify a subset
#that is common to all datasets, for comparability reasons. 
covidSce <- readRDS("../External/Data/Covid/2_normalised.rds")
healthySce <- readRDS("../External/Data/Healthy/3_post_Kotliarov.rds")
influensaSce <- readRDS("../External/Data/Influensa/4_post_monaco.rds")
msRameshCsfSce <- readRDS("../External/Data/MS_Ramesh_CSF/3_post_Kotliarov.rds")
msRameshPbmcSce <- readRDS("../External/Data/MS_Ramesh_PBMC/3_post_Kotliarov.rds")
msShafCsfSce <- readRDS("../External/Data/MS_Schafflick_CSF/3_post_Kotliarov.rds")
msShafPbmcSce <- readRDS("../External/Data/MS_Schafflick_PBMC/3_post_Kotliarov.rds")

sceList <- list("AE_CSF" = aeSce,
                "Covid_PBMC" = covidSce, 
                "Healthy_PBMC" = healthySce, 
                "Influensa_PBMC" = influensaSce, 
                "MS_Ramesh_CSF" = msRameshCsfSce, 
                "MS_Ramesh_PBMC" = msRameshPbmcSce, 
                "MS_Shafflick_CSF" = msShafCsfSce,
                "MS_Shafflick_PBMC" = msShafPbmcSce)
rowDatList <- lapply(sceList, rowData)

commonEnsembl <- Reduce(intersect, lapply(rowDatList, function(x) x$ensembl_gene_id))
length(commonEnsembl)
#16543
commonHgnc <- Reduce(intersect, lapply(rowDatList, function(x) x$hgnc_symbol))
length(commonHgnc)
#17544

#Interestingly, the number of shared transcripts is considerably larger with the 
#hgnc symbols. It is however also so, that the hgnc symbols are not necessarily
#unique, but can map to multiple ensembl ids. If we only look at the 
#ensembl ids, however, we lose genes like TNF and IFNGR2, which we would like to include for 
#downstream analyses. We will now check how many genes we would lose if we exlude
#all tuplicated hgnc symbols. 
duplicatedHgnc <- unique(unlist(lapply(rowDatList, function(x) 
    x$hgnc_symbol[which(duplicated(x$hgnc_symbol))])))

length(duplicatedHgnc)
#20
#So that is a clear win. How many of these are shared?
ducplicatedSharedHgnc <- duplicatedHgnc[which(duplicatedHgnc %in% commonHgnc)]
length(ducplicatedSharedHgnc)
#11
#These are: 
#"PINX1"   "ATXN7"   "ABCF2"   "POLR2J3" "GGT1"    "AHRR"    NA        "GOLGA8M"
#""SIGLEC5" "PRPF31"  "TMC4"  

#These will have to be removed, sadly. 

#We will now start by adding the unique hgncs to each file. 
commonHgncUnique <- commonHgnc[-which(commonHgnc %in% duplicatedHgnc)]
length(commonHgncUnique)
#17533


#Now, we reduce all datasets to these genes, change the row names
#and order them. 

sceListRed <- lapply(sceList, function(x){
    xTemp <- x[which(rowData(x)$hgnc_symbol %in% commonHgncUnique),]
    row.names(xTemp) <- rowData(xTemp)$hgnc_symbol
    xTemp <- xTemp[order(row.names(xTemp)),]
})

#Are they now all identical?
all(sapply(sceListRed, function(x) identical(row.names(x), row.names(sceListRed[[1]]))))
#TRUE


saveRDS(sceListRed, "Data/Comp_to_others/All_sce_common_genes_pre_all.rds")