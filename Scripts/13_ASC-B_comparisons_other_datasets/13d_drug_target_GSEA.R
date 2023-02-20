#Here, we are using a targeted gene-set enrichment-based approach from this publication: 
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03929-0
#to assess whether we enrich for any useful drug-related genes in our dataset.
#Here, bash and R code is mixed. All bash code has double hashtags.

##cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/Results
##git clone https://github.com/sxf296/drug_targeting

fileList <- list.files("Data/ASC_vs_B", pattern = "EdgeR_results.rds", 
                       recursive = TRUE, full.names = TRUE)

#Here, we remove the AE file. 
fileList <- fileList[-grep("/AE/", fileList)]

resList <- lapply(fileList, readRDS)

names(resList) <- gsub("Data/ASC_vs_B/|/EdgeR_results.rds", "", fileList)

#Now, we loop what we have previously done for individual datasets. 
lapply(names(resList), function(x){
    locDat <- resList[[which(names(resList) == x)]]
    locDatOrd <- locDat$table
    locDatOrd <- locDatOrd[order(locDatOrd$F, decreasing = TRUE),]
    locDatOrd$t <- sqrt(locDatOrd$F)
    write.csv(locDatOrd, paste0("Data/drug_targeting/ASC_vs_B_", x, ".csv"))
})

##cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/drug_targeting

##python3 dpGSEA.py -tt ASC_vs_B_Covid_PBMC.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_ASC_vs_B_Covid_PBMC.tsv
##python3 dpGSEA.py -tt ASC_vs_B_Healthy_PBMC.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_ASC_vs_B_Healthy_PBMC.tsv
##python3 dpGSEA.py -tt ASC_vs_B_Influensa_PBMC.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_ASC_vs_B_Influensa_PBMC.tsv
##python3 dpGSEA.py -tt ASC_vs_B_MS_Ramesh_CSF.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_ASC_vs_B_MS_Ramesh_CSF.tsv
##python3 dpGSEA.py -tt ASC_vs_B_MS_Ramesh_PBMC.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_ASC_vs_B_MS_Ramesh_PBMC.tsv
##python3 dpGSEA.py -tt ASC_vs_B_MS_Shafflick_CSF.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_ASC_vs_B_MS_Shafflick_CSF.tsv
##python3 dpGSEA.py -tt ASC_vs_B_MS_Shafflick_PBMC.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_ASC_vs_B_MS_Shafflick_PBMC.tsv

#We had to open this file in excel and re-save it as there were some formatting issues with the .tsv version
dtGseaResFileList <- list.files("Data/drug_targeting", full.names = TRUE, 
                                pattern = "L1000_50_results_ASC_vs_B_.+PBMC.csv|L1000_50_results_ASC_vs_B_.+CSF.csv")

allResOrdList <- lapply(dtGseaResFileList, function(x){
    locRes <- read.csv(x)
    locResOrd <- locRes[order(locRes$ES_p),]
    })

sigResList <- lapply(allResOrdList, function(x){
    xOrd <- x[order(x$ES_p),]
    xSig <- xOrd[which(xOrd$ES_p < 0.01 & xOrd$TCS_p < 0.01),]
})

saveRDS(sigResList, "Results/Drug_targeting/L1K_P50_significants_ASC_vs_B_other_datasets.rds")

#Now, what drugs, if any, are unique among the top ones we have identified?
aeSigRes <- read.csv("Results/drug_targeting/L1K_P50_significants_ASC_vs_B_AE.csv")

allNonAEHitsList <- lapply(sigResList, function(x) x$drug)
allNonAEHits <- unique(unlist(allNonAEHitsList))
aeSigRes$drug[which(gsub("|_.+", "", aeSigRes$drug)  %in% 
          gsub("|_.+", "", allNonAEHits))]
#"pazopanib_HELA"
#This is a multi-target tyrosine kinase inhibitor. All the others are unique, which
#is of course suspicious. 

#How much overlap do we see between the different conditions? 
Reduce(intersect, allNonAEHits)
#None are common to all. 

drugHitOverlap <- do.call("rbind", lapply(allNonAEHitsList, function(x){
    unlist(sapply(allNonAEHitsList, function(y){
        length(which(x %in% y))/length(x)
    }))
}))

rownames(drugHitOverlap) <- colnames(drugHitOverlap) <- gsub("Data/drug_targeting/L1000_50_results_ASC_vs_B_|.csv", "",
                                 dtGseaResFileList)
drugHitOverlap
#                  Covid_PBMC Healthy_PBMC Influensa_PBMC MS_Ramesh_CSF MS_Ramesh_PBMC MS_Shafflick_CSF MS_Shafflick_PBMC
#Covid_PBMC                 1   0.00000000      0.0000000     0.0000000      0.0000000                0        0.00000000
#Healthy_PBMC               0   1.00000000      0.2750000     0.0500000      0.0500000                0        0.12500000
#Influensa_PBMC             0   0.21568627      1.0000000     0.1176471      0.1764706                0        0.07843137
#MS_Ramesh_CSF              0   0.14285714      0.4285714     1.0000000      0.3571429                0        0.00000000
#MS_Ramesh_PBMC             0   0.04651163      0.2093023     0.1162791      1.0000000                0        0.00000000
#MS_Shafflick_CSF         NaN          NaN            NaN           NaN            NaN              NaN               NaN
#MS_Shafflick_PBMC          0   0.19230769      0.1538462     0.0000000      0.0000000                0        1.00000000

#This makes me doubt the veracity of the results: only a maximum of 42% of the drugs are shared anywhere, and in most
#cases this number is considerably lower. 

##################
#Association to known drugs
##################
#These drugs are either used or in trials for AE
drugList <- c("methylprednisolone", "dexamethasone", "varlilumab", "daratumumab",
              "plerixafor", "BL-8040", "LY2510924", "ulocuplumab", "CX-01",
              "rituximab", "ofatumumab", "ocrelizumab", "ublituximab", "obinutuzumab", 
              "ocaratuzumab", "veltuzumab", "tositumomab", "tafasitamab", 
              "loncastuximab", "inebiluzimab", "obexelimab",
              "indatuximab", "epratuzumab", "inotuzumab", 
              "belimumab", "ianalumab", "telitacicept", "atacicept",
              "NN 9828", "NNC0114-0005", "tocilizumab", "satriluzimab",
              "siltuximab", "ibrutinib", "tolebrutinib", "evobrutinib",
              "bortezomib", "carfilzomib", "natalizumab", 
              "dacetuzumab", "CDX-1140", "APX005M", "CP-870", "ADC-1013",
              "fingolimod", "siponimod", 
              "baricitinib", "tofacitinib", "upadacitinib")

#Now, which of these are present?
drugAndNothingElse <- sapply(allResOrdList[[1]]$drug, function(x){
    xClean <- gsub("|_.+", "", x)})

drugInData <- sapply(drugAndNothingElse, function(x){
    if(x %in% drugList){
        TRUE
    } else {
        FALSE
    }
})
table(drugInData)
#drugInData
#FALSE  TRUE 
#15415    76 
#So of course identical to the setup with our cells. 8 drugs are present here. 

#Now, for our cells, the list looked like this. So how many have hits that 
#are lower than 439, and which ones are they?

#carfilzomib.carfilzomib_HA1E                baricitinib.baricitinib_PC3 
#439                                        648 
#dexamethasone.dexamethasone_PC3                  ibrutinib.ibrutinib_HUES3 
#841                                       1290 
#tofacitinib.tofacitinib_MCF7                 bortezomib.bortezomib_A375 
#3965                                       4473 
#methylprednisolone.methylprednisolone_A375                 fingolimod.fingolimod_MCF7 
#5257                                       5354 

lapply(allResOrdList, function(x){
    drugAndNothingElse <- sapply(x$drug, function(y){
        xClean <- gsub("|_.+", "", y)})
    drugShortList <- unique(drugAndNothingElse[drugInData])
    #Here, we only include the top hit for each drig. 
    drugTopPlaces <- sapply(drugShortList, function(x){ which(drugAndNothingElse == x)[1]
    })
    print(drugTopPlaces[which(drugTopPlaces <= 439 & names(drugTopPlaces) %in% drugInData)])
    #We also check if any known drugs are significantly involved
    print(x[which(drugAndNothingElse %in% drugInData &
                      x$ES_p < 0.01 & x$TCS_p < 0.01)])
})
#This is completely negative. 





#In other words, none are enriched, which fits with the over-representation
#analysis, not indicating any of these drug targets as being engaged. 
#Therefore, we will not push this further. 












