#Here, we are using a targeted gene-set enrichment-based approach from this publication: 
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03929-0
#to assess whether we enrich for any useful drug-related genes in our dataset.
#Here, bash and R code is mixed. All bash code has double hashtags.

##cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/Results
##git clone https://github.com/sxf296/drug_targeting

res <- readRDS("Results/Specific_vs_not/EdgeR_result.rds")

#First, we rank the values based on the F statistic
resDatOrd <- res$table
resDatOrd <- resDatOrd[order(resDatOrd$F, decreasing = TRUE),]

#Then we create a T statistic based on the F statistic. If I have correctly
#understood this, then the F ant T stastistics are related as follows: 
#F = T^2. Followingly, sqrt(F) = T.
resDatOrd$t <- sqrt(resDatOrd$F)

#So now, we will export this, and use it for the drug screening work. 
write.csv(resDatOrd, "Data/drug_targeting/Spec_vs_non.csv")

##cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Data/drug_targeting

##python3 dpGSEA.py -tt Spec_vs_non.csv -dr pms/L1K_P50.csv -i 1000 -o L1000_50_results_spec_vs_non.tsv
#We had to open this file in excel and re-save it as there were some formatting issues with the .tsv version

allRes <- read.csv("Data/drug_targeting/L1000_50_results_spec_vs_non.csv")
    
allResOrd <- allRes[order(allRes$ES_p),]

allResOrdLow <- allResOrd[which(allResOrd$ES_p < 0.01 & 
                                          allResOrd$TCS_p < 0.01),]

dir.create("Results/Drug_targeting")
write.csv(allResOrdLow, "Results/drug_targeting/L1K_P50_significants_spec_vs_non.csv", row.names = FALSE)

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
drugAndNothingElse <- sapply(allResOrd$drug, function(x){
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
#Now, how many of these are unique, in other words, how many drugs are really present here?

table(drugAndNothingElse[drugInData])
#baricitinib         bortezomib        carfilzomib      dexamethasone         fingolimod 
#7                  7                  7                  9                  7 
#ibrutinib methylprednisolone        tofacitinib 
#24                  7                  8

#So 8 uniwue drugs.

drugShortList <- unique(drugAndNothingElse[drugInData])
#Now, how do they rank?
sapply(drugShortList, function(x){ which(drugAndNothingElse == x)[1]
})
#carfilzomib.carfilzomib_HA1E                tofacitinib.tofacitinib_NPC 
#336                                        535 
#ibrutinib.ibrutinib_CD34 methylprednisolone.methylprednisolone_YAPC 
#827                                       1115 
#baricitinib.baricitinib_PC3                 fingolimod.fingolimod_HA1E 
#2327                                       2353 
#dexamethasone.dexamethasone_A549                 bortezomib.bortezomib_MCF7 
#3565                                       5221

#In other words, none are enriched, which fits with the over-representation
#analysis, not indicating any of these drug targets as being engaged. 
#Therefore, we will not push this further. 












