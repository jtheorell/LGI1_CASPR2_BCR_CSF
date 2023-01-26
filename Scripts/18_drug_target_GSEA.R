#Here, we are using a targeted gene-set enrichment-based approach from this publication: 
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03929-0
#to assess whether we enrich for any useful drug-related genes in our dataset.
#Here, bash and R code is mixed. All bash code has double hashtags.

##cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/Results
##git clone https://github.com/sxf296/drug_targeting

integrDat <- readRDS("Results/Comp_to_others/All_DE_integrated_downsamp.rds")$Median

#First, we rank the values based on the F statistic
integrDatOrd <- integrDat[order(integrDat$F_stat, decreasing = TRUE),]

#Then we create a T statistic based on the F statistic. If I have correctly
#understood this, then the F ant T stastistics are related as follows: 
#F = T^2. Followingly, sqrt(F) = T.
integrDatOrd$t <- sqrt(integrDatOrd$F_stat)

#So now, we will export this, and use it for the drug screening work. 
write.csv(integrDatOrd, "Results/drug_targeting/AE_vs_all.csv")

##cd ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/Results/drug_targeting

##python3 dpGSEA.py -tt AE_vs_all.csv -dr pms/L1K_P50.csv -i 1000 -o Super_deep_results.tsv
#We had to open this file in excel and re-save it as there were some formatting issues with the .tsv version
allResSuperDeep <- read.csv("Results/drug_targeting/Super_deep_results.csv")
    
allSupResOrd <- allResSuperDeep[order(allResSuperDeep$ES_p),]

allSupResOrdLow <- allSupResOrd[which(allSupResOrd$ES_p < 0.01 & 
                                          allSupResOrd$TCS_p < 0.01),]

write.csv(allSupResOrdLow, "Results/drug_targeting/L1K_P50_significants.csv", row.names = FALSE)

#Now, we identify all drug targets that are associated to our genes of interest, by
#simply subsetting the result list based on our genes of interest. 
#First, all drugs that have an association to our genes are identified: 

targetsAE <- c("CD27", "CD38", "CXCR4", "CXCR5", "MS4A1", "CD19", "SDC1", "CD22", "TNFSF13B",
               "TNFRSF13C", "ANP32B", "TNFRSF13B", "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IFNG",
               "IFNGR1","IFNGR2", "TNF", "TNFRSF1A", "TNFRSF1B", "IL10", "IL10RA","IL10RB", 
               "CSF2", "CSF2RA", "CSF2RB", "IL6", "IL6R", "BTK", "ITGA4", "ITGB1", "CD40")

#Now, we are going to go through all the drugs in the result to identify the ones 
#that match at least one of these drugs.
drugMatch <- sapply(allSupResOrd$genes, function(x){
    xClean <- gsub("\\[|\\]| |'", "", x)
    if(grepl(",", xClean)){
        xClean <- strsplit(xClean, ",")
    }
    if(any(xClean %in% targetsAE)){
        TRUE
    } else {
        FALSE
    }
})

#So there are 38 that match. These are: 
allSupResOrd$drug[drugMatch]


#And which positions are these on? 
allSupResOrd$drugMatch <- drugMatch

ggplot(allSupResOrd, aes(x=ES, color=drugMatch, fill=drugMatch)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, 
                   position="identity") +theme_bw()
ggsave("Results/drug_targeting/Hist_of_selected_drugs_ES.pdf")

ggplot(allSupResOrd, aes(x=-log10(ES_p), color=drugMatch, fill=drugMatch)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, 
                   position="identity") +theme_bw()
ggsave("Results/drug_targeting/Hist_of_selected_drugs_ES-p.pdf")

#Now, which one is the highest?
allSupResOrd$drug[drugMatch][1]
#perampanel_HA1E

#Ok, and what place does this drug have? 
which(drugMatch)[1]
#It is place 759, and it is associated to CD138. 


#Now, we will try to create another list based on the genes instead
#CD27: varlilumab (agonist)
#CD38: Daratumumab
#CXCR4: Plerixafor, BL-8040, LY2510924, ulocuplumab, CX-01 (dociparstat sodium)
#CXCR5: can only find CAR-T cells
#MS4A1: rituximab, ofatumumab, ocrelizumab, Ublituximab, Obinutuzumab, Ocaratuzumab, Veltuzumab, Tositumomab
#CD19: tafasitamab, loncastuximab tesirine, Inebiluzimab, Obexelimab
#SDC1: Indatuximab ravtansine
#CD22: epratuzumab, inotuzumab ozogamicin
#BAFF/BAFFR: Belimumab, Ianalumab, Telitacicept, Atacicept (not efficient in MS) (2 last are TACI-Ig-fusion proteins), 
#APRIL/TACI: See BAFF above
#IL-21/IL-21R: One called Ab-01, but seems far from clinical
#IL-6R: tocilizumab, sarilumab
#IL6: siltuximab
#BTK: ibrutinib, Tolebrutinib, Evobrutinib
#proteosome pathways: bortezomib
#a1B4 integrin: natalizumab
#CD40 (agonists): dacetuzumab, CDX-1140, APX005M, CP-870,893, ADC-1013
#SIPR: 
    