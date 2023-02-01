

#Now, we identify all drug targets that are associated to our genes of interest, by
#simply subsetting the result list based on our genes of interest. 
#First, all drugs that have an association to our genes are identified: 

targetsAE <- gsub("Data/STRING_data/|_string_interactions_short.tsv", 
                  "", list.files("Data/STRING_data", recursive = TRUE, full.names = TRUE))


#The names that might be hard to get from the list are: 
#MS4A1 = CD20
#SDC1 = CD138
#TNFSF13B = BAFF
#TNFRSF13C = BAFFR
#ANP32B = APRIL
#TNFRSF13B = TACI
#TNFRSF1A = TNF receptor 1a
#TNFRSF1B = TNF receptor 1b
#CSF2 = GM-CSF
#CSF2RA = GM-CSF-RA
#CSF2RB = GM-CSF-RB

#Now, we are going to go through all the drugs in the result to identify the ones 
#that match at least one of these drugs.
drugMatch <- sapply(allResOrd$genes, function(x){
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

#So there are 20 that match. These are: 
allResOrd$drug[drugMatch]

#And which positions are these on? 
allResOrd$drugMatch <- drugMatch

ggplot(allResOrd, aes(x=ES, color=drugMatch, fill=drugMatch)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, 
                   position="identity") +theme_bw()
ggsave("Results/drug_targeting/Hist_of_selected_drugs_ES_ASC_vs_B_AE.pdf")

ggplot(allResOrd, aes(x=-log10(ES_p), color=drugMatch, fill=drugMatch)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, 
                   position="identity") +theme_bw()
ggsave("Results/drug_targeting/Hist_of_selected_drugs_ES-p_ASC_vs_B_AE.pdf")

#Now, which one is the highest?
allResOrd$drug[drugMatch][1]
#anecortave-acetate_HT29

#Ok, and what place does this drug have? 
which(drugMatch)[1]
#It is place 2597 and it is associated to CD38. 



