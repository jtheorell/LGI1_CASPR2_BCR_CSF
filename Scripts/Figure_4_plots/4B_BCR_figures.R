library(ggnetwork)
library(viridis)
library(scatterpie)
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

BCR_extras <- read.csv("Data/BCR_auxiliaries/CDR3_length_etc.csv")
#Before doing anything else, we will exclude all cells with an unknown specificity, 
#as they are not useful in this analysis. 
BCR_spec <- BCR_all[-which(BCR_all$Specific == "Not_tested"),]

#Now, we will generate a V-gene family thing. 
BCR_spec$V_gene_fam <- gsub("|-.+", "", BCR_spec$V_CALL)

#And now before going to the heavy chains, we will add the light chain
#V-family usage to all chains. 
BCR_spec$V_gene_fam_L <- sapply(BCR_spec$CELL, function(x){
    BCR_spec$V_gene_fam[which(BCR_spec$LOCUS != "H" &
                                  BCR_spec$CELL == x)][1]
})
BCR_spec$V_gene_fam_L <- gsub("D", "", BCR_spec$V_gene_fam_L)
#For the visualisations, we do not need the light chains, as
#the heavy chains are already categorised according to light chain usage, and the
#light chain distances, that will be used for collapsing BCRs, are also present 
#for the heavy chains. 
#THis step is unnecessary here, as the light chains are not IgG. 
BCR_H <- BCR_spec[which(BCR_spec$LOCUS == "H"),]

#Here, we integrate the extra columns. 
BCR_extra_H <- BCR_extras[which(BCR_extras$SEQUENCE_ID %in% BCR_H$SEQUENCE_ID),]

BCR_extraOrd_H <- BCR_extra_H[match(BCR_H$SEQUENCE_ID,
                                    BCR_extra_H$SEQUENCE_ID),]

identical(BCR_H$SEQUENCE_ID, BCR_extraOrd_H$SEQUENCE_ID) #TRUE

BCR_H$CDR3_length <- BCR_extraOrd_H$CDR3_Length
BCR_H$Charge <- BCR_extraOrd_H$CDR3_Net_Charge_at_pH_7.0
BCR_H$Isoel <- BCR_extraOrd_H$CDR3_Isoelectric_Point

#Now, we are going to start with anyway taking in information from the light chains: 
#we will add a number of mutational scores for each cell. 
BCR_H$All_mutations_H_L <- sapply(BCR_H$CELL, function(x){
    sum(BCR_spec$All_mutations[which(BCR_spec$CELL == x)])
})

BCR_H$Non_silent_mutations_H_L <- sapply(BCR_H$CELL, function(x){
    sum(BCR_spec$Non_silent_mutations[which(BCR_spec$CELL == x)])
})

BCR_H$CDR_mutations_H_L <- sapply(BCR_H$CELL, function(x){
    sum(BCR_spec$CDR_mutations[which(BCR_spec$CELL == x)])
})

BCR_H$Non_silent_CDR_mutations_H_L <- sapply(BCR_H$CELL, function(x){
    sum(BCR_spec$Non_silent_CDR_mutations[which(BCR_spec$CELL == x)])
})

#Here, we are going to add a dummy donor will all the info from all donors, 
#that will be displayed in the paper. 
BCR_H2 <- BCR_H
BCR_H2$Sample <- "All"
BCR_H_double <- rbind(BCR_H, BCR_H2)
#Now over to the plotting. 

#Now, we will create simple stacked bar graphs for the above parameters. 
barGraphDat <- BCR_H_double[,c("Sample","Clonal", "light_type", "ISOTYPE", "CDR3_length",
                        "Specific", "Specific_antigen", "Isoel", "Charge", "V_gene_fam", 
                        "V_gene_fam_L",
                        "All_mutations_H_L",
                        "Non_silent_mutations_H_L",
                        "CDR_mutations_H_L",
                        "Non_silent_CDR_mutations_H_L",
                        "All_mutations",
                        "Non_silent_mutations",
                        "CDR_mutations",
                        "Non_silent_CDR_mutations")]
barGraphDat$light_type <- factor(barGraphDat$light_type, levels = c("K", "L"))
barGraphDat$ISOTYPE <- factor(barGraphDat$ISOTYPE, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4"))
barGraphDat$V_gene_fam <- factor(barGraphDat$V_gene_fam, levels = c("IGHV1", "IGHV2", "IGHV3", "IGHV4", "IGHV5"))
barGraphDat$V_gene_fam_L <- factor(barGraphDat$V_gene_fam_L, levels = 
                                       c("IGKV1","IGKV2","IGKV3","IGKV4",
                                         "IGLV1","IGLV2","IGLV3","IGLV6"))
barGraphDat$Specific_antigen <- factor(barGraphDat$Specific_antigen, levels = c("FALSE", "LGI1", "CASPR2"))
colsAndNames <- list(list("V_gene_fam", gray.colors(12, start =1, end = 0)[3:8]),
                     list("V_gene_fam_L", gray.colors(12, start =1, end = 0)[4:12]),
                     list("light_type", c("#D40041", "#2E00B7")),
                     list("ISOTYPE", c("#B1CFFF", "#9E4532", "#541352", "#238A8D"))
                     )
dir.create("Results/Figure_4_plots/4B_CloneBarGraphs", recursive = TRUE)
for(i in colsAndNames){
    locDat <- as.data.frame(table(barGraphDat[,i[[1]]], barGraphDat$Sample, barGraphDat$Specific_antigen))
    colnames(locDat) <- c(i[[1]], "Sample","Specific_antigen", "value")
    p <- ggplot(locDat,                         
           aes(x = Specific_antigen,
               y = value,
               fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        facet_grid(~ Sample) + scale_fill_manual(values=i[[2]], drop = FALSE)
    p
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 9, height = 6)
    
}

#And now, we are doing the same but for the ones that need a continuous scale


for(i in colnames(barGraphDat)[grep("mutations|Charge|Isoel|CDR3", colnames(barGraphDat))]){
    locDat <- as.data.frame(table(barGraphDat[,i], barGraphDat$Sample, barGraphDat$Specific_antigen))
    locDat[,1] <- as.numeric(as.character(locDat[,1]))
    colnames(locDat) <- c(i, "Sample","Specific_antigen", "value")
    if(grepl("mutations", i)){
        colorOption <- "G"
        scaleLims <- c(0, max(locDat[,1]))
    } else if(grepl("Charge", i)){
        colorOption <- "E"
        scaleLims <- c(min(locDat[,1]), max(locDat[,1]))
    } else if(grepl("Isoel", i)){
        colorOption <- "F"
        scaleLims <- c(min(locDat[,1]), max(locDat[,1]))
    } else if(grepl("CDR3", i)){
        colorOption <- "A"
        scaleLims <- c(min(locDat[,1]), max(locDat[,1]))
    }
    
    
    p <- ggplot(locDat,                         
                aes(x = Specific_antigen,
                    y = value,
                    fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        facet_grid(~ Sample) + scale_fill_viridis(limits = scaleLims, 
                                                  discrete = FALSE, 
                                                  option = colorOption)
    p
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 9, height = 6)
}

#Now, as I am lazy and we want plots without the third, empty category for the
#supplementary, I will create a simplified version here of all the plots. 
dir.create("Results/Figure_4_plots/4B_CloneBarGraphs_for_supplement")

for(i in colsAndNames){
    locDat <- as.data.frame(table(barGraphDat[,i[[1]]], barGraphDat$Sample, barGraphDat$Specific))
    colnames(locDat) <- c(i[[1]], "Sample","Specific", "value")
    p <- ggplot(locDat,                         
                aes(x = Specific,
                    y = value,
                    fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        facet_grid(~ Sample) + scale_fill_manual(values=i[[2]], drop = FALSE)
    p
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs_for_supplement/", i[[1]], "_no_legend.pdf"),
           width = 8, height = 6)
    
}

#And now, we are doing the same but for the ones that need a continuous scale


for(i in colnames(barGraphDat)[grep("mutations|Charge|Isoel|CDR3", colnames(barGraphDat))]){
    locDat <- as.data.frame(table(barGraphDat[,i], barGraphDat$Sample, barGraphDat$Specific))
    locDat[,1] <- as.numeric(as.character(locDat[,1]))
    colnames(locDat) <- c(i, "Sample","Specific", "value")
    if(grepl("mutations", i)){
        colorOption <- "G"
        scaleLims <- c(0, max(locDat[,1]))
    } else if(grepl("Charge", i)){
        colorOption <- "E"
        scaleLims <- c(min(locDat[,1]), max(locDat[,1]))
    } else if(grepl("Isoel", i)){
        colorOption <- "F"
        scaleLims <- c(min(locDat[,1]), max(locDat[,1]))
    } else if(grepl("CDR3", i)){
        colorOption <- "A"
        scaleLims <- c(min(locDat[,1]), max(locDat[,1]))
    }
    
    
    p <- ggplot(locDat,                         
                aes(x = Specific,
                    y = value,
                    fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        facet_grid(~ Sample) + scale_fill_viridis(limits = scaleLims, 
                                                  discrete = FALSE, 
                                                  option = colorOption)
    p
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Figure_4_plots/4B_CloneBarGraphs_for_supplement/", i[[1]], "_no_legend.pdf"),
           width = 8, height = 6)
}


#FOR FIGURE: 
table(substr(BCR_H$CELL, 1,4))
#1166 1227 1284 
#98   16   52

#Some statistics. 
#The V gene family question needs to be dealt with by ANOVA. 
#TO do all these, we add another category, namely the LGI1 patients together
newDat <- barGraphDat[which(barGraphDat$Sample %in% c("1166_1", "1227_1")),]
newDat$Sample <- "LGI1"
barGraphDat2 <- rbind(barGraphDat, newDat)

datSplit <- split(barGraphDat2, barGraphDat2$Sample)
unlist(lapply(datSplit, function(x){
    datDf <- as.data.frame(table(x$Specific, x$V_gene_fam))
    summary(aov(Freq ~ Var1, data = datDf))[[1]]$'Pr(>F)'[1]
}))
#   1166_1    1227_1    1284_1       All      LGI1 
#0.1139950 0.1822570 0.2211275 0.1238940 0.1153242

unlist(lapply(datSplit, function(x){
    datDf <- as.data.frame(table(x$Specific, x$V_gene_fam_L))
    summary(aov(Freq ~ Var1, data = datDf))[[1]]$'Pr(>F)'[1]
}))
#    1166_1     1227_1     1284_1        All       LGI1 
#0.16228226 0.06528795 0.11278756 0.09837277 0.13516677

#AS can be seen, the Pr values are all clearly insignificant, with the
#all analysis giving 0.124. 

#For Isotype, we divide into IgG4 and not
unlist(lapply(datSplit, function(x){
    x$IgG4 <- "Other"
    x$IgG4[which(x$ISOTYPE == "IGHG4")] <- "IgG4"
    datMat <- as.matrix(table(x$Specific, x$IgG4))
    fisher.test(datMat)$p.value
}))
#      1166_1       1227_1       1284_1          All         LGI1 
#0.5666465422 0.1357142857 0.0006715205 0.0037115955 0.2824445403 

#This is significant for all togheter, but when separating into donors, it is driven
#mainly by JR1284, and 1166 is not showing any signs of the same trend.

#Odds ratio for the combined: 
oddsDat <- datSplit[[4]]
oddsDat$IgG4 <- "Other"
oddsDat$IgG4[which(oddsDat$ISOTYPE == "IGHG4")] <- "IgG4"
oddsratio(table(oddsDat$IgG4, oddsDat$Specific), rev = "rows",
          method = "fisher")$measure
#       odds ratio with 95% C.I.
#       estimate     lower     upper
#Other 1.000000       NA       NA
#IgG4  3.214267 1.452963 7.262386

#Light type
unlist(lapply(datSplit, function(x){
    datMat <- as.matrix(table(x$Specific, x$light_type))
    fisher.test(datMat)$p.value
}))
#    1166_1     1227_1     1284_1        All       LGI1 
#0.77325039 0.26153846 0.08660473 0.06269744 0.30473277 

oddsratio(table(oddsDat$light_type, oddsDat$Specific))$measure
#   odds ratio with 95% C.I.
#   estimate     lower    upper
#K 1.0000000        NA       NA
#L 0.4590353 0.2070451 1.014351

focNames <- colnames(datSplit[[1]])[grep("mutations|Charge|Isoel|CDR3", colnames(datSplit[[1]]))]
do.call("rbind", lapply(datSplit, function(x){
    
    locRes <- unlist(lapply(focNames, function(y){
        specGroup <- x[which(x$Specific == "TRUE"),which(colnames(x) == y)]
        nonSpecGroup <- x[which(x$Specific == "FALSE"),which(colnames(x) == y)]
        wilcox.test(specGroup, nonSpecGroup, exact = FALSE)$p.value
    }))
    names(locRes) <- focNames
    locRes
}))

#       CDR3_length      Isoel     Charge All_mutations_H_L Non_silent_mutations_H_L
#1166_1 0.051753850 0.17861656 0.11398432        0.02102681               0.02050094
#1227_1 0.490067552 0.24367848 0.57703051        0.39031778               0.35880660
#1284_1 0.105742080 0.00225387 0.00587588        0.71863067               0.34041269
#All    0.006062449 0.01933258 0.01662936        0.14464783               0.21095821
#LGI1   0.020987701 0.39753469 0.29263663        0.03573621               0.03388839
#       CDR_mutations_H_L Non_silent_CDR_mutations_H_L All_mutations Non_silent_mutations
#1166_1        0.01545805                    0.1474349     0.2881002            0.9595261
#1227_1        0.26757532                    0.5807834     0.4271364            0.2141407
#1284_1        0.27933354                    0.1657060     0.8798416            0.5061666
#All           0.57918490                    0.9786757     0.4532041            0.8533768
#LGI1          0.10004473                    0.4137104     0.3514799            0.9769716
#       CDR_mutations Non_silent_CDR_mutations
#1166_1     0.6335600               0.28779656
#1227_1     0.1760614               0.15406949
#1284_1     0.5062490               0.32091814
#All        0.7483219               0.09745026
#LGI1       0.9106638               0.16799061

#And here, the medians: 
sigFocNames <- colnames(datSplit[[1]])[grep("Charge|Isoel|CDR3", colnames(datSplit[[1]]))]

mediansOfInterest <- do.call("rbind", lapply(sigFocNames, function(y){
    specGroupMed <- median(oddsDat[which(oddsDat$Specific == "TRUE"),which(colnames(oddsDat) == y)])
    nonSpecGroupMed <- median(oddsDat[which(oddsDat$Specific == "FALSE"),which(colnames(oddsDat) == y)])
    c(specGroupMed, nonSpecGroupMed)
}))
rownames(mediansOfInterest ) <- sigFocNames
colnames(mediansOfInterest ) <- c("Median_spec", "Median_non_spec")
mediansOfInterest 

#            Median_spec Median_non_spec
#CDR3_length        13.0              15
#Isoel               4.3               7
#Charge             -1.0               0


