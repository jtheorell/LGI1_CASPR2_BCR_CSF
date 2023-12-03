#Here, we remake the plot for 2B, as Ruby clearly shows something else. 
library(ggplot2)
library(DepecheR)
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

BCR_h_all <- BCR_all[which(BCR_all$LOCUS == "H"),]
BCR_h_all$ClonalSize <- 1
BCR_h_all$ClonalSize[which(BCR_h_all$Clonal)] <- sapply(BCR_h_all$CELL[which(BCR_h_all$Clonal)], function(x){
    locCloneSize <- length(which(BCR_h_all$HamClone == BCR_h_all$HamClone[which(BCR_h_all$CELL == x)]))
    if(locCloneSize >= 4){
        "≥4"
    } else if(locCloneSize >= 2){
        "2-3"
    } else {
        as.character(locCloneSize)   
        }
})
BCR_h_all$ClonalSize <- factor(BCR_h_all$ClonalSize, levels = c("1", "2-3", "≥4"))

#There are also a few cells that bind tissue but that do not bind either CASPR2
#or LGI1. These are: 
BCR_h_all$SpecTiss <- BCR_h_all$Specific_antigen
BCR_h_all$SpecTiss[which(BCR_h_all$CELL %in% c("1166_1_2_H12",
                                       "1166_1_2_A15",
                                       "1227_1_2_J08", 
                                       "1227_1_2_O05",
                                       "1166_1_2_J05", 
                                       "1166_1_2_L15",
                                       "1166_1_2_C11", 
                                       "1166_1_2_H08", 
                                       "1166_1_2_N07"))] <- "Tissue"

BCR_h_all$SpecTiss <- factor(BCR_h_all$SpecTiss, levels = c("Tissue", "FALSE", "Not_tested","CASPR2", "LGI1"))
BCR_h_all$Specific_antigen <- factor(BCR_h_all$Specific_antigen, levels = c("FALSE", "Not_tested", "CASPR2", "LGI1"))

#Here, the analysis forks, briefly, focusing on the tested cells only. 
BCR_h <- BCR_h_all[which(BCR_h_all$Specific_antigen != "Not_tested"),]

clonSpecMat <- as.matrix(table(BCR_h$ClonalSize, BCR_h$Specific_antigen))
clonSpecDf <- as.data.frame(clonSpecMat)
colnames(clonSpecDf) <- c("CloneSize", "Specific_antigen", "Freq")
ggplot(clonSpecDf, aes(x = CloneSize, y = Freq, fill = Specific_antigen)) + 
    geom_bar(stat="identity") + theme_bw() + 
    scale_fill_manual(values = c("black", "grey", "#FF6633", "orange")) + 
    scale_y_continuous(limits = c(0, 85), expand = c(0, 0, 0, 0))
ggsave("Results/Figure_2_plots/2B_clonality_vs_specificity.pdf", width = 4.5, height = 4)

#And the percentages: 
round(100*t(apply(clonSpecMat, 1, function(x) x/sum(x))))
#      FALSE Not_tested CASPR2 LGI1
#  1      38          0      5   57
#  2-3    21          0     17   62
#  ≥4      0          0     56   44

#Now, a version splitting the donors
sizeDon <- paste0(BCR_h$ClonalSize, "_", BCR_h$Sample)
clonSpecMat <- as.matrix(table(sizeDon, BCR_h$Specific_antigen))
clonSpecDf <- as.data.frame(clonSpecMat)
colnames(clonSpecDf) <- c("CloneSize", "Specific_antigen", "Freq")
clonSpecDf$CloneSize <- factor(clonSpecDf$CloneSize, 
                               levels = c("1_1166_1", "2-3_1166_1", "≥4_1166_1",
                                          "1_1227_1", "2-3_1227_1", "≥4_1227_1",
                                          "1_1284_1", "2-3_1284_1", "≥4_1284_1"))


ggplot(clonSpecDf, aes(x = CloneSize, y = Freq, fill = Specific_antigen)) + 
    geom_bar(stat="identity") + theme_bw() + 
    scale_fill_manual(values = c("black", "grey", "#FF6633", "orange")) + 
    scale_y_continuous(limits = c(0, 60), expand = c(0, 0, 0, 0))  +
    scale_x_discrete(drop=FALSE)
ggsave("Results/Figure_2_plots/2B_clonality_vs_specificity_per_donor.pdf", 
       width = 6, height = 4)

round(100*t(apply(clonSpecMat, 1, function(x) x/sum(x))))

#sizeDon      FALSE Not_tested TRUE
#  ≥4_1166_1      0          0  100
#  ≥4_1284_1      0          0  100
#  1_1166_1      27          0   73
#  1_1227_1      29          0   71
#  1_1284_1      78          0   22
#  2-3_1166_1    22          0   78
#  2-3_1227_1    22          0   78
#  2-3_1284_1    19          0   81

#And we also create one for the supplement with the percentages of other neuronal 
#tissue-specific cells in it. 

clonSpecMat <- as.matrix(table(BCR_h$ClonalSize, BCR_h$SpecTiss))
clonSpecDf <- as.data.frame(clonSpecMat)
colnames(clonSpecDf) <- c("CloneSize", "Specific_incl_tissue", "Freq")
ggplot(clonSpecDf, aes(x = CloneSize, y = Freq, fill = Specific_incl_tissue)) + 
    geom_bar(stat="identity") + theme_bw() + 
    scale_fill_manual(values = c("yellow","black", "grey", "#FF6633", "orange")) + 
    scale_y_continuous(limits = c(0, 80), expand = c(0, 0, 0.05, 0))
ggsave("Results/Figure_2_plots/2B_clonality_vs_specificity_with_tissue.pdf", width = 4.5, height = 4)

#And the percentages: 
round(100*t(apply(clonSpecMat, 1, function(x) x/sum(x))))
#      Tissue FALSE Not_tested CASPR2 LGI1
#  1        5    33          0      5   57
#  2-3      9    12          0     17   62
#  ≥4       0     0          0     56   44

#And the same but for individual donors. 
clonSpecMat <- as.matrix(table(sizeDon, BCR_h$SpecTiss))
clonSpecDf <- as.data.frame(clonSpecMat)
colnames(clonSpecDf) <- c("CloneSize", "specTiss", "Freq")
clonSpecDf$CloneSize <- factor(clonSpecDf$CloneSize, 
                               levels = c("1_1166_1", "2-3_1166_1", "≥4_1166_1",
                                          "1_1227_1", "2-3_1227_1", "≥4_1227_1",
                                          "1_1284_1", "2-3_1284_1", "≥4_1284_1"))

ggplot(clonSpecDf, aes(x = CloneSize, y = Freq, fill = specTiss)) + 
    geom_bar(stat="identity") + theme_bw() + 
    scale_fill_manual(values = c("yellow", "black", "grey", "#FF6633", "orange")) + 
    scale_y_continuous(limits = c(0, 60), expand = c(0, 0, 0, 0))  +
    scale_x_discrete(drop=FALSE)
ggsave("Results/Figure_2_plots/2B_clonality_vs_specificity_with_tissue_per_donor.pdf", 
       width = 6, height = 4)

round(100*t(apply(clonSpecMat, 1, function(x) x/sum(x))))
#sizeDon      Tissue FALSE Not_tested CASPR2 LGI1
#  ≥4_1166_1       0     0          0      0  100
#  ≥4_1284_1       0     0          0    100    0
#  1_1166_1        8    19          0      0   73
#  1_1227_1        0    29          0      0   71
#  1_1284_1        0    78          0     22    0
#  2-3_1166_1     10    12          0      0   78
#  2-3_1227_1     22     0          0      0   78
#  2-3_1284_1      0    19          0     81    0

#Now, we will plot the percentages of observed non-binders, untested and binders. 
specTab <- table(BCR_h_all$Specific_antigen, BCR_h_all$Sample)
ggDat <- data.frame(specTab)
colnames(ggDat) <- c("Specific_antigen", "Sample", "Freq")
ggplot(ggDat, aes(x=Sample, y=Freq, fill=Specific_antigen)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("black", "grey", "#FF6633", "orange")) +
    theme_void()
ggsave("Results/Figure_2_plots/Percentage_specific.pdf")

#Now, for the percentages
tSpecTab <- t(specTab)
tSpecTab/rowSums(tSpecTab)
#              FALSE Not_tested     CASPR2       LGI1
#  1166_1 0.08490566 0.53773585 0.00000000 0.37735849
#  1227_1 0.06153846 0.75384615 0.00000000 0.18461538
#  1284_1 0.09615385 0.50000000 0.40384615 0.00000000

specTissTab <- table(BCR_h_all$SpecTiss, BCR_h_all$Sample)
ggDat <- data.frame(specTissTab)
colnames(ggDat) <- c("SpecTiss", "Sample", "Freq")
ggplot(ggDat, aes(x=Sample, y=Freq, fill=SpecTiss)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("yellow", "black", "grey", "#FF6633",  "orange")) +
    theme_void()
ggsave("Results/Figure_2_plots/Percentage_specific_plus_tissue.pdf")

#Now, for the percentages
tSpecTissTab <- t(specTissTab)
tSpecTissTab/rowSums(tSpecTissTab)
#             Tissue      FALSE Not_tested     CASPR2       LGI1
#  1166_1 0.03301887 0.05188679 0.53773585 0.00000000 0.37735849
#  1227_1 0.03076923 0.03076923 0.75384615 0.00000000 0.18461538
#  1284_1 0.00000000 0.09615385 0.50000000 0.40384615 0.00000000

#Now, as there turns out to be such a skew and as we are skewing the percentages
#of specific cells towards the populations with 100% specificity, we will redo 
#the calculations, based on the total number of found sequences, including the 
#non-tested ones. We will do this per donor, per clonal size. 
#To simplify things, we will make a version of this, with TRUE instead of LGI1
#and CASPR2. 
BCR_h_all$specTissSimple <- as.character(BCR_h_all$SpecTiss)
BCR_h_all$specTissSimple[which(BCR_h_all$SpecTiss %in% c("CASPR2", "LGI1"))] <- "TRUE"
estSpecTiss <- as.table(do.call("rbind", lapply(split(BCR_h_all, BCR_h_all$Sample), function(x){
    colSums(do.call("rbind", lapply(split(x, x$ClonalSize), function(y){
        if(nrow(y) == 0){
            NULL
        } else {
            if(any(y$SpecTiss == "Not_tested")){
                locKnown <- y[-which(y$specTissSimple == "Not_tested"),]
            } else {
                locKnown <- y
            }
            knownTab <- table(locKnown$specTissSimple)
            knownPerc <- knownTab/sum(knownTab)
            nonTestSize <- length(which(y$specTissSimple == "Not_tested"))
            approxSpec <- unname(nonTestSize*knownPerc["TRUE"])
            if(is.na(knownPerc["FALSE"]) == FALSE){
                approxNonSpec <- unname(round(nonTestSize*knownPerc["FALSE"]))
            } else {
                approxNonSpec <- 0
            }
            
            c("EstTissue" = round(nonTestSize-(approxSpec+approxNonSpec)),
              "EstFALSE" = round(approxNonSpec),
              "EstTRUE" = round(approxSpec))
        }
        
    })))
})))

#Now, we bring the CASPR2/LGI1 division back. 
estSpecTissDf <- as.data.frame.matrix(estSpecTiss)
estSpecTissDf$EstLGI1 <- estSpecTissDf$EstCASPR2 <- estSpecTissDf$EstTRUE
estSpecTissDf$EstTRUE <- NULL
estSpecTissDf$EstCASPR2[-which(row.names(estSpecTissDf) == "1284_1")] <- 0
estSpecTissDf$EstLGI1[which(row.names(estSpecTissDf) == "1284_1")] <- 0
estSpecTissDf
#       EstTissue EstFALSE EstCASPR2 EstLGI1
#1166_1        10       21         0      83
#1227_1         1       13         0      35
#1284_1         0       39        13       0

estSpec <- estSpecTissDf[,2:4]
estSpec[,1] <- estSpec[,1]+estSpecTissDf[,1]
estSpec
#       EstFALSE EstCASPR2 EstLGI1
#1166_1       31         0      83
#1227_1       14         0      35
#1284_1       39        13       0


#And now, a figure. 
ggDatEstSpec <- reshape2::melt(rbind(specTab[-2,], t(estSpec)))
colnames(ggDatEstSpec) <- c("Specific", "Sample", "Freq")
ggDatEstSpec$Specific <- factor(ggDatEstSpec$Specific, levels = c("FALSE", "EstFALSE",
                                                         "EstCASPR2", "CASPR2", "EstLGI1", "LGI1"))
#Colors are identified for the estimates: 
dColorVector(1:3, colorScale = c("black", "grey"))
#5F5F5F
dColorVector(1:3, colorScale = c("#FF6633", "grey"))
#DE9278
dColorVector(1:3, colorScale = c("orange", "grey"))
#DEB15F

ggplot(ggDatEstSpec, aes(x=Sample, y=Freq, fill=Specific)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("black", "#5F5F5F","#DE9278", "#FF6633", "#DEB15F", "orange")) +
    theme_void()
ggsave("Results/Figure_2_plots/Percentage_specific_with_est.pdf")

#The tissue version
ggDatEstTiss <- reshape2::melt(rbind(specTissTab[-3,],t(estSpecTissDf)))
colnames(ggDatEstTiss) <- c("Specific", "Sample", "Freq")
ggDatEstTiss$Specific <- factor(ggDatEstTiss$Specific, levels = c("EstTissue", "Tissue",
                                                                  "FALSE", "EstFALSE",
                                                                  "EstCASPR2", "CASPR2", 
                                                                  "EstLGI1", "LGI1"))
#Colors are identified for the estimates: 
dColorVector(1:3, colorScale = c("yellow", "grey"))
#DEDE5F

ggplot(ggDatEstTiss, aes(x=Sample, y=Freq, fill=Specific)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("#DEDE5F", "yellow", 
                                 "black", "#5F5F5F", 
                                 "#DE9278", "#FF6633", 
                                 "#DEB15F", "orange")) +
    theme_void()
ggsave("Results/Figure_2_plots/Percentage_specific_plus_tissue_with_est.pdf")

#And finally, percentages. Here, we will look at the total, observed plus estimated, 
#only
totalSpecTab <- specTab[-2,] + t(estSpec)

apply(totalSpecTab, 2, function(x) x/sum(x))
#            1166_1    1227_1    1284_1
#  FALSE  0.2311321 0.2769231 0.4711538
#  CASPR2 0.0000000 0.0000000 0.5288462
#  LGI1   0.7688679 0.7230769 0.0000000

totalTissTab <- specTissTab[-3,]+t(estSpecTissDf)

apply(totalTissTab, 2, function(x) x/sum(x))
#             1166_1     1227_1    1284_1
#  Tissue 0.08018868 0.04615385 0.0000000
#  FALSE  0.15094340 0.23076923 0.4711538
#  CASPR2 0.00000000 0.00000000 0.5288462
#  LGI1   0.76886792 0.72307692 0.0000000

#Now, finally, we make a simplified version of the graphs, excluding the non-tested. 
#Now, we will plot the percentages of observed non-binders, untested and binders.
tested <- BCR_h_all[which(BCR_h_all$Specific != "Not_tested"),]
specTab <- table(tested$Specific_antigen, tested$Sample)
ggDat <- data.frame(specTab)
colnames(ggDat) <- c("Specific_antigen", "Sample", "Freq")
ggplot(ggDat, aes(x=Sample, y=Freq, fill=Specific_antigen)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("black", "grey", "#FF6633",  "orange")) +
    theme_void()
ggsave("Results/Figure_2_plots/Percentage_specific_no_untested.pdf")

apply(specTab, 2, function(x) x/sum(x))
#                1166_1 1227_1    1284_1
#  FALSE      0.1836735   0.25 0.1923077
#  Not_tested 0.0000000   0.00 0.0000000
#  CASPR2     0.0000000   0.00 0.8076923
#  LGI1       0.8163265   0.75 0.0000000

specTissTab <- table(tested$SpecTiss, tested$Sample)
ggDat <- data.frame(specTissTab)
colnames(ggDat) <- c("SpecTiss", "Sample", "Freq")
ggplot(ggDat, aes(x=Sample, y=Freq, fill=SpecTiss)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("yellow", "black", "grey", "#FF6633",  "orange")) +
    theme_void()
ggsave("Results/Figure_2_plots/Percentage_specific_plus_tissue_no_untested.pdf")

apply(specTissTab, 2, function(x) x/sum(x))

#                 1166_1 1227_1    1284_1
#  Tissue     0.07142857  0.125 0.0000000
#  FALSE      0.11224490  0.125 0.1923077
#  Not_tested 0.00000000  0.000 0.0000000
#  CASPR2     0.00000000  0.000 0.8076923
#  LGI1       0.81632653  0.750 0.0000000

