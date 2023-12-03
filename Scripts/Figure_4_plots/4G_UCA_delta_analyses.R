#Here, we are going to investigate distances in mutations, EC50 and isotype switch, 
#for the first to within clones, from UCA to first member and from first member
#to last, and for the final one, number of switches within clones. 
library(ggforce)
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

#counts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "XBP1"),
#              which(aeSce$Specific == TRUE & aeSce$cellType == "B" &
#                        aeSce$Clonal == FALSE)]
BCR_dat <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")
#First observation that needs to be stressed is that four of the five specific B
#cells are clonal and in their clones there are also ASC in all cases. 

#We start with the EPDs. We have the EPD distribution known for eight clones.

EPDs <- read.csv("Data/BCR_auxiliaries/End-point_dilutions.csv")
#We know that the clones with more than one EPD value are: 
epdMoreThanOne <- c(5, 83, 84, 158, 202, 212, 214, 291)
#So now, we calculate the internal range of the log transformed final concentrations
#WE SHOULD NOT CALL IT END-POINT DILUTIONS BUT END-POINT CONCENTRATIONS
EPDs$logEPD <- log10(EPDs$EPD)

clonalEpdRanges <- unlist(lapply(epdMoreThanOne, function(x){
    locCells <- BCR_dat$CELL[which(BCR_dat$HamClone == x)]
    locLogEpds <- EPDs$logEPD[which(EPDs$CELL %in% locCells)]
    max(locLogEpds)-min(locLogEpds)
}))
names(clonalEpdRanges) <- paste0("Clone_", epdMoreThanOne)
clonalEpdRanges
#  Clone_5  Clone_83  Clone_84 Clone_158 Clone_202 Clone_212 Clone_214 Clone_291 
#0.0000000 0.0000000 0.3010300 0.6020600 0.4771213 0.6020601 0.6020600 0.0000000 

#What about the full range?
max(EPDs$logEPD)-min(EPDs$logEPD) #3.91339

EPD_UCA <- read.csv("Data/BCR_auxiliaries/UCA_EPD_and_mut.csv")
#The unique clones with an UCA are: 
unique(EPD_UCA$HamClone) #111 429 333 318 203 256  NA   5  96  20  44 461 325

#It is very much worth noting, too, that all the clones with more than one member
#have been tested here, and therefore,they can all be included at a high value as 20 is the highest value they were tested at. 
BCR_H <- BCR_dat[which(BCR_dat$LOCUS == "H"),]
BCR_foc <- BCR_H[which(BCR_H$HamClone %in% 
                           epdMoreThanOne[-which(epdMoreThanOne %in% unique(EPD_UCA$HamClone))]),]
extra_UCA <- data.frame("CELL" = BCR_foc$CELL, 
                        "HamClone" = BCR_foc$HamClone, 
                        "EPD_.ug.mL._UCA_BS" = 100,
                        "All_mutations" = BCR_foc$All_mutations)
EPD_UCA_full <- rbind(EPD_UCA, extra_UCA)

#Of the ones that are incomplete, all but 111 have the first value in the clone tested, after looking at the
#trees. For all of them
#all the cells are either unique, or it is only the closest one that has been tested. 
UCAPlusFirst <- unique(EPD_UCA_full$HamClone)
UCAPlusFirst <- UCAPlusFirst[-which(UCAPlusFirst == 111 | is.na(UCAPlusFirst))]

EPD_UCA_full$logEPD <- log10(EPD_UCA_full$EPD_.ug.mL._UCA_BS)

#So now, we calculate the EPD distance for this latter group, as we. have already done so
#for the first. 
ucaEpdRanges <- unlist(lapply(UCAPlusFirst, function(x){
    locCells <- BCR_dat$CELL[which(BCR_dat$HamClone == x)]
    locLogEpds <- EPDs$logEPD[which(EPDs$CELL %in% locCells)][1]
    
    locLogUcaEpds <- EPD_UCA_full$logEPD[which(EPD_UCA_full$HamClone == x)][1]
    locLogUcaEpds -locLogEpds
}))
names(ucaEpdRanges) <- paste0("Clone_", UCAPlusFirst)
ucaEpdRanges

# Clone_429  Clone_333  Clone_318  Clone_203  Clone_256    Clone_5   Clone_96   Clone_20 
#-0.3010300 -0.6020601  3.3113299  1.8061800  2.7092699  0.6020599  0.9030899  1.2041199 
#  Clone_44  Clone_461  Clone_325  Clone_202  Clone_212   Clone_83  Clone_214  Clone_158 
# 1.8061799  2.4082400  3.3113300  3.7092700  4.3113299  4.0103000  3.1072100  1.6020600 
#  Clone_84  Clone_291 
# 3.7092700  2.8061800 
clonalEpdPlotDf <- data.frame("Clone" = c(names(ucaEpdRanges), names(clonalEpdRanges)),
                          "variable" = c(rep("ucaEpdRange",length(ucaEpdRanges)), 
                                         rep("cloneEpdRange", length(clonalEpdRanges))), 
                          "EPD" = c(ucaEpdRanges, clonalEpdRanges))
    

#We need to make it even clearer that we are dealing with negative UCAs for some of them
negClones <- names(clonalEpdRanges)[-1]
#As a later addition, we need to bring in the specificity here. 
clonalEpdPlotDf$Specificity <- "LGI1"
clonalEpdPlotDf$Specificity[which(substring(clonalEpdPlotDf$Clone, 7) %in% 
                                      EPD_UCA_full$HamClone[which(substr(EPD_UCA_full$CELL, 1, 4) == "1284")])] <- "CASPR2"
clonalEpdPlotDf$SpecCol <- sapply(clonalEpdPlotDf$Specificity, switch, 
                                   "FALSE" = "black", 
                                   "LGI1" = "orange", 
                                   "CASPR2" = "#FF6633")

clonalEpdPlotDf$EPD[which(clonalEpdPlotDf$variable == "ucaEpdRange" &
                                clonalEpdPlotDf$Clone %in% negClones)] <- 4
clonalEpdPlotDf$variable <- factor(clonalEpdPlotDf$variable, levels = c("ucaEpdRange", "cloneEpdRange"))
set.seed(10)
plotPoss <- ggplot(clonalEpdPlotDf, aes(x = variable, y = EPD)) + 
    geom_sina(jitter_y = FALSE,maxwidth = 0.4,scale = "width")

x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y


dir.create("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots")
set.seed(10)
p <- ggplot(clonalEpdPlotDf, aes(x = variable, y = EPD)) + 
    geom_path(x = x_coords, aes(group = interaction(Clone), color = SpecCol)) + 
    geom_sina(jitter_y = FALSE, scale = "width", maxwidth = 0.4, size = 2.5, color = clonalEpdPlotDf$SpecCol) + 
    geom_violin(alpha = 0, scale = "width", width = 0.4, linewidth = 1) +
    theme_bw() + scale_color_manual(values = c("#FF6633", "orange"))
p
set.seed(10)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/EPD_left_is_diff_from_UCA_to_first_member_right_is_max_epd_dist_in_clone.pdf",
       width = 3, height = 3.5)

set.seed(10)
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
set.seed(10)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/EPD_left_is_diff_from_UCA_to_first_member_right_is_max_epd_dist_in_clone_no_text.pdf",
       width = 3, height = 3.5)

#Statistics: 
wilcox.test(clonalEpdPlotDf$EPD[which(clonalEpdPlotDf$variable == "ucaEpdRange")],
            clonalEpdPlotDf$EPD[which(clonalEpdPlotDf$variable == "cloneEpdRange")],
            paired = FALSE, exact = FALSE)
#p-value = 0.003841
median(clonalEpdPlotDf$EPD[which(clonalEpdPlotDf$variable == "ucaEpdRange")])
#3.0103
median(clonalEpdPlotDf$EPD[which(clonalEpdPlotDf$variable == "cloneEpdRange")])
#0.5395906

#Paired stats
clonalEpdParied <- clonalEpdPlotDf[which(clonalEpdPlotDf$Clone %in% 
                                             clonalEpdPlotDf$Clone[which(clonalEpdPlotDf$variable == "cloneEpdRange")]),]
clonalEpdPariedOrd <- clonalEpdParied[order(clonalEpdParied$Clone),]
median(clonalEpdPariedOrd$EPD[which(clonalEpdPariedOrd$variable == "ucaEpdRange")])
#4
median(clonalEpdPariedOrd$EPD[which(clonalEpdPariedOrd$variable == "cloneEpdRange")])
#0.5395906
wilcox.test(clonalEpdPariedOrd$EPD[which(clonalEpdPariedOrd$variable == "ucaEpdRange")],
            clonalEpdPariedOrd$EPD[which(clonalEpdPariedOrd$variable == "cloneEpdRange")],
            paired = TRUE, exact = FALSE)
#P-value 0.01403

#Now, we are going to investigate the intraclonal distances and compare these to
#the UCA to founder distances. 
#In an earlier iteration, it was clear that the small clones dominate the dataset by number,
#and as they per definition have a smaller intraclonal range, given that all clonal
#members only need to pass below a certain threshold of distance to one of the other
#cells in the same clone to be defined as clonal. Therefore, we divide the ones that 
#are larger than four from the rest. 

treeDists <- read.csv("Results/Figure_4_plots/4E_Trees/treeDists_founder_to_furthest.csv")

largeClones <- gsub("...._1_.+_|", "", treeDists$cloneNum[which(treeDists$size > 3)])
largeCloneDat <- clonalEpdPlotDf[which(gsub("Clone_|", "", clonalEpdPlotDf$Clone) %in% largeClones),]
set.seed(10)
plotPoss <- ggplot(largeCloneDat, aes(x = variable, y = EPD)) + geom_sina(jitter_y = FALSE,
                                                                            maxwidth = 0.4,
                                                                            scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

set.seed(10)
p <- ggplot(largeCloneDat, aes(x = variable, y = EPD)) + 
    geom_path(x = x_coords, aes(group = interaction(Clone), color = SpecCol)) + 
    geom_sina(jitter_y = FALSE, scale = "width", maxwidth = 0.4, size = 2.5, color = largeCloneDat$SpecCol) + 
    geom_violin(alpha = 0, scale = "width", width = 0.4, linewidth = 1) +
    theme_bw() + scale_color_manual(values = c("orange", "#FF6633"))
p
set.seed(10)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/Large_clones_EPD_left_is_diff_from_UCA_to_first_member_right_is_max_epd_dist_in_clone.pdf",
       width = 3, height = 3.5)

smallClones <- gsub("...._1_.+_|", "", treeDists$cloneNum[which(treeDists$size <= 3)])
smallCloneDat <- clonalEpdPlotDf[which(gsub("Clone_|", "", clonalEpdPlotDf$Clone) %in% smallClones),]
set.seed(10)
plotPoss <- ggplot(smallCloneDat, aes(x = variable, y = EPD)) + geom_sina(jitter_y = FALSE,
                                                                          maxwidth = 0.4,
                                                                          scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

set.seed(10)
p <- ggplot(smallCloneDat, aes(x = variable, y = EPD)) + 
    geom_path(x = x_coords, aes(group = interaction(Clone), color = SpecCol)) + 
    geom_sina(jitter_y = FALSE, scale = "width", maxwidth = 0.4, size = 2.5, color = smallCloneDat$SpecCol) + 
    geom_violin(alpha = 0, scale = "width", width = 0.4, linewidth = 1) +
    theme_bw()  + scale_color_manual(values = c("#FF6633", "orange"))
p
set.seed(10)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/Small_clones_EPD_left_is_diff_from_UCA_to_first_member_right_is_max_epd_dist_in_clone.pdf",
       width = 3, height = 3.5)

#One question here is: are the large clones less likely to have a positive UCA?
allUCAs <- read.csv("Data/BCR_auxiliaries/All_with_UCA.csv")
allUCAs$largeClone <- "no"
allUCAs$largeClone[which(allUCAs$HamClone %in% largeClones)] <- "yes"
table(allUCAs$largeClone, allUCAs$Specific_UCA_BS)
#     FALSE TRUE
#no    108   67
#yes    80    0

#So that looks convincing. However, what about if we do this per clone? 
allUCAsRed <- allUCAs[-which(duplicated(allUCAs$HamClone, incomparables = NA)),]
table(allUCAsRed$largeClone, allUCAsRed$Specific_UCA_BS)
#    FALSE TRUE
#no     47   29
#yes     7    0
fisher.test(table(allUCAsRed$largeClone, allUCAsRed$Specific_UCA_BS))
#p-value = 0.09014

treeDists$Specific <- sapply(gsub("...._._.+_|", "", treeDists$cloneNum), function(x){
    BCR_dat$Specific[which(BCR_dat$HamClone == x)][1]
})

#And the non-tested clones are excluded, as well as the ones where the data is incomplete. 
treeDists <- treeDists[-which(treeDists$Specific == "Not_tested"),]

##We start with the small
treeDistsSmall <- treeDists[which(treeDists$size %in% 2:3),]
wilcox.test(treeDistsSmall$germToFounder, 
            treeDistsSmall$maxToFounder, 
            paired = TRUE, exact = FALSE)
#p-value = 4.998e-07
summary(treeDistsSmall$germToFounder)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.02561 0.09854 0.14906 0.15633 0.18169 0.35721 
summary(treeDistsSmall$maxToFounder)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   0.000   0.000   0.020   0.017   0.145 

clonalDistPlotDf <- reshape2::melt(treeDistsSmall[,-which(colnames(treeDistsSmall) == "size")], id.vars = c("cloneNum", "Specific"))

#Here again, late in the day, we change the specificity, to include both LGI1 and CASPR2. 
clonalDistPlotDf$Specific[which(clonalDistPlotDf$Specific == "TRUE" &
                                    substr(clonalDistPlotDf$cloneNum, 1, 4) != "1284")] <- "LGI1"
clonalDistPlotDf$Specific[which(clonalDistPlotDf$Specific == "TRUE" &
                                    substr(clonalDistPlotDf$cloneNum, 1, 4) == "1284")] <- "CASPR2"

#And we order this to make the false come on top. 
clonalDistPlotDf <- clonalDistPlotDf[order(clonalDistPlotDf$variable, factor(clonalDistPlotDf$Specific), decreasing = c(FALSE, TRUE)),]
clonalDistPlotDf$SpecCol <- sapply(clonalDistPlotDf$Specific, switch, 
                                   "FALSE" = "black", 
                                   "LGI1" = "orange", 
                                   "CASPR2" = "#FF6633")
set.seed(2)
plotPoss <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) + geom_sina(maxwidth = 0.4,
                                                                              jitter_y = FALSE,
                                                                              scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

set.seed(2)
p <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) +
    geom_path(x = x_coords, aes(group = interaction(cloneNum)), color = clonalDistPlotDf$SpecCol) + 
    geom_sina(maxwidth = 0.4, jitter_y = FALSE, scale = "width", size = 2.5, color = clonalDistPlotDf$SpecCol) +
    geom_violin(alpha = 0, width = 0.4, scale = "width", linewidth = 1) + theme_bw() + 
    scale_y_continuous(limits = c(0,0.4))
p
set.seed(2)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_small_clones.pdf",
       width = 3, height = 3)
set.seed(2)
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
set.seed(2)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_small_clones_no_text.pdf",
       width = 3, height = 3.5)


#But if we exclude the small clones, the whole effect goes. 
treeDistsLarge <- treeDists[-which(treeDists$size %in% 2:3),]
wilcox.test(treeDistsLarge$germToFounder, 
            treeDistsLarge$maxToFounder, 
            paired = TRUE, exact = FALSE)
#p-value = 0.2361

#So now, we make a new figure version including only the big ones.
clonalDistPlotDf <- reshape2::melt(treeDistsLarge[,-which(colnames(treeDistsLarge) == "size")], id.vars = c("cloneNum", "Specific"))
#Here again, late in the day, we change the specificity, to include both LGI1 and CASPR2. 
clonalDistPlotDf$Specific[which(clonalDistPlotDf$Specific == "TRUE" &
                                    substr(clonalDistPlotDf$cloneNum, 1, 4) != "1284")] <- "LGI1"
clonalDistPlotDf$Specific[which(clonalDistPlotDf$Specific == "TRUE" &
                                    substr(clonalDistPlotDf$cloneNum, 1, 4) == "1284")] <- "CASPR2"

#And we order this to make the false come on top. 
clonalDistPlotDf <- clonalDistPlotDf[order(clonalDistPlotDf$variable, factor(clonalDistPlotDf$Specific), decreasing = c(FALSE, TRUE)),]
clonalDistPlotDf$SpecCol <- sapply(clonalDistPlotDf$Specific, switch, 
                                   "FALSE" = "black", 
                                   "LGI1" = "orange", 
                                   "CASPR2" = "#FF6633")

set.seed(2)
plotPoss <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) + geom_sina(maxwidth = 0.4,
                                                                               jitter_y = FALSE,
                                                                               scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

set.seed(2)
p <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) +
    geom_path(x = x_coords, aes(group = interaction(cloneNum), color = clonalDistPlotDf$SpecCol)) + 
    geom_sina(maxwidth = 0.4, jitter_y = FALSE, size = 2.5, scale = "width", color = clonalDistPlotDf$SpecCol) + 
    geom_violin(alpha = 0, width = 0.4, scale = "width", linewidth = 1) + theme_bw() + 
    scale_y_continuous(limits = c(0,0.4)) + scale_color_manual(values = c("#FF6633", "orange"))
p
set.seed(2)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_large_clones.pdf",
       width = 3, height = 3)
set.seed(2)
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
set.seed(2)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_large_clones_no_text.pdf",
       width = 3, height = 3.5)

