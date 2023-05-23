#Here, we are going to investigate distances in mutations, EC50 and isotype switch, 
#for the first to within clones, from UCA to first member and from first member
#to last, and for the final one, number of switches within clones. 
library(ggforce)
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

#counts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "XBP1"),
#              which(aeSce$Specific == TRUE & aeSce$cellType == "B" &
#                        aeSce$Clonal == FALSE)]
BCR_dat <- read.csv("Data/BCR_database_versions/8_new_mutations_added.csv")
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
unique(EPD_UCA$HamClone) #5  20  44  96 111 203 256 318 333 429  NA

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

#Clone_5   Clone_20   Clone_44   Clone_96  Clone_203  Clone_256  Clone_318  Clone_333 
#0.6020599  1.2041199  1.8061799  0.9030899  1.8061800  2.7092699  3.3113299 -0.6020601 
#Clone_429  Clone_202  Clone_212   Clone_83  Clone_214  Clone_158   Clone_84  Clone_291 
#-0.3010300  3.0314893  3.6335492  3.3325193  2.4294293  0.9242793  3.0314893  2.1283993
clonalEpdPlotDf <- data.frame("Clone" = c(names(ucaEpdRanges), names(clonalEpdRanges)),
                          "variable" = c(rep("ucaEpdRange",length(ucaEpdRanges)), 
                                         rep("cloneEpdRange", length(clonalEpdRanges))), 
                          "EPD" = c(ucaEpdRanges, clonalEpdRanges))
    

#We need to make it even clearer that we are dealing with negative UCAs for some of them
negClones <- names(clonalEpdRanges)[-1]

clonalEpdPlotDf$EPD[which(clonalEpdPlotDf$variable == "ucaEpdRange" &
                                clonalEpdPlotDf$Clone %in% negClones)] <- 4
clonalEpdPlotDf$variable <- factor(clonalEpdPlotDf$variable, levels = c("ucaEpdRange", "cloneEpdRange"))
set.seed(10)
plotPoss <- ggplot(clonalEpdPlotDf, aes(x = variable, y = EPD)) + geom_sina(jitter_y = FALSE,
                                                                            maxwidth = 0.4,
                                                                            scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

dir.create("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots")
set.seed(10)
p <- ggplot(clonalEpdPlotDf, aes(x = variable, y = EPD)) + 
    geom_path(x = x_coords, aes(group = interaction(Clone)), color = "orange") + 
    geom_sina(jitter_y = FALSE, scale = "width", maxwidth = 0.4, color = "orange", size = 2.5) + 
    geom_violin(alpha = 0, scale = "width", width = 0.4, linewidth = 1) +
    theme_bw() 
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

#We also check if the large clones stick out in any way compared to the small. 
treeDists <- read.csv("Results/Figure_4_plots/Trees/treeDists.csv")

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
    geom_path(x = x_coords, aes(group = interaction(Clone)), color = "orange") + 
    geom_sina(jitter_y = FALSE, scale = "width", maxwidth = 0.4, color = "orange", size = 2.5) + 
    geom_violin(alpha = 0, scale = "width", width = 0.4, linewidth = 1) +
    theme_bw() 
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
    geom_path(x = x_coords, aes(group = interaction(Clone)), color = "orange") + 
    geom_sina(jitter_y = FALSE, scale = "width", maxwidth = 0.4, color = "orange", size = 2.5) + 
    geom_violin(alpha = 0, scale = "width", width = 0.4, linewidth = 1) +
    theme_bw() 
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

#Now, we do the same for the IgPhyML distances. 
treeDists <- read.csv("Results/Figure_4_plots/Trees/treeDists.csv")

#In an earlier iteration, it is clear that the small clones dominate the dataset by number,
#and as they per definition have a smaller intraclonal range, given that all clonal
#members only need to pass below a certain threshold of distance to one of the other
#cells in the same clone to be defined as clonal. Therefore, we divide the ones that 
#are larger than four from the rest. 
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

#And we order this to make the false come on top. 
clonalDistPlotDf <- clonalDistPlotDf[order(clonalDistPlotDf$variable, factor(clonalDistPlotDf$Specific), decreasing = c(FALSE, TRUE)),]
set.seed(2)
plotPoss <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) + geom_sina(maxwidth = 0.4,
                                                                              jitter_y = FALSE,
                                                                              scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

clonalDistPlotDf$SpecCol <- sapply(clonalDistPlotDf$Specific, switch, "TRUE" = "orange", "FALSE" = "black")

set.seed(2)
p <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) +
    geom_path(x = x_coords, aes(group = interaction(cloneNum), color = Specific)) + 
    geom_sina(maxwidth = 0.4, jitter_y = FALSE, scale = "width", size = 2.5, color = clonalDistPlotDf$SpecCol) +
    geom_violin(alpha = 0, width = 0.4, scale = "width", linewidth = 1) + theme_bw() + 
    scale_color_manual(values = c("black", "orange")) + 
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
            treeDistsLarge$maxIntraClonal, 
            paired = TRUE, exact = FALSE)
#p-value = 0.9057

#So now, we make a new figure version including only the big ones.
clonalDistPlotDf <- reshape2::melt(treeDistsLarge[,-which(colnames(treeDistsLarge) == "size")], id.vars = c("cloneNum", "Specific"))
set.seed(2)
plotPoss <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) + geom_sina(maxwidth = 0.4,
                                                                               jitter_y = FALSE,
                                                                               scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

set.seed(2)
p <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) +
    geom_path(x = x_coords, aes(group = interaction(cloneNum), color = "orange")) + 
    geom_sina(maxwidth = 0.4, jitter_y = FALSE, size = 2.5, scale = "width", color = "orange") + 
    geom_violin(alpha = 0, width = 0.4, scale = "width", linewidth = 1) + theme_bw() + 
    scale_color_manual(values = c("orange")) + 
    scale_y_continuous(limits = c(0,0.4))
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

########################
#######################
#Here, we instead look at founder to furthest instead of max intraclonal. 
treeDists <- read.csv("Results/Figure_4_plots/Trees/treeDists_founder_to_furthest.csv")
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

clonalDistPlotDf <- reshape2::melt(treeDistsSmall[,-which(colnames(treeDistsSmall) == "size")], id.vars = c("cloneNum", "Specific"))
set.seed(2)
plotPoss <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) + geom_sina(maxwidth = 0.4,
                                                                               jitter_y = FALSE,
                                                                               scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

clonalDistPlotDf$SpecCol <- sapply(clonalDistPlotDf$Specific, switch, "TRUE" = "orange", "FALSE" = "black")

set.seed(2)
p <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) +
    geom_path(x = x_coords, aes(group = interaction(cloneNum), color = Specific)) + 
    geom_sina(maxwidth = 0.4, jitter_y = FALSE, scale = "width", size = 2.5, color = clonalDistPlotDf$SpecCol) +
    geom_violin(alpha = 0, width = 0.4, scale = "width", linewidth = 1) + theme_bw() + 
    scale_color_manual(values = c("black", "orange")) + 
    scale_y_continuous(limits = c(0,0.4))
p
set.seed(2)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_small_clones_max_to_founder.pdf",
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
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_small_clones_max_to_founder_no_text.pdf",
       width = 3, height = 3.5)


#But if we exclude the small clones, the whole effect goes. 
treeDistsLarge <- treeDists[-which(treeDists$size %in% 2:3),]
wilcox.test(treeDistsLarge$germToFounder, 
            treeDistsLarge$maxToFounder, 
            paired = TRUE, exact = FALSE)
#p-value = 0.2361

#So now, we make a new figure version including only the big ones.
clonalDistPlotDf <- reshape2::melt(treeDistsLarge[,-which(colnames(treeDistsLarge) == "size")], id.vars = c("cloneNum", "Specific"))
set.seed(2)
plotPoss <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) + geom_sina(maxwidth = 0.4,
                                                                               jitter_y = FALSE,
                                                                               scale = "width")
x_coords <- ggplot_build(plotPoss)$data[[1]]$x
y_coords <- ggplot_build(plotPoss)$data[[1]]$y

set.seed(2)
p <- ggplot(clonalDistPlotDf, aes(x = variable, y = value)) +
    geom_path(x = x_coords, aes(group = interaction(cloneNum), color = "orange")) + 
    geom_sina(maxwidth = 0.4, jitter_y = FALSE, size = 2.5, scale = "width", color = "orange") + 
    geom_violin(alpha = 0, width = 0.4, scale = "width", linewidth = 1) + theme_bw() + 
    scale_color_manual(values = c("orange")) + 
    scale_y_continuous(limits = c(0,0.4))
p
set.seed(2)
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_large_clones_max_to_founder.pdf",
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
ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/IgPhyML_delta_UCA_vs_intraclonal_large_clones_max_to_founder_no_text.pdf",
       width = 3, height = 3.5)

#DEPRECATED
#Now for the mutations. 
#mutDistData <- read.csv("Results/Figure_4_plots/Mutation_analysis/mutStatsPerCell.csv")
#clonalMutDist <- data.frame(do.call("rbind", lapply(unique(mutDistData$HamClone), function(x){
#    locClone <- mutDistData[which(mutDistData$HamClone == x),]
#    c(min(locClone$Germline_all_mut), max(locClone$maxAllMutDist))
#})))
#
#colnames(clonalMutDist) <- c("Min_germ_dist", "Max_intraclone_dist")
#clonalMutDist$Clone <- unique(mutDistData$HamClone)
#clonalMutDist$Specific <- sapply(clonalMutDist$Clone, function(x){
#    aeSce$Specific[which(aeSce$HamClone == x)][1]
#})
##Now, we remove rows with inf values. 
#clonalMutDist <- clonalMutDist[-which(is.infinite(clonalMutDist$Max_intraclone_dist)),]
#
##And the ones with unknown specificity
#clonalMutDist <- clonalMutDist[-which(clonalMutDist$Specific == "Not_tested" | 
#                                          is.na(clonalMutDist$Specific)),]

#clonalMutPlotDf <- reshape2::melt(clonalMutDist, id.vars = c("Clone", "Specific"))
#set.seed(1)
#plotPoss <- ggplot(clonalMutPlotDf, aes(x = variable, y = value)) + geom_sina(maxwidth = 0.4,
#                                                                              jitter_y = FALSE,
#                                                                              scale = "width")
#x_coords <- ggplot_build(plotPoss)$data[[1]]$x
#y_coords <- ggplot_build(plotPoss)$data[[1]]$y
#
#set.seed(1)
#p <- ggplot(clonalMutPlotDf, aes(x = variable, y = value)) +
#    geom_sina(maxwidth = 0.4, jitter_y = FALSE, scale = "width") + 
#    geom_path(x = x_coords, aes(group = interaction(Clone), color = Specific)) + 
#    geom_violin(alpha = 0, width = 0.4, scale = "width") +
#    theme_bw() +scale_y_continuous(limits = c(0,60), expand = c(0,0)) + 
#    scale_color_manual(values = c("black", "orange"))
#p
#set.seed(1)
#ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/Mutational_delta_UCA_vs_intraclonal.pdf",
#       width = 3, height = 3)
#set.seed(1)
#p + theme(legend.position = "none", 
#            panel.grid.major = element_blank(), 
#            panel.grid.minor = element_blank(),
#            axis.text.x = element_blank(),
#            axis.text.y = element_blank(),
#            strip.text.x = element_blank(),
#            axis.title.x = element_blank(),
#            axis.title.y = element_blank())
#set.seed(1)
#ggsave("Results/Figure_4_plots/4F_UCA_Mutations_EPD_plots/Mutational_delta_UCA_vs_intraclonal_no_text.pdf",
#       width = 3, height = 3)
#
##And some statistics. 
#wilcox.test(clonalMutDist$Min_germ_dist, 
#            clonalMutDist$Max_intraclone_dist, 
#            paired = TRUE, exact = FALSE)
##p-value = 3.302e-08