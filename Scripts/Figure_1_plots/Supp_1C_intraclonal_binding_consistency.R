#Here, we are going to make a chart showing the intraclonal binding consistency, 
#motivating why we later define all cells from each clone as either binders or
#non-bnders depending on if one cell has shown to be a binder or not. 
#first_clonal_selection <- read.csv("../External/data/20220426 BCR_selected_db_with exclusions_light chains-used_for_intraclonal_consistency.csv")

original_selection <- read.csv("Data/BCR_auxiliaries/Originally_selected_BCRs.csv")[,1]
BCR_all <- read.csv( "Data/BCR_database_versions/7_Subclones_included.csv")

#Now we zoom in
clonal1166Sel <- BCR_all[which(BCR_all$Sample == "1166_1" & 
                                BCR_all$Clonal &
                                   BCR_all$CELL %in% original_selection &
                                   BCR_all$LOCUS == "H"),]

#Here, we zoom in even further, only considering the clones with more than 
#one member:
isCloneLargeTested <- unlist(lapply(split(clonal1166Sel$HamClone, clonal1166Sel$HamClone), 
                                    function(x){
                                        if(length(x) > 1){
                                            TRUE
                                        } else {
                                            FALSE
                                        }
                                    }))

largeTested <- names(isCloneLargeTested)[which(isCloneLargeTested)]

BCRLargeTested <- clonal1166Sel[which(clonal1166Sel$HamClone %in% largeTested),]

#After creating this, it turns out that we lose clone 354, that was originally tested. 
#unclear why it is not registered as such. 
BCRLargeTested <- rbind(BCRLargeTested, BCR_all[which(BCR_all$HamClone == "354" &
                                                          BCR_all$LOCUS == "H"),])
table(BCRLargeTested$HamClone, BCRLargeTested$Specific)
BCRLargeTested$Specific <- factor(BCRLargeTested$Specific, levels = c("FALSE", "TRUE"))

testTab <- data.frame(table(BCRLargeTested$HamClone, BCRLargeTested$Specific))
colnames(testTab) <- c("Clone", "Binder", "Freq")
library(ggplot2)
ggplot(testTab, aes(fill=Binder, y=Freq, x=Clone)) + 
    geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("black", "orange"))+
    theme_bw()
ggsave("Results/CloneBarGraphs/Specificity_per_clone.pdf", width = 6, height = 4)


