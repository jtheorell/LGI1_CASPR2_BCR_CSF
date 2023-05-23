#Here, we are identifying the number of clonally expanded CD4 T cells. We base these
#calculations on Pappalardo et al Sci Imm 2020:
#https://pubmed.ncbi.nlm.nih.gov/32948672/
#This study provide the following numbers: 
#Total number of T cells: 7199
#Total number of unexpanded CD4 T cells: 4883
#Total number of unexpanded CD8 T cells: 751
#Total number of duplicated T cells: 792 
#Total number of duplicated CD4 T cells: 
#Total number of duplicated CD8 T cells: 
#Total number of highly expanded T cells: 708
#Total number of highly expanded CD4 T cells: 344
#Total number of highly expanded CD8 T cells: 312

#Furthermore, in supplementary table 5, tab "Healthy Clonal T clusters", they
#show that the CD8 T cells contribute about 30% to the duplicated pool. This 
#does not fully square with the fact that figure 4D shows that the CD8T cells seem
#to make up sligthly less than 25% of the duplicated pool. We therefore use the
#more conservative estimate seen from the CD4 T cell side, namely that 30% of the cells
#are not belonging to this pool. with that, it would mean that we have 751*0.7≈525 
#duplicated CD4 T cells. This gives us a total whooping percentage of clonality of 
344+512 #856
(344+512)/(4883+344+512) #15%. 

#After these nice deliberations, I was provided with a number from a more suitable 
#age range, namely a healthy control group from Mb Alzheimer patient study. Adam
#Handel has calculated in this dataset that 23% of the CD4 T cells are clonal in the CSF. 

#We then make our data. 
#The CD4 clonality data can be found at the end of the T cell clonality analysis page:
#          1166_1     1227_1    1284_1
#FALSE 0.97520661 0.97927461 0.8921283
#TRUE  0.02479339 0.02072539 0.1078717

#The B cell clonality information is not shown straight on anywhere, so that will
#Have to be retrieved. 
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
clonTab <- table(aeSce$donor, aeSce$Clonal)
clonTab/rowSums(clonTab)
#           FALSE      TRUE
#  1166 0.6453202 0.3546798
#  1227 0.8076923 0.1923077
#  1284 0.5833333 0.4166667
library(ggplot2)
#So, those are our datapoints then: 
ggDat <- data.frame("Donor" = c("#1", "#2", "#3"),
                    "T_clonal" = 100*c(0.02479339, 0.02072539, 0.1078717),
                    "B_clonal" = 100*c(0.3546798, 0.1923077, 0.4166667),
                    "Disease" = c("LGI1", "LGI1", "CASPR2"))

ggplot(ggDat, aes(x = T_clonal, y = B_clonal, fill = Disease)) + 
    geom_point(shape = 23, size = 4) + 
    scale_x_continuous(limits = c(0,50), expand = c(0,0)) + 
    scale_y_continuous(limits = c(-1,50), expand = c(0,0)) + 
    theme_bw() + geom_vline(xintercept = 23, linetype="dotted", 
                            linewidth=1) +
    geom_hline(yintercept = 0, linetype="dotted", 
                 size=1) + coord_fixed() + scale_fill_manual(values = c("white", "black"))
ggsave("Results/Figure_4_plots/T_vs_B_clonality.pdf", height = 5, width = 5)
