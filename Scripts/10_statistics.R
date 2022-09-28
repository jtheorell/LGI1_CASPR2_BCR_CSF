#This starts like the figure 2 data. 

BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")
BCR_H_tested <- BCR_all[which(BCR_all$LOCUS == "H" & BCR_all$Specific != "Not_tested"),]

#In this version, we will collapse the BCRs within the same subclone, as defined
#by the hierarchical tree analysis. 
BCR_reduced <- BCR_H_tested[,c("Sample","HamClone", "Mutations", "light_type", "ISOTYPE",
                        "Specific")]
#BCR_reduced$n_identical <- 1

#Now, we separate information about the cell type, the isotype and the light type
#into separate columns
for(i in c("light_type", "ISOTYPE")){
    newCols <- unique(BCR_reduced[,i])
    newDf <- as.data.frame(matrix(data = 0, nrow = nrow(BCR_reduced), 
                                  ncol = length(newCols),
                                  dimnames= list(seq(1, nrow(BCR_reduced)),newCols)))
    for(j in newCols){
        newDf[,j][which(BCR_reduced[,i] == j)] <- 1
    }
    BCR_reduced <- cbind(BCR_reduced[,-which(colnames(BCR_reduced) == i)], newDf)
}

BCR_reduced$Specific01 <- 0
BCR_reduced$Specific01[which(BCR_reduced$Specific == TRUE)] <- 1
BCR_reduced$Clonal <- 0
BCR_reduced$Clonal[which(is.na(BCR_reduced$HamClone))] <- 1
BCR_reduced$Mutations <- BCR_reduced$Mutations/max(BCR_reduced$Mutations)
#Given that we only have three independent samples, this cannot be 
#considered valid from a statistical point of view. However, it is almost
#entirely negative, and thus might still give an indication. Only kappa
#comes out as potentially significant. This will be tested also with a WIlcoxon below

binMod <- glm(Specific01~Clonal+Mutations+K+L+double+IGHG1+IGHG2+IGHG4,
              family=binomial,data=BCR_reduced) 
summary(binMod) 
#Call:
#    glm(formula = Specific01 ~ Clonal + Mutations + K + L + double + 
#            IGHG1 + IGHG2 + IGHG4, family = binomial, data = BCR_reduced)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.1671  -1.1062   0.5476   0.7963   1.5230  
#
#Coefficients: (2 not defined because of singularities)
#Estimate Std. Error z value Pr(>|z|)  
#(Intercept)   0.9036     0.7590   1.190   0.2339  
#Clonal       -0.6002     0.4025  -1.491   0.1359  
#Mutations    -1.4093     1.0950  -1.287   0.1981  
#K             1.4096     0.6642   2.122   0.0338 *
#    L             0.3053     0.6798   0.449   0.6533  
#double            NA         NA      NA       NA  
#IGHG1        -0.3853     0.4393  -0.877   0.3805  
#IGHG2        -0.8024     0.5781  -1.388   0.1651  
#IGHG4             NA         NA      NA       NA  
#---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 197.83  on 165  degrees of freedom
#Residual deviance: 179.23  on 159  degrees of freedom
#AIC: 193.23
#
#Number of Fisher Scoring iterations: 4

BCR_reducedSpecific <- BCR_reduced[which(BCR_reduced$Specific == "TRUE"),]
BCR_reducedNonSpecific <- BCR_reduced[which(BCR_reduced$Specific == "FALSE"),]

table(BCR_reducedSpecific$K, BCR_reducedSpecific$Sample)
#1166_1 1227_1 1284_1
#0     20      4     10
#1     58      8     19

table(BCR_reducedNonSpecific$K, BCR_reducedNonSpecific$Sample)
#1166_1 1227_1 1284_1
#0      8      3     15
#1     12      1      8
#WHen running a Fisher test on this, i.e. on the table: 
#        K dominant    L dominant
#Spec     3             0
#NonSpec  1             2

#We get a p-value of 0.4. 
#Even in the ideal scenario, i.e. 3,0,0,3, we would only get a p-value of 0.1
#Due to the small sample size. 

#When running a Wilcoxon: 
specKTab <- as.data.frame.matrix(table(BCR_reducedSpecific$K, BCR_reducedSpecific$Sample))
nonSpecKTab <- as.data.frame.matrix(table(BCR_reducedNonSpecific$K, BCR_reducedNonSpecific$Sample))

specKTabFreq <- apply(specKTab, 2, function(x) x/sum(x))

nonSpecKTabFreq <- apply(nonSpecKTab, 2, function(x) x/sum(x))

wilcox.test(specKTabFreq[2,], nonSpecKTabFreq[2,], paired = TRUE)
#Gives a p-value of 0.25. 

#So essentially we see no differences that can be statistically verified, 
#simply due to the low numbers. 

#Another thing we want to investigate is the Levenshtein distance between
#the non-mutated chains in the specific and non-specific compartments
germDists <- adist(BCR_H_tested$GERMLINE_IMGT[which(BCR_H_tested$Specific == "TRUE")],
      BCR_H_tested$GERMLINE_IMGT[which(BCR_H_tested$Specific == "FALSE")])

gerDistMins <- apply(germDists, 2, min)
summary(gerDistMins)

#And then we compare this to the non-tested cells, to get a general idea. 
BCR_allHGermUntested <- BCR_all$GERMLINE_IMGT[which(BCR_all$LOCUS == "H" & BCR_all$Specific == "Not_tested" &
                                                        BCR_all$Clonal == FALSE)]
germDistsUntested <- adist(BCR_allHGermUntested,BCR_allHGermUntested)
diag(germDistsUntested) <- 1000
germDistsUntestedMins <- 
    apply(germDistsUntested, 2, min)
summary(germDistsUntestedMins)

#And this was generally very unimpressive too, with some of the 
#distance between the singleton sequences in the non-tested group 
#being closer than the non-specific and the specific.  




