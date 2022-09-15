

#Now we run a logistic regression model, more specifically a binominal such. 
binModDat <- BCR_reduced2[-which(BCR_reduced2$Specific == "Not_tested"),]
binModDat$Specific <- sapply(binModDat$Specific, switch, "no" = 0, "yes" = 1)
binModDat$Clonal <- sapply(as.character(binModDat$Clonal), switch, "FALSE" = 0, "TRUE" = 1)
binModDat$light_type <- sapply(as.character(binModDat$light_type), switch, 
                               "double" = 0.5, 
                               "K" = 0,
                               "L" = 1)
binModDat$ISOTYPE <- sapply(as.character(binModDat$ISOTYPE), 
                            switch, "IGHG1" = 0, "IGHG2" = 0.33, 
                            "IGHG3" = 0.66, "IGHG4" = 1)

binMod <- glm(Specific~Clonal+light_type+Mutations+ISOTYPE,
              family=binomial,data=binModDat) 
summary(binMod) 
#Gives: 
#Call:
#    glm(formula = Specific ~ Clonal + light_type + Mutations + ISOTYPE, 
#        family = binomial, data = binModDat)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.1897  -1.0244   0.5219   0.8633   1.8232  
#
#Coefficients:
#    Estimate Std. Error z value Pr(>|z|)   
#(Intercept)  1.88535    0.93058   2.026  0.04277 * 
#    Clonal      -0.16397    0.54235  -0.302  0.76240   
#light_type  -1.00427    0.65687  -1.529  0.12629   
#Mutations   -0.10765    0.03936  -2.735  0.00624 **
#    ISOTYPE      1.54958    0.59062   2.624  0.00870 **
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 108.267  on 83  degrees of freedom
#Residual deviance:  91.427  on 79  degrees of freedom
#AIC: 101.43
#
#Number of Fisher Scoring iterations: 4
#