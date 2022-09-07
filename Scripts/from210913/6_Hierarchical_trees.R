#I am here heavily reliant on 
#https://www.antibodysociety.org/wordpress/wp-content/uploads/2021/11/Immcantation-webinar-slides.pdf
BCR_all <- read.csv("Complete_BCR_db_with_specific.csv")

BCR_clonal <- BCR_all[which(BCR_all$Clonal == TRUE),]
BCR_clonal$seq_id <- paste0(BCR_clonal$CELL, BCR_clonal$LOCUS)
clones = formatClones(BCR_clonal,id = "seq_id",
                      seq = "SEQUENCE_IMGT", 
                      germ = "GERMLINE_IMGT",
                      v_call = "V_CALL", j_call = "J_CALL", 
                      junc_len = "JUNCTION_LENGTH",
                      traits = c("Donor", "ISOTYPE", "Specific", "Cell_type"),
                      clone = "HamClone", heavy = "H", cell = "CELL", 
                      locus = "LOCUS", add_count = TRUE)

#This is now using maximum parsimony. 
trees = getTrees(clones,nproc=2)
plots = plotTrees(trees, tips="ISOTYPE",
                  tipsize=2)
plots[[1]]

treesToPDF(plots, file="Results/final_data_trees.pdf",
           nrow=2, ncol=2)