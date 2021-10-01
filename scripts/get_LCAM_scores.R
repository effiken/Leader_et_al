

###   GET_LCAM_SCORES.R

###   ANDREW M. LEADER
###   7/16/2020
###   ANDREW.LEADER@ICAHN.MSSM.EDU

###   DESCRIPTION:
###   Computes LCAM scores (as defined in Leader, et al. Biorxiv 2020) over a matrix of expression values from bulkRNAseq values.
###   This analysis is meant to be performed on resection specimens or biopsies of primary NSCLC tissue (such as in TCGA).

###   Usage:
###   LCAMscore_mat <-  get_LCAM_scores(expr_mat)

###   input: EXPR_MAT is a GENES x SAMPLES matrix of bulk RNA gene expression values

###   output: a SAMPLES x 3 dimension matrix.
###       Column 1 is an LCAMhi score per sample
###       Column 2 is an LCAMlo score per sample
###       Column 3 is equal to the difference (LCAMhi score - LCAMlo score); this is the main quantity that we correlate to.

library(matrixStats)

LCAMhi_subtypes <- c("IgG_plasma","SPP1_mac","T_activated")
LCAMlo_subtypes <- c("B","AM","cDC2","AZU1_mac","cDC1")
genes <- strsplit("IGHG3,IGHG4,IGHG1,MZB1,FAM92B,SPAG4,DERL3,JSRP1,TBCEL,LINC01485,SLC2A5,SPP1,CCL7,HAMP,CXCL13,ZBED2,GNG4,KRT86,TNFRSF9,LAYN,GZMB,KIR2DL4,GAPT,CTB-133G6.1,AKR1C3,RSPO3,MCEMP1,RND3,HP,FOLR3,GPD1,GS1-600G8.5,PCOLCE2,CAMP,CCL17,CD1B,CD1C,CD1E,PLD4,TRPC6,CLEC10A,FCER1A,PLA1A,CFP,CEACAM8,AZU1,PKP2,RETN,LYZ,ELANE,ATP6AP1L,RP6-159A1.4,THBS1,CFD,P2RY14,DNASE1L3,PTGER3,CST3,BATF3,CPVL,RGS18,SERPINF2,LGALS2",",")[[1]]
subtypes <- rep(c(LCAMhi_subtypes,LCAMlo_subtypes),times=c(10,4,8,2,10,10,10,9))
names(subtypes) <- genes




get_LCAM_scores <- function(expr_mat){
  
  message("computing LCAM scores on bulkRNA expression matrix")
  # identify matching gene symbols
  genes_matching <- intersect(genes,rownames(expr_mat))
  subtypes_matching <- subtypes[genes_matching]
  missing_genes <- setdiff(genes,genes_matching)
  if(!is.null(missing_genes)){
    message(paste("inconsistent gene symbols:",paste(missing_genes,collapse=", ")))
  }
  
  # normalize and scale expression values
  expr_mat <- t(t(expr_mat)/rowSums(t(expr_mat)))
  zmat <- t(scale(t(log10(1e-6+expr_mat[genes_matching,]))))
  
  #subset expression matrix to only the genes needed for the scores
  #zmat <- zmat[genes_matching,]
  
  # compute subtype scores
  score_mat <- matrix(0,nrow=length(genes_matching),ncol=length(unique(subtypes_matching)),dimnames=list(genes_matching,unique(subtypes_matching)))
  for(subtype_iter in subtypes_matching){
    score_mat[genes_matching[subtypes_matching==subtype_iter],subtype_iter] <- 1
  }
  
  subtype_scores <- t(scale(t(t(score_mat)%*%zmat)))
  
  #summarize the subtype scores
  res <- cbind(colMeans(subtype_scores[LCAMhi_subtypes,]),colMeans(subtype_scores[LCAMlo_subtypes,]))
  res <- cbind(res,res[,1]-res[,2])
  colnames(res) <- c("LCAMhi","LCAMlo","difference")
  return(res)

}
