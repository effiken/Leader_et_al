library(matrixStats)
library(Matrix)
library(Matrix.utils)
library(seriation)
library(scales)
library(data.table)
library(tidyr)
library(viridis)
library(CePa)

figure_6ac_s6abcdefjklmnopqr_7abcd_s7ab <- function(){

rm(list=ls())

scClustering_dir <- "scripts/scClustering/"
tcga_dir <- "data"


output_dir <- "output/figures/"
source(file.path(scClustering_dir,"DE.r"))

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)

annots_list$norm_group[annots_list$lineage=="MNP"] <- "MNP"

load("data/lung_ldm.rd")


immune_ep_l2fc <- read.csv("input_tables/immune_vs_ep_de.csv",r=1,h=1,stringsAsFactors = F)

score1_pats <- strsplit("408,403,714,522,371,570",",")[[1]]
score2_pats <- strsplit("571,596,393,593,626,378",",")[[1]]

deg_pairs <- list()
deg_pairs$Bplasma <- list()
deg_pairs$Bplasma[[1]] <- "IgG"
deg_pairs$Bplasma[[2]] <- "B"
deg_pairs$Mac <- list()
deg_pairs$Mac[[1]] <- c("MoMac-II")
deg_pairs$Mac[[2]] <- c("AM","cDC2","AZU1_mac","cDC1")
deg_pairs$Tcells <- list()
deg_pairs$Tcells[[1]] <- c("T_activated")
deg_pairs$Tcells[[2]] <- c("Tcm/naive_II")

tissue <- sample_annots[lung_ldm$dataset$cell_to_sample,"tissue"]
patient <- sample_annots[lung_ldm$dataset$cell_to_sample,"patient_ID"]
prime <- sample_annots[lung_ldm$dataset$cell_to_sample,"prime"]
lineage <- annots_list[as.character(lung_ldm$dataset$cell_to_cluster),"lineage"]
sub_lineage <- annots_list[as.character(lung_ldm$dataset$cell_to_cluster),"sub_lineage"]
#total exprs

bg_mask <- tissue=="Tumor" & prime=="3" & patient%in%score2_pats & lineage!="epi_endo_fibro_doublet"
fg_mask <- tissue=="Tumor" & prime=="3" & patient%in%score1_pats & lineage!="epi_endo_fibro_doublet"


#if we want to subsample evenly:
names(patient) <- names(lung_ldm$dataset$cell_to_sample)
bg_list <- split(names(patient)[bg_mask],factor(patient[bg_mask]))
bg_list <- lapply(bg_list,sample,1409)
bg_mask <- unlist(bg_list)
fg_list <- split(names(patient)[fg_mask],factor(patient[fg_mask]))
fg_list <- lapply(fg_list,sample,1409)
fg_mask <- unlist(fg_list)

# bg_mask <- colnames(lung_ldm$dataset$umitab)[sample(which(bg_mask),5000)]
# fg_mask <- colnames(lung_ldm$dataset$umitab)[sample(which(fg_mask),5000)]

names(patient) <- names(lung_ldm$dataset$cell_to_sample)
#table(patient[bg_mask])
#table(patient[fg_mask])

# DE_total <- DE_between_two_sets(ldm = lung_ldm,mask_bg = bg_mask,mask_fg = fg_mask,nmin_umi_thresh = 0,nmin_cells_with_min_umi = 3,
#                                 reg=1e-6,nchunks = 50,n_per_chunk = 1e3,noise_correction = F)
#save(DE_total,file="input_tables/DE_LCAMhi_vs_LCAMlo_pseudobulk.rd")
load("input_tables/DE_LCAMhi_vs_LCAMlo_pseudobulk.rd")

# png("output/figures/fig6_bulk_analysis/fig_S6a_revised.png",height=2,width=2,units="in",pointsize=6,res=300)
# plot(log10(rowMeans(DE_total[,c("freq_fg","freq_bg")])),DE_total$log2_FC,
#      xlab="Log10 Expression",
#      ylab="Log2FC",
#      bty="n")
# sig_up_mask <- DE_total$adj.p.value<1e-3 & DE_total$log2_FC>1 & DE_total$freq_fg>1e-6
# sig_down_mask <- DE_total$adj.p.value<1e-3 & DE_total$log2_FC< -1 & DE_total$freq_bg > 1e-6
# points(log10(rowMeans(DE_total[sig_up_mask,c("freq_fg","freq_bg")])),DE_total$log2_FC[sig_up_mask],col="purple")
# points(log10(rowMeans(DE_total[sig_down_mask,c("freq_fg","freq_bg")])),DE_total$log2_FC[sig_down_mask],col="green")
# dev.off()

# png("output/figures/fig6_bulk_analysis/fig_S6a_volcanoplot_revised.png",height=2,width=2,units="in",pointsize=6,res=300)
# plot(DE_total$log2_FC,-log10(1e-4+DE_total$adj.p.value),xlab="Log2FC LCAMhi vs. LCAMlo",ylab="-Log10(Padj)",
#      pch=".")
# dev.off()


fg_genes <- rownames(DE_total)[DE_total$log2_FC > 1 & DE_total$freq_fg > 1e-6 & DE_total$adj.p.value < 1e-3]
bg_genes <- rownames(DE_total)[DE_total$log2_FC < -1 & DE_total$freq_bg > 1e-6 & DE_total$adj.p.value < 1e-3]

#to exclude genes that are not related to epithelial
fg_genes <- fg_genes[immune_ep_l2fc[fg_genes,]$l2fc > -1]
bg_genes <- bg_genes[immune_ep_l2fc[bg_genes,]$l2fc > -1]


module1_types <- unlist(lapply(deg_pairs,function(x){x[[1]]}))
module2_types <- unlist(lapply(deg_pairs,function(x){x[[2]]}))

exprs_by_celltype <- matrix(NA,nrow=nrow(lung_ldm$dataset$umitab),ncol=length(c(module1_types,module2_types)),
                            dimnames=list(rownames(lung_ldm$dataset$umitab),c(module1_types,module2_types)))
message("Identifying unique genes to each LCAMhi/lo cell type")
for(col in colnames(exprs_by_celltype)){
  message(col)
  mask <- sub_lineage==col
  rs <- rowSums(lung_ldm$dataset$umitab[,mask])
  exprs_by_celltype[,col] <- rs/sum(rs)
}


m1 <- rowMaxs(exprs_by_celltype[,module1_types])
m2 <- rowMaxs(exprs_by_celltype[,module2_types])
m_l2fc <- log2((m1+1e-7)/(m2+1e-7))
names(m_l2fc) <- rownames(exprs_by_celltype)

###  for plotting sup 5B

mat <- exprs_by_celltype[c(fg_genes,bg_genes),]
mat <- mat[,rev(c("MoMac-II","IgG","T_activated","cDC1","Tcm/naive_II","AZU1_mac","cDC2","AM","B"))]

m_l2fc <- m_l2fc[rownames(mat)]

mat1 <- mat[m_l2fc>3,]
mat1 <- mat1[rev(order(apply(mat1,1,which.max),apply(mat1,1,max))),]

mat2 <- mat[m_l2fc< -3,]
mat2 <- mat2[rev(order(apply(mat2,1,which.max),apply(mat2,1,max))),]

mat3 <- mat[m_l2fc > -3 & m_l2fc < 3,]
mat3 <-  mat3[rev(order(apply(mat3,1,which.max),apply(mat3,1,max))),]

new_mat <- rbind(mat1,mat3,mat2)
new_mat <- log2((1e-6+new_mat)/(1e-6+rowMeans(new_mat)))

thresh <- 3
new_mat[new_mat > thresh] <- thresh
new_mat[new_mat < -thresh] <- -thresh

h=seq(-.5/9,9.5/9,10/9/9)[7]
v <- cumsum(c(nrow(mat1),nrow(mat3)))/nrow(mat)

# png("output/figures/fig6_bulk_analysis/figS6B_revised.png",height=2.27,width=4.55,units="in",res=300,pointsize=6,bg="transparent")
# par(oma=c(5,5,5,5),mar=c(2,2,2,2))
# layout(matrix(1:2,ncol=2),widths=c(20,3))
# 
# image(new_mat,xaxt="n",yaxt="n",col=bluered(50))
# mtext(colnames(mat),side=2,at=seq(0,1,1/(ncol(mat)-1)),las=2,line=0.25)
# abline(h=h,lwd=2,col="green")
# abline(v=v,lwd=2,lty=2)
# mtext("Celltype expression",cex=2,line=1)
# mtext("Differentially expressed genes",line=0.5,side=1)
# box()
# 
# par(pin=c(0.125,1))
# image(t(1:100),col=bluered(100),xaxt="n",yaxt="n")
# mtext(paste(c("<",">"),thresh),at=c(0,1),side=4,las=2,line=0.5)
# mtext("Norm. expression",line=1)
# box()
# dev.off()
####
m_l2fc <- log2((m1+1e-7)/(m2+1e-7))
names(m_l2fc) <- rownames(exprs_by_celltype)

fg_genes <- fg_genes[m_l2fc[fg_genes]>3]
bg_genes <- bg_genes[m_l2fc[bg_genes]< -3]



n_genes_per <- 10


blacklist_genes <- c(lung_ldm$model$params$genes_excluded,lung_ldm$model$params$genes_excluded_from_seeding,
                     lung_ldm$model$params$insilico_gating$EP$genes,
                     grep("^HSP",rownames(lung_ldm$dataset$umitab),v=T),
                     grep("^SFT",rownames(lung_ldm$dataset$umitab),v=T),
                     grep("^SCG",rownames(lung_ldm$dataset$umitab),v=T),
                     grep("^IGL",rownames(lung_ldm$dataset$umitab),v=T),
                     grep("^IGHV",rownames(lung_ldm$dataset$umitab),v=T),
                     grep("^HIST",rownames(lung_ldm$dataset$umitab),v=T),
                     grep("^IGHJ",rownames(lung_ldm$dataset$umitab),v=T),
                     grep("^IGK",rownames(lung_ldm$dataset$umitab),v=T))

fg_genes <- setdiff(fg_genes,blacklist_genes)
bg_genes <- setdiff(bg_genes,blacklist_genes)

max_col <- apply(exprs_by_celltype[c(fg_genes,bg_genes),],1,which.max)

gene_type_list <- list()
for(iter in 1:ncol(exprs_by_celltype)){
  genes <- names(max_col)[max_col==iter]
  fc <- abs(DE_total[genes,]$log2_FC)
  genes <- genes[order(fc,decreasing=T)[1:pmin(n_genes_per,length(genes))]]
  gene_type_list[[colnames(exprs_by_celltype)[iter]]] <- genes
}

genes1 <- unlist(gene_type_list[module1_types])
genes2 <- unlist(gene_type_list[module2_types])

genes1 <- genes1[!is.na(genes1)]
genes2 <- genes2[!is.na(genes2)]

#save(genes1,genes2,max_col,exprs_by_celltype,file="intermediates/bulk_genes_revised.rd")
#load("intermediates/bulk_genes_revised.rd")

score_mat <- matrix(0,nrow=length(max_col),ncol=ncol(exprs_by_celltype),dimnames=list(names(max_col),colnames(exprs_by_celltype)))
for(col in seq(ncol(score_mat))){
  score_mat[names(max_col)[which(max_col==col)],col] <- 1
}
score_mat <- score_mat[c(genes1,genes2),]

score_mat <- score_mat[,colSums(score_mat)>0]

###########
# start TCGA analysis
#############


if(!file.exists("data/LUAD_exprs.rd")){
  message("Downloading TCGA LUAD RNA expression data")
  data_url <- "https://www.dropbox.com/s/xxktkx1jduid1n6/LUAD_exprs.rd?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = file.path("data/LUAD_exprs.rd"),mode="wb")
  }else{
    download.file(url = data_url,destfile = file.path("data/LUAD_exprs.rd"))
  }
}

load(file.path(tcga_dir,"LUAD_exprs.rd"))

tcga_tissue <- factor(substr(colnames(expr_mat),14,15))

expr_mat <- expr_mat[,tcga_tissue=="01"]

colnames(expr_mat) <- substr(colnames(expr_mat),1,15)
expr_mat <- expr_mat[,unique(colnames(expr_mat))]


expr_mat <- t(t(expr_mat)/rowSums(t(expr_mat)))

#GEP_score <- colMeans(expr_mat[GEP,])


z_mat <- t(scale(t(log10(expr_mat+1e-6))))


genes <- rownames(score_mat)
genes <- intersect(genes,rownames(z_mat))
genes <- genes[!is.nan(rowSums(z_mat[genes,]))]
score_mat <- score_mat[genes,]


score <- function(zmat,score_mat,module_types){
  zmat <- zmat[rownames(score_mat),]
  pat_scores <- t(scale(t(t(score_mat)%*%zmat)))
  score1 <- colMeans(pat_scores[rownames(pat_scores)%in%module_types,,drop=F])
  return(score1)
}

ep_scores <- read.table("input_tables/TCGA_ESTIMATE.txt",r=1,h=1,stringsAsFactors = F)

z_mat <- z_mat[,intersect(colnames(z_mat),rownames(ep_scores))]
ep_score <- ep_scores[colnames(z_mat),"Immune_score"]


 
 
ep_group <- cut(ep_score,quantile(ep_score,seq(0,1,.1)))



tumor_split <- split(colnames(z_mat),ep_group)

score1 <- score(z_mat,score_mat,module1_types)
score2 <- score(z_mat,score_mat,module2_types)

save(score1,score2,z_mat,file = "output/statistics/TCGA_scores.rd")

# layout(matrix(1:12,nrow=3,ncol=4,byrow=T))
# for(iter in 1:10){
#   plot(score1[tumor_split[[iter]]],
#        score2[tumor_split[[iter]]],main=paste("immune decile ",iter),
#        xlim=range(score1),ylim=range(score2))
#   grid()
# }

png("output/figures/figure_s6abc.png",height=1.35,width=3*1.28,pointsize=8,units="in",res=300,bg="transparent")
layout(matrix(1:3,nrow=1))
par(mgp=c(2,1,0))
plot(ep_score,score1,col=alpha(1,0.5),xlab="Immune ESTIMATE",ylab="LCAMhi score",cex=0.5,pch=16)
#mtext(paste("cor=",round(cor(ep_score,score1),2),sep=""))
plot(ep_score,score2,col=alpha(1,0.5),xlab="Immune ESTIMATE",ylab="LCAMlo score",cex=0.5,pch=16)
#mtext(paste("cor=",round(cor(ep_score,score2),2),sep=""))
plot(ep_score,score1-score2,pch=16,col=alpha(1,0.5),xlab="Immune ESTIMATE",ylab="LCAMhi-LCAMlo score",cex=0.5)
#mtext(paste("cor=",round(cor(ep_score,score1-score2),2),sep=""))

dev.off()


#tumor_split <- lapply(tumor_split,function(x){t(scale(t(log10(expr_mat[,x]+1e-6))))})
split_cor <- lapply(tumor_split,function(x){
  return(cor.test(score1[x],score2[x],method="spearman"))})
#split_cor <- unlist(split_cor)

cor_list <- unlist(lapply(split_cor,function(x){x$estimate}))
lower_bound <- unlist(lapply(split_cor,function(x){x$conf.int[[1]]}))
upper_bound <- unlist(lapply(split_cor,function(x){x$conf.int[[2]]}))


#with spearman correlation:

sig <- 1.96/sqrt(unlist(lapply(tumor_split,length))-3)
lower_bound <- tanh(atanh(unlist(lapply(split_cor,function(x){x$estimate})))-sig)
upper_bound <- tanh(atanh(unlist(lapply(split_cor,function(x){x$estimate})))+sig)


png("output/figures/figure_s6d.png",height=1.4,width=1.7,units="in",res=300,pointsize=6,bg="transparent")
par(mgp=c(2,1,0))
plot(1:10,cor_list,cex=0,ylim=c(-1,1)*max(abs(c(lower_bound,upper_bound))),xlim=c(0.5,10.5),xaxt="n",xlab="",
     ylab="Spearman rho\n(LCAMhi score, LCAMlo score)")
segments(x0=1:10-0.15,x1=1:10+0.15,y0=cor_list,y1=cor_list,lwd=3)
segments(x0=1:10,x1=1:10,y0=lower_bound,y1=upper_bound)
abline(h=0,col="red",lty=2)
mtext(side=1,at=1:10,1:10)
mtext(side=1,line=2,"Decile of immune infiltrate")
mtext("LCAM-hi/lo summary scores")
dev.off()

png("output/figures/figure_s6e.png",bg="transparent",
    height=2,width=2,units="in",res=300,pointsize=6)
par(mgp=c(2,1,0))
plot(1:10,xlim=c(-2.5,2),ylim=c(-2,3),col="white",xlab="LCAMhi score",ylab="LCAMlo score")
#layout(matrix(1:3,nrow=1))
for(decile in c(1:10)){
  # plot(score1[tumor_split[[decile]]],score2[tumor_split[[decile]]],xlim=range(score1),ylim=c(-2,2.5),pch=16,col=alpha(1,0.5),
  #      xlab="score1",ylab="score2")
  # abline(lm(score2[tumor_split[[decile]]]~score1[tumor_split[[decile]]]))
  # mtext(paste("decile",decile),line=0.25)
  # grid()
  
  # plotting PCA line instead of linear regression line
  x <- score1[tumor_split[[decile]]]
  y <- score2[tumor_split[[decile]]]
  xyNorm <- cbind(x=x-mean(x), y=y-mean(y))
  #plot(xyNorm)
  
  #covariance
  xyCov <- cov(xyNorm)
  eigenValues <- eigen(xyCov)$values
  eigenVectors <- eigen(xyCov)$vectors
  eigenValues
  eigenVectors
  
  #plot(x,y,xlim=range(score1),ylim=c(range(score2)[1],2),pch=16,col=alpha(1,0.5))
  if(decile%in%c(1,3,10)){
  points(x,y,xlim=range(score1),ylim=c(range(score2)[1],2),pch=16,col=alpha(decile,0.5))
  lines(xyNorm[,1]+mean(x), eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y),col=decile)
  text(min(xyNorm[,1])+mean(x),max(eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y)),decile,col=decile,font=2)
  
  #abline( lm((eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y))~(xyNorm[,1]+mean(x))))
  #text(x=min(xyNorm[,1]+mean(x)),max(eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y)),decile,col=decile)
  
  #}else{
 # lines(xyNorm[,1]+mean(x), eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y),col=alpha("grey",0.5))
#   text(x=min(xyNorm[,1]+mean(x)),max(eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y)),decile,col="grey")
  }else if(decile!=9){
    lines(xyNorm[,1]+mean(x), eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y),col=alpha(1,0.3))
    text(min(xyNorm[,1])+mean(x),max(eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y)),decile,col=alpha(1,0.3))
  }else{
    lines(xyNorm[,1]+mean(x), eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y),col=alpha(1,0.3))
    text(max(xyNorm[,1])+mean(x),min(eigenVectors[2,1]/eigenVectors[1,1] * xyNorm[,1]+mean(y)),decile,col=alpha(1,0.3))
  
}

}
dev.off()

#score_var <- unlist(lapply(tumor_split,function(x){var(score1[x]-score2[x])}))

#png("output/figures/fig6_bulk_analysis/fig_sup6h_scorevar_by_decile_revise.png",height=1.4,width=1.7,units="in",res=300,pointsize=6)
#barplot(score_var,names.arg=1:10,col="black",xlab="Decile of immune infiltrate",ylab="var[LCAMhi-LCAMlo score]")
#cor.test(1:10,score_var,method="spearman")
#dev.off()





sliced <- z_mat[,ep_scores[colnames(z_mat),"Immune_score"]>quantile(ep_scores[colnames(z_mat),"Immune_score"],0.1)]
score1 <- score(sliced,score_mat,module1_types)
score2 <- score(sliced,score_mat,module2_types)
#plot(score1,score2,pch=16,col=alpha(1,0.5))



ord <- order(score1-score2)
#sliced <- z_mat
mat <- sliced[rownames(score_mat),ord]

thresh <- 2
mat[mat > thresh] <- thresh
mat[mat < -thresh] <- -thresh

v <- c(-0.5:(nrow(mat)+0.5)/c(nrow(mat)-1))[cumsum(colSums(score_mat))[-ncol(score_mat)]+1]

# png("output/figures/fig6_bulk_analysis/fig6a_revise.png",height=2.75,width=4.7,units="in",res=300,pointsize = 5)
# layout(matrix(1:2,nrow=1),widths=c(20,3))
# 
# image(mat,col=bluered(50),xaxt="n",yaxt="n")
# mtext(side=1,rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),las=2,line=0.25)
# abline(v=v)
# abline(v=v[3],lwd=3)
# box()
# mtext("TCGA: LUAD",cex=2,line=1)
# 
# par(pin=c(0.125,1))
# image(t(1:100),col=bluered(100),xaxt="n",yaxt="n")
# box()
# mtext(paste(c("< -","> "),thresh,sep=""),at=c(0,1),las=2,line=0.25,side=4)
# mtext("Z-score",line=0.5)
# 
# dev.off()

png("output/figures/figure_6a.png",height=3.6,width=5.66,units="in",res=300,pointsize = 5)
layout(matrix(1:4,nrow=2),widths=c(20,3),heights=c(1.5,20))
par(mar=c(0,3,2,2))
image(as.matrix(rep(c(1,2),times=c(length(genes1),length(genes2)))),col=c(rgb(t(col2rgb("purple")*.7/255)),rgb(0,.7,0)),xaxt="n",yaxt="n")
abline(v=v[3],lwd=3)
box()

par(mar=c(10,3,0,2))

image(mat,col=bluered(50),xaxt="n",yaxt="n")
mtext(side=1,rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),las=2,line=0.25,font=3)
abline(v=v)
abline(v=v[3],lwd=3)
box()
# 
# plot.new()
# par(pin=c(0.125,1))
# 
# image(t(1:100),col=bluered(100),xaxt="n",yaxt="n")
# box()
# mtext(paste(c("< -","> "),thresh,sep=""),at=c(0,1),las=2,line=0.25,side=4)
# mtext("Z-score",line=0.5)

dev.off()

## mutation analyses

if(!file.exists("data/TCGA_LUAD_mutect_somatic.maf.gz")){
  message("Downloading TCGA LUAD mutation data")
  data_url <- "https://www.dropbox.com/s/dg75ne7j3t1382p/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf.gz?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = "data/TCGA_LUAD_mutect_somatic.maf.gz",mode="wb")
  }else{
    download.file(url = data_url,destfile = "data/TCGA_LUAD_mutect_somatic.maf.gz")
  }
}


mut_data <- fread("data/TCGA_LUAD_mutect_somatic.maf.gz")
mut_data <- mut_data[mut_data$Variant_Classification%in%c("Missense_Mutation","Nonsense_Mutation")]
mut_tab <- table(mut_data$Hugo_Symbol,mut_data$Tumor_Sample_Barcode)
colnames(mut_tab) <- substr(colnames(mut_tab),1,15)
n_mut <- colSums(mut_tab)
pats <- intersect(colnames(sliced),colnames(mut_tab))

res <- cor.test(score1[pats]-score2[pats],log10(1+n_mut[pats]))


GEP <- "CCL5, CD27, CD274, CD276, CD8A, CMKLR1, CXCL9, CXCR6, HLA-DQA1, HLA-DRB1, HLA-E, IDO1, LAG3, NKG7, PDCD1LG2, PSMB10, STAT1, TIGIT"
GEP <- strsplit(GEP,", ")[[1]]
GEP_score <- colMeans(z_mat[GEP,])

png("output/figures/figure_s6lm.png",height=4,width=8,units="in",res=300,bg="transparent")
layout(matrix(1:2,nrow=1))
par(mgp=c(2,1,0))
plot(log10((1+n_mut[pats])/48.2),ep_scores[pats,"Immune_score"],xlab=expression(paste("Log"["10"],"[TMB/Mb]",sep="")),ylab="Immune ESTIMATE",
     pch=16,col=alpha(1,0.5))
res <- cor.test(log10(1+n_mut[pats]),ep_scores[pats,"Immune_score"])
abline(lm(ep_scores[pats,"Immune_score"]~log10(1+n_mut[pats])))
mtext(paste("Immune content: rho=",round(res$estimate,2)),cex=1.5)

plot(log10((1+n_mut[pats])/48.2),GEP_score[pats],xlab=expression(paste("Log"["10"],"[TMB/Mb]",sep="")),ylab="GEP score",pch=16,col=alpha(1,0.5))
res <- cor.test(log10(1+n_mut[pats]),GEP_score[pats])
abline(lm(GEP_score[pats]~log10(1+n_mut[pats])))
mtext(paste("GEP score: rho=",round(res$estimate,2)),cex=1.5)
dev.off()


# do mutational signatures analysis here:
########################################
source("scripts/mutation_signatures.R")






mut_sig <- fread("input_tables/mSignatureDB_profiles.txt")

mut_sig_mat <- pivot_wider(data=mut_sig,id_cols = c("Signature","Tumor_Sample_Barcode","Contribution"),names_from="Tumor_Sample_Barcode",values_from="Contribution")
mut_sig_mat <- as.matrix(mut_sig_mat)
rownames(mut_sig_mat) <- mut_sig_mat[,1]
mut_sig_mat <- mut_sig_mat[,-1]
mut_sig_mat <- matrix(as.numeric(mut_sig_mat),nrow=nrow(mut_sig_mat),dimnames=dimnames(mut_sig_mat))

colnames(mut_sig_mat) <- sub("\\.","-",colnames(mut_sig_mat))
colnames(mut_sig_mat) <- sub("\\.","-",colnames(mut_sig_mat))
colnames(mut_sig_mat) <- paste(colnames(mut_sig_mat),"-01",sep="")

sig_pats <- intersect(names(score1),colnames(mut_sig_mat))

mut_sig_mat <- mut_sig_mat[,sig_pats]

mut_sig_mat <- t(t(mut_sig_mat)/rowSums(t(mut_sig_mat)))

mut_sig_mat <- log10(1e-2+mut_sig_mat)

#set.seed(910430)
# set.seed(043091)
# k <- kmeans(t(mut_sig_mat),9)
# 
# clust_ord <- hclust(as.dist(1-cor(t(k$centers))))$order
# sig_ord <- hclust(as.dist(1-cor(t(mut_sig_mat),method="spearman")))$order
# pat_order <- order(factor(k$cluster,clust_ord))
# mat <- mut_sig_mat[sig_ord,pat_order]
# 
# rownames(mat) <- unlist(lapply(strsplit(rownames(mat),"\\."),function(x){x[[2]]}))
# 
# bet_clusts <- which(diff(k$cluster[pat_order])!=0)/length(pat_order)
# clust_name_x <- (c(0,bet_clusts)+c(bet_clusts,1))/2
# 
# cols <- c(2:8,rgb(t(col2rgb(1:8)*0.7)/255))
# 
# png(file.path(output_dir,"fig_6d_mutational_sig_heatmap.png"),height=4,width=4,units="in",res=1000,pointsize=8)
# image(t(mat),col=viridis(100),xaxt="n",yaxt="n")
# abline(v=bet_clusts,col="white")
# mtext(at=clust_name_x,clust_ord,line=0.25,col=cols[clust_ord],font=2,cex=1.5)
# mtext(rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),side=2,las=2,line=0.15)
# dev.off()

# layout(matrix(1:4,nrow=2))
# for(iter in 1:max(clust_ord)){
#   mask <- names(k$cluster)[k$cluster==iter]
#   plot(log10(1+n_mut[mask]),score1[mask]-score2[mask],xlim=range(log10(1+n_mut)),ylim=range(score1-score2))
#   grid()
#   res <- cor.test(log10(1+n_mut[mask]),score1[mask]-score2[mask])
#   mtext(paste("clust ",iter,": r=",signif(res$estimate,2),"; p=",signif(res$p.value,2),sep=""),line=0.25)
#   abline(lm((score1[mask]-score2[mask])~log10(1+n_mut[mask])))
# }
# png("output/figures/fig6_bulk_analysis/fig_6e_mutational_sig_group_scatters.png",height=4,width=4,units="in",res=1000,pointsize=12)
# 
# plot(log10((1+n_mut[sig_pats])/48.2),score1[sig_pats]-score2[sig_pats],col=alpha(cols[k$cluster[sig_pats]],0.5),pch=16,
#      xlab=expression(paste("Log"["10"],"[TMB/Mb]",sep="")),
#      ylab="LCAMhi-LCAMlo score",bty="n")
# for(iter in 1:max(clust_ord)){
#   message(iter)
#   mask <- names(k$cluster)[k$cluster==iter]
#   abline(lm((score1[mask]-score2[mask])~log10((1+n_mut[mask])/48.2)),col=cols[iter],lwd=3)
#   print(length(mask))
#   print(cor.test((score1[mask]-score2[mask]),log10(1+n_mut[mask]))$estimate)
#   print(cor.test((score1[mask]-score2[mask]),log10(1+n_mut[mask]))$p.value)
# }
# dev.off()

mask <- mut_sig_mat["Signature.4",]==-2
#mask <- colnames(mut_sig_mat)[mask]
#mask <- names(k$cluster)[mask]

mask1 <- mut_sig_mat["Signature.4",]== -2

png("output/figures/figure_6c.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
plot(log10((1+n_mut[sig_pats])/48.2),score1[sig_pats]-score2[sig_pats],col=alpha(1+as.integer(mask),0.8),pch=16,
     xlab=expression(paste("Log"["10"],"[TMB/Mb]",sep="")),
     ylab="LCAM score",cex=1)
abline(lm((score1[sig_pats[mask]]-score2[sig_pats[mask]])~log10((1+n_mut[sig_pats[mask]])/48.2)),lwd=3,col="red")
abline(lm((score1[sig_pats[!mask]]-score2[sig_pats[!mask]])~log10((1+n_mut[sig_pats[!mask]])/48.2)),lwd=3,col="black")

cor.test(log10(1+n_mut[sig_pats[mask]]),score1[sig_pats[mask]]-score2[sig_pats[mask]])
cor.test(log10(1+n_mut[sig_pats[!mask]]),score1[sig_pats[!mask]]-score2[sig_pats[!mask]])
mtext("Smoking signature present",cex=1.5)
#cor.test(log10(1+n_mut[sig_pats]),score1[sig_pats]-score2[sig_pats])
dev.off()



##########################################


clin <- fread("input_tables/TCGA-LUAD_clinical.csv",sep=",",stringsAsFactors = F)
clin <- as.matrix(clin)
rownames(clin) <- clin[,1]; clin <- clin[,-1]
clin <- clin[!grepl("_1",rownames(clin)),]
rownames(clin) <- paste(rownames(clin),"-01",sep="")
clin <- clin[names(score1),]

#feature <- "ajcc_pathologic_t"
#png("output/figures/fig6_bulk_analysis/stage_vs_LCAM_revise.png",height=4,width=6,units="in",res=300)
#par(mar=c(10,5,2,2))
#boxplot((score1-score2)~clin[,feature],las=2,range=0,ylab="LCAMhi-LCAMlo score")
#points(jitter(as.integer(factor(clin[,feature]))),score1-score2,pch=16,col=alpha(1,0.5),cex=0.5)
#summary(aov((score1-score2)~as.factor(clin[,feature])))
#dev.off()

# feature <- "site_of_resection_or_biopsy" 
# png("output/figures/fig6_bulk_analysis/site_vs_LCAM_revise.png",height=4,width=6,units="in",res=300)
# par(mar=c(10,5,2,2))
# boxplot((score1-score2)~clin[,feature],las=2,range=0,ylab="LCAMhi-LCAMlo score")
# points(jitter(as.integer(factor(clin[,feature]))),score1-score2,pch=16,col=alpha(1,0.5),cex=0.5)
# summary(aov((score1-score2)~as.factor(clin[,feature])))
# dev.off()

# plot(as.numeric(clin[,"pack_years_smoked"]),score1-score2,col=alpha(1,0.4),pch=16)
# cor.test(score1-score2,log10(1e-2+as.numeric(clin[,"pack_years_smoked"])),method="spearman")
# abline(lm((score1-score2)~as.numeric(clin[,"age_at_index"])))
# cor.test(as.numeric(clin[,"age_at_index"]),score1-score2,method="spearman")
# text(paste("rho="))
# cor.test(as.numeric(clin[,"age_at_index"]),as.numeric(clin[,"pack_years_smoked"]),method="spearman")
# 
# cor.test(as.numeric(clin[,"age_at_index"]),log2(1+n_mut))

# 
# mut_type_tab <- table(mut_data$Variant_Classification,mut_data$Tumor_Sample_Barcode)
# colnames(mut_type_tab) <- substr(colnames(mut_type_tab),1,15)
# 
# mut_type_cor <- cor(log10(1e-4+t(mut_type_tab)),method="pearson")
# type_ord <- hclust(as.dist(1-mut_type_cor))$order
# pimage(mut_type_cor[type_ord,type_ord],axes="both")

# plot(sort(rowSums(mut_type_tab),d=T),ylab="# of mutations",xaxt="n",log="y",pch=16,xlab="")
# lines(sort(rowSums(mut_type_tab),d=T),lty=2)
# mtext(side=1,las=2,at=seq(nrow(mut_type_tab)),names(sort(rowSums(mut_type_tab),d=T)),line=0.25)
# mtext("Frequency of mutation types")

# fit <- lm((score1[pats]-score2[pats])~log10(1+colSums((mut_type_tab[!rownames(mut_type_tab)%in%c("Silent"),pats]))))
# cor.test(log10(1+colSums(mut_type_tab[c("Silent"),pats,drop=F])),fit$residuals,method="pearson")
# 
# fit <- lm((score1[pats]-score2[pats])~log10(1+colSums((mut_type_tab[!rownames(mut_type_tab)%in%c("Frame_Shift_Del","Frame_Shift_Ins"),pats]))))
# cor.test(log10(1+colSums(mut_type_tab[c("Frame_Shift_Del","Frame_Shift_Ins"),pats])),fit$residuals,method="pearson")
# 
# fit <- lm((score1[pats]-score2[pats])~log10(1+colSums((mut_type_tab[!rownames(mut_type_tab)%in%c("Missense_Mutation"),pats]))))
# cor.test(log10(1+colSums(mut_type_tab[c("Missense_Mutation"),pats,drop=F])),fit$residuals,method="pearson")

# mask <- mut_tab["TP53",pats]>0
# #mask <- mut_tab["PMAIP1",pats]
# points(log10(1+n_mut[pats[mask]]),score1[pats[mask]]-score2[pats[mask]],pch=16,col="red",cex=0.5)
# 
# 
# res <- cor.test(score1_unadj[pats]-score2_unadj[pats],log10(1+n_mut[pats]),method="spearman")
# plot(score1_unadj[pats]-score2_unadj[pats],log10(1+n_mut[pats]),
#      xlab="Score1_unadj - Score2_unadj",
#      ylab="Log10(1+#Mutations)")
# abline(lm(log10(1+n_mut[pats])~score1_unadj[pats]-score2_unadj[pats]))
# text(x=1,y=0.75,paste("rho=",round(res$estimate,2),"; p=1e",round(log10(res$p.value)),sep=""))



# single gene correlation

gene_cor_with_mut <- cor(t(z_mat[rownames(score_mat),pats]),n_mut[pats],method="spearman",use="p")[,1]
png("output/figures/figure_s6n.png",height=1.14,width=4.34,units="in",res=300,pointsize=4)
b <- barplot(gene_cor_with_mut,las=2,xaxs="i",main="Correlation with TMB",col=c("purple","green")[1+as.integer(names(gene_cor_with_mut)%in%genes2)])
abline(v=b[cumsum(colSums(score_mat))[-ncol(score_mat)]]+0.5)
#abline(v=b[cumsum(colSums(score_mat))[3]]+0.5,col="green",lwd=2)
dev.off()


# single gene effects

#mask <- mut_tab["TP53",pats]>0

mut_data <- fread("data/TCGA_LUAD_mutect_somatic.maf.gz")

mut_data <- mut_data[!mut_data$Variant_Classification%in%c("RNA","Silent"),]
#mut_data <- mut_data[mut_data$Variant_Classification%in%c("Missense_Mutation","Nonsense_Mutation")]
mut_tab <- table(mut_data$Hugo_Symbol,mut_data$Tumor_Sample_Barcode)
colnames(mut_tab) <- substr(colnames(mut_tab),1,15)

class <- array(NA,ncol(mut_tab),dimnames=list(colnames(mut_tab)))
class[colSums(mut_tab[c("TP53","KRAS","STK11","EGFR"),])==0] <- 1
class[colSums(mut_tab[c("TP53","KRAS","STK11"),])==0 & mut_tab["EGFR",]>0] <- 2
class[colSums(mut_tab[c("TP53","KRAS","EGFR"),])==0 & mut_tab["STK11",]>0] <- 3
class[colSums(mut_tab[c("TP53","STK11","EGFR"),])==0 & mut_tab["KRAS",]>0] <- 4
class[colSums(mut_tab[c("TP53","EGFR"),])==0 & mut_tab["KRAS",]>0 & mut_tab["STK11",]>0] <- 5
class[colSums(mut_tab[c("KRAS","STK11"),])==0 & mut_tab["TP53",]>0 & mut_tab["EGFR",]>0] <- 6
class[colSums(mut_tab[c("KRAS","STK11","EGFR"),])==0 & mut_tab["TP53",]>0] <- 7
class[colSums(mut_tab[c("KRAS","EGFR"),])==0 & mut_tab["TP53",]>0 & mut_tab["STK11",]>0] <- 8
class[colSums(mut_tab[c("STK11","EGFR"),])==0 & mut_tab["TP53",]>0 & mut_tab["KRAS",]>0] <- 9
class[colSums(mut_tab[c("EGFR"),,drop=F])==0 & mut_tab["TP53",]>0 & mut_tab["KRAS",]>0 & mut_tab["STK11",]>0] <- 10

class_mat <- matrix(c(0,0,0,0,
                      0,0,0,1,
                      0,0,1,0,
                      0,1,0,0,
                      0,1,1,0,
                      1,0,0,1,
                      1,0,0,0,
                      1,1,0,0,
                      1,0,1,0,
                      1,1,1,0),nrow=4,dimnames=list(c("TP53","KRAS","STK11","EGFR"),NULL))
class_mat <- class_mat[4:1,]
nclass <- ncol(class_mat)

png("output/figures/figure_7a.png",height=1.43,width=2,units="in",res=300,pointsize=6)
layout(matrix(1:2,nrow=2),heights=c(2,1))
par(mar=c(0,4,2,2))
b <- boxplot((score1[pats]-score2[pats])~class[pats],range=0,xaxt="n",xaxs="i",xlim=c(0.75,10.25),
             ylab="LCAM score",axes=F)
box()
axis(side=2,at=c(-2,-0,2),las=2)
points(jitter(class[pats]),score1[pats]-score2[pats],pch=16,col=alpha(1,0.5),cex=0.5,
       ylab="LCAMhi-LCAMlo score")
par(mar=c(1,4,0.5,2))
image(t(class_mat),col=c(0,1),xaxt="n",yaxt="n"); box()
mtext(rownames(class_mat),side=2,las=2,at=seq(0,1,1/(nrow(class_mat)-1)),font=3,las=2,line=0.25)
abline(v=c(seq(-1/2/(nclass-1),1+1/2/(nclass-1),(1+1/(nclass-1))/nclass)))
abline(h=c(seq(-1/2/(4-1),1+1/2/(4-1),(1+1/(4-1))/4)))
dev.off()


png("output/figures/figure_7b.png",height=1.43,width=2,units="in",res=300,pointsize=6)
layout(matrix(1:2,nrow=2),heights=c(2,1))
par(mar=c(0,4,2,2))
b <- boxplot(n_mut[pats]/48.2~class[pats],log="y",range=0,xaxt="n",xaxs="i",xlim=c(0.75,10.25),
             ylab="TMB/Mb",axes="F")
box()
axis(side=2,at=c(0.05,1,20),labels=F)
mtext(side=2,at=c(0.05,1,20),c("0.05","1","20"),las=2,line=0.75)
points(jitter(class[pats]),n_mut[pats]/48.2,pch=16,col=alpha(1,0.5),cex=0.5,
       ylab="TMB/Mb")
par(mar=c(1,4,0.5,2))
image(t(class_mat),col=c(0,1),xaxt="n",yaxt="n"); box()
mtext(rownames(class_mat),side=2,las=2,at=seq(0,1,1/(nrow(class_mat)-1)),font=3,las=2,line=0.25)
abline(v=c(seq(-1/2/(nclass-1),1+1/2/(nclass-1),(1+1/(nclass-1))/nclass)))
abline(h=c(seq(-1/2/(4-1),1+1/2/(4-1),(1+1/(4-1))/4)))
dev.off()


fit <- lm((score1[pats]-score2[pats])~1+log10((n_mut[pats]+1)/48.2))

t.test(fit$residuals~mut_tab["TP53",pats]>0)
t.test(fit$residuals~mut_tab["KRAS",pats]>0)
t.test(fit$residuals~mut_tab["STK11",pats]>0)
t.test(fit$residuals~mut_tab["EGFR",pats]>0)

genes <- c("TP53","KRAS","STK11","EGFR")

png("output/figures/figure_7cd.png",height=2.59,width=3.58,units="in",res=300,pointsize=6)
layout(matrix(1:4,nrow=2,byrow=T))
par(mgp=c(2,1,0),mar=c(3,4,2,1))
xlim <- range(fit$residuals)
for(gene_iter in genes){
  res <- t.test(fit$residuals[pats]~mut_tab[gene_iter,pats]>0)
  d_wt <- density(fit$residuals[pats[mut_tab[gene_iter,pats]==0]])
  d_mut <- density(fit$residuals[pats[mut_tab[gene_iter,pats]>0]])
  plot(d_wt,ylim=c(0,.7),xlab="Residuals",main=NA,xlim=xlim)
  mtext(gene_iter,line=0.25,cex=1.3)
  lines(d_mut,col="red")
  legend("topleft",legend=c(paste("WT (N=",sum(mut_tab[gene_iter,pats]==0),")",sep=""),
                            paste("mut (N=",sum(mut_tab[gene_iter,pats]>0),")",sep="")),lty=1,col=c(1,2))
  text(x=1.5,y=0.6,paste("p=",signif(res$p.value,2),sep=""))
}
dev.off()

# png("output/figures/fig6_bulk_analysis/residual_boxplots_revise.png",height=2.59,width=3.58,units="in",res=300,pointsize=6)
# layout(matrix(1:4,nrow=2,byrow=T))
# par(mgp=c(2,1,0),mar=c(3,4,2,1))
# xlim <- range(fit$residuals)
# for(gene_iter in genes){
#   res <- t.test(fit$residuals[pats]~mut_tab[gene_iter,pats]>0)
#   boxplot(fit$residuals[pats]~mut_tab[gene_iter,pats]>0,
#           names=paste(c("WT (","mut ("),
#                       c(sum(mut_tab[gene_iter,pats]==0),sum(mut_tab[gene_iter,pats]>0)),")",sep=""),
#           range=0,ylim=c(min(fit$residuals),max(fit$residuals+1)),ylab="residuals")
#   points(jitter(as.integer(factor(mut_tab[gene_iter,pats]>0))),fit$residuals[pats],pch=16,col=alpha(1,0.4))
#   text(x=1.5,y=3.5,paste("p=",signif(res$p.value,2),sep=""))
#   mtext(gene_iter,cex=1.5)
# }
# dev.off()

#Stage T1a vs. others
# png("output/figures/fig6_bulk_analysis/residual_boxplot_T1a_revise.png",height=4,width=4,units="in",res=300)
# par(mgp=c(2,1,0))
# res <- t.test(fit$residuals~factor(clin[names(fit$residuals),"ajcc_pathologic_t"]!="T1a"))
# xlim <- range(fit$residuals)
# boxplot(fit$residuals~factor(clin[names(fit$residuals),"ajcc_pathologic_t"]!="T1a"),range=0,
#         names=paste(c("T1a (",">T1a ("),
#                       c(sum(clin[names(fit$residuals),"ajcc_pathologic_t"]=="T1a"),sum(clin[names(fit$residuals),"ajcc_pathologic_t"]!="T1a")),
#                     ")",sep=""),
#                     ylim=c(xlim[1],xlim[2]+1),
#         ylab="residuals")
# points(jitter(as.integer(factor(clin[names(fit$residuals),"ajcc_pathologic_t"]!="T1a"))),fit$residuals,
#        pch=16,col=alpha(1,0.4))
# text(x=1.5,y=3.5,paste("p=",signif(res$p.value,3)))
# mtext("T-stage",cex=1.5)
# dev.off()

# png("output/figures/fig6_bulk_analysis/LCAM_boxplot_T1a_revise.png",height=4,width=4,units="in",res=300)
# par(mgp=c(2,1,0))
# res <- t.test((score1[pats]-score2[pats])~factor(clin[names(fit$residuals),"ajcc_pathologic_t"]!="T1a"))
# xlim <- range(fit$residuals)
# boxplot((score1[pats]-score2[pats])~factor(clin[pats,"ajcc_pathologic_t"]!="T1a"),range=0,
#         names=paste(c("T1a (",">T1a ("),
#                     c(sum(clin[names(fit$residuals),"ajcc_pathologic_t"]=="T1a"),sum(clin[names(fit$residuals),"ajcc_pathologic_t"]!="T1a")),
#                     ")",sep=""),
#         ylim=c(xlim[1],xlim[2]+1),
#         ylab="LCAM")
# points(jitter(as.integer(factor(clin[names(fit$residuals),"ajcc_pathologic_t"]!="T1a"))),fit$residuals,
#        pch=16,col=alpha(1,0.4))
# text(x=1.5,y=3.5,paste("p=",signif(res$p.value,3)))
# mtext("T-stage",cex=1.5)
# dev.off()

mask <- clin[pats,"ajcc_pathologic_t"]=="T1a"

Tstage_simp <- clin[,"ajcc_pathologic_t"]
Tstage_simp <- substr(Tstage_simp,1,2)
Tstage_simp[Tstage_simp=="TX"] <- NA
res <- cor.test((score1[pats]-score2[pats]),as.integer(factor(Tstage_simp[pats])),method="spearman")

png("output/figures/figure_s6j.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
boxplot((score1[pats]-score2[pats])~Tstage_simp[pats],range=0,ylab="LCAM score",
        border=1+1:4)
points(jitter(as.integer(factor(Tstage_simp[pats]))),score1[pats]-score2[pats],
       col=alpha(1,0.4),pch=16)
#text(x=3.5,y=-3.25,
#     paste("rho=",signif(res$estimate,2),"\np=",signif(res$p.value,2),sep=""))
mtext("Tstage",cex=1.5)
dev.off()

#t.test((score1[pats]-score2[pats])~Tstage_simp[pats]!="T1")

png("output/figures/figure_s6o.png",height=3,width=12,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
layout(matrix(1:4,nrow=1))
for(Tstage in 1:4){
  mask <- as.integer(factor(Tstage_simp))==Tstage
  plot(log10((1+n_mut[pats])/48.2)[mask],(score1[pats]-score2[pats])[mask],
       col=alpha(1,0.4),lwd=2,ylab="LCAM score",
       xlab=expression(paste("Log"["10"],"[TMB/Mb]",sep="")),
       xlim=range(log10((1+n_mut[pats])/48.2)),ylim=range(score1[pats]-score2[pats]))
  abline(lm((score1[pats]-score2[pats])[mask]~
              log10((1+n_mut[pats])/48.2)[mask]),
         col=rgb(t(col2rgb(Tstage+1)*.7/255)),lwd=2)
  mtext(paste("T",Tstage,sep=""),cex=1.5)
}
dev.off()


thorrsen <- read.csv("input_tables/Thorsson_et_al_TCGA_tumor_immunology_metrics.csv",
                     r=1,h=1,stringsAsFactors = F)

thorrsen <- thorrsen[substr(pats,1,12),]
pats <- pats[rownames(thorrsen)!="NA" & !is.na(thorrsen$CTA.Score)]
thorrsen <- thorrsen[substr(pats,1,12),]


#plot(thorrsen$CTA.Score,score1-score2,ylab="Intratumor heterogeneity",xlab="LCAM score")
#cor.test(log10(1+n_mut[pats]),log10(1+thorrsen$CTA.Score))

fit <- lm((score1[pats]-score2[pats])~log10(1+n_mut[pats]))
#res <- cor.test(log10(1+thorrsen$CTA.Score),fit$residuals[pats])
#cor.test(log10(1+thorrsen$Indel.Neoantigens),fit$residuals)



#plot(log10(1+thorrsen$SNV.Neoantigens),fit$residuals)

#   
# plot(log10(1+thorrsen$CTA.Score),fit$residuals[pats],pch=16,col=alpha(1,0.5),xlab="logCTA score",ylab="residuals")
# abline(lm(fit$residuals[pats]~log10(1+thorrsen$CTA.Score)))

fit <- lm((score1[pats]-score2[pats])~
            log10(1+n_mut[pats]))

# par(mar=c(15,5,2,2))
# barplot(sort(cor(thorrsen[,4:ncol(thorrsen)],fit$residuals,method="spearman",use="p")[,1]),las=2)

png(file.path(output_dir,"figure_s6r.png"),units="in",height=4,width=4,res=300,bg="transparent")
par(mgp=c(2,1,0))
res <- cor.test(log10(1+thorrsen$CTA.Score),fit$residuals[pats])
plot(log10(1+thorrsen$CTA.Score),fit$residuals[pats],
     ylab="residuals",xlab="CTA score",col=alpha(1,0.5),pch=16)
abline(lm(fit$residuals[pats]~log10(1+thorrsen$CTA.Score)))
text(x=0.7,y=-2.5,paste("r=",signif(res$estimate,2),"; p=",signif(res$p.value,2),sep=""))
mtext("Cancer testis antigens",cex=1.5)
dev.off()

png(file.path(output_dir,"figure_s6p.png"),units="in",height=4,width=4,res=300,bg="transparent")
par(mgp=c(2,1,0))
res <- cor.test(log10(1+thorrsen$Indel.Neoantigens),fit$residuals[pats])
plot(log10(1+thorrsen$Indel.Neoantigens),fit$residuals[pats],
     ylab="residuals",xlab="Log#Indel NeoAntigens",col=alpha(1,0.5),pch=16)
abline(lm(fit$residuals[pats]~log10(1+thorrsen$Indel.Neoantigens)))
text(x=0.7,y=-2.5,paste("r=",signif(res$estimate,2),"; p=",signif(res$p.value,2),sep=""))
mtext("Indel Neoantigens",cex=1.5)
dev.off()
png(file.path(output_dir,"figure_s6q.png"),units="in",height=4,width=4,res=300,bg="transparent")
par(mgp=c(2,1,0))
res <- cor.test(log10(1+thorrsen$SNV.Neoantigens),fit$residuals[pats])
plot(log10(1+thorrsen$SNV.Neoantigens),fit$residuals[pats],
     ylab="residuals",xlab="Log#SNV NeoAntigens",col=alpha(1,0.5),pch=16)
abline(lm(fit$residuals[pats]~log10(1+thorrsen$SNV.Neoantigens)))
text(x=0.7,y=-2.5,paste("r=",signif(res$estimate,2),"; p=",signif(res$p.value,2),sep=""))
mtext("SNV Neoantigens",cex=1.5)
dev.off()

#######################

# CPTAC validation
###################


if(!file.exists("data/CPTAC_LUAD_exprs.gct")){
  message("Downloading CPTAC LUAD RNA expression data")
  data_url <- "https://www.dropbox.com/s/5aa6yp564sr0tsr/luad-v3.0-rnaseq-prot-uq-rpkm-log2-NArm.gct?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = file.path("data/CPTAC_LUAD_exprs.gct"),mode="wb")
  }else{
    download.file(url = data_url,destfile = file.path("data/CPTAC_LUAD_exprs.gct"))
  }
}

gct <- read.gct("data/CPTAC_LUAD_exprs.gct")

expr_gct <- gct

if(!file.exists("data/CPTAC_LUAD_mut.gct")){
  message("Downloading CPTAC LUAD mutation data")
  data_url <- "https://www.dropbox.com/s/lyd91sfkzarns20/luad-v3.0-any-somatic-mutation-freq-by-gene.gct?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = file.path("data/CPTAC_LUAD_mut.gct"),mode="wb")
  }else{
    download.file(url = data_url,destfile = file.path("data/CPTAC_LUAD_mut.gct"))
  }
}

gct <- read.gct("data/CPTAC_LUAD_mut.gct")

mut_gct <- gct
rm(gct)
expr_gct <- expr_gct[,expr_gct["Type",]=="Tumor"]
mut_gct <- mut_gct[,mut_gct["Type",]=="Tumor"]

pats <- intersect(colnames(expr_gct),colnames(mut_gct))

expr_gct <- expr_gct[,pats]
mut_gct <- mut_gct[,pats]

mut_gct <- mut_gct[79:nrow(mut_gct),]
expr_gct <- expr_gct[86:nrow(expr_gct),]

score_mat <- score_mat[intersect(rownames(score_mat),rownames(expr_gct)),]

expr_mat <- 2^matrix(as.numeric(expr_gct),nrow=nrow(expr_gct),ncol=ncol(expr_gct),dimnames=dimnames(expr_gct))
expr_mat[is.na(expr_mat)] <- 0

expr_mat <- t(t(expr_mat)/rowSums(t(expr_mat)))

z_mat <- t(scale(log10(t(expr_mat+1e-6))))

source("scripts/get_LCAM_scores.R")

lcam_cptac <- get_LCAM_scores(expr_mat)

score1 <- lcam_cptac[,1]
score2 <- lcam_cptac[,2]

plot(score1,score2)
pat_ord <- order(score1-score2)

mat <- z_mat[rownames(score_mat),pat_ord]


thresh <- 2
mat[mat > thresh] <- thresh
mat[mat < -thresh] <- -thresh

v <- c(-0.5:(nrow(mat)+0.5)/c(nrow(mat)-1))[cumsum(colSums(score_mat))[-ncol(score_mat)]+1]

png("output/figures/figure_s6f.png",height=2.75,width=4.16,units="in",res=300,pointsize = 5)
layout(matrix(1:4,nrow=2),widths=c(20,3),heights=c(1.5,20))
par(mar=c(0,3,2,2))
image(as.matrix(rep(c(1,2),times=c(18,nrow(score_mat)-18))),col=c(rgb(t(col2rgb("purple")*.7/255)),rgb(0,.7,0)),xaxt="n",yaxt="n")
abline(v=v[3],lwd=3)
box()

par(mar=c(10,3,0,2))


image(mat,col=bluered(50),xaxt="n",yaxt="n")
mtext(side=1,rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),las=2,line=0.25,font=3)
abline(v=v)
abline(v=v[3],lwd=3)
box()
#mtext("CPTAC: LUAD",cex=2,line=1)
dev.off()

# 
# png("output/figures/fig6_bulk_analysis/outfig_sup6i_revise.png",height=2.75,width=4.7,units="in",res=300,pointsize = 5)
# layout(matrix(1:2,nrow=1),widths=c(20,3))
# 
# image(mat,col=bluered(50),xaxt="n",yaxt="n")
# mtext(side=1,rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),las=2,line=0.25)
# abline(v=v)
# abline(v=v[3],col="green",lwd=2)
# box()
# mtext("CPTAC: LUAD",cex=2,line=1)
# 
# par(pin=c(0.125,1))
# image(t(1:100),col=bluered(100),xaxt="n",yaxt="n")
# box()
# mtext(paste(c("< -","> "),thresh,sep=""),at=c(0,1),las=2,line=0.25,side=4)
# mtext("Z-score",line=0.5)
# 
# dev.off()

mut_tab <- matrix(as.numeric(mut_gct),nrow=nrow(mut_gct),dimnames = dimnames(mut_gct))

n_mut <- colSums(mut_tab)
pats <- intersect(names(n_mut),names(score1))

res <- cor.test(score1[pats]-score2[pats],log10(1+n_mut[pats]),method="spearman")

png("output/figures/figure_s6k.png",height=4,width=4,units="in",res=300)
par(mgp=c(2,1,0))
plot(log10((1+n_mut[pats])/48.2),score1[pats]-score2[pats],
     xlab=expression(paste("Log"["10"],"[TMB/Mb]",sep="")),
     ylab="LCAM score",pch=16,col=alpha(1,0.4))
abline(lm((score1[pats]-score2[pats])~log10((1+n_mut[pats])/48.2)))
mtext("CPTAC: LUAD",cex=1.5)
dev.off()
#cor.test(log10(1+n_mut[pats]),score1[pats]-score2[pats],method="spearman")
# 
# mask <- mut_gct["Stage",pats]=="1A"
# plot(log10(1+n_mut[pats[!mask]]),score1[pats[!mask]]-score2[pats[!mask]],
#      xlab="LogTMB",
#      ylab="LCAM score",xlim=c(0.5,3.5),ylim=c(-3,2.5),pch=16,col=alpha(1,0.4))
# abline(lm((score1[pats[!mask]]-score2[pats[!mask]])~log10(1+n_mut[pats[!mask]])))
# points(log10(1+n_mut[pats[mask]]),score1[pats[mask]]-score2[pats[mask]],
#      xlab="LogTMB",
#      ylab="LCAM score",xlim=c(0.5,3.5),ylim=c(-3,2.5),pch=16,col=alpha(2,0.4))
# abline(lm((score1[pats[mask]]-score2[pats[mask]])~log10(1+n_mut[pats[mask]])),col="red")
# mtext("CPTAC: LUAD",cex=1.5)
# 
# boxplot((score1-score2)~factor(mut_gct["Stage",names(score1)]))
# points(factor(mut_gct["Stage",names(score1)]),score1-score2)
# TMB / single gene boxplots
class <- array(NA,ncol(mut_tab),dimnames=list(colnames(mut_tab)))
class[colSums(mut_tab[c("TP53","KRAS","STK11","EGFR"),])==0] <- 1
class[colSums(mut_tab[c("TP53","KRAS","STK11"),])==0 & mut_tab["EGFR",]>0] <- 2
class[colSums(mut_tab[c("TP53","KRAS","EGFR"),])==0 & mut_tab["STK11",]>0] <- 3
class[colSums(mut_tab[c("TP53","STK11","EGFR"),])==0 & mut_tab["KRAS",]>0] <- 4
class[colSums(mut_tab[c("TP53","EGFR"),])==0 & mut_tab["KRAS",]>0 & mut_tab["STK11",]>0] <- 5
class[colSums(mut_tab[c("KRAS","STK11"),])==0 & mut_tab["TP53",]>0 & mut_tab["EGFR",]>0] <- 6
class[colSums(mut_tab[c("KRAS","STK11","EGFR"),])==0 & mut_tab["TP53",]>0] <- 7
class[colSums(mut_tab[c("KRAS","EGFR"),])==0 & mut_tab["TP53",]>0 & mut_tab["STK11",]>0] <- 8
class[colSums(mut_tab[c("STK11","EGFR"),])==0 & mut_tab["TP53",]>0 & mut_tab["KRAS",]>0] <- 9
class[colSums(mut_tab[c("EGFR"),,drop=F])==0 & mut_tab["TP53",]>0 & mut_tab["KRAS",]>0 & mut_tab["STK11",]>0] <- 10

class_mat <- matrix(c(0,0,0,0,
                      0,0,0,1,
                      0,0,1,0,
                      0,1,0,0,
                      0,1,1,0,
                      1,0,0,1,
                      1,0,0,0,
                      1,1,0,0,
                      1,0,1,0,
                      1,1,1,0),nrow=4,dimnames=list(c("TP53","KRAS","STK11","EGFR"),NULL))
class_mat <- class_mat[4:1,]
nclass <- ncol(class_mat)

png("output/figures/figure_s7a.png",height=1.43,width=2,units="in",res=300,pointsize=6)
layout(matrix(1:2,nrow=2),heibbbbghts=c(2,1))
par(mar=c(0,4,2,2))
b <- boxplot((score1[pats]-score2[pats])~class[pats],range=0,xaxt="n",xaxs="i",xlim=c(0.75,10.25),
             ylab="LCAM score",axes=F); box(); axis(side=2,at=c(-2,0,2),las=2)
points(jitter(class[pats]),score1[pats]-score2[pats],pch=16,col=alpha(1,0.5),cex=0.5,
       ylab="LCAMhi-LCAMlo score")
par(mar=c(1,4,0.5,2))
image(t(class_mat),col=c(0,1),xaxt="n",yaxt="n"); box()
mtext(rownames(class_mat),side=2,las=2,at=seq(0,1,1/(nrow(class_mat)-1)),font=3,las=2,line=0.25)
abline(v=c(seq(-1/2/(nclass-1),1+1/2/(nclass-1),(1+1/(nclass-1))/nclass)))
abline(h=c(seq(-1/2/(4-1),1+1/2/(4-1),(1+1/(4-1))/4)))
dev.off()


png("output/figures/figure_s7b.png",height=1.43,width=2,units="in",res=300,pointsize=6)
layout(matrix(1:2,nrow=2),heights=c(2,1))
par(mar=c(0,4,2,2))
b <- boxplot(n_mut[pats]/48.2~class[pats],log="y",range=0,xaxt="n",xaxs="i",xlim=c(0.75,10.25),
             ylab="TMB/Mb",axes=F); box(); axis(side=2,at=c(0.02,0.5,20),labels=F)
mtext(side=2,at=c(0.02,0.5,20),c("0.02","0.5","20"),line=0.75,las=2)
points(jitter(class[pats]),n_mut[pats]/48.2,pch=16,col=alpha(1,0.5),cex=0.5,
       ylab="TMB")
par(mar=c(1,4,0.5,2))
image(t(class_mat),col=c(0,1),xaxt="n",yaxt="n"); box()
mtext(rownames(class_mat),side=2,las=2,at=seq(0,1,1/(nrow(class_mat)-1)),font=3,las=2,line=0.25)
abline(v=c(seq(-1/2/(nclass-1),1+1/2/(nclass-1),(1+1/(nclass-1))/nclass)))
abline(h=c(seq(-1/2/(4-1),1+1/2/(4-1),(1+1/(4-1))/4)))
dev.off()

}



