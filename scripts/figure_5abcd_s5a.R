
# Figure 5 part 1: correlation of normalized celltypes across patients, 
# heatmaps of patient frequencies ordered by LCAM score

library(matrixStats)
library(seriation)
library(RColorBrewer)

figure_5abcd_s5a <- function(){

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)
annots_list$norm_group[annots_list$lineage=="MNP"] <- "MNP"

## This part of the code figures out the right cluster order, and plots the correlation matrix, using only the V2 beads
#####

if(!exists("lung_ldm")){
  message("Loading Mt. Sinai scRNAseq data into R")
  load("data/lung_ldm.rd")
}

cell_mask <- names(lung_ldm$dataset$cell_to_sample)[
  sample_annots[lung_ldm$dataset$cell_to_sample,"library_chemistry"]=="V2" &
    sample_annots[lung_ldm$dataset$cell_to_sample,"prep"]=="beads" &
    sample_annots[lung_ldm$dataset$cell_to_sample,"tissue"]=="Tumor"
  ]

tab <- table(sample_annots[lung_ldm$dataset$cell_to_sample[cell_mask],"patient_ID"],
             annots_list[lung_ldm$dataset$cell_to_cluster[cell_mask],"sub_lineage"])

tab <- tab[,-1]

tab <- tab/rowSums(tab)

for(norm_group in c("T","B&plasma","MNP","lin_neg")){
  tab_tmp <- tab[,annots_list$norm_group[match(colnames(tab),annots_list$sub_lineage)]==norm_group]
  tab_tmp <- tab_tmp / rowSums(tab_tmp)
  tab[,colnames(tab_tmp)] <- tab_tmp
}

clust_cor <- cor(log10(1e-2+tab),method="spearman")


LCAM <- c("T_activated","IgG","MoMac-II")
LCAM_score <- rowSums(log(tab[,LCAM]+1e-2))

#clust_ord <- order(cor(tab,LCAM_score[rownames(tab)],method="spearman")[,1])

clust_ord <- get_order(seriate(cor(tab,method="spearman")),method="OLO")


mat <- clust_cor[clust_ord,clust_ord]

thresh <- 0.7
mat[mat > thresh] <- thresh
mat[mat < -thresh] <- -thresh
mat <- mat/thresh

mat <- round((mat+1)/2*49)+1
#diag(mat) <- max(mat)+1
mat <- mat[,ncol(mat):1]

pal_cor <- colorRampPalette(rev(brewer.pal(11,"PiYG")))
col <- c(pal_cor(50)[min(mat):max(mat)])

#plotting 5a

png(file.path("output/figures/figure_5a.png"),height=3,width=3.5,pointsize=5,res=300,units="in",bg="transparent")

layout(matrix(1:2,nrow=1),widths=c(10,2))

par(pin=c(1.5,1.5),mar=c(2,2,2,2),oma=c(5,5,5,5))

image(mat,col=col,xaxt="n",yaxt="n")
box()
mtext(side=1,at=seq(0,1,1/(nrow(mat)-1)),rownames(mat),las=2,line=0.25)
mtext(side=2,at=seq(0,1,1/(nrow(mat)-1)),rev(rownames(mat)),las=2,line=0.25)
mtext("Celltype frequency spearman cor.",line=0.5,cex=2)

par(pin=c(0.125,1))
image(t(1:100),col=pal_cor(100),xaxt="n",yaxt="n")
mtext(paste(c("<",">"),thresh),at=c(0,1),side=4,las=2,line=0.5)
box()
dev.off()



###### Figure 5B, C, D and S5A


tab <- table(apply(sample_annots[lung_ldm$dataset$cell_to_sample,c("patient_ID","tissue")],1,paste,collapse=" "),
             annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"])

rm("lung_ldm")
message("loading Lambrechts et al. data into R")
load("data/lambrechts_ldm.rd")
tab <- rbind(tab,table(apply(sample_annots[lambrechts_ldm$dataset$cell_to_sample,c("patient_ID","tissue")],1,paste,collapse=" "),
                       annots_list[lambrechts_ldm$dataset$cell_to_cluster,"sub_lineage"]))
rm("lambrechts_ldm")
message("loading Zilionis et al. data into R")
load("data/zilionis_ldm.rd")
tab <- rbind(tab,table(apply(sample_annots[lung_ldm_zili$dataset$cell_to_sample,c("patient_ID","tissue")],1,paste,collapse=" "),
                       annots_list[lung_ldm_zili$dataset$cell_to_cluster,"sub_lineage"]))
rm("lung_ldm_zili")
tab <- tab[,-1]

tab <- tab/rowSums(tab)

tab_raw <- tab
for(norm_group in c("T","B&plasma","MNP","lin_neg")){
  tab_tmp <- tab[,annots_list$norm_group[match(colnames(tab),annots_list$sub_lineage)]==norm_group]
  tab_tmp <- tab_tmp / rowSums(tab_tmp)
  tab[,colnames(tab_tmp)] <- tab_tmp
}


LCAM_score <- rowSums(log(tab[,LCAM]+1e-2))
resting_clusts <- c("B","AM","cDC2","AZU1_mac","Tcm/naive_II","cDC1")
resting_score <- rowSums(log(tab[,resting_clusts]+1e-2))

pat_ord_tumor <- order((LCAM_score-resting_score)[grep("Tumor",rownames(tab),v=T)])
pat_ord_normal <- order((LCAM_score-resting_score)[grep("Normal",rownames(tab),v=T)])

# plot highlighted clusters
tab_tumor <- tab[grep("Tumor",rownames(tab)),]
tab_normal <- tab[grep("Normal",rownames(tab)),]

mat_tumor <- t(tab_tumor[pat_ord_tumor,rev(c(resting_clusts,LCAM))])

norm_ord <- match(unlist(lapply(strsplit(colnames(mat_tumor)," "),function(x){x[1]})),
                  unlist(lapply(strsplit(rownames(tab_normal)," "),function(x){x[1]})))

mat_normal <- t(tab_normal[norm_ord,rev(c(resting_clusts,LCAM))])
mat_normal <- mat_normal[,!is.na(colnames(mat_normal))]

clust.means <- rowMeans(cbind(mat_normal,mat_tumor))
mat_normal <- log2((1e-2+mat_normal)/(1e-2+clust.means))
mat_tumor <- log2((1e-2+mat_tumor)/(1e-2+clust.means))

thresh <- 2
mat_normal[mat_normal < -thresh] <- -thresh
mat_normal[mat_normal > thresh] <- thresh
mat_tumor[mat_tumor < -thresh] <- -thresh
mat_tumor[mat_tumor > thresh] <- thresh

mat_normal <- mat_normal/thresh
mat_tumor <- mat_tumor/thresh

mat_normal <- round((mat_normal + 1)/2*49)+1
mat_tumor <- round((mat_tumor+1)/2*49)+1

col_normal <- bluered(50)[min(mat_normal):max(mat_normal)]
col_tumor <- bluered(50)[min(mat_normal):max(mat_normal)]

h <- seq(-0.5,8.5)[4]/8

col_dataset_normal <- array(1,ncol(mat_normal)) + as.numeric(grepl("Lambrechts",colnames(mat_normal)))
col_dataset_tumor <- array(1,ncol(mat_tumor)) + as.numeric(grepl("Lambrechts",colnames(mat_tumor)))+2*as.numeric(grepl("zilionis",colnames(mat_tumor)))

png(file.path("output/figures/figure_5b.png"),height=2.7,width=4.76,pointsize=5,res=300,units="in",bg="transparent")

par(oma=c(5,5,5,5))

layout(matrix(1:6,nrow=2,ncol=3,byrow = T),widths = c(10,10,3),heights=c(6,2.5))

image(t(mat_normal),col=col_normal,xaxt="n",yaxt="n")
mtext(rownames(mat_normal),side=2,at=seq(0,1,1/(nrow(mat_normal)-1)),las=2,line=0.25)
mtext(lapply(strsplit(colnames(mat_normal)," "),function(x){x[1]}),side=1,at=seq(0,1,1/(ncol(mat_normal)-1)),las=2,line=0.25,cex=0.7)
mtext("nLung samples",cex=2,line=0.5)
abline(h=h,col="green",lwd=2)
box()

segments(y0=c(-0.5,2.5,8.5)/8,y1=c(-0.5,2.5,8.5)/8,x0=23.5/23,x1=1.25,xpd=NA,lty=2)

image(t(mat_tumor),col=col_tumor,xaxt="n",yaxt="n")
mtext(lapply(strsplit(colnames(mat_tumor)," "),function(x){x[1]}),side=1,at=seq(0,1,1/(ncol(mat_tumor)-1)),las=2,line=0.25,cex=0.7)
mtext("Tumor samples",cex=2,line=0.5)
abline(h=h,col="green",lwd=2)
box()

#par(pin=c(0.125,0.5))
image(t(1:100),col=bluered(100),xaxt="n",yaxt="n"); box()
mtext("Norm. frequency",line=0.5)
mtext(side=4,at=c(0,1),paste(c("< -",">"),thresh,sep=""),line=0.25,las=2)

par(pin=c())
image(as.matrix(col_dataset_normal),col=c(1,2),xaxt="n",yaxt="n")
abline(v=seq(-.5/length(col_dataset_normal),1+.5/(length(col_dataset_normal)+1),(1+2*.5/length(col_dataset_normal))/length(col_dataset_normal)),col="white",lwd=0.2)
box()

image(as.matrix(col_dataset_tumor),col=c(1,2,3),xaxt="n",yaxt="n")
abline(v=seq(-.5/length(col_dataset_tumor),1+.5/(length(col_dataset_tumor)+1),(1+2*.5/length(col_dataset_tumor))/length(col_dataset_tumor)),col="white",lwd=0.2)
box()

dev.off()

###  show all clusters in supp:


tab_tumor <- tab[grep("Tumor",rownames(tab)),]
tab_normal <- tab[grep("Normal",rownames(tab)),]

mat_tumor <- t(tab_tumor[pat_ord_tumor,rev(clust_ord)])

norm_ord <- match(unlist(lapply(strsplit(colnames(mat_tumor)," "),function(x){x[1]})),
                  unlist(lapply(strsplit(rownames(tab_normal)," "),function(x){x[1]})))

mat_normal <- t(tab_normal[norm_ord,rev(clust_ord)])
mat_normal <- mat_normal[,!is.na(colnames(mat_normal))]

clust.means <- rowMeans(cbind(mat_normal,mat_tumor),na.rm=T)
mat_normal <- log2((1e-2+mat_normal)/(1e-2+clust.means))
mat_tumor <- log2((1e-2+mat_tumor)/(1e-2+clust.means))

thresh <- 2
mat_normal[mat_normal < -thresh] <- -thresh
mat_normal[mat_normal > thresh] <- thresh
mat_tumor[mat_tumor < -thresh] <- -thresh
mat_tumor[mat_tumor > thresh] <- thresh

mat_normal <- mat_normal/thresh
mat_tumor <- mat_tumor/thresh

mat_normal <- round((mat_normal + 1)/2*49)+1
mat_tumor <- round((mat_tumor+1)/2*49)+1

mat_normal[is.na(mat_normal)] <- 25
mat_tumor[is.na(mat_tumor)] <- 25

col_normal <- bluered(50)[min(mat_normal):max(mat_normal)]
col_tumor <- bluered(50)[min(mat_normal):max(mat_normal)]

h <- seq(-0.5,8.5)[4]/8


plot_me <- function(){
  par(oma=c(5,5,5,5))
  
  layout(matrix(1:3,nrow=1),widths = c(10,10,3))
  
  image(t(mat_normal),col=col_normal,xaxt="n",yaxt="n")
  mtext(rownames(mat_normal),side=2,at=seq(0,1,1/(nrow(mat_normal)-1)),las=2,line=0.25,cex=0.5)
  mtext(lapply(strsplit(colnames(mat_normal)," "),function(x){x[1]}),side=1,at=seq(0,1,1/(ncol(mat_normal)-1)),las=2,line=0.25,cex=0.6)
  mtext("nLung samples",cex=2,line=0.5)
  box()
  
  image(t(mat_tumor),col=col_tumor,xaxt="n",yaxt="n")
  mtext(lapply(strsplit(colnames(mat_tumor)," "),function(x){x[1]}),side=1,at=seq(0,1,1/(ncol(mat_tumor)-1)),las=2,line=0.25,cex=0.6)
  mtext("Tumor samples",cex=2,line=0.5)
  box()
  
  par(pin=c(0.125,0.5))
  image(t(1:100),col=bluered(100),xaxt="n",yaxt="n"); box()
  mtext("Norm. frequency",line=0.5)
  mtext(side=4,at=c(0,1),paste(c("< -",">"),thresh,sep=""),line=0.25,las=2)
}

png(file.path("output/figures/figure_s5a.png"),height=2,width=4.76,pointsize=5,res=500,units="in",bg="transparent")
plot_me()
dev.off()

# Plot heatmap by lineage

lin_ord <- rev(c("NK","T","MNP","pDC","B&plasma","mast"))
tab <- tab_raw
lin_tab <- matrix(NA,nrow=nrow(tab),ncol=length(lin_ord),dimnames=list(rownames(tab),lin_ord))

for(lin in lin_ord){
  lin_tab[,lin] <- rowSums(tab[,annots_list$lineage[match(colnames(tab),annots_list$sub_lineage)]==lin,drop=F])
}

tab <- lin_tab

tab_tumor <- tab[grep("Tumor",rownames(tab)),]
tab_normal <- tab[grep("Normal",rownames(tab)),]

mat_tumor <- t(tab_tumor[pat_ord_tumor,])

norm_ord <- match(unlist(lapply(strsplit(colnames(mat_tumor)," "),function(x){x[1]})),
                  unlist(lapply(strsplit(rownames(tab_normal)," "),function(x){x[1]})))

mat_normal <- t(tab_normal[norm_ord,])
mat_normal <- mat_normal[,!is.na(colnames(mat_normal))]



clust.means <- rowMeans(cbind(mat_normal,mat_tumor),na.rm=T)
mat_normal <- log2((1e-2+mat_normal)/(1e-2+clust.means))
mat_tumor <- log2((1e-2+mat_tumor)/(1e-2+clust.means))

thresh <- 2
mat_normal[mat_normal < -thresh] <- -thresh
mat_normal[mat_normal > thresh] <- thresh
mat_tumor[mat_tumor < -thresh] <- -thresh
mat_tumor[mat_tumor > thresh] <- thresh

mat_normal <- mat_normal/thresh
mat_tumor <- mat_tumor/thresh

mat_normal <- round((mat_normal + 1)/2*49)+1
mat_tumor <- round((mat_tumor+1)/2*49)+1

mat_normal[is.na(mat_normal)] <- 25
mat_tumor[is.na(mat_tumor)] <- 25

col_normal <- bluered(50)[min(mat_normal):max(mat_normal)]
col_tumor <- bluered(50)[min(mat_normal):max(mat_normal)]

h <- seq(-0.5,8.5)[4]/8



plot_me <- function(){
  par(oma=c(5,5,5,5))
  
  layout(matrix(1:3,nrow=1),widths = c(10,10,3))
  
  image(t(mat_normal),col=col_normal,xaxt="n",yaxt="n")
  mtext(rownames(mat_normal),side=2,at=seq(0,1,1/(nrow(mat_normal)-1)),las=2,line=0.25,cex=1)
  mtext(lapply(strsplit(colnames(mat_normal)," "),function(x){x[1]}),side=1,at=seq(0,1,1/(ncol(mat_normal)-1)),las=2,line=0.25,cex=0.6)
  #mtext("nLung samples",cex=2,line=0.5)
  box()
  
  image(t(mat_tumor),col=col_tumor,xaxt="n",yaxt="n")
  mtext(lapply(strsplit(colnames(mat_tumor)," "),function(x){x[1]}),side=1,at=seq(0,1,1/(ncol(mat_tumor)-1)),las=2,line=0.25,cex=0.6)
  #mtext("Tumor samples",cex=2,line=0.5)
  box()
  
  par(pin=c(0.125,0.5))
  image(t(1:100),col=bluered(100),xaxt="n",yaxt="n"); box()
  mtext("Norm. frequency",line=0.5)
  mtext(side=4,at=c(0,1),paste(c("< -",">"),thresh,sep=""),line=0.25,las=2)
}

png("output/figures/figure_5c.png",height=1.5,width=4.76,pointsize=5,res=500,units="in",bg="transparent")
plot_me()
dev.off()


tab <- tab_raw

for(norm_group in c("T","B&plasma","MNP","lin_neg")){
  tab_tmp <- tab[,annots_list$norm_group[match(colnames(tab),annots_list$sub_lineage)]==norm_group]
  tab_tmp <- tab_tmp / rowSums(tab_tmp)
  tab[,colnames(tab_tmp)] <- tab_tmp
}

LCAM_score <- rowSums(log10(tab[,LCAM]+1e-2))
resting_clusts <- c("B","AM","cDC2","AZU1_mac","Tcm/naive_II","cDC1")
resting_score <- rowSums(log10(tab[,resting_clusts]+1e-2))

tumor_names <- names(LCAM_score)[grepl("Tumor",names(LCAM_score))]
col <- 1+as.numeric(grepl("Lambrechts",tumor_names))+2*as.numeric(grepl("zilionis",tumor_names))

png("output/figures/figure_5d.png",height=2,width=1.84,units="in",res=300,pointsize=5)
par(mar=c(5,5,5,1))
plot(10^LCAM_score[grepl("Tumor",names(LCAM_score))],10^resting_score[grepl("Tumor",names(LCAM_score))],log="xy",
     col=col,pch=16,cex=1.5,xaxt="n",yaxt="n",xlab="",ylab="")
mtext("LCAMhi vs. LCAMlo scores",cex=1.7,line=0.1)
legend(c("Dataset:","Mt. Sinai","Lambrechts","Zilionis"),x="bottomleft",pch=16,col=c(0,1:3))
mtext(side=1,"LCAMhi score",line=3,cex=1.5)
mtext(side=2,"LCAMlo score",line=3,cex=1.5)
axis(side=1,at=c(1e-5,1e-3,1e-1))
axis(side=2,at=c(1e-9,1e-7,1e-5))
dev.off()


}