library(seriation)
library(Matrix)

figure_s1hi <- function(){



sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)

if(!exists("lambrechts_ldm")){
load("data/lambrechts_ldm.rd")
}
cell_to_annot <- annots_list[lambrechts_ldm$dataset$cell_to_cluster,]$sub_lineage

patient_site <- apply(sample_annots[lambrechts_ldm$dataset$cell_to_sample,c("patient_ID","biopsy_site")],1,paste,collapse="_")

tab <- table(cell_to_annot,patient_site)

tab <- tab[,!grepl("NA",colnames(tab))]
tab <- tab[-1,]

tab <- t(t(tab)/rowSums(t(tab)))

samp_h <- hclust(as.dist(1-cor(tab,method="spearman")))
samp_ord <- samp_h$order
clust_ord <- hclust(as.dist(1-cor(t(tab),method="spearman")))$order

mat <- log2((1e-3+tab)/rowMeans(1e-3+tab))
mat[mat > 3] <- 3
mat[mat < -3] <- -3

png("output/figures/figure_s1h.png",height=8,width=8,units="in",res=300,bg="transparent")
layout(matrix(1:2,nrow=2),heights=c(1,3))
par(mar=c(0,6,3,3))
plot(as.dendrogram(samp_h),horiz=F,xlab="sqrt((1-cor)/2)",type="triangle",xaxs="i",axes=F,leaflab="none")

par(mar=c(10,6,0,3))
image(t(mat[clust_ord,samp_ord]),col=bluered(50),las=2,xaxt="n",yaxt="n")
mtext(side=2,rownames(mat)[clust_ord],at=seq(0,1,1/(nrow(mat)-1)),las=2,line=0.25)
mtext(side=1,colnames(mat)[samp_ord],at=seq(0,1,1/(ncol(mat)-1)),las=2,line=0.25,
      col=rgb(0.7/255*t(col2rgb(1+as.integer(factor(unlist(lapply(strsplit(colnames(mat)[samp_ord],"_"),function(x){x[2]}))))))),
      font=2)
box()
dev.off()

patient <- unlist(lapply(strsplit(colnames(tab),"_"),function(x){x[2]}))

mat <- log10(1e-3+tab)
dists <- matrix(NA,nrow=ncol(tab),ncol=ncol(tab),dimnames=list(colnames(tab),colnames(tab)))
for(row in colnames(tab)){
  for(col in colnames(tab)){
    dists[row,col] <- sqrt(sum((mat[,row]-mat[,col])^2))
  }
}
diag(dists) <- NA

same_mask <- list()
for(pat in 1:8){
  same_mask[[pat]] <- matrix(1,3,3)
}
same_mask <- bdiag(same_mask)

same_dists <- dists*same_mask*lower.tri(same_mask)
same_dists[same_dists==0] <- NA

diff_mask <- 1-same_mask
diff_dists <- dists*diff_mask*lower.tri(diff_mask)
diff_dists[diff_dists==0] <- NA

dists <- c(array(same_dists),array(diff_dists))
type <- rep(c("same","diff"),times=c(length(same_dists),length(diff_dists)))

png("output/figures/figure_s1i.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
boxplot(dists~type,range=0,names=c("Patient-\nunmatched","Patient-\nmatched"),las=2,ylab="Euclidean distance",ylim=c(min(dists,na.rm=T),max(dists,na.rm=T)+1))
points(jitter(as.integer(factor(type))),dists,pch=16,cex=0.5,col=alpha(1,0.6))
segments(x0=1,x1=2,y0=5.5,y1=5.5)
text("***",x=1.5,y=5.8,cex=1.5)
mtext("Inter- vs. intra-patient\ndiversity",cex=1.5)
dev.off()
}