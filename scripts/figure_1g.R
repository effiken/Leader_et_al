
library(matrixStats)

figure_1g <- function(){

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",
                          r=1,h=1,stringsAsFactors = F)
annot_lists <- read.csv("input_tables/annots_list.csv",
                        r=1,h=1,stringsAsFactors=F)

if(!exists("lung_ldm")){
  message("Loading Sinai scRNA data into R")
  load("data/lung_ldm.rd")
}

c2cluster_sinai <- lung_ldm$dataset$cell_to_cluster
c2sample_sinai <- lung_ldm$dataset$cell_to_sample

rm(lung_ldm)

if(!exists("lambrechts_ldm")){
  message("Loading Lambrechts et al. data into R")
  load("data/lambrechts_ldm.rd")
}
c2cluster_lam <- lambrechts_ldm$dataset$cell_to_cluster
c2sample_lam <- lambrechts_ldm$dataset$cell_to_sample
rm(lambrechts_ldm)

patient_tissue_to_cluster_sinai <- table(
  apply(sample_annots[as.character(c2sample_sinai),c("patient_ID","tissue")],1,paste,collapse="_"),
  annot_lists[as.character(c2cluster_sinai),"sub_lineage"]
)

patient_tissue_to_cluster_lam <- table(
  apply(sample_annots[as.character(c2sample_lam),c("patient_ID","tissue")],1,paste,collapse="_"),
  annot_lists[as.character(c2cluster_lam),"sub_lineage"]
)

patient_tissue_to_cluster_sinai <- patient_tissue_to_cluster_sinai[,-1]
patient_tissue_to_cluster_lam <- patient_tissue_to_cluster_lam[,-1]

patient_tissue_to_cluster_sinai <- patient_tissue_to_cluster_sinai/rowSums(patient_tissue_to_cluster_sinai)
patient_tissue_to_cluster_lam <- patient_tissue_to_cluster_lam/rowSums(patient_tissue_to_cluster_lam)


normal_sinai <- unlist(lapply(strsplit(grep("Normal",rownames(patient_tissue_to_cluster_sinai),v=T),"_"),function(x){x[1]}))
tumor_sinai <- unlist(lapply(strsplit(grep("Tumor",rownames(patient_tissue_to_cluster_sinai),v=T),"_"),function(x){x[1]}))

normal_lam <- unlist(lapply(strsplit(grep("Normal",rownames(patient_tissue_to_cluster_lam),v=T),"_"),function(x){paste(x[1:2],collapse="_")}))
tumor_lam <- unlist(lapply(strsplit(grep("Tumor",rownames(patient_tissue_to_cluster_lam),v=T),"_"),function(x){paste(x[1:2],collapse="_")}))

sinai_pats <- intersect(normal_sinai,tumor_sinai)
lam_pats <- intersect(normal_lam,tumor_lam)

reg <- 1e-3
delta.sinai <- log2((reg+patient_tissue_to_cluster_sinai[paste(sinai_pats,"Tumor",sep="_"),])/
                      (reg+patient_tissue_to_cluster_sinai[paste(sinai_pats,"Normal",sep="_"),]))
delta.lam <- log2((reg+patient_tissue_to_cluster_lam[paste(lam_pats,"Tumor",sep="_"),])/
                      (reg+patient_tissue_to_cluster_lam[paste(lam_pats,"Normal",sep="_"),]))


## plot
m.sinai <- colMeans(delta.sinai)
sem.sinai <- colSds(delta.sinai)/sqrt(nrow(delta.sinai))
m.lam <- colMeans(delta.lam)
sem.lam <- colSds(delta.lam)/sqrt(nrow(delta.lam))


xlim <- range(c(m.sinai-sem.sinai,m.sinai+sem.sinai))
ylim <- range(c(m.lam-sem.lam,m.lam+sem.lam))
xlim <- c(-5,5)
ylim <- c(-5,5)

r2 <- cor(m.lam,m.sinai)^2

png("output/figures/figure_1g.png",height=4,width=6,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
plot(m.sinai,m.lam,xlim=xlim,ylim=ylim,bty="n",col="white",
     xlab="Mt. Sinai",
     ylab="Lambrechts")
segments(x0=m.sinai-sem.sinai,x1=m.sinai+sem.sinai,y0=m.lam,y1=m.lam,col="grey")
segments(x0=m.sinai,x1=m.sinai,y0=m.lam-sem.lam,y1=m.lam+sem.lam,col="grey")
abline(h=0,col="red")
abline(v=0,col="red")
abline(c(0,1),lty=2)

clust_mask <- names(m.sinai)%in%c("T_activated","MoMac-II","B","AM","NK")
text(x=m.sinai[clust_mask],y=m.lam[clust_mask],names(m.sinai)[clust_mask])
text(x=2,y=-2,expression(paste("R"^"2","=",0.78,sep="")))
mtext("Log2FC: Tumor vs. nLung",cex=1.5)
dev.off()

}
