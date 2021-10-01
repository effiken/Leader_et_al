library(scales)
library(Matrix)
library(viridis)
library(seriation)
library(matrixStats)
library(MASS)

figure_3def_s3i <- function(){
  
  
source(file.path(scClustering_dir,"DE.r"))

  
#    DE_between_two_sets(ldm = lung_ldm,mask_bg = bg_mask,mask_fg = fg_mask,nmin_umi_thresh = 0,nmin_cells_with_min_umi = 3,
#                      #                                 reg=1e-6,nchunks = 50,n_per_chunk = 1e3,noise_correction = F)
  
  


# 1. load the data
output_dir <- "output/figures/"

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors=F)
annots_list <- read.csv("input_tables/annots_list.csv",
                        r=1,h=1,stringsAsFactors = F)


v2_tissue_mask <- sample_annots[lung_ldm$dataset$cell_to_sample,]$library_chemistry=="V2" #& sample_annots[lung_ldm$dataset$cell_to_sample,]$tissue=="Tumor"

cell_to_annot <- annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"]

s <- split(colnames(lung_ldm$dataset$umitab),cell_to_annot)
s <- s[c("CD14 mono","cDC2",grep("MoMac",names(s),v=T),"cDC1","DC3","AM","AZU1_mac","IM","CD16 mono")]

expr <- lapply(s,function(x){rs <- rowSums(lung_ldm$dataset$umitab[,x]); return(rs/sum(rs))})
expr <- do.call(cbind,expr)

mono_expr <- expr[,"CD14 mono"]
dc_expr <- expr[,"cDC2"]
mac_expr <- rowMaxs(expr[,grep("MoMac",colnames(expr),v=T)])
moDC_expr <-  expr[,"DC3"]
AM_expr <- expr[,"AM"]

reg <- 1e-6

AM_v_mono <- log2((reg+AM_expr)/(reg+mono_expr))
AM_v_momac <- log2((reg+AM_expr)/(reg+mac_expr))
momac_v_mono <- log2((reg+mac_expr)/(reg+mono_expr))

max_cells <- 10000
set.seed(910340)

if(!file.exists("output/statistics/figure_s3i_DE.rd")){
  message("Performing differential expression between AMs and monocytes")
DE_am_vs_mono <- DE_between_two_sets(ldm=lung_ldm,
                                     mask_fg = sample(s$AM,pmin(length(s$AM),max_cells)),
                                     mask_bg = sample(s$`CD14 mono`,pmin(length(s$`CD14 mono`),max_cells)),
                                     nmin_cells_with_min_umi = 1,
                                     nmin_umi_thresh = 1,
                                     nchunks = 20,
                                     n_per_chunk = 100,
                                     reg = 1e-6,
                                     noise_correction = F)
message("Performing differential expression between AMs and MoMacs")
DE_am_vs_momac <- DE_between_two_sets(ldm=lung_ldm,
                                      mask_fg = sample(s$AM,pmin(length(s$AM),max_cells)),
                                      mask_bg = sample(unlist(s[grep("MoMac",names(s))]),pmin(length(unlist(s[grep("MoMac",names(s))])),max_cells)),
                                      nmin_cells_with_min_umi = 1,
                                      nmin_umi_thresh = 1,
                                      nchunks = 20,
                                      n_per_chunk = 100,
                                      reg = 1e-6,
                                      noise_correction = F)
message("Performing differential expression between MoMacs and monocytes")
DE_momac_vs_mono <- DE_between_two_sets(ldm=lung_ldm,
                                        mask_fg = sample(unlist(s[grep("MoMac",names(s))]),pmin(length(unlist(s[grep("MoMac",names(s))])),max_cells)),
                                        mask_bg = sample(s$`CD14 mono`,pmin(length(s$`CD14 mono`),max_cells)),
                                        nmin_cells_with_min_umi = 1,
                                        nmin_umi_thresh = 1,
                                        nchunks = 20,
                                        n_per_chunk = 100,
                                        reg = 1e-6,
                                        noise_correction = F)
save(DE_am_vs_mono,DE_am_vs_momac,DE_momac_vs_mono,file="output/statistics/figure_s3i_DE.rd")
}else{
  load("output/statistics/figure_s3i_DE.rd")
}

am_genes <- intersect(rownames(DE_am_vs_mono),rownames(DE_am_vs_momac))
am_genes_new <- am_genes[AM_v_mono[am_genes] > 0.75 & AM_v_momac[am_genes] > 1 &
                       DE_am_vs_mono[am_genes,]$adj.p.value < 1e-3 &
                       DE_am_vs_momac[am_genes,]$adj.p.value < 1e-3]
  
momac_genes <- intersect(rownames(DE_am_vs_momac),rownames(DE_momac_vs_mono))
momac_genes_new <- momac_genes[AM_v_momac[momac_genes] < -1 & momac_v_mono[momac_genes] > 1 &
                             DE_am_vs_momac[momac_genes,]$adj.p.value < 1e-3 &
                             DE_momac_vs_mono[momac_genes,]$adj.p.value < 1e-3]

mono_genes <- intersect(rownames(DE_am_vs_mono),rownames(DE_momac_vs_mono))
mono_genes_new <- mono_genes[AM_v_mono[mono_genes] < -1 & momac_v_mono[mono_genes] < -1 &
                           DE_am_vs_mono[mono_genes,]$adj.p.value < 1e-3 &
                           DE_momac_vs_mono[mono_genes,]$adj.p.value < 1e-3]



mono_genes <- names(mono_expr)[AM_v_mono < -1 & momac_v_mono < -1]
am_genes <- names(mono_expr)[AM_v_mono > 0.75 & AM_v_momac > 1 ]
momac_genes <- names(mono_expr)[AM_v_momac < -1 & momac_v_mono > 1]

other_genes <- setdiff(names(momac_v_mono),c(mono_genes,am_genes,momac_genes))
png(file.path(output_dir,"figure_s3i.png"),
    height=1.09,width=4.33,pointsize = 6,res=1000,units="in")
layout(matrix(1:3,nrow=1))
par(bty="n",mgp=c(1.5,0.75,0),mar=c(2.3,2.3,0,2.1))

m1 <- log10(rowMeans(cbind(mono_expr,mac_expr))+1e-6)
plot(m1[other_genes],momac_v_mono[other_genes],xlab="",ylab="",pch=16,cex=0.35,col=alpha("grey",.4),ylim=c(-6,6))
points(m1[mono_genes],momac_v_mono[mono_genes],pch=16,cex=0.35,col=rgb(t(col2rgb("orange")*.7/255)))
points(m1[momac_genes],momac_v_mono[momac_genes],col=alpha(rgb(0,.7,0),.9),pch=16,cex=0.35)
points(m1[am_genes],momac_v_mono[am_genes],col=alpha("purple",.6),pch=16,cex=0.35)

m1 <- log10(rowMeans(cbind(mono_expr,AM_expr))+1e-6)
plot(m1[other_genes],AM_v_mono[other_genes],xlab="",ylab="",pch=16,cex=0.35,col=alpha("grey",.4),ylim=c(-6,6))
points(m1[mono_genes],AM_v_mono[mono_genes],pch=16,cex=0.35,col=rgb(t(col2rgb("orange")*.7/255)))
points(m1[momac_genes],AM_v_mono[momac_genes],col=alpha(rgb(0,.7,0),.9),pch=16,cex=0.35)
points(m1[am_genes],AM_v_mono[am_genes],col=alpha("purple",.6),pch=16,cex=0.35)

m1 <- log10(rowMeans(cbind(mac_expr,AM_expr))+1e-6)
plot(m1[other_genes],AM_v_momac[other_genes],xlab="",ylab="",pch=16,cex=0.35,col=alpha("grey",.4),ylim=c(-6,6))
points(m1[mono_genes],AM_v_momac[mono_genes],pch=16,cex=0.35,col=rgb(t(col2rgb("orange")*.7/255)))
points(m1[momac_genes],AM_v_momac[momac_genes],col=alpha(rgb(0,.7,0),.9),pch=16,cex=0.35)
points(m1[am_genes],AM_v_momac[am_genes],col=alpha("purple",.6),pch=16,cex=0.35)
dev.off()









annots_ord <- c("AM","MoMac-IV","MoMac-III","MoMac-II","MoMac-I","CD14 mono")

cells <- annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"]%in%annots_ord

s <- split(colnames(lung_ldm$dataset$umitab[,cells]),annots_list[lung_ldm$dataset$cell_to_cluster[cells],"sub_lineage"])
s <- lapply(s,function(x){sample(x,min(2000,length(x)))})
cells <- unlist(s)

numi <- colSums(lung_ldm$dataset$umitab[,cells])
xmin <- 10^-4
ymin <- 10^-4
AM_score <- log10(xmin+colSums(lung_ldm$dataset$umitab[am_genes,cells])/numi)
#IM_score <- log10(ymin+colSums(lung_ldm$dataset$umitab[IM_genes,cells])/numi)
mono_score <- log10(1e-4+colSums(lung_ldm$dataset$umitab[mono_genes,cells])/numi)
momac_score <- log10(1e-4+colSums(lung_ldm$dataset$umitab[momac_genes,cells])/numi)

ord <- sample(1:length(AM_score))

genes <- c("PPARG","SERPINA1","MME","FOLR3","TREM1","TREM2","SPP1","MERTK","FOLR2")

xlim <- c(-3.3,-1.4); ylim<-c(-2.3,-0.8)

annot_to_cols <- rgb(t(col2rgb(c("yellow",6,1,2,3,4)))*.7/255)
names(annot_to_cols) <-  c("CD14 mono","AM","MoMac-I","MoMac-II","MoMac-III","MoMac-IV")



png(file.path(output_dir,"figure_3de.png"),height=1,width=3,units="in",res=1000,pointsize=6)
layout(matrix(1:3,nrow=1))
par(pin=c(0.7,0.7),mar=c(1,1,2,1))

plot(AM_score[ord],momac_score[ord],
     col=alpha(annot_to_cols[annots_list[lung_ldm$dataset$cell_to_cluster[cells],"sub_lineage"][ord]],0.2),
     pch=16,cex=0.2,xlim=xlim,ylim=ylim,
     ylab="",xlab="",xaxt="n",yaxt="n")

plot(AM_score[ord],momac_score[ord],
     col="white",
     pch=16,cex=0.2,xlim=xlim,ylim=ylim,
     ylab="",xlab="",xaxt="n",yaxt="n")

for(annot in unique(names(annot_to_cols))){
  mask <- annots_list[lung_ldm$dataset$cell_to_cluster[cells],"sub_lineage"][ord]==annot
  z <- kde2d(x=AM_score[ord[mask]],y=momac_score[ord[mask]],n=50)
  contour(x=z$x,y=z$y,z=z$z,add=T,col=annot_to_cols[annot],nlevels=3,drawlabels=F,lwd=0.7)
}


val <- mono_score
val[val < -1.9] <- -1.9
val[val > -1] <- -1
val <- val-min(val)
val <- round(val/range(val)*99+1)
plot(AM_score[cells],momac_score[cells],
     col=rev(alpha(viridis(100),0.8))[val],
     pch=16,cex=0.2,xlim=xlim,ylim=ylim,
     xaxt="n",yaxt="n",ylab="",xlab="")

dev.off()


png(file.path(output_dir,"figure_3f.png"),height=3,width=3,units="in",res=1000,pointsize=6)
layout(matrix(1:9,nrow=3))
par(mar=c(1,1,2,1))
for(gene in genes){
  val <- log10(1e-4+lung_ldm$dataset$umitab[gene,cells]/numi)
  val <- val-min(val)
  val <- round(val/range(val)*99+1)
  plot(AM_score[cells],momac_score[cells],
       col=rev(alpha(magma(120)[21:120],0.8))[val],
       pch=16,cex=0.2,xlim=xlim,ylim=ylim,
       xaxt="n",yaxt="n",ylab="",xlab="")
  mtext(gene,font=3)
}
dev.off()
}

