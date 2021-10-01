library(scales)
library(Matrix)

figure_2g_s2de <- function(){
  
  
source(file.path(scClustering_dir,"DE.r"))

# 1. load the data
output_dir <- "output/figures/"

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors=F)
annots_list <- read.csv("input_tables/annots_list.csv",
                        r=1,h=1,stringsAsFactors = F)



v2_tissue_mask <- sample_annots[lung_ldm$dataset$cell_to_sample,]$library_chemistry=="V2"
mono_mask <- v2_tissue_mask & annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"]=="CD14 mono"
momac_mask <- v2_tissue_mask & grepl("MoMac",annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"])
dc2_mask <- v2_tissue_mask & annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"]=="cDC2"
moDC_mask <- v2_tissue_mask & annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"]=="DC3"
dc1_mask <- v2_tissue_mask & annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"]=="cDC1"


#Note: here "DC3" is called "moDC" in the code

if(!file.exists("output/statistics/figure_s2d_DE.rd")){

max_cells <- 10000
set.seed(910340)
message("Performing differential expression analysis between MoMacs and monocytes")
DE_momac_vs_mono <- DE_between_two_sets(ldm=lung_ldm,
                                        mask_fg = sample(colnames(lung_ldm$dataset$umitab)[momac_mask],pmin(sum(momac_mask),max_cells)),
                                        mask_bg = sample(colnames(lung_ldm$dataset$umitab)[mono_mask],pmin(sum(mono_mask),max_cells)),
                                        nmin_cells_with_min_umi = 1,
                                        nmin_umi_thresh = 1,
                                        nchunks = 20,
                                        n_per_chunk = 100,
                                        reg = 1e-6,
                                        noise_correction = F)
message("Performing differential expression analysis between cDC2 and MoMacs")
DE_dc2_vs_momac <- DE_between_two_sets(ldm=lung_ldm,
                                        mask_fg = sample(colnames(lung_ldm$dataset$umitab)[dc2_mask],pmin(sum(dc2_mask),max_cells)),
                                        mask_bg = sample(colnames(lung_ldm$dataset$umitab)[momac_mask],pmin(sum(momac_mask),max_cells)),
                                        nmin_cells_with_min_umi = 1,
                                        nmin_umi_thresh = 1,
                                        nchunks = 20,
                                        n_per_chunk = 100,
                                        reg = 1e-6,
                                        noise_correction = F)
message("Performing differential expression analysis between cDC2 and monocytes")
DE_dc2_vs_mono <- DE_between_two_sets(ldm=lung_ldm,
                                        mask_fg = sample(colnames(lung_ldm$dataset$umitab)[dc2_mask],pmin(sum(dc2_mask),max_cells)),
                                        mask_bg = sample(colnames(lung_ldm$dataset$umitab)[mono_mask],pmin(sum(mono_mask),max_cells)),
                                        nmin_cells_with_min_umi = 1,
                                        nmin_umi_thresh = 1,
                                        nchunks = 20,
                                        n_per_chunk = 100,
                                        reg = 1e-6,
                                        noise_correction = F)
save(DE_momac_vs_mono,DE_dc2_vs_momac,DE_dc2_vs_mono,file="output/statistics/figure_s2d_DE.rd")
}else{
  message("Loading myeloid sup-population differential expression results")
  load("output/statistics/figure_s2d_DE.rd")
}


mono_expr <- rowSums(lung_ldm$dataset$umitab[,mono_mask]); mono_expr <- mono_expr/sum(mono_expr)
momac_expr <- rowSums(lung_ldm$dataset$umitab[,momac_mask]); momac_expr <- momac_expr/sum(momac_expr)
dc2_expr <- rowSums(lung_ldm$dataset$umitab[,dc2_mask]); dc2_expr <- dc2_expr/sum(dc2_expr)
moDC_expr <- rowSums(lung_ldm$dataset$umitab[,moDC_mask]); moDC_expr <- moDC_expr/sum(moDC_expr)

reg <- 1e-6

momac_v_mono <- log2((reg+momac_expr)/(reg+mono_expr))
dc2_v_momac <- log2((reg+dc2_expr)/(reg+momac_expr))
dc2_v_mono <- log2((reg+dc2_expr)/(reg+mono_expr))
mono_v_dc <- -dc2_v_mono
mono_v_mac <- -momac_v_mono
dc_v_mac <- dc2_v_momac



#mono_genes <- names(mono_expr)[momac_v_mono < -1 & dc2_v_mono < -1 & mono_expr > 10^-6]
mono_genes <- intersect(rownames(DE_dc2_vs_mono),rownames(DE_momac_vs_mono))
mono_genes <- mono_genes[momac_v_mono[mono_genes] < -1 & dc2_v_mono[mono_genes] < -1 & mono_expr[mono_genes] > 1e-6 &
                           DE_dc2_vs_mono[mono_genes,]$adj.p.value<1e-3 & DE_momac_vs_mono[mono_genes,]$adj.p.value<1e-3]
#dc2_genes <- names(mono_expr)[dc2_v_momac > 1 & dc2_v_mono > 1 & dc2_expr > 10^-6]
dc2_genes <- intersect(rownames(DE_dc2_vs_mono),rownames(DE_dc2_vs_momac))
dc2_genes <- dc2_genes[dc2_v_momac[dc2_genes] > 1 & dc2_v_mono[dc2_genes]>1 & dc2_expr[dc2_genes] > 1e-6 &
                         DE_dc2_vs_momac[dc2_genes,]$adj.p.value < 1e-3 & 
                         DE_dc2_vs_mono[dc2_genes,]$adj.p.value < 1e-3]

#momac_genes <- names(mono_expr)[dc2_v_momac < -1 & momac_v_mono > 1 & momac_expr > 10^-6]
momac_genes <- intersect(rownames(DE_dc2_vs_momac),rownames(DE_momac_vs_mono))
momac_genes <- momac_genes[dc2_v_momac[momac_genes] < -1 & momac_v_mono[momac_genes] > 1 & momac_expr[momac_genes] > 1e-6 &
                             DE_dc2_vs_momac[momac_genes,]$adj.p.value<1e-3 & DE_momac_vs_mono[momac_genes,]$adj.p.value<1e-3]

other_genes <- setdiff(names(mono_expr),c(mono_genes,dc2_genes,momac_genes))


png(file.path(output_dir,"figure_s2d.png"),
    height=1.09,width=4.33,pointsize = 6,res=1000,units="in")
layout(matrix(1:3,nrow=1))
par(bty="n",mgp=c(1.5,0.75,0),mar=c(2.3,2.3,0,2.1))

m1 <- log10(rowMeans(cbind(mono_expr,dc2_expr))+1e-6)
plot(m1[other_genes],mono_v_dc[other_genes],xlab="",ylab="",pch=16,cex=0.35,col=alpha("grey",.4))
points(m1[mono_genes],mono_v_dc[mono_genes],pch=16,cex=0.35,col=rgb(t(col2rgb("orange")*.7/255)))
points(m1[dc2_genes],mono_v_dc[dc2_genes],col=alpha(rgb(0,.7,0),.9),pch=16,cex=0.35)
points(m1[momac_genes],mono_v_dc[momac_genes],col=alpha("purple",.6),pch=16,cex=0.35)

m1 <- log10(1e-6+rowMeans(cbind(mono_expr,momac_expr)))
plot(m1[other_genes],mono_v_mac[other_genes],xlab="",ylab="",pch=16,cex=0.35,col=alpha("grey",.4))
points(m1[mono_genes],mono_v_mac[mono_genes],pch=16,cex=0.35,col=rgb(t(col2rgb("orange")*.7/255)))
points(m1[dc2_genes],mono_v_mac[dc2_genes],col=alpha(rgb(0,.7,0),.9),pch=16,cex=0.35)
points(m1[momac_genes],mono_v_mac[momac_genes],col=alpha("purple",.6),pch=16,cex=0.35)

m1 <- log10(1e-6+rowMeans(cbind(dc2_expr,momac_expr)))
plot(m1[other_genes],dc_v_mac[other_genes],xlab="",ylab="",pch=16,cex=0.35,col=alpha("grey",.4))
points(m1[mono_genes],dc_v_mac[mono_genes],pch=16,cex=0.35,col=rgb(t(col2rgb("orange")*.7/255)))
points(m1[dc2_genes],dc_v_mac[dc2_genes],col=alpha(rgb(0,.7,0),.9),pch=16,cex=0.35)
points(m1[momac_genes],dc_v_mac[momac_genes],col=alpha("purple",.6),pch=16,cex=0.35)
dev.off()






numi <- colSums(lung_ldm$dataset$umitab)
mono_score <- colSums(lung_ldm$dataset$umitab[mono_genes,])/numi
dc_score <- colSums(lung_ldm$dataset$umitab[dc2_genes,])/numi
mac_score <- colSums(lung_ldm$dataset$umitab[momac_genes,])/numi
names(mono_score)=names(dc_score)=names(mac_score) <- colnames(lung_ldm$dataset$umitab)




suppressWarnings(cor.test(mono_score[moDC_mask],dc_score[moDC_mask],method="spearman"))
suppressWarnings(cor.test(mono_score[dc2_mask],dc_score[dc2_mask],method="spearman"))

mono_genes_plot <- names(mono_v_dc[mono_expr>10^-3.5])[order(rank(mono_v_dc[mono_expr>10^-3.5])+rank(mono_v_mac[mono_expr>10^-3.5]),decreasing =T)[1:20]]
dc_genes_plot <- names(mono_v_dc[dc2_expr>10^-3.5])[order(rank(-mono_v_dc[dc2_expr>10^-3.5])+rank(dc_v_mac[dc2_expr>10^-3.5]),decreasing =T)[1:20]]
mac_genes_plot <- names(mono_v_dc[momac_expr>10^-3.5])[order(rank(-mono_v_mac[momac_expr>10^-3.5])+rank(-dc_v_mac[momac_expr>10^-3.5]),decreasing =T)[1:20]]

cells <- colnames(lung_ldm$dataset$umitab)[mono_mask | momac_mask | moDC_mask | dc2_mask | dc1_mask]

cells <- intersect(cells,colnames(lung_ldm$dataset$ds[[1]]))


mono_genes_plot <- names(sort((mono_v_dc+mono_v_mac)[moDC_expr>1e-4 & dc2_expr < 1e-4],d=T))[1:20]



dc_genes_plot <- rev(names(sort((dc_v_mac-mono_v_dc)[mono_expr<1e-4 & dc2_expr > 1e-4],d=T))[1:20])

mac_genes_plot <- rev(names(sort((dc_v_mac+mono_v_mac)[moDC_expr>1e-4]))[1:20])


reg <- 1e-4

mono_score.n <- log2((reg+mono_score)/mean(reg+mono_score))
dc_score.n <- log2((reg+dc_score)/mean(reg+dc_score))
mac_score.n <- log2((reg+mac_score)/mean(reg+mac_score))

nper <- 1000
set.seed(910430)

cells_list <- list()
cells_list$momac <- sample(intersect(cells,colnames(lung_ldm$dataset$umitab)[momac_mask]),nper)
cells_list$mono  <- sample(intersect(cells,colnames(lung_ldm$dataset$umitab)[mono_mask ]),nper)
cells_list$moDC <- sample(intersect(cells,colnames(lung_ldm$dataset$umitab)[moDC_mask]),nper)
cells_list$dc2 <- sample(intersect(cells,colnames(lung_ldm$dataset$umitab)[dc2_mask]),nper)
cells_list$dc1 <- sample(intersect(cells,colnames(lung_ldm$dataset$umitab)[dc1_mask]),nper)

cells_list <- lapply(cells_list,function(x){
  x <- x[order(dc_score.n[x]-mono_score.n[x])]
 
  return(x)
})
cells <- unlist(cells_list)
mat <- lung_ldm$dataset$ds[[1]][c(mono_genes_plot,dc_genes_plot,mac_genes_plot),rev(cells)]
mat <- log2(1+mat)
mat[mat > 3] <- 3

mat2 <- rbind(mono_score[colnames(mat)],dc_score[colnames(mat)],mac_score[colnames(mat)])
mat2 <- log10(1e-4+mat2)


colgrad_rna <- colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100)


png(file.path(output_dir,"figure_2g.png"),
    height=3.5,width=5,units="in",res=1000,pointsize=4)
par(oma=c(0,5,0,0))
layout(matrix(1:2,ncol=2),widths=c(10,2))
image(as.matrix(mat),xaxt="n",yaxt="n",col=colgrad_rna)
mtext(rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),las=2,line=0.25,side=1)
abline(h=seq(0,4000,1000)[-1]/5000,col="black")
abline(v=cumsum(c(length(mono_genes_plot),length(dc_genes_plot)))/nrow(mat),col=rgb(t(col2rgb("green")*0.7/255)))
mtext(side=2,c("cDC1","cDC2","DC3","CD14 mono",expression(paste("MoM",Phi,sep=""))),at=c(1/10,3/10,5/10,7/10,9/10),line=0.5,cex=2,las=2)
image(as.matrix(mat2),col=colgrad_rna,xaxt="n",yaxt="n")
abline(h=seq(0,4000,1000)[-1]/5000,col="black")

dev.off()


###########################

density_lists <- list()
density_lists$cDC2_score <- list()
density_lists$cDC2_score$mono <- dc_score[mono_mask]
density_lists$cDC2_score$moDC <- dc_score[moDC_mask]
density_lists$cDC2_score$cDC2 <- dc_score[dc2_mask]
density_lists$cDC2_score$momac <- dc_score[momac_mask]
density_lists$cDC2_score$cDC1 <- dc_score[dc1_mask]


density_lists$mac_score <- list()
density_lists$mac_score$mono <- mac_score[mono_mask]
density_lists$mac_score$moDC <- mac_score[moDC_mask]
density_lists$mac_score$cDC2 <- mac_score[dc2_mask]
density_lists$mac_score$momac <- mac_score[momac_mask]
density_lists$mac_score$cDC1 <- mac_score[dc1_mask]


subtype_2_col <- rgb(t(col2rgb(c("blue","yellow","red","green","orange"))*0.7/255))
names(subtype_2_col) <- c("momac","mono","moDC","cDC2","cDC1")


png("output/figures/figure_s2e.png",height=1.3,width=2.24,units="in",bg="transparent",pointsize=6,res=1000)
layout(matrix(1:2,nrow=1))
par(mgp=c(2,1,0),mar=c(3,1,2,1))
for(list_iter in names(density_lists)){
  xlim <- log10(10^-2.8+range(unlist(density_lists[[list_iter]])))
  plot(1:10,col="white",xlim=xlim,ylim=c(0,7),ylab="",yaxt="n",bty="L",xlab="")
  for(iter in rev(seq(length(density_lists[[list_iter]])))){
    d <- density(log10(10^-2.8+density_lists[[list_iter]][[c("cDC1","cDC2","moDC","mono","momac")[iter]]]))
    polygon(d$x,iter-1+d$y,col=subtype_2_col[c("cDC1","cDC2","moDC","mono","momac")[iter]])
    abline(h=iter-1)
  }
}
dev.off()


}

