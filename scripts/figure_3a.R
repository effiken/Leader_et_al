
library(Matrix)
library(viridis)
library(seriation)

figure_3a <- function(lung_ldm){
  
  
  sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
  annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)
  cluster_order_list <- read.table("input_tables/cluster_order_list.txt",sep="\t",r=1,h=1,stringsAsFactors = F)
  genes_in_figures <- read.csv("input_tables/table_s5_genes_used_in_figures.csv",r=1,h=1,stringsAsFactors = F)
  
  cluster_order <- as.character(strsplit(cluster_order_list["Mono_mac",],",")[[1]])
  genes <- strsplit(genes_in_figures["Figure 3A",2],",")[[1]]
  
  
  lineage_col_array <- c("red","green",rgb(t(col2rgb(c(6,
                                                       4,rgb(0,.5,1),rgb(0,.8,1),5,"yellow","orange","red","green"
  ))*0.7/255)))
  names(lineage_col_array) <- unique(annots_list[cluster_order,"sub_lineage"])
  
  
  s <- split(colnames(lung_ldm$dataset$umitab),lung_ldm$dataset$cell_to_cluster)

  clust_avgs <- do.call(cbind,lapply(s[cluster_order],function(x){rs <- rowSums(lung_ldm$dataset$umitab[,x]); rs <- rs[genes]/sum(rs); return(rs)}))


cols <- lineage_col_array[annots_list[cluster_order,]$sub_lineage]

mat <- log2((1e-6+clust_avgs)/(1e-6+rowMeans(clust_avgs)))
mat[mat < -2] <- -2
mat[mat > 2] <- 2
mat <- t(mat)

v <- which(diff(as.numeric(factor(cols)))!=0)
v <- seq(-1/length(cluster_order)/2,1+1/length(cluster_order)/2,(1+1/length(cluster_order))/length(cluster_order))[v+1]

png("output/figures/figure_3a.png",height=2.35,width=4.54,units="in",res=300,pointsize=5,bg="transparent")

layout(matrix(1:4,nrow=2),heights=c(1,10),widths=c(10,4))
par(mar=c(0,0,0,0),oma=c(5,8,5,3))
image(as.matrix(1:length(cols)),col=cols,xaxt="n",yaxt="n")
abline(v=v,col="white")

box()

image(mat,col=bluered(50),xaxt="n",yaxt="n")
mtext(colnames(mat),at=seq(0,1,1/(length(genes)-1)),side=2,las=2,line=.25)
#mtext(ord,at=seq(0,1,1/(length(ord)-1)),side=2,las=2,line=.25)
abline(v=v,col="white")

box()

plot.new()

par(pin=c(.1,.5))
image(t(as.matrix(1:50)),col=bluered(50),xaxt="n",yaxt="n")
box()
mtext("Relative expr.",line=.25)
mtext(c("-2", "0", "2"),at=c(0,.5,1),side=4,las=2,line=.25,cex=0.75)

dev.off()


}
