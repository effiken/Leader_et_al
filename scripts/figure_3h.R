

library(Matrix)
library(viridis)
library(seriation)
library(matrixStats)

figure_3h <- function(){

plotting_dir <- "output/figures/"

if(!exists("lung_ldm")){
  message("Loading scRNAseq data into R")
  load("data/lung_ldm.rd")
}

annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)
annots_ord <- c("IM","AM","AZU1_mac","MoMac-IV","MoMac-III","MoMac-II","MoMac-I","DC3","CD14 mono","cDC2","mregDC","cDC1","CD16 mono")

s <- split(colnames(lung_ldm$dataset$umitab),factor(annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"],annots_ord))

clust_avgs <- do.call(cbind,lapply(s,function(x){rs <- rowSums(lung_ldm$dataset$umitab[,x]); rs <- rs/sum(rs); return(rs)}))


lig_rec <- read.csv("input_tables/ligand_receptor_ramilowski2015.csv",r=1,h=1,stringsAsFactors = F)
ligs <- unique(lig_rec$Ligand.ApprovedSymbol)
ligs["IL8"] <- "CXCL8"
ligs <- intersect(ligs,rownames(lung_ldm$dataset$umitab))
ligs <- ligs[which(rowSums(clust_avgs[ligs,]>1e-5)>0)]

exclude_prefixes <- c("ADA","ANX","B2M","C1Q","CAL","CD1","CD2","CD4","CD5","CFH","CSF","LAM","TFP","LYPD","HSP","GRE","SPO","SFT","SCG","SLP","VWF","CTG","CYR","GZM","IFN",
                      "PIG","NUC","PON","HMG","PPB","MFG","SAA","BST","PF4","SEL","HLA","FAR","LTB","SER","TCT","MDK","POM","GNAS","HEB","ARF","PDA","MLL","QDP","RTN","GPI",
                      "HRA","LIN","PSE","GNA","APP","GST","UBA","RPS","MFN","LTB","FLT","NAM","TGF","ZG1","SPI","ICA","NRG","S10","VAS")
exclude <- unlist(lapply(as.list(exclude_prefixes),function(x){grep(x,ligs,v=T)}))
exclude <- c(exclude,"IL16","ICAM3","TNFSF8","CXCL13","YARS","ITGB3BP","PTMA","PTDSS1","DUSP18","PROK2","SEMA4D","SYTL3","CCL5","TNFSF14","IGF1","ICAM2","GAS6","HBEGF","COL4A1","RNASE2","CLEC11A","CLCF1","COL9A2",
             "TGM2","LRPAP1")
ligs <- setdiff(ligs,exclude)


ligs <- ligs[rowMaxs(clust_avgs[ligs,annots_ord])>10^-4]

ligs <- c(ligs,"TNFSF13","CXCL11","IL6","IL1A","IL1B","IL23A")
ligs <- unique(ligs)

mat <- clust_avgs[ligs,rev(annots_ord)]


mat <- log2((1e-6+mat)/(1e-6+rowMeans(mat)))
thresh <- 2
mat[mat < -thresh] <- -thresh
mat[mat > thresh] <- thresh

# 
#d <- as.dist(1-cor((mat),method="spearman"))
# hc <- hclust(d)
#d1 <- as.dist(1-cor(t(mat),method="spearman"))
#hc1 <- hclust(d1)
 # 
# col_ord <- reorder(hc,d,method="OLO")$order
# mat <- mat[,col_ord]
lig_ord <- rev(order(apply(mat,1,which.max)))
#mat <- mat[hc1$order,hc$order]
mat <- mat[lig_ord,]

png(file.path(plotting_dir,"figure_3h.png"),height=2,width=4.02,units="in",res=300,pointsize=4)
layout(matrix(1:2,nrow=2),heights=c(3,1))
par(oma=c(0,1,0,0),mar=c(5,4,1,1))
image(mat,col=bluered(50),xaxt="n",yaxt="n")
mtext(rownames(mat),at=seq(0,1,1/(nrow(mat)-1)),side=1,las=2,line=.25)
mtext(colnames(mat),at=seq(0,1,1/(ncol(mat)-1)),side=2,las=2,line=.25)
box()
par(pin=c(0.8,0.15))
image(matrix(1:100,100),col=bluered(100),xaxt="n",yaxt="n")
mtext(at=c(0,0.5,1),c(-2,0,2),side=1,line=0.25)
mtext("Normalized expression")
box()
dev.off()
}


