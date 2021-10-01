

figure_s3b <- function(){
  
  if(!exists("lung_ldm")){
    load("data/lung_ldm.rd")
  }
  
genelist1 <- read.csv("input_tables/chakarov_Table_S4.csv",stringsAsFactors=F,r=1,h=1)[,2]
genelist2 <- read.csv("input_tables/chakarov_Table_S5.csv",stringsAsFactors=F,r=1,h=1)[,2]

cells <- lung_ldm$dataset$cell_to_cluster=="8"


genelist1 <- genelist1[1:50]
genelist2 <- genelist2[1:50]

genelist2[genelist2=="C4A,C4B"] <- c("C4A")
genelist2 <- c(genelist2,"C4B")
genelist2[genelist2=="SELENOP"] <- "SEPP1"
genelist2[genelist2=="PRXL2B"] <- "FAM213B"
genelist2 <- unique(genelist2)
# 
# genelist1 <-  genelist1[!grepl("ICOSLG|FCGR3A|FCGR3B|GRK3|Sep-09|H3F3A|H3F3B|RACK1",genelist1)]
# genelist1 <- c(genelist1,"FCGR3A","FCGR3B","ICOSLG","ADRBK2","SEPT9","H3F3A","H3F3B","GNB2L1")
# 
# genelist2 <- intersect(genelist2,rownames(lung_ldm$dataset$umitab))
# genelist2 <- c(genelist2,"C4A","C4B","SEPP1","FAM213B","HBA1","HBA2","HSPA1A","HSPA1B","GRAMD3","SEP15","SEPT8","ALOX12","SHFM1","KIAA0196","NBL1")





score1 <- colSums(lung_ldm$dataset$umitab[genelist1,cells])/colSums(lung_ldm$dataset$umitab[,cells])
score2 <- colSums(lung_ldm$dataset$umitab[genelist2,which(cells)])/colSums(lung_ldm$dataset$umitab[,which(cells)])

#plot(log10(1e-4+score1),log10(1e-4+score2))
message("correlation analysis between LYVE1+ and CX3CR1+ IM gene scores")
suppressWarnings(cor.test(score1,score2,method="spearman"))

ds <- lung_ldm$dataset$ds[[1]][c(genelist1,rev(genelist2)),
                               lung_ldm$dataset$cell_to_cluster[colnames(lung_ldm$dataset$ds[[1]])]==8]
ds <- ds[,order(score1[colnames(ds)]-score2[colnames(ds)])]

mat <- log2(1+ds)
mat[mat > 3] <- 3
#image(as.matrix(t(mat)),xaxt="n",yaxt="n")
#mtext(side=2,rownames(ds),at=seq(0,1,1/(nrow(ds)-1)),las=2)

colgrad_rna <- colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100)

genes_to_plot <- c(genelist1[1:15],genelist2[1:15])

png("output/figures/figure_s3b.png",height=2.6,width=2.41,units="in",res=1000,pointsize=5)
par(mar=c(5.1,6,4.1,2.1))
image(as.matrix(t(mat[genes_to_plot,])),xaxt="n",yaxt="n",col=colgrad_rna)
mtext(side=2,genes_to_plot,at=seq(0,1,1/(length(genes_to_plot)-1)),las=2,font=3,cex=0.7,line=0.25)
abline(h=0.5,col="black")
box()
mtext(c(expression(paste("LYVE1+ IM",Phi,sep="")),expression(paste("CX3CR1+ IM",Phi,sep=""))),at=c(0.15,0.85),side=3,cex=1.2)
mtext(c(expression(paste("top CX3CR1+ IM",Phi," genes",sep="")),expression(paste("top LYVE1+ IM",Phi," genes",sep=""))),side=2,line=4,cex=0.7,at=c(0.25,0.75))

dev.off()

}