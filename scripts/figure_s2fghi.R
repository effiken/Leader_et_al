library(viridis)
library(seriation)
library(scDissector)
library(matrixStats)
library(R.devices)

parse_modules_leader <- function (ldm, cormat, ds_version, modules_version = "", nmods = 50, 
                                  reg = 1e-06, mod_size = 4, min_mod_cor = 0.1, zlim = c(-0.9, 
                                                                                         0.9),
                                  tables_path = "") 
{
  gene_mask2 = names(which(apply(cormat, 1, quantile, 1 - mod_size/ncol(cormat), 
                                 na.rm = T) >= min_mod_cor))
  ds = ldm$dataset$ds[[match(ds_version, ldm$dataset$ds_numis)]]
  mods = cutree(hclust(as.dist(1 - cormat[gene_mask2, gene_mask2])), 
                nmods)
  modsl = split(names(mods), mods)
  
  write.table(file = paste(tables_path, modules_version, "_modules.txt", 
                           sep = ""), sapply(modsl, paste, collapse = ","), row.names = T, 
              col.names = F, quote = F)
  return(modsl)
}

figure_s2fghi <- function(){

  
source(file.path(scClustering_dir,"clustering3.r"))
  
cluster_order_list <- read.table("input_tables/cluster_order_list.txt",sep="\t",r=1,h=1,stringsAsFactors = F)
clusts <- strsplit(cluster_order_list["cDC",],",")[[1]]

cells <- names(lung_ldm$dataset$cell_to_cluster)[lung_ldm$dataset$cell_to_cluster%in%clusts]
cells <- intersect(cells,colnames(lung_ldm$dataset$ds[[1]]))

set.seed(910430)
cells <- sample(cells,10000)

message("Identifying highly variable genes")

nulldev()
hvg <- get_highly_variable_genes(ds=lung_ldm$dataset$ds[[1]][,cells],
                                 seeding_genes = rownames(lung_ldm$dataset$ds[[1]]),min_umis_per_var_gene = 70,varmean_quantile = .8)
dev.off()
#Excluding antibody genes from analysis--  these are due to likely noise from batch effect
excluded_genes <- grep("^IGH|^IGK|^IGL|JCHAIN",hvg,v=T)
hvg <- setdiff(hvg,excluded_genes)

message("Computing gene-gene correlations per sample and averaging")
cmat <- get_avg_gene_to_gene_cor(lung_ldm$dataset$ds[[1]][hvg,cells],
                                 cell_to_sample = lung_ldm$dataset$cell_to_sample[cells])

message("parsing and saving modules")
modules <- parse_modules_leader(ldm=lung_ldm,cormat = cmat,ds_version = "2000", modules_version = "DC",
                                      tables_path = "output/modules/",nmods = 50,min_mod_cor = 0.06)

#modules <- read.table("tables/fig_3_modules.txt",sep="\t",stringsAsFactors = F)$V1
modules <- read.table("output/modules/DC_modules.txt",sep="\t",stringsAsFactors=F)$V1

modules <- lapply(modules,function(x){strsplit(x," ")[[1]][[2]]})

genes <- unlist(unlist(lapply(modules,strsplit,",")))

mod_mat <- matrix(0,nrow=length(modules),ncol=length(genes),dimnames=list(1:length(modules),genes))
for(iter in 1:length(modules)){
  mod_mat[iter,strsplit(modules[[iter]],",")[[1]]] <- 1
}



cells <- names(lung_ldm$dataset$cell_to_cluster)[lung_ldm$dataset$cell_to_cluster%in%clusts]

mod_scores <- t(t(mod_mat%*%lung_ldm$dataset$umitab[genes,cells])/colSums(lung_ldm$dataset$umitab[,cells]))

nmods <- length(modules)
names(modules) <- seq(nmods)
mods <- lapply(modules,function(x){strsplit(x,",")[[1]]})

#mod_exclude <- c(1,2,18)

# message("computing gene-gene correlations of highly variable genes")
# cmat <- scDissector::get_avg_gene_to_gene_cor(
#   ds = lung_ldm$dataset$ds[[1]][unlist(mods),intersect(colnames(lung_ldm$dataset$ds[[1]]),cells)],
#   cell_to_sample = lung_ldm$dataset$cell_to_sample[intersect(colnames(lung_ldm$dataset$ds[[1]]),cells)])

cmat <- cmat[unlist(mods),unlist(mods)]

#avg_cor <- array(NA,length(modules[-mod_exclude]),dimnames=list(as.character(seq(nmods)[-mod_exclude])))
avg_cor <- array(NA,length(modules),dimnames=list(as.character(seq(nmods))))
#mod_iter <- 3
for(mod_iter in names(avg_cor)){
  genes <- mods[[mod_iter]]
  mat.iter <- cmat[genes,genes]
  diag(mat.iter) <- NA
  avg_cor[mod_iter] <- mean(mat.iter,na.rm=T)
}

for(mod_iter in names(mods)){
  genes <- mods[[mod_iter]]
  mat.iter <- cmat[genes,genes]
  diag(mat.iter) <- NA
  mods[[mod_iter]] <- genes[order(rowSums(mat.iter,na.rm=T),decreasing=T)]
}

rs <- rowSums(lung_ldm$dataset$umitab[,cells])
expr <- rs[unlist(mods)]/sum(rs)
#expr <- rowMeans(lung_ldm$models[unlist(mods),])
expr[expr > 7/2000] <- 15/2000
expr <- log10(expr)
expr <- round((expr-min(expr))/diff(range(expr))*99)+1

#mod_cor <- cor(as.matrix(t(log10(1e-6+mod_scores[-mod_exclude,]))),method="spearman")
mod_cor <- cor(as.matrix(t(log10(1e-6+mod_scores))),method="spearman")
mod_ord <- hclust(as.dist(1-mod_cor))$order
mod_ord <- get_order(seriate(as.dist(1-mod_cor),method="OLO"))


mod_cor[mod_cor > -min(mod_cor)] <- -min(mod_cor)

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)

cell_to_patient <- sample_annots[lung_ldm$dataset$cell_to_sample,]$patient_ID
names(cell_to_patient) <- names(lung_ldm$dataset$cell_to_sample)
cell_to_tissue <- sample_annots[lung_ldm$dataset$cell_to_sample,"tissue"]
names(cell_to_tissue) <- names(lung_ldm$dataset$cell_to_sample)
s <- split(colnames(mod_scores),
           list(
             factor(cell_to_patient[colnames(mod_scores)]),
             factor(cell_to_tissue[colnames(mod_scores)])))

tmp <- do.call(cbind,lapply(s,function(x){rowMeans(mod_scores[,x])}))
tmp_n <- tmp[,grep("Normal",colnames(tmp))]
tmp_t <- tmp[,grep("Tumor",colnames(tmp))]
tmp_t <- tmp_t[,paste(substr(colnames(tmp_n),1,3),".Tumor",sep="")]

l2fc <- log2((tmp_t+10^-(4.5))/(tmp_n+10^(-4.5)))
#l2fc <- l2fc[-mod_exclude,]



text_spots <- seq(0,length(mod_ord)+.5)/(length(mod_ord)+1)
ngenes <- 5
at <- seq(0,1,1/(length(mod_ord)-1))

colgrad_rna <- colorRampPalette(c(colors()[378],"orange", "tomato","mediumorchid4"))(100)


message("Performing and saving statistical tests")
t.res <- matrix(NA,nrow=nrow(l2fc),ncol=1,dimnames=list(rownames(l2fc),c("p.value")))
for(row_iter in 1:nrow(t.res)){
  t.iter<- wilcox.test(l2fc[row_iter,])
  t.res[row_iter,"p.value"] <- c(t.iter$p.value)
}
t.res <- cbind(as.numeric(rownames(t.res)),t.res)

adj.p.value <- p.adjust(t.res[,"p.value"],method="bonferroni")
t.res <- as.data.frame(cbind(t.res,adj.p.value))
colnames(t.res)[1] <- "module"

write.csv(t.res,"output/statistics/figure_s2h.csv")

png("output/figures/figure_s2fgh.png",
    height=3.17,width=7.32,units="in",pointsize=6,res=300)


layout(matrix(1:4,nrow=1),widths=c(.3,3.2,2,1.6))
par(pin=c(3.1,3.1/47),mar=c(2,2,2,1))
image(t(as.matrix(log10(avg_cor[rownames(mod_cor)[rev(mod_ord)]]))),col=rev(viridis(100)),xaxt="n",yaxt="n")
mtext(side=2,rownames(mod_cor)[rev(mod_ord)],at=at,las=2,line=.5)



par(pin=c(3.1,3.1),mar=c(2,2,2,2))

image(t(mod_cor[rev(mod_ord),mod_ord]),xaxt="n",yaxt="n",col=bluered(50))
mtext(side=4,rownames(mod_cor)[rev(mod_ord)],at=at,las=2,line=.5)

image(t(mod_cor[mod_ord,rev(mod_ord)]),col="white",xaxt="n",yaxt="n",bty="n")


for(mod_iter in rownames(mod_cor)[mod_ord]){
  genes <- mods[[as.numeric(mod_iter)]][1:5]
  expr_genes <- expr[genes]
  text(x=(seq(length(genes))-1)/ngenes,y=rev(at)[match(mod_iter,rownames(mod_cor)[mod_ord])],genes,
       adj=0,xpd=F,font=2,
       col=colgrad_rna[expr_genes])#rgb(1-cex.genes,1-cex.genes,1-cex.genes))
}
h <- (at[-1]+at[-length(at)])/2
abline(h=h,col="grey",lwd=.5)

plot(NULL,xlim=c(-6,6),ylim=c(.5,nrow(l2fc)+.5),xaxs="i",yaxs="i",bty="n",xaxt="n",yaxt="n",ylab="",xlab="Log2(T vs. N)")
grid(ny=0,nx=NULL)
boxplot(t(l2fc[rev(mod_ord),]),horizontal=T,range=0,xaxs="i",yaxs="i",col="grey",xlim=c(1,nrow(l2fc)),add=T,at=seq(nrow(l2fc)),las=2)
abline(v=0)
boxplot(t(l2fc[rev(mod_ord),]),horizontal=T,range=0,col="grey",yaxt="n",xaxt="n",add=T,lty=1)
points(x=l2fc[rev(mod_ord),],y=jitter(matrix(1:nrow(l2fc),nrow=nrow(l2fc),ncol=ncol(l2fc)),.3),cex=.5)

text("*",
     x=array(-5,sum(as.numeric(t.res$adj.p.value)<5e-2 & as.numeric(t.res$adj.p.value) >=1e-2)),
     y=which(t.res[rev(mod_ord),]$adj.p.value<5e-2 & t.res[rev(mod_ord),]$adj.p.value>=1e-2))
text("**",
     x=array(-5,sum(t.res$adj.p.value<1e-2 & t.res$adj.p.value >=1e-3)),
     y=which(t.res[rev(mod_ord),]$adj.p.value<1e-2 & t.res[rev(mod_ord),]$adj.p.value>=1e-3))

text("***",
     x=array(-5,sum(t.res$adj.p.value<1e-3)),
     y=which(t.res[rev(mod_ord),]$adj.p.value<1e-3))

dev.off()

rm(list=c("cmat","modules"))

annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)

lineage_col_array <- c(rgb(t(col2rgb(c(2,3,6,"orange"))*0.7/255)))
names(lineage_col_array) <- unique(annots_list[clusts,"sub_lineage"])
s <- split(cells,lung_ldm$dataset$cell_to_cluster[cells])
mod_avgs <- do.call(cbind,lapply(s[clusts],function(x){rowMeans(mod_scores[mod_ord,x])}))

n_mod_avgs <- log2((1e-4+mod_avgs)/rowMeans(1e-4+mod_avgs))
n_mod_avgs[n_mod_avgs < -2] <- -2
n_mod_avgs[n_mod_avgs > 2] <- 2

cols <- lineage_col_array[annots_list[clusts,]$sub_lineage]


v <- which(diff(as.numeric(factor(cols)))!=0)
v <- seq(-1/length(clusts)/2,1+1/length(clusts)/2,(1+1/length(clusts))/length(clusts))[v+1]

png("output/figures/figure_s2i.png",height=1.39,width=4.49,units="in",res=300,pointsize=5,bg="transparent")
par(mar=c(0,1,0,1.5),oma=c(5,8,5,3))

layout(matrix(1:3,nrow=1),widths=c(1,10,2))
image(t(as.matrix(1:length(cols))),col=cols,xaxt="n",yaxt="n")
abline(h=v,col="white",lwd=1.5)

box()

image(n_mod_avgs,col=bluered(50),xaxt="n",yaxt="n")
mtext(colnames(n_mod_avgs),at=seq(0,1,1/(ncol(n_mod_avgs)-1)),side=2,las=2,line=.25)
mtext(rownames(n_mod_avgs),at=seq(0,1,1/(nrow(n_mod_avgs)-1)),side=1,las=2,line=0.25)
#mtext(ord,at=seq(0,1,1/(length(ord)-1)),side=2,las=2,line=.25)
abline(h=v,col="white",lwd=1.5)

box()
par(mar=c(4,5,4,3))

image(t(as.matrix(1:50)),col=bluered(50),xaxt="n",yaxt="n")
box()
mtext("Relative expr.",line=.25)
mtext(c("-2", "0", "2"),at=c(0,.5,1),side=4,las=2,line=.25,cex=0.75)

dev.off()



}






         