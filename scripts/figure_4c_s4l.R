library(scales)


figure_4c_s4l <- function(lung_ldm){
  

# set appropriate directories

output_dir <- "output/figures/"

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)

patient_tissue <- apply(sample_annots[lung_ldm$dataset$cell_to_sample,c("patient_ID","tissue")],1,paste,collapse="_")

sub_lineage <- annots_list[lung_ldm$dataset$cell_to_cluster,]$sub_lineage

tab <- table(sub_lineage,patient_tissue)
#tab <- tab[annots_list[match(rownames(tab),annots_list$sub_lineage),]$lineage=="MNP",]

tab <- tab[c("NK","T_Nklike","T_GZMK","CD8 Trm","T_activated","CD4 Trm","Tcm/naive_II","Tcm/naive_I","Treg","Tcycle"),]
tab <- t(t(tab)/rowSums(t(tab)))
tab_normal <- tab[,grepl("Normal",colnames(tab))]
tab_tumor <- tab[,grep("Tumor",colnames(tab))]
colnames(tab_normal) <- substr(colnames(tab_normal),1,3)
colnames(tab_tumor) <- substr(colnames(tab_tumor),1,3)
tab_normal <- tab_normal[,intersect(colnames(tab_normal),colnames(tab_tumor))]
tab_tumor <- tab_tumor[,intersect(colnames(tab_normal),colnames(tab_tumor))]

delta <- log2((1e-2+tab_tumor)/(1e-2+tab_normal))
mat <- delta[c("NK","T_Nklike","T_GZMK","CD8 Trm","T_activated","CD4 Trm","Tcm/naive_II","Tcm/naive_I","Treg","Tcycle"),]

res <- apply(mat,1,function(x){suppressWarnings(wilcox.test(x))})

statistic <- unlist(lapply(res,function(x){x$statistic}))
pval <- unlist(lapply(res,function(x){x$p.value}))
adj_pval <- p.adjust(pval,method="bonferroni")
res_mat <- cbind(statistic,pval,adj_pval)
rownames(res_mat) <- rownames(mat)
write.csv(res_mat,"output/statistics/figure_4c.csv")

sig_vec <- p.adjust(lapply(res,function(x){x$p.value}),method="bonferroni")

cols <- c(rgb(t(col2rgb(c(2,3,4,5,6,7,8
))*0.7/255)),1,2,3)

png(file.path(output_dir,"figure_4c.png"),height=4,width=7.59,units="in",res=300,pointsize=12)

boxplot.matrix(t(mat),range=0,col=cols,ylim=c(min(mat),max(mat)+1),yaxt="n",xaxt="n")
axis(side=2,at=seq(-3,round(max(mat)+1),2),las=2)
x <- jitter(matrix(1:10,nrow=10,ncol=ncol(delta)))
points(x,mat,col=alpha(1,0.5),pch=16)
abline(h=0,col="red")
if(sum(sig_vec<5e-2 & sig_vec>1e-2)>0){
  text(x=(1:10)[sig_vec<5e-2 & sig_vec>1e-2],y=array(max(mat)+0.5,sum(sig_vec<5e-2 & sig_vec>1e-2)),"*",cex=2)
}
if(sum(sig_vec<1e-2 & sig_vec>1e-3)>0){
  text(x=(1:10)[sig_vec<1e-2 & sig_vec>1e-3],y=array(max(mat)+0.5,sum(sig_vec<1e-2 & sig_vec>1e-3)),"**",cex=2)
}
if(sum(sig_vec<1e-3)>0){
  text(x=(1:10)[sig_vec<1e-3],y=array(max(mat)+0.5,sum(sig_vec<1e-3)),"***",cex=2)
}
dev.off()


tab <- table(sub_lineage,patient_tissue)
#tab <- tab[annots_list[match(rownames(tab),annots_list$sub_lineage),]$lineage=="MNP",]

tab <- tab[c("B","IgD","IgM","IgA","IgG"),]
tab <- t(t(tab)/rowSums(t(tab)))
tab_normal <- tab[,grepl("Normal",colnames(tab))]
tab_tumor <- tab[,grep("Tumor",colnames(tab))]
colnames(tab_normal) <- substr(colnames(tab_normal),1,3)
colnames(tab_tumor) <- substr(colnames(tab_tumor),1,3)
tab_normal <- tab_normal[,intersect(colnames(tab_normal),colnames(tab_tumor))]
tab_tumor <- tab_tumor[,intersect(colnames(tab_normal),colnames(tab_tumor))]

delta <- log2((1e-2+tab_tumor)/(1e-2+tab_normal))
mat <- delta[c("B","IgD","IgM","IgA","IgG"),]

res <- apply(mat,1,function(x){suppressWarnings(wilcox.test(x))})

statistic <- unlist(lapply(res,function(x){x$statistic}))
pval <- unlist(lapply(res,function(x){x$p.value}))
adj_pval <- p.adjust(pval,method="bonferroni")
res_mat <- cbind(statistic,pval,adj_pval)
rownames(res_mat) <- rownames(mat)
write.csv(res_mat,"output/statistics/figure_s4l.csv")



sig_vec <- p.adjust(lapply(res,function(x){x$p.value}),method="bonferroni")

cols <- rgb(t(col2rgb(c(2,3,4,5,6
))*0.7/255))

png(file.path(output_dir,"figure_s4l.png"),height=4,width=7.59,units="in",res=300,pointsize=12)

boxplot.matrix(t(mat),range=0,col=cols,ylim=c(min(mat),max(mat)+1),yaxt="n",xaxt="n")
axis(side=2,at=seq(-3,round(max(mat)+1),2),las=2)
x <- jitter(matrix(1:5,nrow=5,ncol=ncol(delta)))
points(x,mat,col=alpha(1,0.5),pch=16)
abline(h=0,col="red")
if(sum(sig_vec<5e-2 & sig_vec>1e-2)>0){
  text(x=(1:5)[sig_vec<5e-2 & sig_vec>1e-2],y=array(max(mat)+0.5,sum(sig_vec<5e-2 & sig_vec>1e-2)),"*",cex=2)
}
if(sum(sig_vec<1e-2 & sig_vec>1e-3)>0){
  text(x=(1:5)[sig_vec<1e-2 & sig_vec>1e-3],y=array(max(mat)+0.5,sum(sig_vec<1e-2 & sig_vec>1e-3)),"**",cex=2)
}
if(sum(sig_vec<1e-3)>0){
  text(x=(1:5)[sig_vec<1e-3],y=array(max(mat)+0.5,sum(sig_vec<1e-3)),"***",cex=2)
}
dev.off()
}
