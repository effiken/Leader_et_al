library(matrixStats)
library(skmeans)
library(scales)

figure_s4fgh <- function(){

rm(list=ls())

output_dir <- "output/figures/"

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)

load("data/lung_ldm.rd")

v2_mask <- sample_annots[lung_ldm$dataset$cell_to_sample,]$library_chemistry=="V2"
cell_to_annot <- annots_list[lung_ldm$dataset$cell_to_cluster,"sub_lineage"]
names(cell_to_annot) <- names(lung_ldm$dataset$cell_to_cluster)

#compute the l2fc of each cluster to the max of the other clusters

s <- split(names(cell_to_annot)[v2_mask],cell_to_annot[v2_mask])
s <- s[unique(annots_list[annots_list$lineage=="T",]$sub_lineage)]
annot_avgs <- lapply(s,function(x){rs <- rowSums(lung_ldm$dataset$umitab[,x]); return(rs/sum(rs))})
t_avg <- rowSums(lung_ldm$dataset$umitab[,unlist(s)]); t_avg <- t_avg/sum(t_avg)
annot_avgs <- do.call(cbind,annot_avgs)

#gene lists are based on genes that are higher in one cluster than in any of the others
l2fc <- matrix(NA,nrow=nrow(annot_avgs),ncol=ncol(annot_avgs),dimnames=dimnames(annot_avgs))
for(col in colnames(annot_avgs)){
  l2fc[,col] <- log2((1e-6+annot_avgs[,col])/rowMaxs(1e-6+annot_avgs[,setdiff(colnames(annot_avgs),col)]))
}

# generate gene signatures for each cluster based on those
gene_lists <- list()
for(iter in colnames(l2fc)){

  gene_lists[[iter]] <- rownames(l2fc)[l2fc[,iter]> 0.25 & annot_avgs[,iter]>1e-5]
}

exclude <- c(grep("^MT",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^IGH",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^RP",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^SFT",rownames(lung_ldm$dataset$umitab),v=T),"XIST",
             grep("^IGL",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^IGK",rownames(lung_ldm$dataset$umitab),v=T),"JCHAIN",
             grep("^TRAV",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^TRAD",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^TRAJ",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^TRBV",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^TRBD",rownames(lung_ldm$dataset$umitab),v=T),
             grep("^TRBJ",rownames(lung_ldm$dataset$umitab),v=T))

gene_lists <- lapply(gene_lists,function(x){x <- setdiff(x,exclude)})

#estimate the # of UMI that each cell has that is not related to the cycling score
numi <- colSums(lung_ldm$dataset$umitab)
cycling_score <- colSums(lung_ldm$dataset$umitab[gene_lists[["Tcycle"]],])/numi
numi_no_cycle <- numi - colSums(lung_ldm$dataset$umitab[gene_lists[["Tcycle"]],]) - colSums(lung_ldm$dataset$umitab[exclude,])

#generate gene signatures that are completely non-overlapping with the cycling signature
gene_lists <- lapply(gene_lists,function(x){x <- setdiff(x,gene_lists[["Tcycle"]])})
gene_lists$Tnaive <- unique(gene_lists$`Tcm/naive_I`,gene_lists$`Tcm/naive_II`)
gene_lists <- gene_lists[setdiff(names(gene_lists),c("Tcm/naive_I","Tcm/naive_II"))]
# score all the T cells using those gene scores
gene_list_mat <- matrix(0,nrow=length(gene_lists),ncol=length(unique(unlist(gene_lists))),
                        dimnames=list(names(gene_lists),unique(unlist(gene_lists))))
for(iter in names(gene_lists)){
  gene_list_mat[iter,gene_lists[[iter]]] <- 1
}

t_cell_mask <- annots_list[lung_ldm$dataset$cell_to_cluster,]$lineage=="T"
gene_list_scores <- t(t(gene_list_mat%*%lung_ldm$dataset$umitab[colnames(gene_list_mat),t_cell_mask])/numi_no_cycle[t_cell_mask])

set.seed(910430)

mask <- sample(colnames(gene_list_scores),20000)

# png(file.path(output_dir,"cycling_scores_by_subtype.png"),height=5,width=5,pointsize=6,units="in",res=300)
# layout(matrix(1:9,nrow=3))
# for(iter in names(gene_lists)){
# 
#   plot(1e-4+gene_list_scores[iter,mask],1e-4+cycling_score[mask],log="xy",pch=".",col="grey",main=iter)
#   points(1e-4+gene_list_scores[iter,mask[cell_to_annot[mask]==iter]],1e-4+cycling_score[mask[cell_to_annot[mask]==iter]],
#          pch=".")
#   points(1e-4+gene_list_scores[iter,mask[cell_to_annot[mask]=="Tcycle"]],1e-4+cycling_score[mask[cell_to_annot[mask]=="Tcycle"]],
#          pch=".",col="red")
# }
# dev.off()

cell_to_tissue <- sample_annots[lung_ldm$dataset$cell_to_sample,"tissue"]
cell_to_patient <- sample_annots[lung_ldm$dataset$cell_to_sample,"patient_ID"]
names(cell_to_tissue)=names(cell_to_patient) <- names(lung_ldm$dataset$cell_to_sample)

mask <- cell_to_annot=="Tcycle"
cells <- names(cell_to_annot)[mask]

dat <- as.matrix(gene_list_scores[,cells])
dat <- t(t(dat)/rowSums(t(dat)))
dat <- dat[setdiff(rownames(dat),"Tcycle"),]

k <- skmeans(t(dat),7)

k_ord <- c(4,1,5,2,3,7,6)
row_ord <- rev(c("T_Nklike","T_GZMK","CD8 Trm","T_activated","CD4 Trm","Tnaive","Treg"))

mat <- log2((1e-2+dat)/rowMeans((1e-2+dat)))

mat[mat > 2] <- 2
mat[mat < -2] <- -2

clust_lines <- cumsum(table(factor(k$cluster,k_ord)))/length(k$cluster)
clust_lines <- clust_lines[-7]
clust_ctrs <- 1/2*(c(0,clust_lines)+c(clust_lines,1))

png(file.path(output_dir,"figure_s4f.png"),height=4,width=6,units="in",res=1000)
layout(matrix(1:2,ncol=2),widths = c(10,2))
par(oma=c(0,2,0,0),mar=c(1,5,2,2))
image(t(mat[row_ord,order(factor(k$cluster,k_ord))]),col=bluered(50),xaxt="n",yaxt="n")
abline(v=clust_lines)
mtext(side=2,at=seq(0,1,1/(nrow(mat)-1)),row_ord,las=2,line=0.25)
mtext(at=clust_ctrs,seq(max(k_ord)),line=0.25)
box()
par(mar=c(1,1,1,1),pin=c(0.25,1))
image(t(1:50),col=bluered(100),xaxt="n",yaxt="n"); box()
mtext(side=4,c(-2,0,2),at=c(0,.5,1),las=2,line=0.25)
dev.off()

# 
# 
#barplot(table(cell_to_tissue[cells],factor(k$cluster[cells],k_ord)),beside=T,col=c("blue","brown"),las=2)

patient_tissue <- paste(cell_to_patient,cell_to_tissue,sep="_")
names(patient_tissue) <- names(cell_to_patient)

cycle_tab <- table(patient_tissue[cells],k$cluster[cells])

tcells <- names(lung_ldm$dataset$cell_to_cluster)[annots_list[lung_ldm$dataset$cell_to_cluster,"lineage"]=="T"]
#tcells <- setdiff(tcells,cells)

non_cycle_tab <- table(patient_tissue[tcells],annots_list[lung_ldm$dataset$cell_to_cluster[tcells],"sub_lineage"])

cycle_annots <- c("nklike_gzmk","cd8_trm","activated","cd4_trm","cm_naive","treg")

#cycle_tab <- cycle_tab[,c(5,1,4,6,2,3)]
cycle_tab <- cycle_tab[,-5]
colnames(cycle_tab) <- cycle_annots

#wrangle non-cycle_tab into similar matrix
non_cycle_tab_annots <- matrix(NA,nrow=nrow(cycle_tab),ncol=ncol(cycle_tab),dimnames=dimnames(cycle_tab))
non_cycle_tab_annots[,"nklike_gzmk"] <- rowSums(non_cycle_tab[rownames(non_cycle_tab_annots),c("T_GZMK","T_Nklike")])
non_cycle_tab_annots[,"cd8_trm"] <- rowSums(non_cycle_tab[rownames(non_cycle_tab_annots),c("CD8 Trm"),drop=F])
non_cycle_tab_annots[,"activated"] <- rowSums(non_cycle_tab[rownames(non_cycle_tab_annots),c("T_activated"),drop=F])
non_cycle_tab_annots[,"cd4_trm"] <- rowSums(non_cycle_tab[rownames(non_cycle_tab_annots),c("CD4 Trm"),drop=F])
non_cycle_tab_annots[,"cm_naive"] <- rowSums(non_cycle_tab[rownames(non_cycle_tab_annots),c("Tcm/naive_I","Tcm/naive_II"),drop=F])
non_cycle_tab_annots[,"treg"] <- rowSums(non_cycle_tab[rownames(non_cycle_tab_annots),c("Treg"),drop=F])

cycle_tab <- cycle_tab/rowSums(non_cycle_tab[rownames(cycle_tab),])


non_cycle_tab_annots_normed <- non_cycle_tab_annots/rowSums(non_cycle_tab_annots[rownames(cycle_tab),])
non_cycle_tab_annots_normed[non_cycle_tab_annots<50] <- NA


delta <- 100*((1e-3+cycle_tab)/(1e-3+cycle_tab+non_cycle_tab_annots_normed))

delta_tumor <- delta[grepl("Tumor",rownames(delta)),]
delta_normal <- delta[grepl("Normal",rownames(delta)),]
x_normal <- c(1+3*(0:(ncol(delta)-1)))
x_tumor <- c(2+3*(0:(ncol(delta)-1)))

png(file.path(output_dir,"figure_s4h.png"),height=5,width=5,units="in",res=300)
par(oma=c(3,0,0,5), bty="L")

boxplot.matrix(delta_normal,at=x_normal,las=2,range=0,xlim=c(0.5,max(x_tumor)+0.5),col=alpha("blue",0.4),xaxt="n",log="y",ylim=range(delta,na.rm=T),yaxt="n")
boxplot.matrix(delta_tumor,at=x_tumor,las=2,range=0,add=T,col=alpha("brown",0.4),xaxt="n",yaxt="n")
axis(side=2,at=c(0.1,0.3,1,5,20),las=2)
abline(h=c(0.1,0.3,1,5,20),col="grey",lty=2)

points(jitter(matrix(x_normal,byrow=T,ncol=ncol(delta_normal),nrow=nrow(delta_normal)),0.5),delta_normal,pch=16,cex=0.5)
points(jitter(matrix(x_tumor,byrow=T,ncol=ncol(delta_tumor),nrow=nrow(delta_tumor)),0.5),delta_tumor,pch=16,cex=0.5)

mtext(side=1,at=1/2*(x_normal+x_tumor),c("NK-like/GZMK","CD8 Trm","Activated","CD4 Trm","Tcm/Naive-like","Treg"),las=2,line=0.5)
mtext(side=2,line=3,"% cycling of celltype")

legend(x=max(x_tumor)+1,y=max(delta,na.rm=T),legend=c("nLung","Tumor"),fill=alpha(c("blue","brown"),0.4),xpd=NA,bty="n")

dev.off()
###################### continue revising from here


# 
cell_to_patient <- sample_annots[lung_ldm$dataset$cell_to_sample,]$patient_ID
names(cell_to_patient) <- names(lung_ldm$dataset$cell_to_sample)
# table(cell_to_patient[cells],k$cluster)



#############
#settings for Total cell # barplot
tissue_clust_numbers <- list()
tissue_ord <- c("Normal","Tumor")
clusts_ord <- k_ord
tissue_iter <- "Normal"
for(tissue_iter in tissue_ord){
  tissue_mask <- cell_to_tissue[cells]==tissue_iter
  tab <- table(cell_to_patient[cells][tissue_mask],k$cluster[cells][tissue_mask])
  tissue_clust_numbers[[tissue_iter]] <- matrix(0,nrow=length(unique(cell_to_patient)),ncol=length(clusts_ord),
                                                dimnames=list(unique(cell_to_patient),clusts_ord))
  tissue_clust_numbers[[tissue_iter]][rownames(tab),colnames(tab)] <- tab
}

#settings for cell frequency boxplot
tissue_clust_freqs <- lapply(tissue_clust_numbers,function(x){x/rowSums(x)})


png(file.path(output_dir,"figure_s4g.png"),height=5,width=5,units="in",res=300)
par(oma=c(3,0,0,5), bty="L")

n_tissues <- length(tissue_ord)
n_clusts <- length(clusts_ord)
ymax <- ceiling(max(unlist(lapply(tissue_clust_numbers,colSums)))/100)*100
tissue_cols <- alpha(c("red","blue","brown"),.4)
names(tissue_cols) <- c("PBMC","Normal","Tumor")
xmax <- (1+n_tissues)*n_clusts
for(tissue_iter in tissue_ord){
  add <- tissue_iter!=tissue_ord[1]
  space <- c(.5+1*(match(tissue_iter,tissue_ord)-1),array(n_tissues,n_clusts-1))
  barplot(tissue_clust_numbers[[tissue_iter]],
          space=space,add=add,ylim=c(0,ymax*1.1),xlim=c(0,xmax),
          col=tissue_cols[tissue_iter],
          horiz=F,names.arg=NULL,xaxt="n",
          xaxs="i",xaxs="i",yaxs="i")
}
#axis(side=2,at=c(0,5000,10000,15000,20000),las=2)
box(lwd=1)
#abline(v=clust_breaks*xmax,col="grey",lwd=.1)
mtext(expression("Total cells observed"),side=2,line=3)
#axis(side=1,at=1:30)
mtext(at=c(1.5+3*0:6),c("NK-like/GZMK","CD8 Trm","Activated","CD4 Trm","mix","Tcm/Naive-like","Treg"),side=1,las=2,line=0.75)
dev.off()




}
