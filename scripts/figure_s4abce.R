

library(scDissector)
library(scTools)
library(RColorBrewer)
library(scales)

figure_s4abce <- function(){

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)

if(!exists("lung_ldm")){
  load("data/lung_ldm.rd")
}

adt_list <- lung_ldm$dataset$adt_by_sample[unlist(lapply(lung_ldm$dataset$adt_by_sample,length))!=0]
adtmat <- adt_list_to_matrix(adt_list)
adtmat <- adtmat[,lung_ldm$dataset$cell_to_cluster[colnames(adtmat)]=="44"]
adtmat <- adtmat[,!is.na(adtmat["CD4",])]

tfh <- log2((5+adtmat["CD8",])/(5+adtmat["CD4",])) < -1

texh <- log2((5+adtmat["CD8",])/(5+adtmat["CD4",])) > 1

cd8_v_4_score <- log2((5+adtmat["CD8",])/(5+adtmat["CD4",]))

tfh <- colnames(adtmat)[tfh]
texh <- colnames(adtmat)[texh]

tfh_train <- tfh[sample_annots[lung_ldm$dataset$cell_to_sample[tfh],"patient_ID"]%in%c("706","695")]
texh_train <- texh[sample_annots[lung_ldm$dataset$cell_to_sample[texh],"patient_ID"]%in%c("706","695")]

tfh_test <- setdiff(tfh,tfh_train)
texh_test <- setdiff(texh,texh_train)

tfh_exprs <- rowSums(lung_ldm$dataset$umitab[,tfh_train])
tfh_exprs <- tfh_exprs/sum(tfh_exprs)
texh_exprs <- rowSums(lung_ldm$dataset$umitab[,texh_train])
texh_exprs <- texh_exprs/sum(texh_exprs)

m <- (texh_exprs+tfh_exprs)/2
l2fc <- log2((1e-6+tfh_exprs)/(1e-6+texh_exprs))

png("output/figures/figure_s4a.png",height=4,width=4,units="in",res=300)
plot(log10(1e-6+m),l2fc,pch=".")
points(x=log10(1e-6+m[c("CD4","CD40LG","BCL6","CXCL13","ENTPD1","ITGA1","CD8A")]),y=l2fc[c("CD4","CD40LG","BCL6","CXCL13","ENTPD1","ITGA1","CD8A")],
       pch=16,col="red")
text(c("CD4","CD40LG","BCL6","CXCL13","ENTPD1","ITGA1","CD8A"),x=log10(1e-6+m[c("CD4","CD40LG","BCL6","CXCL13","ENTPD1","ITGA1","CD8A")]),y=l2fc[c("CD4","CD40LG","BCL6","CXCL13","ENTPD1","ITGA1","CD8A")],
     adj=1)
dev.off()

tfh_genes <- names(m)[m > 1e-4 & l2fc > 1]
texh_genes <- names(m)[m > 1e-4 & l2fc < -1]

genes_exclude <- c("STMN1",grep("^AC0",names(m),v=T),grep("^LINC",names(m),v=T),grep("^HSP",names(m),v=T),
                   grep("^IG",names(m),v=T),"JCHAIN",
                   grep("^MT",names(m),v=T),grep("^RP",names(m),v=T),grep("^SCG",names(m),v=T),"XIST",grep("^HIST",names(m),v=T),grep("^CTD-",names(m),v=T),
                   grep("^SCG",names(m),v=T))


tfh_genes <- setdiff(tfh_genes,genes_exclude)
texh_genes <- setdiff(texh_genes,genes_exclude)

clust_mat <- lung_ldm$dataset$umitab[,lung_ldm$dataset$cell_to_cluster=="44"]

tfh_score <- colSums(clust_mat[tfh_genes,])/colSums(clust_mat)
texh_score <- colSums(clust_mat[texh_genes,])/colSums(clust_mat)

accuracy <- function(x,cells){
  cite_array <- cd8_v_4_score[cells]
  gene_score <- log2((1e-3+tfh_score[cells])/(1e-3+texh_score[cells]))
  mask <- abs(cite_array) > 1
  cite_array <- cite_array[mask]
  gene_score <- gene_score[mask]
  
  tab <- table(cite_array > 0,gene_score < x)
  return((tab[1,1] + tab[2,2])/sum(tab))
}
cutoff <- seq(-4,4,0.1) 
acc_vec <- array(NA,length(cutoff))
for(cut in cutoff){
  acc_vec[which(cutoff==cut)] <- accuracy(cut,cells=c(tfh_train,texh_train))
}

#plot(cutoff,acc_vec,xlab="Score ratio Cutoff",ylab="Accuracy")

cut <- cutoff[which.max(acc_vec)]
train_acc <- accuracy(cut,cells=c(tfh_train,texh_train))
test_acc <- accuracy(cut,cells=c(tfh_test,texh_test))

my_cells <- c(tfh,texh)

png("output/figures/figure_s4b.png",height=4,width=4,units="in",res=300)
plot(log2((tfh_score[my_cells]+1e-3)/(texh_score[my_cells]+1e-3)),-cd8_v_4_score[my_cells],col=alpha(1,.4),pch=c(1,16)[1+as.integer(my_cells%in%c(tfh_test,texh_test))],
     xlab="Gene score ratio",ylab="CITEseq CD4/8 ratio")
abline(v=cut,col="red")
legend("bottomright",
       legend=c(paste(c("training data; accuracy = ","test data; accuracy="),signif(c(train_acc,test_acc),2),sep=""),
                "descriminant line"),pch=c(1,16,NA),lty=c(0,0,1),col=c(alpha(1,.4),alpha(1,.4),2),cex=0.5)
mtext("CD4/8 classifier")
dev.off()


cd4_only <- colnames(lung_ldm$dataset$umitab)[lung_ldm$dataset$umitab["CD4",]>0 &
                                                colSums(lung_ldm$dataset$umitab[c("CD8A","CD8B"),])==0 &
                                                lung_ldm$dataset$cell_to_cluster=="44"]
cd8_only <- colnames(lung_ldm$dataset$umitab)[lung_ldm$dataset$umitab["CD4",]==0 &
                                                colSums(lung_ldm$dataset$umitab[c("CD8A","CD8B"),])>0 &
                                                lung_ldm$dataset$cell_to_cluster=="44"]

clust_mat <- lung_ldm$dataset$umitab[,c(cd4_only,cd8_only)]

tfh_score <- colSums(clust_mat[tfh_genes,])/colSums(clust_mat)
texh_score <- colSums(clust_mat[texh_genes,])/colSums(clust_mat)
score_ratio <- log2((1e-3+tfh_score)/(1e-3+texh_score))
acc <- sum(diag(table(score_ratio>cut,names(score_ratio) %in%cd4_only)))/length(score_ratio)

png("output/figures/figure_s4c.png",height=4,width=4,units="in",res=300)
plot(log10(1e-3+tfh_score),log10(1e-3+texh_score),col=alpha(3+as.integer(names(tfh_score)%in%cd4_only),0.4),pch=16,cex=0.5,
     xlab="CD4 gene score",ylab="CD8 gene score",
     xlim=c(-2.5,-1),ylim=c(-2.3,-1.2))
abline(c(cut*log10(2),1),col="red")
legend("bottomleft",legend=c("CD8 UMIs only","CD4 UMIs only",paste("descriminant line; accuracy=",signif(acc,2),sep="")),col=c(3,4,2),pch=c(16,16,NA),
       lty=c(NA,NA,1),cex=0.5)
mtext("Tactivated cells without CITEseq data")
dev.off()


###########

#load("intermediates/lung_ldm_metadata_200327.rd")


total_cells <- table(apply(sample_annots[lung_ldm$dataset$cell_to_sample,c("patient_ID","tissue")],1,paste,collapse="_"),
                     annots_list[lung_ldm$dataset$cell_to_cluster,]$sub_lineage)
total_cells <- total_cells[,-1]
total_cells <- total_cells[grepl("Tumor",rownames(total_cells)),]
rownames(total_cells) <- unlist(lapply(strsplit(rownames(total_cells),"_"),function(x){x[[1]]}))
total_cells <- rowSums(total_cells)

# 
# 
# metadata$tissue <- sample_annots[as.character(metadata$sample_id),]$tissue
# metadata$patient <- sample_annots[as.character(metadata$sample_id),]$patient_ID
# metadata$lineage <- annots_list[as.character(metadata$cluster_assignment),"lineage"]
# 
# metadata <- metadata[metadata$tissue=="Tumor",]
# metadata <- metadata[metadata$lineage!="epi_endo_fibro_doublet",]

#total_cells <- table(metadata$patient)

clust_mat <- lung_ldm$dataset$umitab[,lung_ldm$dataset$cell_to_cluster=="44" & 
                                       sample_annots[lung_ldm$dataset$cell_to_sample,]$tissue=="Tumor"]

tfh_genes <- setdiff(tfh_genes,"CD4")
texh_genes <- setdiff(texh_genes,c("CD8A","CD8B"))
tfh_score <- colSums(clust_mat[tfh_genes,])/colSums(clust_mat)
texh_score <- colSums(clust_mat[texh_genes,])/colSums(clust_mat)
score_ratio <- log2((1e-3+tfh_score)/(1e-3+texh_score))
cell_counts <- table(score_ratio > cut, sample_annots[lung_ldm$dataset$cell_to_sample[colnames(clust_mat)],]$patient_ID)
rownames(cell_counts) <- c("CD8","CD4")
cell_freqs <- as.matrix(cell_counts)/matrix(total_cells[colnames(cell_counts)],nrow=nrow(cell_counts),ncol=ncol(cell_counts),byrow=T)

cell_freqs <- cbind(cell_freqs,c(0,0))
png("output/figures/figure_s4e.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
plot(0.01+cell_freqs[1,]*100,0.01+cell_freqs[2,]*100,xlab="CD8 Tactivated",ylab="CD4 Tactivated",log="xy",pch=16,xaxt="n",yaxt="n")
axis(side=1,at=c(0.01,0.1,1,10))
axis(side=2,at=c(0.01,0.1,1,10))
mtext("% of immune cells",cex=1.5)
dev.off()


cor.test(cell_freqs[1,],cell_freqs[2,],method="spearman")

}
