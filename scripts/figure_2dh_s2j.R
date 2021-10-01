library(scales)

figure_2dh_s2j <- function(){

output_dir <- "output/figures/"

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)

if(!exists("lung_ldm")){
  message("Loading Sinai scRNA data into R")
  load("data/lung_ldm.rd")
}

patient_tissue <- apply(sample_annots[lung_ldm$dataset$cell_to_sample,c("patient_ID","tissue")],1,paste,collapse="_")

sub_lineage <- annots_list[lung_ldm$dataset$cell_to_cluster,]$sub_lineage

tab <- table(sub_lineage,patient_tissue)
#tab <- tab[annots_list[match(rownames(tab),annots_list$sub_lineage),]$lineage=="MNP",]
tab <- tab[c("DC3","cDC2","mregDC","cDC1"),]
tab <- t(t(tab)/rowSums(t(tab)))
tab_normal <- tab[,grepl("Normal",colnames(tab))]
tab_tumor <- tab[,grep("Tumor",colnames(tab))]
colnames(tab_normal) <- substr(colnames(tab_normal),1,3)
colnames(tab_tumor) <- substr(colnames(tab_tumor),1,3)
tab_normal <- tab_normal[,intersect(colnames(tab_normal),colnames(tab_tumor))]
tab_tumor <- tab_tumor[,intersect(colnames(tab_normal),colnames(tab_tumor))]

delta <- log2((1e-2+tab_tumor)/(1e-2+tab_normal))
mat <- delta[c("DC3","cDC2","mregDC","cDC1"),]

############
tab <- table(sub_lineage,patient_tissue)
#tab <- tab[annots_list[match(rownames(tab),annots_list$sub_lineage),]$lineage=="MNP",]
tab <- tab[c("DC3","cDC2","mregDC","cDC1"),]
#tab <- t(t(tab)/rowSums(t(tab)))
tab_normal <- tab[,grepl("Normal",colnames(tab))]
tab_tumor <- tab[,grep("Tumor",colnames(tab))]
colnames(tab_normal) <- substr(colnames(tab_normal),1,3)
colnames(tab_tumor) <- substr(colnames(tab_tumor),1,3)
tab_normal <- tab_normal[,intersect(colnames(tab_normal),colnames(tab_tumor))]
tab_tumor <- tab_tumor[,intersect(colnames(tab_normal),colnames(tab_tumor))]


samp_mask <- colSums(tab_tumor)>30 & 
  colSums(tab_normal)>30

mat <- mat[,samp_mask]
##############

res <- apply(mat,1,wilcox.test)
sig_vec <- p.adjust(lapply(res,function(x){x$p.value}),method="bonferroni")

cols <- rgb(t(col2rgb(c("2","3","6","orange"))*0.7/255))



tab <- table(sub_lineage,patient_tissue)
tab <- tab[annots_list[match(rownames(tab),annots_list$sub_lineage),]$lineage=="MNP",]
#tab <- tab[c("DC3","cDC2","mregDC","cDC1"),]
tab <- t(t(tab)/rowSums(t(tab)))
tab_normal <- tab[,grepl("Normal",colnames(tab))]
tab_tumor <- tab[,grep("Tumor",colnames(tab))]
colnames(tab_normal) <- substr(colnames(tab_normal),1,3)
colnames(tab_tumor) <- substr(colnames(tab_tumor),1,3)
tab_normal <- tab_normal[,intersect(colnames(tab_normal),colnames(tab_tumor))]
tab_tumor <- tab_tumor[,intersect(colnames(tab_normal),colnames(tab_tumor))]

delta <- log2((1e-2+tab_tumor)/(1e-2+tab_normal))
mat <- delta[c("DC3","cDC2","mregDC","cDC1"),]

############
tab <- table(sub_lineage,patient_tissue)
tab <- tab[annots_list[match(rownames(tab),annots_list$sub_lineage),]$lineage=="MNP",]
#tab <- tab[c("DC3","cDC2","mregDC","cDC1"),]
#tab <- t(t(tab)/rowSums(t(tab)))
tab_normal <- tab[,grepl("Normal",colnames(tab))]
tab_tumor <- tab[,grep("Tumor",colnames(tab))]
colnames(tab_normal) <- substr(colnames(tab_normal),1,3)
colnames(tab_tumor) <- substr(colnames(tab_tumor),1,3)
tab_normal <- tab_normal[,intersect(colnames(tab_normal),colnames(tab_tumor))]
tab_tumor <- tab_tumor[,intersect(colnames(tab_normal),colnames(tab_tumor))]


samp_mask <- colSums(tab_tumor)>250 & 
  colSums(tab_normal)>250

mat <- mat[,samp_mask]
##############
res <- apply(mat,1,wilcox.test)
sig_vec <- p.adjust(lapply(res,function(x){x$p.value}),method="bonferroni")

cols <- rgb(t(col2rgb(c("2","3","6","orange"))*0.7/255))


mask <- sub_lineage%in%c("DC3","cDC1","cDC2","mregDC")

cell_split <- split(colnames(lung_ldm$dataset$umitab)[mask],f = list(patient_tissue[mask],sub_lineage[mask]))

avgs <- do.call(cbind,lapply(cell_split,function(x){
  if(length(x)>10){
    rs <- rowSums(lung_ldm$dataset$umitab[,x,drop=F])
    return(rs/sum(rs))
  }else{
    return(NA)
  }
}))

subtype <- factor(unlist(lapply(strsplit(colnames(avgs),"\\."),function(x){x[[2]]})),c("DC3","cDC2","mregDC","cDC1"))
genes <- c("LAMP3","CD274")


###################
# barplots for mregDC gene expression

cell_split <- split(colnames(lung_ldm$dataset$umitab)[mask],f = list(factor(sub_lineage[mask],c("DC3","cDC2","mregDC","cDC1"))))
total_avgs <- do.call(cbind,lapply(cell_split,function(x){
    rs <- rowSums(lung_ldm$dataset$umitab[,x,drop=F])
    return(rs/sum(rs))

}))

genes <- c("LAMP3","CD274")
png(file.path(output_dir,"figure_2d.png"),height=1.58,width=1.17,units="in",res=300,pointsize=5)
layout(matrix(1:2,nrow=2))
par(mar=c(2,3,1,1),oma=c(2,0,1,0),mgp=c(2,1,0))
for(gene in genes){
  b <- barplot((1e-6+total_avgs[gene,]),log="y",range=0,xaxt="n",yaxt="n",col=cols,ylim=c(1e-6,2.5*max(total_avgs[gene,])),
          border=NA)
  axis(side=2,at=c(1e-6,1e-5,1e-4,1e-3),labels = c(-6,-5,-4,-3),las=2)
  mtext(gene,line=0.25,font=3)
  mtext("Expression (Log10)",side=2,line=2)
  if(gene == tail(genes,1)){
    mtext(levels(subtype),las=2,at=b,line=0.25,side=1,col=cols)
  }
}
dev.off()
#############


tissue <- c("Normal","Tumor")[1+as.integer(grepl("Tumor",colnames(avgs)))]
genes <- c("CD1A","CD207","IL22RA2")
tissue_ord <- c("Normal","Tumor")
first_at <- c(1,4,7,10)

png(file.path(output_dir,"figure_2h.png"),height=1.53,width=5.41,units="in",res=300,pointsize=9)
layout(matrix(1:3,nrow=1))
par(mar=c(5,2,2,1),oma=c(0,2,0,0))
for(gene in genes){

for(iter in seq(length(tissue_ord))){
  x <- c(first_at+iter-1)[as.numeric(subtype[tissue==tissue_ord[iter]])]
  if(iter==1){add=F}else{add=T}
  if(iter==1){
    border=cols; col=NULL; lwd=2
  }else{
    border=NULL; col=cols; lwd=1
  }
  
  boxplot((1e-6+avgs[gene,tissue==tissue_ord[iter]])~subtype[tissue==tissue_ord[iter]],at=first_at+iter-1,xlim=c(0.5,11.5),xaxt="n",range=0,yaxt="n",log="y",add=add,
          ylim=range(1e-6+avgs[gene,],na.rm=T),
          border=border,col=col,lwd=lwd)
  points(jitter(x,0.5),(1e-6+avgs[gene,tissue==tissue_ord[iter]]),pch=16,col=alpha(1,0.5),cex=0.5)
  mtext(side=1,at=first_at+iter-1,tissue_ord[iter],las=2,line=0.25,col=cols)
}
axis(side=2,at=c(1e-6,1e-5,1e-4,1e-3),labels = c(-6,-5,-4,-3),las=2)
if(gene==genes[1]){
  mtext("Expression (Log10)",side=2,line=2)
}
mtext(gene,line=0.25,font=3)
}
dev.off()

rm(lung_ldm)

if(!exists("lambrechts_ldm")){
  message("Loading Lambrechts et al. data into R")
  load("data/lambrechts_ldm.rd")
}



patient_tissue <- apply(sample_annots[lambrechts_ldm$dataset$cell_to_sample,c("patient_ID","tissue")],1,paste,collapse="_")

sub_lineage <- annots_list[lambrechts_ldm$dataset$cell_to_cluster,]$sub_lineage


mask <- sub_lineage%in%c("DC3","cDC1","cDC2","mregDC")

cell_split <- split(colnames(lambrechts_ldm$dataset$umitab)[mask],f = list(patient_tissue[mask],sub_lineage[mask]))

avgs <- do.call(cbind,lapply(cell_split,function(x){
  if(length(x)>10){
    rs <- rowSums(lambrechts_ldm$dataset$umitab[,x,drop=F])
    return(rs/sum(rs))
  }else{
    return(NA)
  }
}))

subtype <- factor(unlist(lapply(strsplit(colnames(avgs),"\\."),function(x){x[[2]]})),c("DC3","cDC2","mregDC","cDC1"))


tissue <- c("Normal","Tumor")[1+as.integer(grepl("Tumor",colnames(avgs)))]
genes <- c("CD1A","CD207","IL22RA2")
tissue_ord <- c("Normal","Tumor")
first_at <- c(1,4,7,10)

png(file.path(output_dir,"figure_s2j.png"),height=1.53,width=5.41,units="in",res=300,pointsize=9)
layout(matrix(1:3,nrow=1))
par(mar=c(5,2,2,1),oma=c(0,2,0,0))
for(gene in genes){
  
  for(iter in seq(length(tissue_ord))){
    x <- c(first_at+iter-1)[as.numeric(subtype[tissue==tissue_ord[iter]])]
    if(iter==1){add=F}else{add=T}
    if(iter==1){
      border=cols; col=NULL; lwd=2
    }else{
      border=NULL; col=cols; lwd=1
    }
    
    boxplot((1e-6+avgs[gene,tissue==tissue_ord[iter]])~subtype[tissue==tissue_ord[iter]],at=first_at+iter-1,xlim=c(0.5,11.5),xaxt="n",range=0,yaxt="n",log="y",add=add,
            ylim=range(1e-6+avgs[gene,],na.rm=T),
            border=border,col=col,lwd=lwd)
    points(jitter(x,0.5),(1e-6+avgs[gene,tissue==tissue_ord[iter]]),pch=16,col=alpha(1,0.5),cex=0.5)
    mtext(side=1,at=first_at+iter-1,tissue_ord[iter],las=2,line=0.25,col=cols)
  }
  axis(side=2,at=c(1e-6,1e-5,1e-4,1e-3),labels = c(-6,-5,-4,-3),las=2)
  if(gene==genes[1]){
    mtext("Expression (Log10)",side=2,line=2)
  }
  mtext(gene,line=0.25,font=3)
}
dev.off()

}
