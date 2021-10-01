
library(scTools)
library(scDissector)
library(matrixStats)
library(scales)

figure_4de_s4ijk <- function(){

rm(list=ls())

sample_annots_path <- "input_tables/table_s1_sample_table.csv"
sample_annots <- read.csv(sample_annots_path,r=1,h=1,stringsAsFactors=F)

sample_sets <- list()
sample_sets$'695' <- c("115","116")
sample_sets$'522' <- c("343","344")
sample_sets$'706' <- c("480","481")

annot_list <- read.csv("input_tables/annots_list.csv",
                       r=1,h=1,stringsAsFactors=F)
clusters.exclude <- rownames(annot_list)[annot_list$lineage!="T"]

load(file="data/lung_ldm_fiveprime.rd")

# THESE LINES BUILD THE OBJECT tcr_sets USING THE FUNCTION tcr_output_to_analysis(), which takes in the raw TCR data output from Cellranger and assembles it into an R structure organized by patient. 
# THE OUTPUT HAS BEEN SAVED AND THIS SCRIPT LOADS THE OUTPUT.
# 
# tcr_sets <- list()
# 
# patient_iter <- names(sample_sets)[1]
# for(patient_iter in names(sample_sets)){
#   print(patient_iter)
#   tcr_sets[[patient_iter]] <- tcr_output_to_analysis(sample_IDs = sample_sets[[patient_iter]],
#                                                       sample_annots_path = sample_annots_path,demux = F,
#                                                       data_path = "/users/andrew leader/google drive/merad/scRNAseq_analysis/sc_data_main/lung_data/")
#   names(tcr_sets[[patient_iter]]$TCRid) <- substr(names(tcr_sets[[patient_iter]]$TCRid),1,nchar(names(tcr_sets[[patient_iter]]$TCRid))-2)
# }

# save(tcr_sets,file="data/tcr_sets_190906.rd")

# For simplicity, instead of running the above script, we will download tcr_sets
if(!file.exists("data/tcr_sets_190906.rd")){
  message("Downloading tcr data")
  data_url <- "https://www.dropbox.com/s/a43ksvxjru3lv6t/tcr_sets_190906.rd?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = file.path("data/tcr_sets_190906.rd"),mode="wb")
  }else{
    download.file(url = data_url,destfile = file.path(data_dir,"data/tcr_sets_190906.rd"))
  }
}
load("data/tcr_sets_190906.rd")

cellnames <- names(lung_ldm_fiveprime$dataset$cell_to_cluster)
cellnames <- substr(cellnames,6,nchar(cellnames))
names(lung_ldm_fiveprime$dataset$cell_to_cluster)=names(lung_ldm_fiveprime$dataset$cell_to_sample)=colnames(lung_ldm_fiveprime$dataset$umitab) <- cellnames

lung_ldm_fiveprime$dataset$samples <- substr(lung_ldm_fiveprime$dataset$samples ,6,nchar(lung_ldm_fiveprime$dataset$samples ))
lung_ldm_fiveprime$dataset$cell_to_sample <- substr(lung_ldm_fiveprime$dataset$cell_to_sample,
                                                    6,nchar(lung_ldm_fiveprime$dataset$cell_to_sample))



cell_to_annot <- annot_list[lung_ldm_fiveprime$dataset$cell_to_cluster,]$sub_lineage
names(cell_to_annot) <- names(lung_ldm_fiveprime$dataset$cell_to_sample)

s <- split(names(lung_ldm_fiveprime$dataset$cell_to_cluster),lung_ldm_fiveprime$dataset$cell_to_sample)

names(tcr_sets$`706`$TCRid) <- paste(names(tcr_sets$`706`$TCRid),"-1",sep="")
names(tcr_sets$`695`$TCRid) <- paste(names(tcr_sets$`695`$TCRid),"-1",sep="")
names(tcr_sets$`522`$TCRid) <- paste(names(tcr_sets$`522`$TCRid),"-1",sep="")


clonality_mat <- matrix(NA,nrow=length(unique(cell_to_annot)),ncol=length(s),
                                       dimnames=list(unique(cell_to_annot),names(s)))

set.seed(910430)
for(samp_iter in names(s)){
  patient <- sample_annots[samp_iter,]$patient_ID
  cells <- intersect(s[[samp_iter]],names(tcr_sets[[patient]]$TCRid))
  
  cells <- cells[!tcr_sets[[patient]]$TCRid[cells]%in%c("multiplet","single_chain")]
  #print(paste(samp_iter,length(cells)))
  #cells <- sample(cells,800)
  
  ss <- split(cells,cell_to_annot[cells])
  ss <- lapply(ss,function(y){sample(y,pmin(30,length(y)))})
  cells <- unlist(ss)
  # 
  tab <- table(cell_to_annot[cells],tcr_sets[[patient]]$TCRid[cells])
  #print(array(rowSums(tab)))
  #tab[rowSums(tab) < 2,] <- NA
  tab[rowSums(tab)<30,] <- NA
  
  #print(summary(rowSums(tab)))
  
  entropy <- apply(tab,1,function(x){
    x <- x[x>0]
    p <- x/sum(x)
    plogp <- p*log2(p)
    return(-sum(plogp)/log2(length(x)))
  })
  
  clonality <- 1-entropy
  
  clonality_mat[,samp_iter] <- clonality[rownames(clonality_mat)]
  
  
}

annot_ord <- c("T_Nklike","T_GZMK","CD8 Trm","T_activated","CD4 Trm","Tcm/naive_II","Tcm/naive_I","Treg","Tcycle")
clonality_mat <- clonality_mat[annot_ord,]
tissue <- sample_annots[colnames(clonality_mat),]$tissue

x_normal <- c(1+3*c(0:(length(annot_ord)-1)))
x_tumor <- c(2+3*c(0:(length(annot_ord)-1)))

png("output/figures/figure_s4i.png",height=4,width=6,units="in",res=300,bg="transparent")

plot(jitter(matrix(x_normal,nrow=nrow(clonality_mat),ncol=sum(tissue=="Normal")),.4),clonality_mat[,tissue=="Normal"],
     ylim=range(clonality_mat[!is.nan(clonality_mat)],na.rm=T),xlim=range(c(x_normal,x_tumor)),pch=16,col="blue",bty="L",
     xaxt="n",ylab="clonality",xlab="")
points(jitter(matrix(x_tumor,nrow=nrow(clonality_mat),ncol=sum(tissue=="Tumor")),.4),clonality_mat[,tissue=="Tumor"],pch=16,col="brown")

segments(x0 = x_normal,x1=x_normal,y0=rowMins(clonality_mat[,tissue=="Normal"],na.rm=T),
         y1=rowMaxs(clonality_mat[,tissue=="Normal"],na.rm=T),col=alpha("blue",.4))
segments(x0 = x_tumor,x1=x_tumor,y0=rowMins(clonality_mat[,tissue=="Tumor"],na.rm=T),
         y1=rowMaxs(clonality_mat[,tissue=="Tumor"],na.rm=T),col=alpha("brown",.4))
segments(x0=x_normal-.5,x1=x_normal+.5,
         y0=rowMedians(clonality_mat[,tissue=="Normal"],na.rm=T),
         y1=rowMedians(clonality_mat[,tissue=="Normal"],na.rm=T),col="blue",lwd=3)
segments(x0=x_tumor-.5,x1=x_tumor+.5,
         y0=rowMedians(clonality_mat[,tissue=="Tumor"],na.rm=T),
         y1=rowMedians(clonality_mat[,tissue=="Tumor"],na.rm=T),col="brown",lwd=3)
abline(v=c(3,6,9,12,15,18,21,24),col="grey",lty=2)

segments(x0=matrix(x_normal,nrow=nrow(clonality_mat),ncol=sum(tissue=="Normal")),
         x1=matrix(x_tumor,nrow=nrow(clonality_mat),ncol=sum(tissue=="Tumor")),
         y0=clonality_mat[,tissue=="Normal"],
         y1=clonality_mat[,tissue=="Tumor"],col="grey")
mtext("Subset clonality",cex=1.5)
dev.off()
#for(samp_iter in names(s)){

#annot_ord <- c("CTL","CD8 Trm","T-Gzmk","Texh","Tcycle","CD4 Trm","Tnaive/Tcm","Treg")
annot_ord <- c("Tcm/naive_II","CD4 Trm","T_Nklike","T_GZMK","CD8 Trm","Tcm/naive_I","Tcycle","Treg","T_activated")


# continue from here
########################################
freqs_list <- list()  

pat_iter <- names(sample_sets)[1]
for(pat_iter in names(sample_sets)){
  
  cells_n <- intersect(s[[sample_sets[[pat_iter]][1]]],names(tcr_sets[[pat_iter]]$TCRid))
  cells_n <- cells_n[!tcr_sets[[pat_iter]]$TCRid[cells_n]%in%c("multiplet","single_chain")]
  cells_t <- intersect(s[[sample_sets[[pat_iter]][2]]],names(tcr_sets[[pat_iter]]$TCRid))
  cells_t <- cells_t[!tcr_sets[[pat_iter]]$TCRid[cells_t]%in%c("multiplet","single_chain")]
  print(length(cells_n))
  print(head(cells_n))
  tab_n <- table(tcr_sets[[pat_iter]]$TCRid[cells_n],cell_to_annot[cells_n])
  tab_t <- table(tcr_sets[[pat_iter]]$TCRid[cells_t],cell_to_annot[cells_t])
  print(sum(tab_t))
  fullmat_n <- matrix(0,nrow=length(unique(c(rownames(tab_n),rownames(tab_t)))),ncol=length(annot_ord),
                      dimnames=list(unique(c(rownames(tab_n),rownames(tab_t))),annot_ord))
  fullmat_t <- fullmat_n
  fullmat_n[rownames(tab_n),] <- tab_n[,colnames(fullmat_n)]
  fullmat_t[rownames(tab_t),] <- tab_t[,colnames(fullmat_t)]
  freqs_list[[pat_iter]] <- list()
  freqs_list[[pat_iter]]$fullmat_n <- fullmat_n
  freqs_list[[pat_iter]]$fullmat_t <- fullmat_t
 
}

#shared_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(x[[1]])>0 & rowSums(x[[2]])>0]})
shared_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(x[[1]])>0 & rowSums(x[[2]])>0 & rowSums(cbind(x[[1]],x[[2]]))>3]})
n_specific_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(x[[1]])>3 & rowSums(x[[2]])==0]})
t_specific_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(x[[1]])==0 & rowSums(x[[2]])>3]})
unique_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(cbind(x[[1]],x[[2]]))==1]})
not_gated_n_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(x[[2]]==0 & rowSums(x[[1]])%in%c(2:3))]})
not_gated_t_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(x[[1]]==0 & rowSums(x[[2]])%in%c(2:3))]})
not_gated_shared_list <- lapply(freqs_list,function(x){rownames(x[[1]])[rowSums(x[[1]])>0 & rowSums(x[[2]])>0 & rowSums(cbind(x[[1]],x[[2]]))<=3]})




######

boxcols <- c("purple","black",rgb(0,.7,0))

x <- jitter(rowSums(freqs_list$`695`$fullmat_n)+1,.5)
y <- jitter(rowSums(freqs_list$`695`$fullmat_t)+1,.5)

col <- array(1,nrow(freqs_list$`695`$fullmat_n));names(col) <- rownames(freqs_list$`695`$fullmat_n)
col[shared_list$`695`] <- rgb(0,.7,0)
col[n_specific_list$`695`] <- "blue"
col[t_specific_list$`695`] <- "brown"
col[unique_list$`695`] <- "grey"

png("output/figures/figure_4d.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
plot(x,y,log="xy",pch=16,bty="L",xlab="1 + # cells in nLung",ylab="1 + # cells in tumor",col=alpha(col,.5))
mtext("Cells per TCR",cex=1.5)
dev.off()
# shared <- rowSums(freqs_list$`695`$fullmat_n) > 0 & rowSums(freqs_list$`695`$fullmat_t) > 0
# n_specific <- rowSums(freqs_list$`695`$fullmat_n) > 1 & rowSums(freqs_list$`695`$fullmat_t) == 0
# t_specific <- rowSums(freqs_list$`695`$fullmat_n) == 0 & rowSums(freqs_list$`695`$fullmat_t) > 1
# unique_in_both <- which(as.logical(1-shared-n_specific-t_specific))
# points(x[shared],y[shared],col=alpha(boxcols[1],.5),pch=1,cex=.5)
# points(x[t_specific],y[t_specific],col=alpha(boxcols[2],.5),pch=1,cex=.5)
# points(x[n_specific],y[n_specific],col=alpha(boxcols[2],.5),pch=1,cex=.5)
# points(x[unique_in_both],y[unique_in_both],col=alpha(boxcols[3],.5),pch=1,cex=.5)




shared_freqs_t <- list()
shared_freqs_n <- list()
n_specific_n <- list()
n_specific_t <- list()
t_specific_t <- list()
t_specific_n <- list()
unique_both_n <- list()
unique_both_t <- list()
not_gated_n <- list()
not_gated_t <- list()


for(pat_iter in names(sample_sets)){
  shared_freqs_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[shared_list[[pat_iter]],])
  shared_freqs_t[[pat_iter]] <- shared_freqs_t[[pat_iter]]/sum(shared_freqs_t[[pat_iter]])
  shared_freqs_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[shared_list[[pat_iter]],])
  shared_freqs_n[[pat_iter]] <- shared_freqs_n[[pat_iter]]/sum(shared_freqs_n[[pat_iter]])
  
  
  n_specific_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[n_specific_list[[pat_iter]],])
  n_specific_t[[pat_iter]] <- n_specific_t[[pat_iter]]/sum(n_specific_t[[pat_iter]])
  n_specific_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[n_specific_list[[pat_iter]],])
  n_specific_n[[pat_iter]] <- n_specific_n[[pat_iter]]/sum(n_specific_n[[pat_iter]])
 

  t_specific_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[t_specific_list[[pat_iter]],])
  t_specific_t[[pat_iter]] <- t_specific_t[[pat_iter]]/sum(t_specific_t[[pat_iter]])
  t_specific_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[t_specific_list[[pat_iter]],])
  t_specific_n[[pat_iter]] <- t_specific_n[[pat_iter]]/sum(t_specific_n[[pat_iter]])

  unique_both_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[unique_list[[pat_iter]],])
  unique_both_t[[pat_iter]] <- unique_both_t[[pat_iter]]/sum(unique_both_t[[pat_iter]])
  unique_both_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[unique_list[[pat_iter]],])
  unique_both_n[[pat_iter]] <- unique_both_n[[pat_iter]]/sum(unique_both_n[[pat_iter]])

}




# layout(matrix(1:4,nrow=4))
# x <- matrix(1:16,ncol=16,nrow=3,byrow=T)
# pch <- matrix(c(1,2,4),ncol=16,nrow=3)
# 
# y <- log10(cbind(do.call(rbind,shared_freqs_n),do.call(rbind,shared_freqs_t))*100+1)
# plot(x,y,col=rep(c("blue","brown"),each=24),xaxt="n",ylim=c(0,2),pch=pch,ylab="",cex=2); box(col="red",lwd=3); abline(v=8.5,lty=2)
# segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=rep(c("blue","brown"),each=8))
# segments(x0=x[1,]-.25,x1=x[1,]+.25,y0=colMeans(y),y1=colMeans(y),col=rep(c("blue","brown"),each=8),lwd=2)
# grid(nx=0,ny=NULL)
# 
# y <- log10(1+cbind(do.call(rbind,n_specific_n),do.call(rbind,n_specific_t))*100)
# plot(x,y,col=rep(c("blue","brown"),each=24),xaxt="n",ylim=c(0,2),pch=pch,ylab="",cex=2); box(col="purple",lwd=3); abline(v=8.5,lty=2)
# segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=rep(c("blue","brown"),each=8))
# segments(x0=x[1,]-.25,x1=x[1,]+.25,y0=colMeans(y),y1=colMeans(y),col=rep(c("blue","brown"),each=8),lwd=2)
# grid(nx=0,ny=NULL)
# 
# y <- log10(1+cbind(do.call(rbind,t_specific_n),do.call(rbind,t_specific_t))*100)
# plot(x,y,col=rep(c("blue","brown"),each=24),xaxt="n",ylim=c(0,2),pch=pch,ylab="",cex=2); box(col="blue",lwd=3); abline(v=8.5,lty=2)
# segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=rep(c("blue","brown"),each=8))
# segments(x0=x[1,]-.25,x1=x[1,]+.25,y0=colMeans(y),y1=colMeans(y),col=rep(c("blue","brown"),each=8),lwd=2)
# grid(nx=0,ny=NULL)
# 
# y <- log10(1+cbind(do.call(rbind,unique_both_n),do.call(rbind,unique_both_t))*100)
# plot(x,y,col=rep(c("blue","brown"),each=24),xaxt="n",ylim=c(0,2),pch=pch,ylab="",cex=2); box(col="orange",lwd=3); abline(v=8.5,lty=2)
# segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=rep(c("blue","brown"),each=8))
# segments(x0=x[1,]-.25,x1=x[1,]+.25,y0=colMeans(y),y1=colMeans(y),col=rep(c("blue","brown"),each=8),lwd=2)
# grid(nx=0,ny=NULL)

#####################
# png("/users/andrew leader/google drive/merad/Leader_et_al_figure_scripts/Fig4/TCR_pop_freqs.png",height=3,width=1.65,units="in",res=300,
#     pointsize=6)
# layout(matrix(1:3,nrow=3))
# par(mar=c(2,5,2,2))
# #x <- matrix(1:16,ncol=16,nrow=3,byrow=T)
# x <- matrix(c(seq(1,15,2),seq(2,16,2)),ncol=16,nrow=3,byrow=T)
# pch <- matrix(c(15:17),ncol=16,nrow=3)
# v <- seq(2.5,14.5,2)
# 
# y <- log10(cbind(do.call(rbind,shared_freqs_n),do.call(rbind,shared_freqs_t))*100+1)
# plot(x,y,col=rep(c("blue","brown"),each=24),xaxt="n",ylim=c(0,2),pch=pch,ylab="",cex=1); box(col=boxcols[1],lwd=1)
# segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=rep(c("blue","brown"),each=8))
# segments(x0=x[1,]-.25,x1=x[1,]+.25,y0=colMeans(y),y1=colMeans(y),col=rep(c("blue","brown"),each=8),lwd=2)
# #grid(nx=0,ny=NULL)
# abline(v=v,lty=2,col="grey")
# 
# 
# y <- log10(1+cbind(do.call(rbind,n_specific_n),do.call(rbind,n_specific_t))*100)
# y1 <- log10(1+cbind(do.call(rbind,t_specific_n),do.call(rbind,t_specific_t))*100)
# y[,9:16] <- y1[,9:16]
# plot(x,y,col=rep(c("blue","brown"),each=24),xaxt="n",ylim=c(0,2),pch=pch,ylab="",cex=1); box(lwd=1,col=boxcols[2])
# segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=rep(c("blue","brown"),each=8))
# segments(x0=x[1,]-.25,x1=x[1,]+.25,y0=colMeans(y),y1=colMeans(y),col=rep(c("blue","brown"),each=8),lwd=2)
# #grid(nx=0,ny=NULL)
# abline(v=v,lty=2,col="grey")
# 
# 
# 
# y <- log10(1+cbind(do.call(rbind,unique_both_n),do.call(rbind,unique_both_t))*100)
# plot(x,y,col=rep(c("blue","brown"),each=24),xaxt="n",ylim=c(0,2),pch=pch,ylab="",cex=1); box(col=boxcols[3],lwd=1)
# segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=rep(c("blue","brown"),each=8))
# segments(x0=x[1,]-.25,x1=x[1,]+.25,y0=colMeans(y),y1=colMeans(y),col=rep(c("blue","brown"),each=8),lwd=2)
# #grid(nx=0,ny=NULL)
# abline(v=v,lty=2,col="grey")
# dev.off()

###############

unique_both_n <- do.call(rbind,unique_both_n)
unique_both_t <- do.call(rbind,unique_both_t)
shared_freqs_n <- do.call(rbind,shared_freqs_n)
shared_freqs_t <- do.call(rbind,shared_freqs_t)
n_specific_n <- do.call(rbind,n_specific_n)
t_specific_t <- do.call(rbind,t_specific_t)
pch <- matrix(c(1,2,4),nrow=3,ncol=24)
x <- matrix(1:27,nrow=3,ncol=27,byrow=T)
v <- seq(3.5,24.5,3)


png("output/figures/figure_4e.png",
    height=1.9,width=3.37,units="in",res=500,pointsize=6)
layout(matrix(1:2,nrow=2),heights=c(10,10))
par(mar=c(2,5,.3,2),oma=c(0,5,2,0))

col <- matrix(c("grey","blue",rgb(0,.6,0)),nrow=3,ncol=24,byrow=T)
y <- (cbind(unique_both_n,n_specific_n,shared_freqs_n)*100+1)
y <- y[,rep(seq(0,18,9),9)+rep(1:9,each=3)]
plot(x,y,pch=pch,col=col,xaxt="n",xlab="",ylab="",ylim=c(1,100),log="y",yaxt="n",yaxs="i")
axis(side=2,at=c(1,5,20,100),las=2)
#box(col="blue",lwd=2)
segments(x0=x[1,],x1=x[1,],y0=colMaxs(y),y1=colMins(y),col=col[1,])
abline(v=v,col="grey")
abline(h=c(5,20),col="grey",lty=2)

col <- matrix(c("grey","brown",rgb(0,.6,0)),nrow=3,ncol=24,byrow=T)
y <- (cbind(unique_both_t,t_specific_t,shared_freqs_t)*100+1)
y <- y[,rep(seq(0,18,9),9)+rep(1:9,each=3)]
plot(x,y,pch=pch,col=col,xaxt="n",xlab="",ylab="",ylim=c(1,100),log="y",yaxt="n",yaxs="i")
#box(col="brown",lwd=2)
segments(x0=x[1,],x1=x[1,],y0=colMaxs(y),y1=colMins(y),col=col[1,])
abline(v=v,col="grey")
abline(h=c(5,20),col="grey",lty=2)
axis(side=2,at=c(1,5,20,100),las=2)


# par(mar=c(0,5,0.3,2))
# plot(1:27,seq(0.5,3.5,3/(27-1)),xaxt="n",yaxt="n",ylab="",col="white")
# text(x=seq(1,25,3),y=rep(3,9),"+",cex=2)
# text(x=seq(2,26,3),y=rep(2,9),"+",cex=2)
# text(x=seq(3,27,3),y=rep(1,9),"+",cex=2)
# abline(v=v,col="grey")
# mtext(side=2,at=3:1,c("Non-clonal","Tissue specific","Shared"),col=c(boxcols[3:1]),las=2,line=.25)
# mtext(side=1,at=seq(2,23,3),las=2,line=.25,annot_ord)
dev.off()

#####################
shared_freqs_t <- list()
shared_freqs_n <- list()
n_specific_n <- list()
n_specific_t <- list()
t_specific_t <- list()
t_specific_n <- list()
unique_both_n <- list()
unique_both_t <- list()
not_gated_n <- list()
not_gated_t <- list()

for(pat_iter in names(sample_sets)){
  shared_freqs_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[shared_list[[pat_iter]],])
  shared_freqs_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[shared_list[[pat_iter]],])
  
  
  n_specific_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[n_specific_list[[pat_iter]],])
  n_specific_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[n_specific_list[[pat_iter]],])
  
  
  t_specific_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[t_specific_list[[pat_iter]],])
  t_specific_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[t_specific_list[[pat_iter]],])
  
  unique_both_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[unique_list[[pat_iter]],])
  unique_both_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[unique_list[[pat_iter]],])
  
  not_gated_t[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_t[c(not_gated_t_list[[pat_iter]],not_gated_shared_list[[pat_iter]]),])
  not_gated_n[[pat_iter]] <- colSums(freqs_list[[pat_iter]]$fullmat_n[c(not_gated_n_list[[pat_iter]],not_gated_shared_list[[pat_iter]]),])
  
}

unique_both_n <- do.call(rbind,unique_both_n)
unique_both_t <- do.call(rbind,unique_both_t)
shared_freqs_n <- do.call(rbind,shared_freqs_n)
shared_freqs_t <- do.call(rbind,shared_freqs_t)
n_specific_n <- do.call(rbind,n_specific_n)
t_specific_t <- do.call(rbind,t_specific_t)
not_gated_n <- do.call(rbind,not_gated_n)
not_gated_t <- do.call(rbind,not_gated_t)

col <- matrix(c(rgb(0,.6,0),"black","purple"),nrow=3,ncol=27,byrow=T)
pch <- matrix(c(1,2,4),nrow=3,ncol=27)
x <- matrix(1:27,nrow=3,ncol=27,byrow=T)
v <- seq(3.5,24.5,3)
# 
# png("output/figures/fig4_t&b/clonal_gating_quant_absolute.png",
#     height=2.34,width=3.37,units="in",res=500,pointsize=6)
# layout(matrix(1:3,nrow=3),heights=c(10,10,5))
# par(mar=c(0.3,5,.3,2),oma=c(10,5,2,0))
# y <- cbind(unique_both_n,n_specific_n,shared_freqs_n)+1
# y <- y[,rep(seq(0,18,9),9)+rep(1:9,each=3)]
# plot(x,y,pch=pch,col=col,xaxt="n",xlab="",ylab="",ylim=c(1,900),log="y",xlim=c(0.5,27.5),xaxs="i")
# box(col="blue",lwd=2)
# segments(x0=x[1,],x1=x[1,],y0=colMaxs(y),y1=colMins(y),col=col[1,])
# abline(v=v,col="grey")
# grid(nx=0,ny=NULL)
# 
# par(mar=c(0.3,5,0.3,2))
# y <- cbind(unique_both_t,t_specific_t,shared_freqs_t)+1
# y <- y[,rep(seq(0,18,9),9)+rep(1:9,each=3)]
# plot(x,y,pch=pch,col=col,xaxt="n",xlab="",ylab="",ylim=c(1,900),log="y",xlim=c(0.5,27.5),xaxs="i")
# box(col="brown",lwd=2)
# segments(x0=x[1,],x1=x[1,],y0=colMaxs(y),y1=colMins(y),col=col[1,])
# abline(v=v,col="grey")
# grid(nx=0,ny=NULL)
# 
# par(mar=c(0,5,0.3,2))
# plot(1:27,seq(0.5,3.5,3/(27-1)),xaxt="n",yaxt="n",ylab="",col="white",xlim=c(0.5,27.5),xaxs="i")
# text(x=seq(1,25,3),y=rep(3,9),"+",cex=2)
# text(x=seq(2,26,3),y=rep(2,9),"+",cex=2)
# text(x=seq(3,27,3),y=rep(1,9),"+",cex=2)
# abline(v=v,col="grey")
# mtext(side=2,at=3:1,c("Non-clonal","Tissue specific","Shared"),col=c(boxcols[3:1]),las=2,line=.25)
# mtext(side=1,at=seq(2,26,3),las=2,line=.25,annot_ord)
# dev.off()

y <- cbind(rowSums(unique_both_n),rowSums(shared_freqs_n),rowSums(n_specific_n),rowSums(not_gated_n),
           rowSums(unique_both_t),rowSums(shared_freqs_t),rowSums(t_specific_t),rowSums(not_gated_t))
x <- matrix(1:8,nrow=3,ncol=8,byrow=T)
pch <- matrix(c(1,2,4),nrow=3,ncol=8)
col <- matrix(c("grey",rgb(0,.7,0),"blue",1,"grey",rgb(0,.7,0),"brown",1),nrow=3,ncol=8,byrow=T)

#par(mar=c(0.3,5,0.3,2),oma=c(0,5,2,0))
png("output/figures/figure_s4j.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
plot(x,y,pch=pch,col=col,log="y",ylim=c(5,5000),xaxt="n",xlab="",yaxt="n",ylab="Cells")
axis(side=2,at=c(5,50,500,5000))
abline(col="grey",v=4.5)
abline(col="grey",h=c(5,50,500,5000),lty=2)
segments(x0=x[1,],x1=x[1,],y0=colMaxs(y),y1=colMins(y),col=col[1,])
dev.off()


## number of unique TCRs
x <- unlist(lapply(unique_both_n,length))

shared_freqs_t <- list()
shared_freqs_n <- list()
n_specific_n <- list()
n_specific_t <- list()
t_specific_t <- list()
t_specific_n <- list()
unique_both_n <- list()
unique_both_t <- list()
not_gated_n <- list()
not_gated_t <- list()
not_gated_shared <- list()

for(pat_iter in names(sample_sets)){
  shared_freqs_t[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_t[shared_list[[pat_iter]],])>0)
  shared_freqs_n[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_n[shared_list[[pat_iter]],])>0)
  
  
  n_specific_t[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_t[n_specific_list[[pat_iter]],])>0)
  n_specific_n[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_n[n_specific_list[[pat_iter]],])>0)
  
  
  t_specific_t[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_t[t_specific_list[[pat_iter]],])>0)
  t_specific_n[[pat_iter]] <- sum(colSums(freqs_list[[pat_iter]]$fullmat_n[t_specific_list[[pat_iter]],])>0)
  
  unique_both_t[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_t[unique_list[[pat_iter]],])>0)
  unique_both_n[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_n[unique_list[[pat_iter]],])>0)
  
  not_gated_t[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_t[not_gated_t_list[[pat_iter]],])>0)
  not_gated_n[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_n[not_gated_n_list[[pat_iter]],])>0)
  not_gated_shared[[pat_iter]] <- sum(rowSums(freqs_list[[pat_iter]]$fullmat_n[not_gated_shared_list[[pat_iter]],])>0)
  
}

y <- cbind(unlist(shared_freqs_t),unlist(n_specific_n),
           unlist(t_specific_t),unlist(unique_both_n),
           unlist(unique_both_t),unlist(not_gated_n),unlist(not_gated_t),unlist(not_gated_shared))
x <- matrix(1:8,nrow=3,ncol=8,byrow=T)
pch=matrix(c(1,2,4),nrow=3,ncol=7)

col=matrix(c(rgb(0,.6,0),"blue","brown","grey","grey","black","black","black"),nrow=3,ncol=8,byrow=T)
png("output/figures/figure_s4k.png",height=4,width=6,units="in",res=300,bg="transparent")
plot(x,y,log="y",xaxt="n",xlab="",ylab="",pch=pch,col=col,yaxt="n",ylim=c(1,1500))
segments(x0=x[1,],x1=x[1,],y0=colMins(y),y1=colMaxs(y),col=col[1,])
axis(side=2,c(1,2,5,10,20,50,100,200,500,1000),las=2)
abline(h=c(1,2,5,10,20,50,100,200,500,1000),col="grey",lty=2)
mtext("# Unique TCRs per group",cex=1.5)
dev.off()
####################
}
