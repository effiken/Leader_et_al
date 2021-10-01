library(scDissector)
setwd("~/GoogleDrive/work/shared/analysis/lung/")
scDissector_datadir="~/GoogleDrive/work/shared/clustering_data/clustering_data_lung4/"
sample_tab=read.csv(paste(scDissector_datadir,"samples.csv",sep="/"),stringsAsFactors = F)

#load("lung_ldm.rd")

annots=read.csv("~/GoogleDrive/work/shared/clustering_data/clustering_data_lung4/metadata/sample_annots.csv",stringsAsFactors = F)
rownames(annots)=annots$sample_ID
samples=lung_ldm$dataset$samples
boxplot(split(lung_ldm$dataset$alpha_noise,annots[match(samples,annots$sample_ID),]$prep),ylab="alpha (Noise) ")

sort(lung_ldm$dataset$alpha_noise[intersect(names(lung_ldm$dataset$alpha_noise),rownames(annots)[annots$prime=="3"&annots$library_chemistry=="V2"])])



if (1==0){

m1=(do.call(cbind,maskl[[1]]))

mv2=m1[,annots[colnames(m1),]$library_chemistry=="V2"&annots[colnames(m1),]$prime=="3"]
m5=m1[,annots[colnames(m1),]$prime=="5"]
msort=m1[,annots[colnames(m1),]$prep=="sort"]
mbeads=m1[,annots[colnames(m1),]$prep=="beads"]
mdead=m1[,annots[colnames(m1),]$prep=="digest/dead cell"]
tail(rownames(lung_ldm$model$models)[order(rowMeans(mdead==2))],100)
g1="IGHG3";g2="IGHG1"
tab_beads_dead=function(g1,g2){
  l=list()
  dsg=lung_ldm$dataset$ds[[1]][c(g1,g2),annots[lung_ldm$dataset$cell_to_sample[colnames(lung_ldm$dataset$ds[[1]])],]$prep=="beads"&annots[lung_ldm$dataset$cell_to_sample[colnames(lung_ldm$dataset$ds[[1]])],]$library_chemistry=="V2"]
  tab=table(floor(log2(1+dsg[1,])),floor(log2(1+dsg[2,])))
  l$beads=round(tab/sum(tab),digits=4)
  dsg=lung_ldm$dataset$ds[[1]][c(g1,g2),annots[lung_ldm$dataset$cell_to_sample[colnames(lung_ldm$dataset$ds[[1]])],]$prep=="digest/dead cell"&annots[lung_ldm$dataset$cell_to_sample[colnames(lung_ldm$dataset$ds[[1]])],]$library_chemistry=="V2"]
  tab=table(floor(log2(1+dsg[1,])),floor(log2(1+dsg[2,])))
  l$dead_digest=round(tab/sum(tab),digits=4)
  print(noise_model_beads[c(g1,g2)])
  print(noise_model_dead[c(g1,g2)])
  l
}

noise_model_beads=rowMeans(lung_ldm$dataset$noise_models[,rownames(annots)[annots$sample_ID%in%samples&annots$prep=="beads"&annots$library_chemistry=="V2"]])
noise_model_dead=rowMeans(lung_ldm$dataset$noise_models[,rownames(annots)[annots$sample_ID%in%samples&annots$prep=="digest/dead cell"&annots$library_chemistry=="V2"]])

}

get_obs_exp_stats=function(ldm,samples,genes){
  obs_l=list()
  exp_l=list()
  obs_numi=list()
  for (samp in samples){
    print(samp)
  
    clust_tab=table(ldm$dataset$cell_to_cluster[ldm$dataset$cell_to_sample==samp])
    clusts=names(clust_tab)[clust_tab>10]
    obs_numi[[samp]]=ldm$dataset$counts[samp,genes,clusts]
    numis_per_clust=colSums(ldm$dataset$counts[samp,genes,clusts])
    obs_l[[samp]]=t(t(ldm$dataset$counts[samp,genes,clusts])/numis_per_clust)
    exp_l[[samp]]=ldm$model$models[genes,clusts]*(1-ldm$dataset$alpha_noise[samp])+ldm$dataset$noise_models[genes,samp]*ldm$dataset$alpha_noise[samp]
  
  }
  return(list(obs=obs_l,exp=exp_l,obs_numi=obs_numi))
}

generate_projection_qc_plots=function(ldm,samples,obs_exp,genes,save_plots=T,path="output/figures/qc_figures/"){

  clustAnnots=paste(rep(names(ldm$cluster_sets),sapply(ldm$cluster_sets,function(z){length(unlist(z))})),unlist(sapply(ldm$cluster_sets,function(x){rep(names(x),sapply(x,length))})))
  names(clustAnnots)=unlist(ldm$cluster_sets)
  maskl=list()
  for (samp in samples){
    print(samp)
    clustv=colnames(obs_exp[["obs"]][[samp]])
    maskl[[samp]]=Matrix(0,length(genes),length(clustv),dimnames = list(genes,clustv))
    for (clust in clustv){
      
      rang=c(5e-6,5e-2)
      
      if (is.null(maskl[[clust]])){
        maskl[[clust]]=list()
      }
      
      obs=obs_exp[["obs"]][[samp]][genes,clust]
      obs_numis=obs_exp[["obs_numi"]][[samp]][genes,clust]
      
      reg=ldm$model$params$reg
      
      if (save_plots){
        dir.create(paste(path,"/",samp,sep=""),showWarnings = F)
        png(paste(path,"/",samp,"/noise_QC_samp_",samp,"_clust_",clust,".png",sep=""),1000,400)
        layout(matrix(1:3,1,3))
        colv=rep(1,length(obs))
        mask_bigger=obs_numis-4*sqrt(sum(obs_numis)*ldm$model$models[genes,clust])>sum(obs_numis)*ldm$model$models[genes,clust] & log10((reg+obs)/(reg+ldm$model$models[genes,clust]))>.5
        colv[mask_bigger]=2
        mask_smaller=obs_numis+4*sqrt(sum(obs_numis)*ldm$model$models[genes,clust])<sum(obs_numis)*ldm$model$models[genes,clust] & log10((reg+obs)/(reg+ldm$model$models[genes,clust]))<(-.5)
        colv[mask_smaller]=4
        cormethod="pearson"
        plot(reg+ldm$dataset$noise_models[genes,samp],reg+ldm$model$models[genes,clust],log="xy",main=paste("alpha =",round(ldm$dataset$alpha_noise[samp],digits=3),sep=""),col=colv,xlab=paste("Noise model - sample",samp),ylab=paste("model cluster",clust),pch=20,cex=ifelse(colv!=1,1,.1),panel.first = grid(lty=1),xlim=rang,ylim=rang)
        abline(0,1)
        mtext(side=3,clustAnnots[clust])
        plot(reg+obs,reg+ldm$model$models[genes,clust],log="xy",col=colv,ylab=paste("model cluster",clust),xlab=paste("cluster",clust,"sample",samp),pch=20,cex=ifelse(colv!=1,1,.1),main=round(cor(log10(reg+obs),log10(reg+ldm$model$models[genes,clust]),method = cormethod),digits=3),panel.first = grid(lty=1),xlim=rang,ylim=rang)
        abline(0,1)
        mtext(side=3,paste("sample ",samp,annots[match(samp,annots$sample_ID),]$prep,annots[match(samp,annots$sample_ID),]$prime,annots[match(samp,annots$sample_ID),]$library_chemistry))
        plot(reg+obs,reg+obs_exp[["exp"]][[samp]][genes,clust],log="xy",col=colv,ylab=paste("predicted cluster",clust),xlab=paste("cluster",clust,"sample",samp),pch=20,cex=ifelse(colv!=1,1,.1),panel.first = grid(lty=1),main=round(cor(log10(reg+obs),log10(reg+ldm$model$models[genes,clust]*(1-ldm$dataset$alpha_noise[samp])+ldm$dataset$noise_models[genes,samp]*ldm$dataset$alpha_noise[samp]),method = cormethod),digits=3),xlim=rang,ylim=rang)
        abline(0,1)
        dev.off()
      }   
      mask_bigger2=obs_numis-4*sqrt(sum(obs_numis)*ldm$model$models[genes,clust])>sum(obs_numis)*obs_exp[["exp"]][[samp]][genes,clust] & log10((reg+obs)/(reg+obs_exp[["exp"]][[samp]][genes,clust]))>.5
      mask_smaller2=obs_numis+4*sqrt(sum(obs_numis)*ldm$model$models[genes,clust])<sum(obs_numis)*obs_exp[["exp"]][[samp]][genes,clust] & log10((reg+obs)/(reg+obs_exp[["exp"]][[samp]][genes,clust]))<(-.5)
      if (sum(mask_bigger2)>0){
        maskl[[samp]][mask_bigger2,clust]=+1
      }
      if (sum(mask_smaller2)>0){
        maskl[[samp]][mask_smaller2,clust]=-1
      }
    }
  }
  return(maskl)
}


obs_exp=get_obs_exp_stats(lung_ldm,lung_ldm$dataset$samples,genes=rownames(lung_ldm$model$models))
maskl=generate_projection_qc_plots(lung_ldm,lung_ldm$dataset$samples,obs_exp,genes=rownames(lung_ldm$model$models),save_plots = T)
maskl=generate_projection_qc_plots(lung_ldm,lung_ldm$dataset$samples,obs_exp,genes=rownames(lung_ldm$model$models),save_plots = F)


obs_exp=get_obs_exp_stats(zilionis_ldm,zilionis_ldm$dataset$samples,genes=genes)
generate_projection_qc_plots(zilionis_ldm,zilionis_ldm$dataset$samples,obs_exp,genes=genes)

e1=c()
e2=c()
for (samp in samples){
  print(samp)
  nclust=ncol(obs_exp[["obs"]][[samp]])
  rs=rowSums(lung_ldm$dataset$counts[samp,,])
  m=rs/sum(rs)
  #gmask=m>1e-6
  gmask=names(m)
  reg1=1e-6
  if (nclust<5){
    next
  }
  c1=sqrt(rowSums((log10(reg1+obs_exp[["obs"]][[samp]][gmask,])-log10(reg1+obs_exp[["exp"]][[samp]][gmask,]))^2))
  c2=sqrt(rowSums((log10(reg1+obs_exp[["obs"]][[samp]][gmask,])-log10(reg1+lung_ldm$model$models[gmask,colnames(obs_exp[["obs"]][[samp]])]^2))))
  e1=cbind(e1,c1)
  e2=cbind(e2,c2)
  png(paste("output/figures/qc_figures/model_errors_",samp,".png",sep=""),400,400)
  plot(c2,c1,pch=".",xlab="error cell type model",ylab="error full model")
  dev.off()
}
colnames(e1)=samples
colnames(e2)=samples
c1=c1[is.finite(c1)]

obs_exp_plot=function(gene,samp){
  layout(matrix(1:2,2,1))
  ylim=c(-7,7)
  xlim=c(-7,-2)
  plot(log10(obs_exp[["exp"]][[samp]][gene,]),log2(obs_exp[["obs"]][[samp]][gene,]/obs_exp[["exp"]][[samp]][gene,]),ylim=ylim,xlim=xlim,panel.first={grid(lty=1);abline(h=0)},ylab="obs/exp",xlab="exp by full model",main=paste(gene,samp))
  plot(log10(lung_ldm$model$models[gene,colnames(obs_exp[["obs"]][[samp]])]),log2(obs_exp[["obs"]][[samp]][gene,]/lung_ldm$model$models[gene,colnames(obs_exp[["obs"]][[samp]])]),panel.first={grid(lty=1);abline(h=0)},ylim=ylim,,xlim=xlim,ylab="obs/exp",xlab="exp by model")
}
obs_exp_plot(gene="CEBPD",samp="120")

