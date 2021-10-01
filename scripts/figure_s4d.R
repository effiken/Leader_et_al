

figure_s4d <- function(){

if(!exists("lung_ldm")){
  load("data/lung_ldm.rd")
}

cells <- names(lung_ldm$dataset$cell_to_cluster)[lung_ldm$dataset$cell_to_cluster=="44"]

cells <- cells[sample_annots[lung_ldm$dataset$cell_to_sample[cells],]$prep=="beads" &
                 sample_annots[lung_ldm$dataset$cell_to_sample[cells],]$library_chemistry=="V2" &
                 sample_annots[lung_ldm$dataset$cell_to_sample[cells],]$tissue=="Tumor"]

sample_to_patient <- sample_annots[lung_ldm$dataset$samples,"patient_ID"]
names(sample_to_patient) <- lung_ldm$dataset$samples

s <- split(cells,sample_to_patient[lung_ldm$dataset$cell_to_sample[cells]])


scores <- unlist(lapply(s,function(x){
  numis <- colSums(lung_ldm$dataset$umitab[,x,drop=F])
  tfh_score <- colSums(lung_ldm$dataset$umitab[tfh_genes,x,drop=F])/numis
  texh_score <- colSums(lung_ldm$dataset$umitab[texh_genes,x,drop=F])/numis
  return(sum(log2((1e-4+tfh_score)/(1e-4+texh_score))>0))
}))

#names(scores) <- names(s)

#barplot(scores,las=2)

tab <- table(lung_ldm$dataset$cell_to_cluster,sample_to_patient[lung_ldm$dataset$cell_to_sample])
#tab <- tab/rowSums(tab)
#plot(tab["44",names(scores)]-scores,scores,log="xy",
#     xlab="#CD8 cells",ylab="#CD4 cells")

#sample_annots$li

#cor.test(tab["44",names(scores)],scores,method="spearman")

cd8_scores <- tab["44",names(scores)]-scores

 scores <- scores/colSums(tab[,names(scores)])
 cd8_scores <- cd8_scores/colSums(tab[,names(scores)])
 
 #plot(scores+1e-4,cd8_scores+1e-4,log="xy",xlab="Tfh fraction of cells",ylab="CD8 Tact fractin of cells",pch=16)

 cells <- intersect(cells,colnames(lung_ldm$dataset$ds[[1]]))
 numis <- colSums(lung_ldm$dataset$umitab[,cells])
 tfh_score <- colSums(lung_ldm$dataset$umitab[tfh_genes,cells,drop=F])/numis
 texh_score <- colSums(lung_ldm$dataset$umitab[texh_genes,cells,drop=F])/numis
 
 cells <- cells[order(log2((1e-4+tfh_score)/(1e-4+texh_score)))]
 
 

 genes <- strsplit("NMB,IL21,CD4,CD200,TCF7,TNFRSF4,G0S2,THADA,LIMS1,NABP1,ITPR1,TRAT1,CXCL13,IL7R,IL6ST,TNFAIP8,TMEM173,NR3C1,AIM1,TMEM243,TIMP1,SMARCA2,CD59,CORO1B,PGM2L1,UPF2,FAM107B,C10orf54,SLC2A3,FOS,MAF,FOSB,VCAM1,XCL2,XCL1,IVNS1ABP,PPM1G,CAPG,GNLY,CD8A,CD8B,BIN1,SRGAP3,CXCR6,HOPX,SPRY1,ITGA1,GZMA,HAVCR2,SERPINB9,HLA-DRB5,CHST12,TRGC2,CTSW,JAML,CRTAM,GATA3,PRF1,ENTPD1,GSTO1,CCND2,CD27,PTMS,LAG3,CLECL1,CLEC2B,GABARAPL1,KLRD1,KLRC1,KRT86,CD63,C12orf75,GZMH,GZMB,GPR65,ITGAE,TOB1,ACP5,HCST,NKG7,IL2RB",",")[[1]]
 genes <- strsplit("NKG7,GZMB,CTSW,GZMA,KLRD1,GNLY,CD63,PRF1,HCST,CD8A,KLRC1,HOPX,TRGC2,GZMH,CD8B,CLEC2B,KRT86,XCL2,VCAM1,TCF7,C10orf54,CD59,NABP1,TRAT1,TMEM243,SMARCA2,FAM107B,ITPR1,TIMP1,IL21,CD200,AIM1,THADA,G0S2,TNFAIP8,IL6ST,CD4,CORO1B,NMB,MAF,FOS,NR3C1,TMEM173,TNFRSF4,IL7R,CXCL13,SH2D1A",",")[[1]]
 
 
mat <- lung_ldm$dataset$ds[[1]][genes,cells] 
mat <- as.matrix(mat)
genes <- genes[order(cor(t(mat),order(log2((1e-4+tfh_score[colnames(mat)])/(1e-4+texh_score[colnames(mat)]))),method="spearman"))]
mat <- mat[genes,]

mat <- log2(1+mat)
thresh <- 2
mat[mat > thresh] <- thresh 
mat <- as.matrix(mat)

colgrad_rna <- colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100)


png("output/figures/figure_s4d.png",height=1.75,width=5,units="in",res=1000,
    pointsize=4)
image(mat,xaxt="n",yaxt="n",col=colgrad_rna)
mtext(genes,side=1,las=2,line=0.25,at=seq(0,1,1/(length(genes)-1)))
abline(h=sum(log2((1e-4+tfh_score)/(1e-4+texh_score))<0)/length(tfh_score))

dev.off()

}