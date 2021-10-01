
figure_s2ab <- function(){

if(!exists("lung_ldm")){
  load("data/lung_ldm.rd")
}

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",h=1,r=1,stringsAsFactors = F)


tab <- table(lung_ldm$dataset$cell_to_cluster)[c("52","29","30","14","41","45")]

png("output/figures/figure_s2a.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
barplot(tab,col=rgb(t(col2rgb(c(2,2,2,3,6,"orange"))*.7/255)),
        names.arg=c("MoDC-52","MoDC-29","MoDC-30","cDC2","mregDC","cDC1"),las=2)
mtext("Total cells per cluster")
dev.off()

cell_to_tissue <- sample_annots[lung_ldm$dataset$cell_to_sample,"tissue"]
cell_to_patient <- sample_annots[lung_ldm$dataset$cell_to_sample,]$patient_ID
names(cell_to_patient)=names(cell_to_tissue) <- names(lung_ldm$dataset$cell_to_sample)
tab <- table(lung_ldm$dataset$cell_to_cluster[cell_to_tissue=="Tumor"],cell_to_patient[cell_to_tissue=="Tumor"])

tab <- tab[c("52","29","30","14","41","45"),]

png("output/figures/figure_s2b.png",height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
boxplot.matrix(t(tab)+1,log="y",range=0,names=c("MoDC-52","MoDC-29","MoDC-30","cDC2","mregDC","cDC1"),
               col=rgb(t(col2rgb(c(2,2,2,3,6,"orange"))*.7/255)) ,las=2)
x <- matrix(1:6,nrow=6,ncol=ncol(tab))
points(jitter(x),tab+1,pch=16)
mtext("Cells per patient per cluster",cex=1.5)
dev.off()

}