
figure_s1g <- function(){

annot_lists <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)
annot_lists$sub_lineage[annot_lists$sub_lineage==""] <- annot_lists$lineage[annot_lists$sub_lineage==""]

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
sample_annots <- sample_annots[grepl("zilionis",sample_annots$patient_ID),]

load("data/zilionis_ldm.rd")

cell_to_batch <- sample_annots[lung_ldm_zili$dataset$cell_to_sample,"metadata_indicator"]
barcode <- unlist(lapply(strsplit(colnames(lung_ldm_zili$dataset$umitab),"_"),function(x){x[2]}))

zili_name <- paste(cell_to_batch,barcode,sep="_")

zilionis_metadata <- read.table(gzfile("input_tables/GSE127465_human_cell_metadata_54773x25.tsv.gz"),h=1,stringsAsFactors = F,sep="\t")

rownames(zilionis_metadata) <- paste(zilionis_metadata$Library,zilionis_metadata$Barcode,sep="_")

tab <- as.matrix(table(annot_lists[lung_ldm_zili$dataset$cell_to_cluster,"sub_lineage"],zilionis_metadata[zili_name,]$Minor.subset))

tab <- tab[!rownames(tab)=="epi_endo_fibro_doublet",]

new_tab <- matrix(NA,nrow=nrow(tab),ncol=ncol(tab),dimnames=dimnames(tab))
new_tab[1:nrow(tab),] <- tab[1:nrow(tab),]
#new_tab <- t(new_tab)


new_tab <- new_tab/rowSums(new_tab)

#m[order(apply(m,1,which.max),]


tmp <- new_tab[order(apply(new_tab,1,which.max)),]
#tmp <- tmp[order(apply(tmp,1,which.max)),]

tmp <- tmp[rev(seq(nrow(tmp))),]
s <- rev(seq(0,1,1/100))

png("output/figures/figure_s1g.png",height=2.73,width=3.83,units="in",res=300,pointsize=6)
layout(matrix(1:2,nrow=1),widths=c(20,5))
par(mar=c(10,5,2,1))
image(as.matrix(t(tmp)),col=rgb(s,s,s),xaxt="n",yaxt="n")
box()

mtext(side=2,at=seq(0,1,1/(nrow(tmp)-1)),rownames(tmp),las=2,line=0.25,cex=0.5)
mtext(side=1,at=seq(0,1,1/(ncol(tmp)-1)),colnames(tmp),las=2,line=0.25,cex=0.5)

par(mar=c(10,5,2,2))
image(t(as.matrix(seq(length(s)))),col=rgb(s,s,s),xaxt="n",yaxt="n")
mtext(side=2,c(0,.5,1),at=c(0,.5,1),line=0.25,las=2,cex=0.75)
mtext(side=4,"Fraction of row",line=0.25)
dev.off()

}