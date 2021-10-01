library(scales)

figure_s1de <- function(){

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)
order_lists <- read.table("input_tables/cluster_order_list.txt",r=1,h=1,stringsAsFactors = F)

if(!exists("lung_ldm")){
  load("data/lung_ldm.rd")
}

alpha <- lung_ldm$dataset$alpha_noise

rm(lung_ldm)
load("data/lung_ldm_fiveprime.rd")
alpha <- c(alpha,lung_ldm_fiveprime$dataset$alpha_noise)

rm(lung_ldm_fiveprime)
load("data/lambrechts_ldm.rd")
alpha <- c(alpha,lambrechts_ldm$dataset$alpha_noise)
rm(lambrechts_ldm)
load("data/zilionis_ldm.rd")
alpha <- c(alpha,lung_ldm_zili$dataset$alpha_noise)

names(alpha)[grepl("filt",names(alpha))] <- unlist(lapply(strsplit(names(alpha)[grepl("filt",names(alpha))],"_"),function(x){x[2]}))

category <- array(NA,length(alpha),dimnames=list(names(alpha)))

category[names(lung_ldm_zili$model$alpha_noise)] <- "clustered V2"

category[!grepl("Lambrechts|zilionis",sample_annots[names(category),]$patient_ID) & sample_annots[names(category),]$library_chemistry=="V2" & is.na(category)] <-
  "projected V2"

category[sample_annots[names(category),]$Project=="Lambrechts"] <- "Lambrechts"
category[sample_annots[names(category),]$Project=="zilionis"] <- "Zilionis"
category[grepl("Lambrechts",sample_annots[names(category),]$patient_ID)] <- "Lambrechts"
category[grepl("zilionis",sample_annots[names(category),]$patient_ID)] <- "Zilionis"

category[sample_annots[names(category),]$prime=="5"] <- "5prime"

category <- category[names(alpha)]

png("output/figures/figure_s1d.png",height=4,width=4,units="in",res=300,bg="transparent")
par(oma=c(2,0,0,0),mgp=c(2,1,0))
boxplot(100*alpha~factor(category,c("clustered V2","projected V2","Lambrechts","5prime","Zilionis")),las=2,range=0,ylab="% noise")
points(jitter(as.integer(factor(category,c("clustered V2","projected V2","Lambrechts","5prime","Zilionis")))),
       100*alpha,pch=16,col=alpha(1,0.4))
mtext("% noise by batch type",cex=1.5)
dev.off()

load("data/lung_ldm.rd")

cluster_order <- c(strsplit(order_lists["lineage","order"],",")[[1]],
                   setdiff(strsplit(order_lists["overall","order"],",")[[1]],strsplit(order_lists["lineage","order"],",")[[1]]))
cluster_order <- cluster_order[!is.na(cluster_order)]

png("output/figures/figure_s1e.png",height=3,width=6,units="in",res=300,bg="transparent",pointsize=8)
par(mgp=c(2,1,0))
barplot(table(lung_ldm$dataset$cell_to_cluster)[cluster_order],col=rgb(0.7/255*t(col2rgb(1+as.integer(factor(annots_list[cluster_order,]$lineage,c("NK","T","MNP","pDC","B&plasma","mast","CD45- & doublets")))))),las=2,ylim=c(0,25000))
mtext("Cells",side=2,line=3)
dev.off()
}
