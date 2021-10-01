library(scales)

figure_1f <- function(){
  
annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)
sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)

# load("intermediates/lung_ldm_metadata_200327.rd")

# metadata <- cbind(metadata,
#                   sample_annots[as.character(metadata$sample_id),c("patient_ID","tissue","prime")])
# metadata <- metadata[metadata$prime=="3",]
# 
# metadata$patient_tissue <- apply(metadata[,c("patient_ID","tissue")],1,paste,collapse="_")
# 
# metadata$sub_lineage <- annots_list[as.character(metadata$cluster_assignment),"sub_lineage"]
# metadata <- metadata[metadata$sub_lineage!="",]

if(!exists("lung_ldm")){
  load("data/lung_ldm.rd")
}

tab <- table(apply(sample_annots[as.character(lung_ldm$dataset$cell_to_sample),c("patient_ID","tissue")],1,paste,collapse="_"),
             annots_list[as.character(lung_ldm$dataset$cell_to_cluster),"sub_lineage"])

#tab <- table(metadata$patient_tissue,metadata$sub_lineage)
tab <- tab[,-1]

tab <- tab/rowSums(tab)

tab <- log10(tab+1e-3)

dists <- matrix(NA,nrow(tab),nrow(tab),dimnames=list(rownames(tab),rownames(tab)))

for(row in rownames(dists)){
  for(col in colnames(dists)){
    dists[row,col] <- sqrt(sum((tab[row,]-tab[col,])^2))
  }
}

diag(dists) <- NA

dists_tumor <- dists[grepl("Tumor",rownames(dists)),grepl("Tumor",colnames(dists))]
dists_tumor <- dists_tumor[upper.tri(dists_tumor,diag=F)]
dists_normal <- dists[grepl("Normal",rownames(dists)),grepl("Normal",colnames(dists))]
dists_normal <- dists_normal[upper.tri(dists_normal,diag=F)]

dists_cross <- dists[grepl("Tumor",rownames(dists)),grepl("Normal",colnames(dists))]
tumors <- unlist(lapply(strsplit(rownames(dists_cross),"_"),function(x){x[1]}))
normals <- unlist(lapply(strsplit(colnames(dists_cross),"_"),function(x){x[1]}))
for(norm_iter in normals){
  dists_cross[which(tumors==norm_iter),which(normals==norm_iter)] <- NA
}
dists_cross <- array(dists_cross)

dists <- c(dists_normal,dists_tumor,dists_cross)
labs <- rep(c("normal","tumor","cross"),times=c(length(dists_normal),length(dists_tumor),length(dists_cross)))

png("output/figures/figure_1f.png",height=4,width=5,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0),oma=c(0,3,0,0))
boxplot(dists~factor(labs,c("cross","tumor","normal")),range=0,col=alpha(c("white","brown","blue"),0.4),yaxt="n",
        xlab="Euclidean distance",ylim=c(min(dists,na.rm=T),max(dists,na.rm=T)+1),horizontal=T)
points(dists,jitter(as.integer(factor(labs,c("cross","tumor","normal")))),pch='.')
mtext("Sample-sample distances",cex=1.5)
segments(y0=1,y1=2,x0=6.5,x1=6.5)
segments(y0=2,y1=3,x0=6,x1=6)
text(y=c(1.5,2.5),x=c(6.7,6.2),c("***","***"),srt=90)
mtext(c("Tumor-nLung\ndistanes","Tumor-Tumor\ndistances","nLung-nLung\ndistances"),side=2,at=1:3,las=2,line=0.25,
      col=alpha(c(1,"brown","blue"),0.8))

dev.off()
}
