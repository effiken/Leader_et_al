##############
# Download all the TCGA data
###############


# Downloads and compiles TCGA gene expression, mutation, and clinical data


library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicDataCommons)
library(data.table)
library(R.utils)
library(scales)


#### USER MUST SET THEIR OWN DOWNLOAD DIRECTORY
## Note: In windows, if this path is too long, it will cause an error for GDCdownload()
working_TCGA_data_dir <- "/data/TCGA/"





source("scripts/get_LCAM_scores.R")

## functions:

TCGAtranslateID = function(file_ids, legacy = FALSE) {
  info = files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

## Start interactive workflow:

tcga_abbreviations <- read.csv("input_tables/TCGA_abbreviations.csv",stringsAsFactors=F,h=F,fileEncoding="UTF-8-BOM")

tcga_abbreviations[tcga_abbreviations=="Brain Lower Grade Glioma"] <- "Low Grade Glioma"
tcga_abbreviations[tcga_abbreviations=="Breast invasive carcinoma"] <- "Breast Cancer"
tcga_abbreviations[tcga_abbreviations=="Colon adenocarcinoma"] <- "colorectal adenocarcinoma"
tcga_abbreviations[tcga_abbreviations=="Kidney Chromophobe"] <- "kidney chromophobe renal cell carcinoma"
tcga_abbreviations[tcga_abbreviations=="Pancreatic adenocarcinoma"] <- "Pancreatic ductal adenocarcinoma"
tcga_abbreviations[tcga_abbreviations=="Pheochromocytoma and Paraganglioma"] <- "Paraganglioma and Pheochromocytoma"
tcga_abbreviations[tcga_abbreviations=="Rectum adenocarcinoma"] <- "colorectal adenocarcinoma"
tcga_abbreviations[tcga_abbreviations=="Cervical squamous cell carcinoma and endocervical adenocarcinoma"] <- "Cervical Carcinoma"
tcga_abbreviations[tcga_abbreviations=="Thyroid carcinoma"] <- "Thyroid papillary carcinoma"

exclude_projects <- c("Controls","FFPE Pilot Phase II",
                      "Miscellaneous","Chronic Myelogenous Leukemia")

tcga_abbreviations <- tcga_abbreviations[!tcga_abbreviations[,2]%in%exclude_projects,]

project_IDs <- tcga_abbreviations[,1]
#project_IDs <- project_IDs[1:3]


if(!dir.exists(working_TCGA_data_dir)){
  dir.create(working_TCGA_data_dir)
}
github_repo_dir <- getwd()
setwd(working_TCGA_data_dir)


# a translation table from ensemble IDs to gene symbols. Taken from 10X workflow so that they are as consistent as possible.
gene_ids <- read.table(file.path(github_repo_dir,"input_tables/ensemble_ids.tsv"),sep="\t",r=1,h=0,stringsAsFactors = F)


# THESE LOOPS DOWNLOADS ALL THE TCGA DATA
# YOU MAY NEED TO RUN IT MULTIPLE TIMES FOR THE DATA TO COMPLETE

for(project in project_IDs){
  if(!file.exists(file.path(working_TCGA_data_dir,paste(project,"_exprs.rd",sep="")))){
  #get transcriptome data
  q <- GDCquery(project = paste("TCGA",project,sep="-"),
                data.category = "Transcriptome Profiling",
                data.type="Gene Expression Quantification",
                workflow.type = "HTSeq - FPKM",experimental.strategy = "RNA-Seq",
                legacy = F)
  tmp <- try({GDCdownload(q,directory = working_TCGA_data_dir,method="api",files.per.chunk = 10)})
  try(GDCdownload(q,directory = working_TCGA_data_dir,method="api",files.per.chunk = 10))
  GDCdownload(q,directory = working_TCGA_data_dir,method="api",files.per.chunk = 10)
  
  uuids <- list.dirs(file.path(working_TCGA_data_dir,paste("TCGA-",project,sep=""),"harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/"))
  uuids <- unlist(lapply(strsplit(uuids,"/"),function(x){tail(x,1)}))
  uuids <- setdiff(uuids,"Gene_Expression_Quantification")
  submitter_ids <- TCGAtranslateID(uuids)
  
  uuids <- list.dirs(file.path(working_TCGA_data_dir,paste("TCGA-",project,sep=""),"harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/"))
  uuids <- uuids[-1]
  fn <- uuids[1]
  fn <- file.path(fn,list.files(fn))
  expr_vec <- read.table(gzfile(fn),sep="\t",stringsAsFactors = F)
  
  expr_mat <- matrix(NA,nrow=nrow(expr_vec),ncol=nrow(submitter_ids),
                     dimnames=list(expr_vec$V1,submitter_ids$file_id))
  message("reading files")
  for(fns in uuids){
    uuid <- strsplit(fns,"/")
    uuid <- tail(uuid[[1]],1)
    fn <- file.path(fns,list.files(fns))
    expr_vec <- read.table(gzfile(fn),sep="\t",stringsAsFactors = F)
    expr_mat[,uuid] <- expr_vec$V2
  }
  colnames(expr_mat) <- submitter_ids[match(colnames(expr_mat),submitter_ids$file_id),]$submitter_id
  
  rownames(expr_mat) <- unlist(lapply(strsplit(rownames(expr_mat),"\\."),function(x){x[1]}))
  
  expr_mat <- expr_mat[intersect(rownames(expr_mat),gene_ids$V2),]
  rownames(expr_mat) <- rownames(gene_ids)[match(rownames(expr_mat),gene_ids$V2)]
  
  message("saving expr")
  save(expr_mat,file=paste(project,"exprs.rd",sep="_"))
  message("removing read files")
  unlink(file.path(working_TCGA_data_dir,paste("TCGA-",project,sep="")),recursive=T)
  }
}

#patient_barcodes <- TCGAtranslateID(uuids)

# get MAF files
for(project in project_IDs){
  if(!dir.exists(file.path(working_TCGA_data_dir,"GDCdata",paste("TCGA-",project,sep="")))){
  try(GDCquery_Maf(tumor=project,pipelines="mutect2"))
    GDCquery_Maf(tumor=project,pipelines="mutect2")
  }
}

# get clinical data
for(project in project_IDs){
  clin <-  GDCquery_clinic(project=paste("TCGA",project,sep="-"),type="clinical",save.csv=TRUE)
  }

###############
#Download all the ep_scores for the TCGA data

not_in_ESTIMATE <- c("Cholangiocarcinoma",
                     "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                     "Mesothelioma",
                     "Sarcoma",
                     "Testicular Germ Cell Tumors",
                     "Thymoma",
                     "Uveal Melanoma")

ESTIMATE_IDs <- tcga_abbreviations[!tcga_abbreviations[,2]%in%not_in_ESTIMATE,1]

for(project_ID in ESTIMATE_IDs){
  project_name <- tcga_abbreviations[match(project_ID,tcga_abbreviations[,1]),2]
  url <- paste("https://bioinformatics.mdanderson.org/estimate/tables/",
               tolower(gsub(" ","_",project_name)),"_RNAseqV2.txt",sep="")
  destfile <- file.path(working_TCGA_data_dir,paste("estimate_",gsub(" ","_",project_name),".txt",sep=""))
  if(!file.exists(destfile)){
    download.file(url,destfile = file.path(working_TCGA_data_dir,paste("estimate_",gsub(" ","_",project_name),".txt",sep="")))

  }
}

# Compute LCAM scores and n_mut

output_list <- list()

for(project_ID in project_IDs){
  message(project_ID)
  expr_fn <- paste(project_ID,"exprs.rd",sep="_")
  load(expr_fn)
  cols <- unique(colnames(expr_mat)[substr(colnames(expr_mat),14,15)%in%c("01","03") &
    substr(colnames(expr_mat),16,16)%in%c("A","-")])
  output_list[[project_ID]] <- get_LCAM_scores(expr_mat[,cols])
  
  estimate_fn <- file.path(working_TCGA_data_dir,paste("estimate_",gsub(" ","_",tcga_abbreviations[match(project_ID,tcga_abbreviations[,1]),2]),".txt",sep=""))
  
  if(file.exists(estimate_fn)){
    message("reading ESTIMATE scores")
    estimate <- read.csv(estimate_fn,stringsAsFactors = F,sep="\t")
    output_list[[project_ID]] <- cbind(output_list[[project_ID]],
                                       estimate[match(substr(rownames(output_list[[project_ID]]),1,15),estimate$ID),"Immune_score"])
    colnames(output_list[[project_ID]])[4] <- "ESTIMATE_Immune"
  }
  
  message("loading mutation data")
  mut_fn <- list.files(file.path(working_TCGA_data_dir,"GDCdata",paste("TCGA-",project_ID,sep="")),r=T)
  mut_fn <- file.path(working_TCGA_data_dir,"GDCdata",paste("TCGA-",project_ID,sep=""),mut_fn)
  file.copy(from=mut_fn,to=file.path(working_TCGA_data_dir,"staged_mut_file.maf.gz"),overwrite=T)
  
  mut_data <- fread(file.path(working_TCGA_data_dir,"staged_mut_file.maf.gz"))
  mut_data <- mut_data[mut_data$Variant_Classification%in%c("Missense_Mutation","Nonsense_Mutation")]
  mut_tab <- table(mut_data$Hugo_Symbol,mut_data$Tumor_Sample_Barcode)
  colnames(mut_tab) <- substr(colnames(mut_tab),1,15)
  n_mut <- colSums(mut_tab)
  
  if(!"n_muts"%in%colnames(output_list)){
  output_list[[project_ID]] <- cbind(output_list[[project_ID]],n_mut[substr(rownames(output_list[[project_ID]]),1,15)])
  colnames(output_list[[project_ID]])[ncol(output_list[[project_ID]])] <- "n_muts"
  }
  
  mut_data <- mut_data[!mut_data$Variant_Classification%in%c("RNA","Silent"),]
  mut_tab <- table(mut_data$Hugo_Symbol,mut_data$Tumor_Sample_Barcode)
  colnames(mut_tab) <- substr(colnames(mut_tab),1,15)
  
  not_in_tab <- !c("TP53","KRAS","EGFR","STK11")%in%rownames(mut_tab)
  driver_mut_tab <- as.matrix(mut_tab[setdiff(c("TP53","KRAS","EGFR","STK11"),c("TP53","KRAS","EGFR","STK11")[not_in_tab]),])
  if(sum(!not_in_tab)==1){
    driver_mut_tab <- t(driver_mut_tab)
    rownames(driver_mut_tab) <- c("TP53","KRAS","EGFR","STK11")[!not_in_tab]
  }
  driver_mut_tab <- rbind(driver_mut_tab,matrix(0,nrow=sum(not_in_tab),ncol=ncol(driver_mut_tab),
                                                dimnames=list(c("TP53","KRAS","EGFR","STK11")[not_in_tab],NULL)))
  no_mut_pats <- setdiff(substr(rownames(output_list[[project_ID]]),1,15),colnames(mut_tab))
  driver_mut_tab <- cbind(driver_mut_tab,matrix(NA,nrow=nrow(driver_mut_tab),ncol=length(no_mut_pats),
                                                dimnames=list(c("TP53","KRAS","EGFR","STK11"),no_mut_pats)))
  
  output_list[[project_ID]] <- cbind(output_list[[project_ID]],t(driver_mut_tab)[substr(rownames(output_list[[project_ID]]),1,15),])
  colnames(output_list[[project_ID]])[(ncol(output_list[[project_ID]])-3):ncol(output_list[[project_ID]])] <- c("TP53","KRAS","EGFR","STK11")
  

}
  
save(output_list,file=file.path(github_repo_dir,"output/statistics/TCGA_results.rd"))
#load("/users/Andrew Leader/Google Drive/merad/scRNAseq_analysis/results/AL/lung_main/output/TCGA_results.rd")
#Make plots

subtype_ord <- order(unlist(lapply(output_list,function(x){cor(x[,"difference"],x[,"n_muts"],method="spearman",use="p")})))

p_value <- unlist(lapply(output_list,function(x){
  cor.test(x[,"difference"],x[,"n_muts"],method="spearman",use="p")$p.value}))
q_value <- p.adjust(p_value,method="bonf")
cor_value <- unlist(lapply(output_list,function(x){cor.test(x[,"difference"],x[,"n_muts"],method="spearman",use="p")$estimate}))

write.csv(cbind(p_value,q_value),file.path(github_repo_dir,"output/statistics/figure_7f.csv"))

#cancer_var <- unlist(lapply(output_list,function(x){var(log10(1+x[,"n_muts"]),na.rm=T)}))

#plot(cancer_var,cor_value,pch=16,xlab="Var( Log10[1+#Mut] )",ylab="Rho( Log10[1+#Mut] , LCAM )",
#     main="Variance of #mutations vs. detected association with LCAM score",col="white")
#text(names(cancer_var),x=cancer_var,y=cor_value,col=1+as.numeric(names(cancer_var)%in%c("SKCM","LUAD","LUSC")),
#     font=1+as.numeric(names(cancer_var)%in%c("SKCM","LUAD","LUSC")))
#text(paste("rho=",round(cor(cancer_var,cor_value,method="spearman"),2),"; p=",signif(cor.test(cancer_var,cor_value,method="spearman")$p.value,2)),x=0.4,y=-0.1)

#plot(cancer_var,-log10(1e-16+q_value),pch=16,xlab="Var( Log10[1+#Mut] )",ylab="-log10(q value) of cor( Log10[1+#Mut] , LCAM )",
#     main="Variance of #mutations vs. detected association with LCAM score",col="white")
#text(names(cancer_var),x=cancer_var,y=-log10(1e-16+q_value),col=1+as.numeric(names(cancer_var)%in%c("SKCM","LUAD","LUSC")),
#     font=1+as.numeric(names(cancer_var)%in%c("SKCM","LUAD","LUSC")))
#text(paste("rho=",round(cor(cancer_var,-log10(q_value+1e-16),method="spearman"),2),"; p=",signif(cor.test(cancer_var,-log10(q_value+1e-16),method="spearman")$p.value,2)),x=0.4,y=14)






###############

# png(file.path(github_repo_dir,"output/figures/figure_7f.png"),
#     units="in",height=4,width=4,res=300,bg="transparent")
# par(mgp=c(2,1,0))
# plot(cor_value,-log10(q_value+1e-16),xlab=expression(rho),xlim=c(-.5,.5),ylab=expression(paste('-log10(p'['adj'],")",sep="")),col=alpha(1,0.5),pch=16)
# 
# abline(lty=2,col="red",h=2)
# text(expression(paste(alpha,"=0.05")),x=-0.4,y=4,col="red")
# mtext("TCGA tumor types",cex=1.5)
# dev.off()

png(file.path(github_repo_dir,"output/figures/figure_7f.png"),
    units="in",height=4,width=4,res=300,bg="transparent")
par(mgp=c(2,1,0))
plot(cor_value,-log10(q_value+1e-16),xlab=expression(rho),xlim=c(-.5,.5),ylab=expression(paste('-log10(p'['adj'],")",sep="")),col=alpha(1,0.5),pch=16)

abline(lty=2,col="red",h=2)
text(expression(paste(alpha,"=0.05")),x=-0.4,y=4,col="red")
#text(x=0.7,y=-2.5,paste("r=",signif(res$estimate,2),"; p=",signif(res$p.value,2),sep=""))
mtext("TCGA tumor types",cex=1.5)
text(cor_value,-log10(q_value+1e-16),names(output_list))

dev.off()



#############





residual_list <- lapply(output_list,function(x){
  mask <- 
  lm(x[,3]~log10(1+x[,"n_muts"]))$residuals
})

#sig_types <- names(q_value)[q_value < 0.05]

sig_types <- names(output_list)

#t_stat_mat <- matrix(NA,nrow=4,ncol=length(sig_types),dimnames=list(c("TP53","KRAS","EGFR","STK11"),sig_types))
t_stat_mat <- matrix(NA,nrow=2,ncol=length(sig_types),dimnames=list(c("TP53","KRAS"),sig_types))

t_pval_mat <- t_stat_mat
frac_mat <- t_stat_mat
count_mat <- t_stat_mat
for(col in colnames(t_stat_mat)){
  for(row in rownames(t_stat_mat)){
    t <- NA
    t <- try(t.test(residual_list[[col]]~output_list[[col]][names(residual_list[[col]]),row]==0))
    if(!is.null(names(t))){
      t_stat_mat[row,col] <- t.test(residual_list[[col]]~output_list[[col]][names(residual_list[[col]]),row]==0)$statistic
      t_pval_mat[row,col] <- t.test(residual_list[[col]]~output_list[[col]][names(residual_list[[col]]),row]==0)$p.value
      frac_mat[row,col] <- sum(output_list[[col]][,row],na.rm=T)/sum(!is.na(output_list[[col]][,"n_muts"]))
      count_mat[row,col] <- sum(output_list[[col]][,row],na.rm=T)
      if(count_mat[row,col] < 10){
        t_stat_mat[row,col] <- NA
        t_pval_mat[row,col] <- NA
      }
      
    }
  }
}

count_mat[is.na(count_mat)] <- 0


t_stat_mat <- t_stat_mat[,setdiff(colnames(t_stat_mat),"LUAD")]
t_pval_mat <- t_pval_mat[,setdiff(colnames(t_pval_mat),"LUAD")]

thorrsen <- read.csv(file.path(github_repo_dir,"input_tables/Thorsson_et_al_TCGA_tumor_immunology_metrics.csv"),
                     r=1,h=1,stringsAsFactors = F)
cta_rho <- array(NA,ncol(t_stat_mat),dimnames=list(colnames(t_stat_mat)))
cta_cor_pval <- cta_rho
for(type in names(cta_rho)){
  cor_res <- try(cor.test(residual_list[[type]],thorrsen[substr(names(residual_list[[type]]),1,12),]$CTA.Score,na.rm=T,method="spearman"))
  cta_rho[type] <-  try(cor_res$estimate)
  cta_cor_pval[type] <- try(cor_res$p.value)
}

res_mat <- rbind(t_stat_mat,as.numeric(cta_rho))
pmat <- rbind(t_pval_mat,as.numeric(cta_cor_pval))

res_vec <- c(pmat[1,],pmat[2,],pmat[3,])
names(res_vec) <- paste(names(res_vec),rep(c("TP53","KRAS","CTA"),each=ncol(pmat)),sep=".")
res_vec <- res_vec[!is.na(res_vec)]
write.csv(cbind(res_vec,p.adjust(res_vec,method="bonf")),file.path(github_repo_dir,"output/statistics/figure_s7cd.csv"))
# 
# png(file.path(github_repo_dir,"output/figures/figure_s7c.png"),
#     units="in",height=4,width=15,res=300,bg="transparent",pointsize=16)
# par(mfrow=c(1,3))
# par(mgp=c(2,1,0))
# 
# plot(res_mat[1,],-log10(1e-4+pmat[1,]*sum(!is.na(res_mat))),ylim=c(-2,4),pch=16,col=alpha(1,.6),xlim=c(-10,10),xlab="t-statistic",
#      ylab=expression(paste('-log10(p'['adj'],")",sep="")))
# abline(h=-log10(1e-4+0.05),col="red",lty=2)
# text(expression(paste(alpha,"=0.05")),x=-8,y=2,col="red")
# mtext("TP53",cex=1.5)
# 
# 
# plot(res_mat[2,],-log10(1e-4+pmat[2,]*sum(!is.na(res_mat))),ylim=c(-2,4),pch=16,col=alpha(1,.6),xlim=c(-10,10),xlab="t-statistic",
#      ylab=expression(paste('-log10(p'['adj'],")",sep="")))
# abline(h=-log10(1e-4+0.05),col="red",lty=2)
# text(expression(paste(alpha,"=0.05")),x=-8,y=2,col="red")
# mtext("KRAS",cex=1.5)
# plot(res_mat[3,],-log10(1e-4+pmat[3,]*sum(!is.na(res_mat))),ylim=c(-2,4),pch=16,col=alpha(1,.6),xlim=c(-.4,.4),xlab=expression(rho),
#      ylab=expression(paste('-log10(p'['adj'],")",sep="")))
# abline(h=-log10(1e-4+0.05),col="red",lty=2)
# text(expression(paste(alpha,"=0.05")),x=-0.3,y=2,col="red")
# mtext("CTA score",cex=1.5)
# 
# dev.off()

png(file.path(github_repo_dir,"output/figures/figure_s7cd.png"),
    units="in",height=4,width=15,res=300,bg="transparent",pointsize=16)
par(mfrow=c(1,3))
par(mgp=c(2,1,0))

plot(res_mat[1,],-log10(1e-4+pmat[1,]*sum(!is.na(res_mat))),ylim=c(-2,4),pch=16,col=alpha(1,.6),xlim=c(-10,10),xlab="t-statistic",
     ylab=expression(paste('-log10(p'['adj'],")",sep="")))
abline(h=-log10(1e-4+0.05),col="red",lty=2)
text(expression(paste(alpha,"=0.05")),x=-8,y=2,col="red")
mtext("TP53",cex=1.5)
text(colnames(res_mat),x=res_mat[1,],y=-log10(1e-4+pmat[1,]*sum(!is.na(res_mat))))

plot(res_mat[2,],-log10(1e-4+pmat[2,]*sum(!is.na(res_mat))),ylim=c(-2,4),pch=16,col=alpha(1,.6),xlim=c(-10,10),xlab="t-statistic",
     ylab=expression(paste('-log10(p'['adj'],")",sep="")))
abline(h=-log10(1e-4+0.05),col="red",lty=2)
text(expression(paste(alpha,"=0.05")),x=-8,y=2,col="red")
mtext("KRAS",cex=1.5)
plot(res_mat[3,],-log10(1e-4+pmat[3,]*sum(!is.na(res_mat))),ylim=c(-2,4),pch=16,col=alpha(1,.6),xlim=c(-.4,.4),xlab=expression(rho),
     ylab=expression(paste('-log10(p'['adj'],")",sep="")))
abline(h=-log10(1e-4+0.05),col="red",lty=2)
text(expression(paste(alpha,"=0.05")),x=-0.3,y=2,col="red")
mtext("CTA score",cex=1.5)
text(colnames(res_mat),x=res_mat[3,],y=-log10(1e-4+pmat[3,]*sum(!is.na(res_mat))))


dev.off()
# 
# #text(x=0.7,y=-2.5,paste("r=",signif(res$estimate,2),"; p=",signif(res$p.value,2),sep=""))
# mtext("TCGA tumor types",cex=1.5)
# text(cor_value,-log10(q_value+1e-16),names(output_list))
# 
# dev.off()


# 
# par(mfrow=c(2,2))
# for(row in rownames(t_stat_mat)[1:2]){
#   frac <- round(unlist(lapply(output_list[sig_types],function(x){sum(x[,row]/sum(!is.na(x[,"n_muts"])),na.rm=T)})),2)
#   labs <- paste(names(frac),frac,"-")
#   labs <- names(frac)
#   # plot(t_stat_mat[row,],-log10(1e-4+t_pval_mat[row,]*sum(!is.na(t_pval_mat))),pch="",main=row,ylim=c(-2,4),xlim=c(-10,10),
#   #      ylab="-log10(q-value)",xlab="t_statistic")
#   x <- t_stat_mat[row,]
#   y <- -log10(1e-4+t_pval_mat[row,]*sum(!is.na(t_stat_mat)))
#   plot(x,y,pch="",main=row,ylim=c(-2,4),xlim=c(-10,10),
#        ylab="-log10(q-value)",xlab="t_statistic")
#   text(x=x,y=y,names(x),xpd=NA)
#   abline(h=-log10(1e-4+0.05),col="red")
#  
# }
# 

thorrsen <- read.csv(file.path(github_repo_dir,"input_tables/Thorsson_et_al_TCGA_tumor_immunology_metrics.csv"),
                     r=1,h=1,stringsAsFactors = F)
cta_rho <- array(NA,ncol(t_stat_mat),dimnames=list(colnames(t_stat_mat)))
cta_cor_pval <- cta_rho
for(type in names(cta_rho)){
  cor_res <- try(cor.test(residual_list[[type]],thorrsen[substr(names(residual_list[[type]]),1,12),]$CTA.Score,na.rm=T,method="spearman"))
  cta_rho[type] <-  try(cor_res$estimate)
  cta_cor_pval[type] <- try(cor_res$p.value)
}

x <- as.numeric(cta_rho)
y <- -log10(as.numeric(cta_cor_pval)*15)
plot(x,y,pch="",xlab="Rho(residuals, CTA score)",ylab="-log10(padj)")
text(x=x,y=y,label=names(cta_rho),xpd=NA)





# 
# genes <- strsplit("IGHG3,IGHG4,IGHG1,MZB1,FAM92B,SPAG4,DERL3,JSRP1,TBCEL,LINC01485,SLC2A5,SPP1,CCL7,HAMP,CXCL13,ZBED2,GNG4,KRT86,TNFRSF9,LAYN,GZMB,KIR2DL4,GAPT,CTB-133G6.1,AKR1C3,RSPO3,MCEMP1,RND3,HP,FOLR3,GPD1,GS1-600G8.5,PCOLCE2,CAMP,CCL17,CD1B,CD1C,CD1E,PLD4,TRPC6,CLEC10A,FCER1A,PLA1A,CFP,CEACAM8,AZU1,PKP2,RETN,LYZ,ELANE,ATP6AP1L,RP6-159A1.4,THBS1,CFD,P2RY14,DNASE1L3,PTGER3,CST3,BATF3,CPVL,RGS18,SERPINF2,LGALS2",",")[[1]]
# 
# 
# figure_path <- "/users/andrew leader/google drive/merad/scRNAseq_analysis/results/AL/lung_main/output/figures/revisions"
# dir.create(figure_path)
# dir.create(file.path(figure_path,"TCGA"))
# 
# for(iter in subtype_ord){
#   
#   #plot heatmap
#   project <- names(output_list)[iter]
#   load(file.path(working_TCGA_data_dir,paste(project,"_exprs.rd",sep="")))
#   patient_ord <- rev(rownames(output_list[[project]])[order(output_list[[project]][,3])])
#   
#   mat <- expr_mat[genes,patient_ord]
#   
#   mat <- t(t(mat)/rowSums(t(mat)))
#   mat <- log10(1e-6+mat)
#   mat <- t(scale(t(mat)))
#   mat[mat > 2] <- 2
#   mat[mat < -2] <- -2
#   pdf(file.path(figure_path,"TCGA",paste(project,"_heatmap.pdf",sep="")),height=3,width=5,pointsize=6)
#   pimage(t(mat),axes="x")
#   dev.off()
#   
#   #plot correlation by group
#   if(colnames(output_list[[project]])[4]=="ESTIMATE_Immune"){
#   ep_group <- cut(output_list[[project]][,4],quantile(output_list[[project]][,4],seq(0,1,.1),na.rm=T))
#   tumor_split <- split(rownames(output_list[[project]]),ep_group)
#   
#   split_cor <- lapply(tumor_split,function(x){
#     return(cor.test(output_list[[project]][x,1],output_list[[project]][x,2],method="spearman"))})
#   #split_cor <- unlist(split_cor)
#   
#   cor_list <- unlist(lapply(split_cor,function(x){x$estimate}))
#   lower_bound <- unlist(lapply(split_cor,function(x){x$conf.int[[1]]}))
#   upper_bound <- unlist(lapply(split_cor,function(x){x$conf.int[[2]]}))
#   
#   
#   #with spearman correlation:
#   
#   sig <- 1.96/sqrt(unlist(lapply(tumor_split,length))-3)
#   lower_bound <- tanh(atanh(unlist(lapply(split_cor,function(x){x$estimate})))-sig)
#   upper_bound <- tanh(atanh(unlist(lapply(split_cor,function(x){x$estimate})))+sig)
#   
#   
#   png(file.path(figure_path,"TCGA",paste(project,"_quantile_cor.png",sep="")),height=1.4,width=1.7,units="in",res=300,pointsize=6,bg="transparent")
#   par(mgp=c(2,1,0))
#   plot(1:10,cor_list,cex=0,ylim=c(-1,1)*max(abs(c(lower_bound,upper_bound))),xlim=c(0.5,10.5),xaxt="n",xlab="",
#        ylab="Spearman rho\n(LCAMhi score, LCAMlo score)")
#   segments(x0=1:10-0.15,x1=1:10+0.15,y0=cor_list,y1=cor_list,lwd=3)
#   segments(x0=1:10,x1=1:10,y0=lower_bound,y1=upper_bound)
#   abline(h=0,col="red",lty=2)
#   mtext(side=1,at=1:10,1:10)
#   mtext(side=1,line=2,"Decile of immune infiltrate")
#   mtext("LCAM-hi/lo summary scores")
#   dev.off()
#   }
#   png(file.path(figure_path,"TCGA",paste(project,"_mut_cor.png",sep="")))
#   plot(log10(1+output_list[[project]][,"n_muts"]),output_list[[project]][,3])
#   dev.off()
# }
# 
# 
# 
# cor_value <- unlist(lapply(output_list,function(x){cor.test(x[,"difference"],x[,"n_muts"],method="spearman",use="p")$estimate}))
# q_value <- unlist(lapply(output_list,function(x){cor.test(x[,"difference"],x[,"n_muts"],method="spearman",use="p")$p.value}))*length(output_list)
# 
# plot(cor_value,-log10(q_value+1e-10),xlab="rho",xlim=c(-.5,.5),ylab="-log10[bonf. adj. p]",pch="")
# text(cor_value,-log10(q_value+1e-10),names(output_list),col=1+as.numeric(tcga_abbreviations$type=="A"))
# abline(lty=2,col="red",h=2)
# #sorting by the above correlation value, draw heatmap, correlation by ep_decile, and logTMB vs. LCAM.
# 
# cor_value <- unlist(lapply(output_list,function(x){cor.test(x[,"LCAMhi"],x[,"n_muts"],method="spearman",use="p")$estimate}))
# q_value <- unlist(lapply(output_list,function(x){cor.test(x[,"LCAMhi"],x[,"n_muts"],method="spearman",use="p")$p.value}))*length(output_list)
# 
# plot(cor_value,-log10(q_value+1e-10),xlab="rho",xlim=c(-.5,.5),ylab="-log10[bonf. adj. p]",pch="")
# text(cor_value,-log10(q_value+1e-10),names(output_list),col=1+as.numeric(tcga_abbreviations$type=="A"))
# abline(lty=2,col="red",h=2)
# 
# 
# cor_value <- unlist(lapply(output_list,function(x){cor.test(x[,"LCAMlo"],x[,"n_muts"],method="spearman",use="p")$estimate}))
# q_value <- unlist(lapply(output_list,function(x){cor.test(x[,"LCAMlo"],x[,"n_muts"],method="spearman",use="p")$p.value}))*length(output_list)
# 
# plot(cor_value,-log10(q_value+1e-10),xlab="rho",xlim=c(-.5,.5),ylab="-log10[bonf. adj. p]",pch="")
# text(cor_value,-log10(q_value+1e-10),names(output_list),col=1+as.numeric(tcga_abbreviations$type=="A"))
# abline(lty=2,col="red",h=2)
# 
# tcga_abbreviations$type <- "NA"
# adenos <- c("ACC","BRCA","COAD","LUAD","OV","PAAD","PRAD","READ","STAD","UCEC")
# tcga_abbreviations$type[tcga_abbreviations$V1%in%adenos] <- "A"
# 
# plot(cor_value,-log10(q_value+1e-10),xlab="rho",xlim=c(-.5,.5),ylab="-log10[bonf. adj. p]",pch="",col=1+as.numeric(tcga_abbreviations$type=="A"))
# text(cor_value,-log10(q_value+1e-10),names(output_list),col=1+as.numeric(tcga_abbreviations$type=="A"))

library(scales)
#plot(log10(1+output_list$LUAD[,"n_muts"]),output_list$LUAD[,3],pch=16,col=alpha(1,.4),xlab="LogTMB",ylab="LCAMhi-LCAMlo")
#points(log10(1+output_list$LUSC[,"n_muts"]),output_list$LUSC[,3],col=alpha("red",.4),pch=16)     


 png(file.path(github_repo_dir,"output/figures/figure_7e.png"),height=4,width=4,units="in",res=300,bg="transparent")
par(mgp=c(2,1,0))
plot(log10((1+output_list$LUAD[,"n_muts"])/48.2),output_list$LUAD[,3],pch=16,col=alpha(1,.4),
      xlab=expression(paste("Log"["10"],"[TMB/Mb]",sep="")),
      ylab="LCAM score",cex=0.5)
points(log10((1+output_list$LUSC[,"n_muts"])/48.2),output_list$LUSC[,3],col=alpha("red",.4),pch=16,cex=0.5)     

dev.off()

