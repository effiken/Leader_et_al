library(scales)

rm(list=ls())

wd <- "/users/andrew leader/google drive/merad/scRNAseq_analysis/results/AL/lung_main/"
setwd(wd)

plot_dir <- "output/figures/fig6_bulk_analysis/"
load("intermediates/TCGA_scores_revise.rd")
load("/users/andrew leader/google drive/merad/Leader_et_al_figure_scripts/Lambrechts_stromal/stromal_gene_lists.rd")

open_plot <- function(fn){
png(fn,height=4,width=4,units="in",res=300,bg="transparent")
  par(mgp=c(2,1,0))
}

stromal_scores <- lapply(gene_list,function(x){colMeans(z_mat[x,])})

stromal_scores <- do.call(cbind,stromal_scores)

cor(stromal_scores[names(score1),],score1-score2)


open_plot(file.path(plot_dir,"CAF_cor.png"))
plot(stromal_scores[names(score1),"fibroblast_tumor"],score1-score2,
     ylab="LCAM score",xlab="CAF score",pch=16,col=alpha(1,0.4))
mtext("CAF",cex=1.5)
res <- cor.test(stromal_scores[names(score1),"fibroblast_tumor"],score1-score2)
legend("bottomright",c(paste("r=",signif(res$estimate,2),sep=""),paste("p=",signif(res$p.value,2),sep="")),bg="white")
dev.off()

open_plot(file.path(plot_dir,"Lymphatics_cor.png"))
plot(stromal_scores[names(score1),"lymphatic"],score1-score2,
     ylab="LCAM score",xlab="Lymphatics score",pch=16,col=alpha(1,0.4))
mtext("Lymphatics",cex=1.5)
res <- cor.test(stromal_scores[names(score1),"lymphatic"],score1-score2)
legend("bottomright",c(paste("r=",signif(res$estimate,2),sep=""),paste("p=",signif(res$p.value,2),sep="")),bg="white")
dev.off()

open_plot(file.path(plot_dir,"Normal BEC_cor.png"))
plot(stromal_scores[names(score1),"blood_endo_normal"],score1-score2,
     ylab="LCAM score",xlab="Normal BEC score",pch=16,col=alpha(1,0.4))
mtext("Normal BEC",cex=1.5)
res <- cor.test(stromal_scores[names(score1),"blood_endo_normal"],score1-score2)
legend("bottomleft",c(paste("r=",signif(res$estimate,2),sep=""),paste("p=",signif(res$p.value,2),sep="")),bg="white")
dev.off()

open_plot(file.path(plot_dir,"Tumor BEC_cor.png"))
plot(stromal_scores[names(score1),"blood_endo_tumor"],score1-score2,
     ylab="LCAM score",xlab="Tumor BEC score",pch=16,col=alpha(1,0.4))
mtext("Tumor BEC",cex=1.5)
res <- cor.test(stromal_scores[names(score1),"blood_endo_tumor"],score1-score2)
legend("bottomleft",c(paste("r=",signif(res$estimate,2),sep=""),paste("p=",signif(res$p.value,2),sep="")),bg="white")
dev.off()

open_plot(file.path(plot_dir,"Normal Fibroblast_cor.png"))
plot(stromal_scores[names(score1),"fibroblast_normal"],score1-score2,
     ylab="LCAM score",xlab="Normal Fibroblast score",pch=16,col=alpha(1,0.4))
mtext("Normal Fibroblast",cex=1.5)
res <- cor.test(stromal_scores[names(score1),"fibroblast_normal"],score1-score2)
legend("bottomleft",c(paste("r=",signif(res$estimate,2),sep=""),paste("p=",signif(res$p.value,2),sep="")),bg="white")
dev.off()

open_plot(file.path(plot_dir,"CAF_cor.png"))
plot(stromal_scores[names(score1),"fibroblast_tumor"],score1-score2,
     ylab="LCAM score",xlab="CAF score",pch=16,col=alpha(1,0.4))
mtext("CAF",cex=1.5)
res <- cor.test(stromal_scores[names(score1),"fibroblast_tumor"],score1-score2)
legend("bottomright",c(paste("r=",signif(res$estimate,2),sep=""),paste("p=",signif(res$p.value,2),sep="")),bg="white")
dev.off()


open_plot(file.path(plot_dir,"CAF-normal_fib_cor.png"))
plot(stromal_scores[names(score1),"fibroblast_tumor"]-stromal_scores[names(score1),"fibroblast_normal"],score1-score2,
     ylab="LCAM score",xlab="CAF - Normal Fibroblast score",pch=16,col=alpha(1,0.4))
mtext("CAF-Normal Fibroblast",cex=1.5)
res <- cor.test(stromal_scores[names(score1),"fibroblast_tumor"]-stromal_scores[names(score1),"fibroblast_normal"],score1-score2)
legend("bottomright",c(paste("r=",signif(res$estimate,2),sep=""),paste("p=",signif(res$p.value,2),sep="")),bg="white")
dev.off()
