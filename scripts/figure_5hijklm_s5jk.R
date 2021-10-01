
library(scales)
library(plotrix)

figure_5hijklm_s5jk <- function(){

figure_path <- "output/figures/"

lcam_hi <- new.env()
load(envir=lcam_hi,"input_tables/lig_rec_lcamhi_vs_normal.rd")
lcam_lo <- new.env()
load(envir=lcam_lo,"input_tables/lig_rec_lcamlo_vs_normal.rd")
hi_vs_lo <- new.env()
load(envir=hi_vs_lo,"input_tables/lig_rec_lcamhi_vs_lcamlo.rd")

highlight_pairs <- list()
highlight_pairs$cDC <- list()
highlight_pairs$cDC$'T' <- c("CCL19_CCR7","CCL17_CCR4","CCL22_CCR4","CXCL9_CXCR3","CXCL10_CXCR3")
highlight_pairs$mac <- list()
highlight_pairs$mac$'T' <- c("IL10_IL10RA","IL10_IL10RB","OSM_IL6ST","CCL3_CCR5","IL6_IL6R",
                             "TNF_TNFRSF1B","CXCL9_CXCR3","CXCL16_CXCR6",
                             "TNF_TNFRSF1B")
highlight_pairs$mac$cDC <- c("IL10_IL10RA","IL10_IL10RB","OSM_IL6ST","VEGFA_FLT1","VEGFA_NRP2","CSF2_CSF2RA","CSF2_CSF2RB","TNF_TNFRSF21",
                             "CCL20_CCR6","TNF_LTBR","CXCL9_CXCR3","TNF_TNFRSF1A","CXCL10_CXCR3",
                             "THBS1_CD36","THBS1_CD47","PROS1_AXL","GAS6_AXL")
highlight_pairs$B <- list()
highlight_pairs$B$'T' <- c("TNFSF9_TNFRSF9")
highlight_pairs$'T' <- list()
highlight_pairs$'T'$B <- c("CXCL13_CXCR5","BTLA_TNFRSF14","IFNG_IFNGR1","TNFSF14_TNFRSF14")
highlight_pairs$'T'$cDC <- c("IFNG_IFNGR1",
                             "CSF2_CSF2RB",
                            "EBI3_IL27RA","LTB_LTBR","BTLA_TNFRSF14","XCL2_XCR1","XCL1_XCR1")
highlight_pairs$'T'$moDC <- c("TNFSF14_TNFRSF14","TNFRSF14_LTBR","BTLA_TNFRSF14","IFNG_IFNGR2","TNFSF4_TNFRSF4")
highlight_pairs$'T'$mac <- c("TNF_TNFRSF1B","TNFSF14_TNFRSF14","TNFRSF14_LTBR","BTLA_TNFRSF14","IFNG_IFNGR2","TNFSF4_TNFRSF4","CSF1_CSF1R",
                             "ANXA1_FPR1","ANXA1_FPR2")

highlight_pairs$moDC <- list()
highlight_pairs$moDC$'T' <- c("CCL20_CCR6","OSM_IL6ST","CCL4_CCR5","CXCL10_CXCR3","CXCL9_CXCR3","CXCL11_CXCR3","CCL17_CCR4")


plot_funk <- function(fn,labels=TRUE,highlight_pairs=NULL){
  
  col_vec <- c("black",rgb(t(col2rgb(c("green","purple","orange"))*0.7/255)))
  
  png(fn,height=25,width=20,units="in",pointsize=9,res=300)
  
  layout(matrix(1:35,nrow=7,ncol=5,byrow=T))
  par(mar=c(10,8,8,4))
  
  all_intersect <- intersect(rownames(lcam_lo$interaction_stats),rownames(lcam_hi$interaction_stats))
  
  for(subtype1 in unique(lcam_lo$interaction_stats$subtype1)){
    for(subtype2 in setdiff(unique(lcam_lo$interaction_stats$subtype2),subtype1)){
      a <- lcam_lo$interaction_stats
      b <- lcam_hi$interaction_stats
      a <- a[!is.na(a$Ligand),]
      b <- b[!is.na(b$Ligand),]
      
      
      a <- a[a$subtype1==subtype1 & a$subtype2==subtype2,]
      b <- b[b$subtype1==subtype1 & b$subtype2==subtype2,]

      #other_rows <- setdiff(c(rownames(a),rownames(b)),rows)
      # 
      # a <- rbind(a,matrix(0,nrow=length(setdiff(other_rows,rownames(a))),ncol=ncol(a),
      #                     dimnames=list(setdiff(other_rows,rownames(a)),colnames(a))))
      # 
      # 
      # b <- rbind(b,matrix(0,nrow=length(setdiff(other_rows,rownames(b))),ncol=ncol(b),
      #                     dimnames=list(setdiff(other_rows,rownames(b)),colnames(b))))
      
      rows <- intersect(rownames(a),rownames(b))
      
      sig_x <- as.integer(a[rows,]$adj.p.value<1e-2 & abs(a[rows,]$log2_ratio_score1_score2)>0.5)
      sig_y <- as.integer(b[rows,]$adj.p.value<1e-2 & abs(b[rows,]$log2_ratio_score1_score2)>0.5)*2
      sig_vec <- 1+sig_x+sig_y
      names(sig_vec) <- rows
      
      plot(lcam_lo$interaction_stats[all_intersect,]$log2_ratio_score1_score2,
           lcam_hi$interaction_stats[all_intersect,]$log2_ratio_score1_score2,pch=16,cex=0.5,col=alpha("grey",0.8),
           xlab="LR ratio LCAMlo",ylab="LR ratio LCAMhi",
           xlim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
                  max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)),
           ylim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
                  max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)))
      
#      points(a[rows,]$log2_ratio_score1_score2,
#           b[rows,]$log2_ratio_score1_score2,
#           pch=16,col=col_vec[sig_vec])
           # points(a[rows,]$log2_ratio_score1_score2,
           #      b[rows,]$log2_ratio_score1_score2,
           #      pch=16,col=1+as.integer(hi_vs_lo$interaction_stats[rows,]$adj.p.value<1e-2))
      points(a[rows,]$log2_ratio_score1_score2,
                   b[rows,]$log2_ratio_score1_score2,
                   pch=16,col=1)
      points(a[rows,]$log)
            abline(c(0,1),lty=2)
      abline(h=0,col=alpha("grey",0.4))
      abline(v=0,col=alpha("grey",0.4))
      # plot(a[rows,]$log2_ratio_score1_score2,
      #      b[rows,]$log2_ratio_score1_score2,
      #      xlab="LR ratio LCAMlo",ylab="LR ratio LCAMhi",
      #      xlim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
      #             max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)),
      #      ylim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
      #             max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)),
      #      pch=16,col=col_vec[sig_vec])
      mtext(paste(subtype1,"\u2b95",subtype2),cex=1.5)
      
      top_a <- rows[order(a[rows,]$log2_ratio_score1_score2,decreasing=T)]
      top_a <- c(head(top_a,3),tail(top_a,3))
      top_b <- rows[order(b[rows,]$log2_ratio_score1_score2,decreasing=T)]
      top_b <- c(head(top_b,3),tail(top_b,3))

      top_a <- top_a[a[top_a,]$adj.p.value<0.01]
      top_b <- top_b[b[top_b,]$adj.p.value<0.01]


      opsign <- rows[a[rows,]$adj.p.value < 0.01 & b[rows,]$adj.p.value<0.01 & a[rows,]$log2_ratio_score1_score2*b[rows,]$log2_ratio_score1_score2 < 0]


      rows <- unique(c(top_a,top_b,opsign))
      
      if(!is.null(labels)){
        par(lwd=0.6)
        
        if(!is.null(highlight_pairs[[subtype1]][[subtype2]])){
          rows <- paste(subtype1,subtype2,highlight_pairs[[subtype1]][[subtype2]],sep="_")
          rows <- intersect(rows,rownames(a))
          rows <- intersect(rows,rownames(b))
          points(a[rows,]$log2_ratio_score1_score2,
                b[rows,]$log2_ratio_score1_score2,
                pch=16,col=2)
        }
        #print(setdiff(rows,rownames(a)))
        if(labels=="spread"){
          if(length(rows)>0){
            spread.labels(a[rows,]$log2_ratio_score1_score2,
                          b[rows,]$log2_ratio_score1_score2,
                          labels=paste(a[rows,]$Ligand,a[rows,]$Receptor,sep="\u2b95"),
                          offsets=c(0.05),ony=T,
                          cex=1,xpd=NA)
            
            # points(a[rows,]$log2_ratio_score1_score2,
            #        b[rows,]$log2_ratio_score1_score2,col=col_vec[sig_vec[rows]],pch=16)
            # points(a[rows,]$log2_ratio_score1_score2,
            #        b[rows,]$log2_ratio_score1_score2,
            #        pch=16,col=1+as.integer(hi_vs_lo$interaction_stats[rows,]$adj.p.value<1e-2))
          }
        }else if(labels=="tiny"){
          text(a[rows,]$log2_ratio_score1_score2,
                        b[rows,]$log2_ratio_score1_score2,
                        labels=paste(a[rows,]$Ligand,a[rows,]$Receptor,sep="\u2b95"),
                        cex=0.7,xpd=NA)        }
              }
      par(lwd=1)
  
    }
  }
  
  # plot.new()
  # #legend("topleft",c("other celltypes","N.S.","adj-p<1e-2 in LCAMlo","adj-p<1e-2 in LCAMhi","adj-p<1e-2 in both"),pch=16,col=c("grey",col_vec),bty="n",cex=1.6)
  # legend("topleft",c("other celltypes","N.S.","FDR<1e-2 in LCAMhi vs lo tumors"),pch=16,col=c("grey",1,2),bty="n",cex=1.6)
  dev.off()
}

#plot_funk(fn=file.path(figure_path,"lig_rec_smallplots_nolabs.png"))
plot_funk(fn=file.path(figure_path,"figure_5ijklm_s5jk.png"),labels = "spread",highlight_pairs = highlight_pairs)
#plot_funk(fn=file.path(figure_path,"lig_rec_smallplots_tiny.png"),labels = "tiny",highlight_pairs = highlight_pairs)

plot_funk_one_interaction <- function(fn,labels=TRUE,highlight_pairs=NULL,subtype1,subtype2){
  
  col_vec <- c("black",rgb(t(col2rgb(c("green","purple","orange"))*0.7/255)))
  
  png(fn,height=2.3,width=3,units="in",res=300,pointsize=7.5,bg="transparent")
  par(pin=c(2.1,1.5),mgp=c(2,1,0))
  
  
  all_intersect <- intersect(rownames(lcam_lo$interaction_stats),rownames(lcam_hi$interaction_stats))
  
      a <- lcam_lo$interaction_stats
      b <- lcam_hi$interaction_stats
      a <- a[!is.na(a$Ligand),]
      b <- b[!is.na(b$Ligand),]
      
      
      a <- a[a$subtype1==subtype1 & a$subtype2==subtype2,]
      b <- b[b$subtype1==subtype1 & b$subtype2==subtype2,]
      
      #other_rows <- setdiff(c(rownames(a),rownames(b)),rows)
      # 
      # a <- rbind(a,matrix(0,nrow=length(setdiff(other_rows,rownames(a))),ncol=ncol(a),
      #                     dimnames=list(setdiff(other_rows,rownames(a)),colnames(a))))
      # 
      # 
      # b <- rbind(b,matrix(0,nrow=length(setdiff(other_rows,rownames(b))),ncol=ncol(b),
      #                     dimnames=list(setdiff(other_rows,rownames(b)),colnames(b))))
      
      rows <- intersect(rownames(a),rownames(b))
      
      sig_x <- as.integer(a[rows,]$adj.p.value<1e-2 & abs(a[rows,]$log2_ratio_score1_score2)>0.5)
      sig_y <- as.integer(b[rows,]$adj.p.value<1e-2 & abs(b[rows,]$log2_ratio_score1_score2)>0.5)*2
      sig_vec <- 1+sig_x+sig_y
      names(sig_vec) <- rows
      
      xlim <- c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
                max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5))
      ylim <- c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
                max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5))
      
      plot(lcam_lo$interaction_stats[all_intersect,]$log2_ratio_score1_score2,
           lcam_hi$interaction_stats[all_intersect,]$log2_ratio_score1_score2,pch=16,cex=0.5,col=alpha("grey",0.8),
           xlab="LR ratio LCAMlo",ylab="LR ratio LCAMhi",
           xlim=xlim,
           ylim=ylim,
           xaxt="n",yaxt="n")
      if(floor(xlim[2])-ceiling(xlim[1])>1){
        axis(side=1,ceiling(xlim[1]):floor(xlim[2]))
        axis(side=2,ceiling(ylim[1]):floor(ylim[2]))
      }else{
        x <- seq(xlim[1],xlim[2],0.5)+0.5-xlim[1]%%0.5
        axis(side=1,x[-length(x)])
        y <- seq(ylim[1],ylim[2],0.5)+0.5-ylim[1]%%0.5
        axis(side=2,y[-length(y)])
        
      }
      
      
      #      points(a[rows,]$log2_ratio_score1_score2,
      #           b[rows,]$log2_ratio_score1_score2,
      #           pch=16,col=col_vec[sig_vec])
      # points(a[rows,]$log2_ratio_score1_score2,
      #      b[rows,]$log2_ratio_score1_score2,
      #      pch=16,col=1+as.integer(hi_vs_lo$interaction_stats[rows,]$adj.p.value<1e-2))
      points(a[rows,]$log2_ratio_score1_score2,
             b[rows,]$log2_ratio_score1_score2,
             pch=16,col=1)
      points(a[rows,]$log)
      abline(c(0,1),lty=2)
      abline(h=0,col=alpha("grey",0.4))
      abline(v=0,col=alpha("grey",0.4))
      # plot(a[rows,]$log2_ratio_score1_score2,
      #      b[rows,]$log2_ratio_score1_score2,
      #      xlab="LR ratio LCAMlo",ylab="LR ratio LCAMhi",
      #      xlim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
      #             max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)),
      #      ylim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
      #             max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)),
      #      pch=16,col=col_vec[sig_vec])
      mtext(paste(subtype1,"\u2b95",subtype2),cex=1.5)
      
      top_a <- rows[order(a[rows,]$log2_ratio_score1_score2,decreasing=T)]
      top_a <- c(head(top_a,3),tail(top_a,3))
      top_b <- rows[order(b[rows,]$log2_ratio_score1_score2,decreasing=T)]
      top_b <- c(head(top_b,3),tail(top_b,3))
      
      top_a <- top_a[a[top_a,]$adj.p.value<0.01]
      top_b <- top_b[b[top_b,]$adj.p.value<0.01]
      
      
      opsign <- rows[a[rows,]$adj.p.value < 0.01 & b[rows,]$adj.p.value<0.01 & a[rows,]$log2_ratio_score1_score2*b[rows,]$log2_ratio_score1_score2 < 0]
      
      
      rows <- unique(c(top_a,top_b,opsign))
      
      if(!is.null(labels)){
        par(lwd=0.6)
        
        if(!is.null(highlight_pairs[[subtype1]][[subtype2]])){
          rows <- paste(subtype1,subtype2,highlight_pairs[[subtype1]][[subtype2]],sep="_")
          rows <- intersect(rows,rownames(a))
          rows <- intersect(rows,rownames(b))
          points(a[rows,]$log2_ratio_score1_score2,
                 b[rows,]$log2_ratio_score1_score2,
                 pch=16,col=2)
        }
        #print(setdiff(rows,rownames(a)))
        if(labels=="spread"){
          if(length(rows)>0){
            spread.labels(a[rows,]$log2_ratio_score1_score2,
                          b[rows,]$log2_ratio_score1_score2,
                          labels=paste(a[rows,]$Ligand,a[rows,]$Receptor,sep="\u2b95"),
                          offsets=c(0.05),ony=T,
                          cex=0.7,xpd=NA)
            
            # points(a[rows,]$log2_ratio_score1_score2,
            #        b[rows,]$log2_ratio_score1_score2,col=col_vec[sig_vec[rows]],pch=16)
            # points(a[rows,]$log2_ratio_score1_score2,
            #        b[rows,]$log2_ratio_score1_score2,
            #        pch=16,col=1+as.integer(hi_vs_lo$interaction_stats[rows,]$adj.p.value<1e-2))
          }
        }else if(labels=="tiny"){
          text(a[rows,]$log2_ratio_score1_score2,
               b[rows,]$log2_ratio_score1_score2,
               labels=paste(a[rows,]$Ligand,a[rows,]$Receptor,sep="\u2b95"),
               cex=0.7,xpd=NA)        }
      }
      par(lwd=1)
      
    
  
  
  # plot.new()
  # #legend("topleft",c("other celltypes","N.S.","adj-p<1e-2 in LCAMlo","adj-p<1e-2 in LCAMhi","adj-p<1e-2 in both"),pch=16,col=c("grey",col_vec),bty="n",cex=1.6)
  # legend("topleft",c("other celltypes","N.S.","FDR<1e-2 in LCAMhi vs lo tumors"),pch=16,col=c("grey",1,2),bty="n",cex=1.6)
  dev.off()
}

# for(subtype1 in names(highlight_pairs)){
#   for(subtype2 in names(highlight_pairs[[subtype1]])){
# plot_funk_one_interaction(fn=file.path(figure_path,paste(subtype1,"_",subtype2,".png",sep="")),labels = "tiny",highlight_pairs = highlight_pairs,
#                           subtype1=subtype1,subtype2=subtype2)
# plot_funk_one_interaction(fn=file.path(figure_path,paste(subtype1,"_",subtype2,"spread.png",sep="")),labels = "spread",highlight_pairs = highlight_pairs,
#                               subtype1=subtype1,subtype2=subtype2)
# plot_funk_one_interaction(fn=file.path(figure_path,paste(subtype1,"_",subtype2,"nolabs.png",sep="")),labels = "NULL",highlight_pairs = highlight_pairs,
#                           subtype1=subtype1,subtype2=subtype2)
#   }
# }




col_vec <- c("black",rgb(t(col2rgb(c("purple","green","orange"))*0.7/255)))


a <- lcam_lo$interaction_stats
b <- lcam_hi$interaction_stats
a <- a[!is.na(a$Ligand),]
b <- b[!is.na(b$Ligand),]

rows <- intersect(rownames(a),rownames(b))

sig_x <- as.integer(a[rows,]$adj.p.value<1e-2 & abs(a[rows,]$log2_ratio_score1_score2)>0.5)
sig_y <- as.integer(b[rows,]$adj.p.value<1e-2 & abs(b[rows,]$log2_ratio_score1_score2)>0.5)*2
sig_vec <- 1+sig_x+sig_y
names(sig_vec) <- rows

png("output/figures/figure_5h.png",height=2.3,width=3,units="in",res=300,pointsize=7.5)
par(pin=c(2.1,1.5),mgp=c(2,1,0))
plot(a[rows,]$log2_ratio_score1_score2,
       b[rows,]$log2_ratio_score1_score2,
       pch=16,col=alpha(col_vec[sig_vec],0.6),xlab="LR ratio LCAMlo",ylab="LR ratio LCAMhi",cex=0.7,
     xlim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
            max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)),
     ylim=c(min(c(min(a[rows,]$log2_ratio_score1_score2)-0.2,min(b[rows,]$log2_ratio_score1_score2)-0.2,-0.5)),
            max(c(max(a[rows,]$log2_ratio_score1_score2)+0.2),max(b[rows,]$log2_ratio_score1_score2)+0.2,0.5)))
abline(c(0,1),lty=2)
abline(h=0,col=alpha("grey",0.6))
abline(v=0,col=alpha("grey",0.6))
mtext("All interactions",cex=1.5)
legend("topleft",c("N.S.","adj-p<1e-2 in LCAMlo","adj-p<1e-2 in LCAMhi","adj-p<1e-2 in both"),pch=16,col=col_vec,bty="n",cex=1)
dev.off()

}