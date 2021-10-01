
library(matrixStats)
library(Matrix)
library(seriation)
library(scales)

figure_4ab <- function(lung_ldm,qn_adtmat){
  
  sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
  annots_list <- read.csv("input_tables/annots_list.csv",r=1,h=1,stringsAsFactors = F)
  cluster_order_list <- read.table("input_tables/cluster_order_list.txt",sep="\t",r=1,h=1,stringsAsFactors = F)
  genes_in_figures <- read.csv("input_tables/table_s5_genes_used_in_figures.csv",r=1,h=1,stringsAsFactors = F)
  markers_in_figures <- read.csv("input_tables/markers_in_figures.csv",r=1,h=1,stringsAsFactors = F)
  
  
  colgrad_rna <- colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100)
  
  genes <- strsplit(genes_in_figures["Figure 4A",2],",")[[1]]
  markers <- strsplit(markers_in_figures["NK_T",],",")[[1]]
  cluster_order <- as.character(strsplit(cluster_order_list["NK_T",],",")[[1]])
  lineage_col_array <- c(rgb(t(col2rgb(c(2,3,4,5,6,7)))*0.7/255),"grey","black","red","green")
  names(lineage_col_array) <- unique(annots_list[cluster_order,]$sub_lineage)
  
  qn_adtmat <- qn_adtmat[markers,!is.na(colSums(qn_adtmat[markers,]))]
  
  qn_adtmat <- qn_adtmat[,lung_ldm$dataset$cell_to_cluster[colnames(qn_adtmat)] %in% cluster_order]
  
  normalized_qn_adtmat <- log2(qn_adtmat/rowMedians(qn_adtmat))
  
  normalized_qn_adtmat <- normalized_qn_adtmat[,intersect(colnames(qn_adtmat),colnames(lung_ldm$dataset$ds[[1]]))]
  rm(qn_adtmat)
  
  n_cells_per_cluster <- min(table(lung_ldm$dataset$cell_to_cluster[colnames(normalized_qn_adtmat)])[cluster_order])
  
  cells_by_cluster <- split(colnames(normalized_qn_adtmat),lung_ldm$dataset$cell_to_cluster[colnames(normalized_qn_adtmat)])[cluster_order]
  set.seed(910430)
  cells_to_plot <- unlist(lapply(cells_by_cluster,sample,n_cells_per_cluster))
  
  rna_mat <- as.matrix(lung_ldm$dataset$ds[[1]][rev(genes),cells_to_plot])
  rna_mat <- log2(1+rna_mat)
  rna_log2_thresh <- 2
  rna_mat[rna_mat > rna_log2_thresh] <- rna_log2_thresh
  
  prot_mat <- normalized_qn_adtmat[rev(markers),cells_to_plot]
  
  prot_thresh <- 5
  prot_mat[prot_mat > prot_thresh] <- prot_thresh
  prot_mat[prot_mat < -prot_thresh] <- -prot_thresh

  cell_to_tissue <- sample_annots[lung_ldm$dataset$cell_to_sample,]$tissue
  cell_to_patient <- sample_annots[lung_ldm$dataset$cell_to_sample,]$patient_ID
  names(cell_to_tissue)=names(cell_to_patient) <- names(lung_ldm$dataset$cell_to_cluster)
  
  cell_tab <- table(
                    cell_to_patient[lung_ldm$dataset$cell_to_cluster%in%cluster_order],
                    lung_ldm$dataset$cell_to_cluster[lung_ldm$dataset$cell_to_cluster%in%cluster_order],
                    cell_to_tissue[lung_ldm$dataset$cell_to_cluster%in%cluster_order]
                    )
  cell_tab_normal <- cell_tab[,,"Normal"]/rowSums(cell_tab[,,"Normal"])
  cell_tab_normal <- cell_tab_normal[,cluster_order]
  cell_tab_tumor <- cell_tab[,,"Tumor"]/rowSums(cell_tab[,,"Tumor"])
  cell_tab_tumor <- cell_tab_tumor[,cluster_order]
  
  at_normal <- seq(1/3,length(cluster_order)-2/3,1)
  at_tumor <- at_normal + 1/3

  lin_labs <- annots_list[cluster_order,]$sub_lineage
  
  clust_breaks <- (1:(length(cluster_order)-1))/length(cluster_order)
  mid_clusts <- 1/2*(c(0,clust_breaks)+c(clust_breaks,1))
  
  message("Saving figure...")
  if(!dir.exists("output/figures")){
    dir.create("output/figures")
  }
  
  output_file <- "output/figures/figure_4ab.png"
  height=4.68
  width=4.72
  res=1000
  pointsize=6
  
  png(output_file, height = height, width = width, units = "in", 
      res = res, pointsize = pointsize,bg="transparent")
  
  layout(matrix(1:6, ncol = 2, byrow = T), heights = c(3, 
                                                        length(genes), 1.5 * length(markers)), 
         widths = c(length(cluster_order), 5))
  par(oma = c(1, 10, 5, 3), mar = c(0, 0, 1, 2))
  
  image(as.matrix(as.numeric(factor(lin_labs,levels=unique(lin_labs)))),col=lineage_col_array,
        xaxt="n",yaxt="n")
  box()
  abline(v=clust_breaks*(1+length(cluster_order))/length(cluster_order)-1/2/length(cluster_order),lwd=0.1)
  mtext(side = 1, at = seq(0,1,1/(length(cluster_order)-1)), cluster_order, las = 1, 
        cex = 0.8, line = 0.1, xpd = NA)
  plot.new()
  
  image(t(rna_mat),col=colgrad_rna,xaxt="n",yaxt="n")
  abline(v = clust_breaks, col = "grey", lwd = 0.5)
  mtext(rownames(rna_mat), side = 2, at = seq(0, 1, 1/(nrow(rna_mat) - 
                                                         1)), las = 2, line = 0.5, cex = 0.7, las = 2,font=3)
  box()
  par(pin = c(0.15, 0.6))
  image(t(as.matrix(1:100)), xaxt = "n", yaxt = "n", col = colgrad_rna)
  mtext("Log2(#UMI\n/2000)", line = 0.25, cex = 1)
  mtext(c(0, paste("\u2265", rna_log2_thresh, sep = "")), at = c(0, 
                                                            1), line = 0.25, cex = 1, side = 4, las = 2)
  box()
  par(mar = c(0, 0, 1, 2))
  image(t(as.matrix(prot_mat)), xaxt = "n", yaxt = "n", col = bluered(50))
  mtext(rev(markers), side = 2, at = seq(0, 1, 1/(length(markers) - 
                                               1)), las = 2, line = 0.5, cex = 1)
  box()
  abline(v = clust_breaks, col = "grey", lwd = 0.5)
  par(pin = c(0.15, 0.6))
  image(t(as.matrix(1:100)), xaxt = "n", yaxt = "n", col = bluered(50))
  mtext("Log2(expression\n/median)", line = 0.25, cex = 1, 
        side = 3)
  mtext(round(range(prot_mat), 1), at = c(0, 1), line = 0.25, 
        cex = 1, side = 4, las = 2)
  box()
  dev.off()
}


