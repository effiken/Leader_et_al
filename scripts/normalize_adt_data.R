
library(scTools)
library(Matrix)




adt_list_to_matrix <- function(adt_list){
  bad_rows <- c("bad_struct_adt", "no_match_adt", "total_reads_adt", 
                "total_reads_adt_pf", "der_bad_struct_adt", "der_no_match_adt", 
                "der_total_reads_adt", "der_total_reads_adt_pf", "der_gex_umi_sum", 
                "der_unique_gene_count", "der_HTO_max1_v_HTO_min", "der_HTO_max2_v_HTO_min")
  for (iter in 1:length(adt_list)) {
    adt_list[[iter]] <- adt_list[[iter]][!rownames(adt_list[[iter]]) %in% 
                                           bad_rows, ]
    rownames(adt_list[[iter]]) <- sub(x = rownames(adt_list[[iter]]), 
                                      pattern = ".IH", replacement = "")
    rownames(adt_list[[iter]]) <- sub(x = rownames(adt_list[[iter]]), 
                                      pattern = "adt_", replacement = "")
    rownames(adt_list[[iter]]) <- gsub(x = rownames(adt_list[[iter]]), 
                                       pattern = "[[:punct:]]", replacement = "")
    if (substr(rownames(adt_list[[iter]])[1], 1, 1) == 5) {
      rownames(adt_list[[iter]]) <- substr(rownames(adt_list[[iter]]), 
                                           2, nchar(rownames(adt_list[[iter]])))
    }
  }
  markers <- unique(unlist(lapply(adt_list, rownames)))
  marker_mat <- matrix(NA, nrow = length(markers), ncol = sum(unlist(lapply(adt_list, 
                                                                            ncol))), dimnames = list(markers, unlist(lapply(adt_list, 
                                                                                                                            colnames))))
  for (iter in 1:length(adt_list)) {
    iter_markers <- intersect(markers, rownames(adt_list[[iter]]))
    marker_mat[iter_markers, colnames(adt_list[[iter]])] <- as.matrix(adt_list[[iter]][iter_markers, 
                                                                                       ])
  }
  return(marker_mat)
}

quantile_normalize_variable_vector_lengths <- function (vectorList, logScale = TRUE, logRegularization = 1) 
{
  res <- 1/(max(unlist(lapply(vectorList, length))) - 1)
  q_list <- lapply(vectorList, quantile, seq(0, 1, res), na.rm = T)
  q_mat <- matrix(unlist(q_list), nrow = length(q_list[[1]]))
  q_mat[, is.na(unlist(lapply(vectorList, function(x) {
    x[[1]][1]
  })))]
  if (logScale == TRUE) {
    ref_dist <- rowMeans(log10(logRegularization + q_mat), 
                         na.rm = T)
  }
  else {
    ref_dist <- rowMeans(q_mat, na.rm = T)
  }
  qn_dist <- list()
  for (dist_iter in 1:length(vectorList)) {
    l <- length(vectorList[[dist_iter]])
    inds <- round((1:l) * length(ref_dist)/l)
    new_dist <- ref_dist[inds]
    if (is.na(vectorList[[dist_iter]][1])) {
      new_dist <- array(NA, length(new_dist))
    }
    qn_dist[[names(vectorList)[dist_iter]]] <- new_dist[match(vectorList[[dist_iter]], 
                                                              vectorList[[dist_iter]][order(vectorList[[dist_iter]])])]
    if (logScale == TRUE) {
      qn_dist[[names(vectorList)[dist_iter]]] <- 10^qn_dist[[names(vectorList)[dist_iter]]] - 
        logRegularization
    }
  }
  return(qn_dist)
}



qn_adttab <- function (adttab, cell_to_batch) 
{
  batch_cells <- split(names(cell_to_batch), cell_to_batch)
  message("Normalizing antibody signal across batches")
  normalized_adttab <- matrix(NA, nrow = nrow(adttab), ncol = ncol(adttab), 
                              dimnames = list(rownames(adttab), colnames(adttab)))
  for (ab_iter in rownames(normalized_adttab)) {
    message(ab_iter)
    dists <- lapply(batch_cells, function(x) {
      adttab[ab_iter, x]
    })
    qn_dists <- quantile_normalize_variable_vector_lengths(dists)
    for (batch_iter in names(batch_cells)) {
      normalized_adttab[ab_iter, batch_cells[[batch_iter]]] <- qn_dists[[batch_iter]]
    }
  }
  return(normalized_adttab)
}

#message("loading expression data into R")
#load("data/lung_ldm.rd")
sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
adt_list <- lung_ldm$dataset$adt_by_sample
adt_mat <- adt_list_to_matrix(adt_list)
adt_mat <- adt_mat[,intersect(colnames(adt_mat),colnames(lung_ldm$dataset$umitab))]
cell_to_batch <- sample_annots[lung_ldm$dataset$cell_to_sample[colnames(adt_mat)],]$amp_batch_ID
  
names(cell_to_batch) <- colnames(adt_mat)
qn_adtmat <- qn_adttab(adt_mat,cell_to_batch)

if(!dir.exists("intermediates")){
  dir.create("intermediates")
}
save(qn_adtmat,file="intermediates/qn_adtmat.rd")
rm(list=c("sample_annots","adt_list","adt_mat","cell_to_batch"))
