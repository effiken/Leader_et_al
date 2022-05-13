
#' tcr_output_to_analysis
#' Take cellranger VDJ output for T cells & convert to a nicer form, allowing analysis across multiple samples
#' Note: Since currently, cellranger VDJ is run on the output of each amp_batch irrespective of CITEseq hashing,
#' the output file is kept with the raw data independent of the dehashed sample_IDs.
#' 
#' KEY ASSUMPTIONS:
#' 1. All relevent information ends up in the filtered_contig_annotations.csv output file from cellranger
#' 2. All relevent TCRs are marked as PRODUCTIVE==TRUE
#' 3. TCR chains are uniquely define by: 1. precise CDR3 AA sequence, precise V, D, J genes, constant chain type (alpha vs. beta). Also allowing rare cases where the constant chain is characterized as "MULTI".
#' 4. Multiplets can be defined as any cells with >3 observed chains of any chain type (alpha or beta)
#' 5. If the chains comprising cells with 3 observed chains exist as a subset of larger multiplets that are observed, they are also considered multiplets
#' 6. If the chains comprising cells with 2 observed chains exist as a subset of remaining cells with 3 observed chains, they are handled as follows:
#' => if it's a unique clone with 3 chains, group together as clonotype with the 2-chain cells
#' => else if there are multiple 3 chained TCR that contain these same 2 chains, don't group.
#' 7. For all clonotypes with single observed chains: find matches as follows:
#         -if !=1 match in the other remaining clonotypes, leave alone
#         -if 1 unique match, group with the other clonotype.
#
#' @param sample_IDs
#' @param sample_annots_path
#' @param demux NOT YET SUPPORTED. For the case where samples are hashed in one 10X lane using CITEseq, and cell to sampleID mapping is ambiguous based on amp_batch alone. FALSE by default.
#' @param data_path The project cellranger output folder. Ex: scRNAseq_analysis/sc_data_main/liver_data/
#' 
#' @return TCR_struct A datastructure with slots cell_ID, TCRid, unique_TCRs, unique_chains, orig_output
#' 
#' @export


tcr_output_to_analysis <- function(sample_IDs,sample_annots_path,demux=F,data_path){

  warning('Note that this function is deprecated and does not handle hashed data - use tcr_preprocessing.R instead')
  
  if(demux!=FALSE){
    stop("LEAVE DEMUX AS DEFAULT! CITEseq dehashing not yet supported.")
  }
  sample_IDs <- as.character(sample_IDs)
  
sample_annots <- read.csv(sample_annots_path,r=1,h=1,stringsAsFactors = F)
amp_batches <- unique(sample_annots[sample_IDs,]$amp_batch_ID)
contig_annot_fns <- paste(data_path,amp_batches,"/vdj_t/filtered_contig_annotations.csv",sep="")

filtered_annots <- lapply(contig_annot_fns,FUN=function(x){read.csv(x,h=1,stringsAsFactors = F)})
filtered_annots <- lapply(filtered_annots,FUN=function(x){x <- x[x$productive=="True",]})

# 1. rename cells & clonotypes with sample ID

if(demux==F){
  for(iter in 1:length(sample_IDs)){
    filtered_annots[[iter]]$barcode <- paste(sample_IDs[iter],filtered_annots[[iter]]$barcode,sep="_")
    filtered_annots[[iter]]$raw_consensus_id <- paste(sample_IDs[iter],filtered_annots[[iter]]$raw_consensus_id,sep="_")
    filtered_annots[[iter]]$raw_clonotype_id <- paste(sample_IDs[iter],filtered_annots[[iter]]$raw_clonotype_id,sep="_")
  }
}

# 2. concatenate across samples
filtered_annots <- do.call(rbind,filtered_annots)

# find unique TCR chains that are observed
unique_chains <- unique(filtered_annots[,c('cdr3','v_gene','d_gene','j_gene','chain')])
rownames(unique_chains) <- seq(nrow(unique_chains))

#map the "consensus sequences" to the unique set of observed TCR chains
consensus_to_chain <- match(apply(filtered_annots[,c('cdr3','v_gene','d_gene','j_gene','chain')],1,paste,collapse="_"),
             apply(unique_chains,1,paste,collapse="_"))
names(consensus_to_chain) <- filtered_annots$raw_consensus_id

#Determine the unique TCR chains existing within each clonotype annotation
tcr_chain_collection <- lapply(split(filtered_annots$raw_consensus_id,filtered_annots$raw_clonotype_id),
              function(x){unique(consensus_to_chain[x])})

max_chains <- max(unlist(lapply(tcr_chain_collection,length)))

tcr_chain_mat <- matrix(NA,nrow=length(tcr_chain_collection),ncol=max_chains)
rownames(tcr_chain_mat) <- names(tcr_chain_collection)

for(iter in 1:length(tcr_chain_collection)){
  tcr_chain_mat[iter,1:length(tcr_chain_collection[[iter]])] <- sort(tcr_chain_collection[[iter]])
}

message("dereplicate clonotypes")

#     - generate matrix of clonotypes (rows) & unique chain ID numbers (row values, sorted)
unique_TCRs <- unique(tcr_chain_mat)
message(nrow(tcr_chain_mat)-nrow(unique_TCRs)," of ",nrow(tcr_chain_mat)," are exact replicates")
rownames(unique_TCRs) <- seq(nrow(unique_TCRs))

#     - group matching clonotypes
clonotype_to_unique_TCR <- match(apply(tcr_chain_mat,1,paste,collapse="_"),apply(unique_TCRs,1,paste,collapse="_"))
names(clonotype_to_unique_TCR) <- rownames(tcr_chain_mat)


message("separate multiplets")
#     - separate all clonotypes with > 3 observed chains (doublets)
multiplet_flag <- !is.na(unique_TCRs[,4])
multiplets <- rownames(unique_TCRs)[multiplet_flag]
message(length(multiplets)," TCRs have >3 unique chains and are therefore excluded multiplets")

#     - For all cells with 3 chains: exist within multiplets matrix?
three_chains_flag <- !is.na(unique_TCRs)[,3] & is.na(unique_TCRs[,4])
three_chains <- rownames(unique_TCRs)[three_chains_flag]

for(row_fg in seq(length(three_chains))){
  for(row_bg in seq(length(multiplets))){
    if(sum(unique_TCRs[three_chains[row_fg],1:3]%in%unique_TCRs[multiplets[row_bg],])==3){
      multiplets <- c(multiplets,three_chains[row_fg])
      three_chains[row_fg] <- NA
      break
    }
  } 
}
message(sum(is.na(three_chains))," cells with 3 unique chains were similar to multiplets and were also excluded")
three_chains <- three_chains[!is.na(three_chains)]


#     - For for all clonotypes with 2 chains:
#         -search for each within cells with 3 chains.

message("Now searching for similarities between double chain & triple chain TCRs")

two_chains_flag <- !is.na(unique_TCRs)[,2] & is.na(unique_TCRs[,3])
two_chains <- rownames(unique_TCRs)[two_chains_flag]
message(length(two_chains)," TCRs found with 2 unique chains")

in_3chainmat <- apply(unique_TCRs[two_chains,1:2],2,
                      function(x){x%in%unique_TCRs[three_chains,1:3]})
two_chains_test <- two_chains[rowSums(in_3chainmat)==2]

message(length(two_chains_test)," 2-chain TCRs exist within 3-chain collection")

group_count <- 0
multimatch_count <- 0
for(iter in two_chains_test){
  mat3 <- unique_TCRs[three_chains,1:3]
  iter_tcr <- unique_TCRs[iter,1:2]
  match_flag <- rowSums(apply(mat3,2,function(x){x%in%iter_tcr}))==2
  if(sum(match_flag)==1){
    clonotype_to_unique_TCR[clonotype_to_unique_TCR==iter] <- three_chains[match_flag]
    group_count <- group_count+1
    two_chains <- setdiff(two_chains,iter)
    #three_chains <- c(three_chains,iter)
  } else if(sum(match_flag)>1){
    multimatch_count <- multimatch_count + sum(match_flag)
    multiplets <- c(multiplets,three_chains[match_flag])
    three_chains <- setdiff(three_chains,three_chains[match_flag])
  }
}
message(group_count, " 2-chain TCRs match uniquely with existing 3-chain TCRs and are being grouped")
message("2-chain TCRs matched ambiguously to ",multimatch_count," 3-chain TCRs; these 3-chain TCRs were discarded as multiplets")

message("Now seeking single-chain TCRs that match uniquely to only one 2-chain TCR")
one_chains_flag <- !is.na(unique_TCRs)[,1] & is.na(unique_TCRs[,2])
one_chains <- rownames(unique_TCRs)[one_chains_flag]
message(length(one_chains)," TCRs found with only 1 chain")

in_2chainmat <- unique_TCRs[one_chains,1]%in%unique_TCRs[c(two_chains,three_chains),1:3]
one_chains_test <- one_chains[in_2chainmat]
message(sum(!in_2chainmat)," 1-chain TCRs excluded because they don't match any 2- or 3-chain TCRs")
clonotype_to_unique_TCR[clonotype_to_unique_TCR%in%one_chains[!in_2chainmat]] <- "single_chain"

group_count <- 0
multimatch_count <- 0
for(iter in one_chains_test){
  mat3 <- unique_TCRs[c(two_chains,three_chains),1:3]
  mat3[is.na(mat3)] <- 0
  iter_tcr <- unique_TCRs[iter,1]
  match_flag <- rowSums(apply(mat3,2,function(x){x==iter_tcr}))==1
  if(sum(match_flag)==1){
    new_tcr <- names(match_flag)[which(match_flag)]
    clonotype_to_unique_TCR[clonotype_to_unique_TCR==iter] <- new_tcr
    group_count <- group_count+1
    one_chains <- setdiff(one_chains,iter)
    # if(new_tcr%in%two_chains){
    #   two_chains <- c(two_chains,iter)
    # }else if(new_tcr%in%three_chains){
    #   three_chains <- c(three_chains,iter)
    # }
  }else if(sum(match_flag)>1){
    multimatch_count <- multimatch_count + 1
    clonotype_to_unique_TCR[clonotype_to_unique_TCR==iter] <- "single_chain"
  }
}
message(group_count, " 1-chain TCRs match uniquely with existing 2- or 3-chain TCRs and are being grouped")
message(multimatch_count, " 1-chain TCRs matched ambiguously to 2- or 3-chain TCRs; these 1-chain TCRs were discarded as ambiguous single-chain events")

clonotype_to_unique_TCR[clonotype_to_unique_TCR%in%multiplets] <- "multiplet"



cell_to_clonotype <- unlist(lapply(split(filtered_annots$raw_clonotype_id,filtered_annots$barcode),unique))

cell_to_unique_TCR <- clonotype_to_unique_TCR[cell_to_clonotype]
names(cell_to_unique_TCR) <- names(cell_to_clonotype)


# output:
# cell ID

TCR_struct <- list()
TCR_struct$cell_ID <- unique(filtered_annots$barcode)
# TCRid
TCR_struct$TCRid <- cell_to_unique_TCR
rownames(unique_TCRs) <- seq(nrow(unique_TCRs))
TCR_struct$unique_TCRs <- unique_TCRs
TCR_struct$unique_chains <- unique_chains
TCR_struct$orig_output <- filtered_annots

return(TCR_struct)
}
