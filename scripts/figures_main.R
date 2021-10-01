
# LEADER_ET_AL/SCRIPTS/MAIN.R
# Andrew M. Leader
# 9/28/2021

# This script reproduces figures from Leader, et al. Cancer Cell, 2021


# Change wd to the path to the downloaded github repo

rm(list=ls())

############################
#Specify a working directory. This is the path where the data will download to, and where figures will be produced.
# It should also be the path to the scripts from the github.
wd <- "/users/andrew leader/Downloads/Leader_et_al"
############################

# set up the directory
if(!dir.exists(wd)){
  dir.create(wd)
}
setwd(wd)
data_dir <- file.path(wd,"data")
if(!dir.exists(data_dir)){
  dir.create(data_dir)
}

if(!file.exists("data/lung_ldm.rd")){
message("Downloading Mt. Sinai scRNAseq and CITEseq data")
data_url <- "https://www.dropbox.com/s/vjbide8ro5iwrfh/lung_ldm.rd?dl=1"
if(Sys.info()["sysname"]=="Windows"){
  download.file(url = data_url,destfile = file.path(data_dir,"lung_ldm.rd"),mode="wb")
}else{
  download.file(url = data_url,destfile = file.path(data_dir,"lung_ldm.rd"))
}
message("Loading Mt. Sinai scRNAseq and CITEseq data")
download_check <- try(load(file.path(data_dir,"lung_ldm.rd")))
if(download_check!="lung_ldm"){stop("Error using download.file(), due to operating system?-- please manually download dataset from https://www.dropbox.com/s/vjbide8ro5iwrfh/lung_ldm.rd?dl=1 and place in /data/ and remove this chunk of code")}
}else{
  message("Loading Mt. Sinai scRNAseq and CITEseq data into R")
  load("data/lung_ldm.rd")
}


# Figures requiring data from Lambrechts, et al.

if(!file.exists("data/lambrechts_ldm.rd")){
  message("Downloading Lambrechts, et al. scRNAseq data")
  data_url <- "https://www.dropbox.com/s/7ksics3xe4igsu8/lambrechts_ldm_200519.rd?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = file.path(data_dir,"lambrechts_ldm.rd"),mode="wb")
  }else{
    download.file(url = data_url,destfile = file.path(data_dir,"lambrechts_ldm.rd"))
  }
}


# Figures requiring data from Zilionis, et al.

if(!file.exists("data/zilionis_ldm.rd")){
  message("Downloading Zilionis, et al. scRNAseq data")
  data_url <- "https://www.dropbox.com/s/puv3ds5b9ty0rn4/zilionis_ldm.rd?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = file.path(data_dir,"zilionis_ldm.rd"),mode="wb")
  }else{
    download.file(url = data_url,destfile = file.path(data_dir,"zilionis_ldm.rd"))
  }
}

# 5' scRNAseq data
if(!file.exists("data/lung_ldm_fiveprime.rd")){
  message("Downloading 5' scRNAseq data")
  data_url <- "https://www.dropbox.com/s/xw9atefiv6jam9u/lung_ldm_fiveprime.rd?dl=1"
  if(Sys.info()["sysname"]=="Windows"){
    download.file(url = data_url,destfile = file.path(data_dir,"lung_ldm_fiveprime.rd"),mode="wb")
  }else{
    download.file(url = data_url,destfile = file.path(data_dir,"lung_ldm_fiveprime.rd"))
  }
}

if(!file.exists("intermediates/qn_adtmat.rd")){
source("scripts/normalize_adt_data.r")
}else{
  message("Loading normalized CITEseq marker expression")
  load("intermediates/qn_adtmat.r")
}

#Joint scRNA-CITEseq heatmaps
source("scripts/figure_1cde.R")
figure_1cde(lung_ldm,qn_adtmat)
source("scripts/figure_2ab.R")
figure_2ab(lung_ldm,qn_adtmat)
source("scripts/figure_4ab.r")
figure_4ab(lung_ldm,qn_adtmat)
source("scripts/figure_4fg.R")
figure_4fg(lung_ldm,qn_adtmat)

#Other scRNAseq and CITEseq expression heatmaps
source("scripts/figure_3a.R")
figure_3a(lung_ldm)
source("scripts/figure_3b.R")
figure_3b(lung_ldm,qn_adtmat)
source("scripts/figure_3h.R")
figure_3h()
source("scripts/figure_s3a.R")
figure_s3a()
source("scripts/figure_s3b.R")
figure_s3b()

#frequency boxplot figures
source("scripts/figure_2c_s2c.R")
figure_2c_s2c(lung_ldm)
source("scripts/figure_3g.R")
figure_3g(lung_ldm)
source("scripts/figure_4c_s4l.R")
figure_4c_s4l(lung_ldm)

#Myeloid sub-population differential expression and gene-score analyses
scClustering_dir <- "scripts/scClustering/"
source("scripts/figure_2g_s2de.R")
figure_2g_s2de()
source("scripts/figure_3def_s3i.R")
figure_3def_s3i()

#Gene Module analyses
source("scripts/figure_s2fghi.R")
figure_s2fghi()
source("scripts/figure_3c_s3efgh.R")
figure_3c_s3efgh()

# Miscellaneous plots
source("scripts/figure_1f.r")
figure_1f()
source("scripts/figure_1g.R")
figure_1g()
source("scripts/figure_s1de.R")
figure_s1de()
source("scripts/figure_s1g.R")
figure_s1g()
source("scripts/figure_s1hi.R")
figure_s1hi()
source("scripts/figure_2dh_s2j.R")
figure_2dh_s2j()
source("scripts/figure_s2ab.R")
figure_s2ab()
source("scripts/figure_s4abce.R")
figure_s4abce()
source("scripts/figure_s4d.R")
figure_s4d()
source("scripts/figure_4de_s4ijk.R")
figure_4de_s4ijk()
source("scripts/figure_s4fgh.R")
figure_s4fgh()

# Ligand-receptor plots
source("scripts/figure_5hijklm_s5jk.R")
figure_5hijklm_s5jk()

# Celltype frequency correlations and development of the LCAM cellular modules
source("scripts/figure_5abcd_s5a.R")
figure_5abcd_s5a()

#TCGA and CPTAC bulk analyses
source("scripts/figure_6ac_s6abcdefjklmnopqr_7abcd_s7ab.R")
figure_6ac_s6abcdefjklmnopqr_7abcd_s7ab()
source("scripts/figure_6b_s6ghi.R")
figure_6b_s6ghi()