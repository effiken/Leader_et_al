
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
figure_dir <- file.path(wd,"figures_out")
if(!dir.exists(figure_dir)){
  dir.create(figure_dir)
}

# to download & load mouse DC data
data_url <- "https://www.dropbox.com/s/vjbide8ro5iwrfh/lung_ldm.rd?dl=1"
if(Sys.info()["sysname"]=="Windows"){
  download.file(url = data_url,destfile = file.path(data_dir,"lung_ldm.rd"),mode="wb")
}else{
  download.file(url = data_url,destfile = file.path(data_dir,"lung_ldm.rd"))
}
download_check <- try(load(file.path(data_dir,"lung_ldm.rd")))
if(download_check!="ldm"){stop("Error using download.file(), due to operating system?-- please manually download human DC dataset from https://www.dropbox.com/s/3ffq3z75af37rgr/mouse_dc.rd?dl=1 and place in working directory and remove this chunk of code")}

