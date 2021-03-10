## this script sets the parameters for all other scripts; EDIT THIS


## directory ####
# root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
setwd(root)


## packages ####
source("src/helpers.R")
if (!"flowLearn"%in%rownames(installed.packages()))
  devtools::install_github("mlux86/flowLearn")
libr(c(
  "plyr", "purrr", "furrr", "doMC", #"rslurm", # parallelization
  "data.table", "stringr", # generic functions
  "flowCore", "flowDensity", "flowLearn", # for 00_data + flowCore for plot_dens function
  "RColorBrewer", "ggplot2",# for plotting scaterplots
  "Rfast", # for calculating distance in good_samples
  "kmed",  # for calculating kmed in good_samples
  "combinat" # for calculating possible combos in scoring gigasom
), FALSE)


## parallelization ####
no_cores <- ifelse(exists("no_cores"), no_cores, 1)
doMC::registerDoMC(no_cores) # for plyr::llply
# future::plan(future::multisession) # for furrr::future_map

## folder_names ####
raw_dir <- paste0(root, "/raw")
data_dir <- paste0(root, "/data")
results_dir <- paste0(root, "/results")
scores_dir <- paste0(root, "/scores")
# plots_dir <- paste0(root, "/plots")

xn_dir <- paste0(raw_dir,"/nD/x")
yn_dir <- paste0(raw_dir,"/nD/y")
plotn_dir <- paste0(raw_dir,"/nD/scatterplots")
x2_dir <- paste0(raw_dir,"/2D/x")
y2_dir <- paste0(raw_dir,"/2D/y")
thres_dir <- paste0(raw_dir,"/2D/thresholds")
filt_dir <- paste0(raw_dir,"/2D/filters")

# replace folder name
gs_xr <- function(flnm,f3,f1="data",single=TRUE) {
  xd <- stringr::str_extract(flnm,"[2n]D")
  fld <- stringr::str_extract(flnm,"[a-zA-Z_]+/[2n]D/[a-zA-Z_+-]+")
  if (single)
    return(gsub(fld[1],paste0(f1,"/",xd[1],"/",f3),flnm))
  return(sapply(seq_len(length(fld)), function(i)
    gsub(fld[i],paste0(f1,"/",xd[i],"/",f3),flnm[i])))
} 

