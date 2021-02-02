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
  "RColorBrewer", # for plotting scaterplots
  "Rfast", # for calculating distance in good_samples
  "kmed",  # for calculating kmed in good_samples
  "combinat" # for calculating possible combos in scoring gigasom
))


## parallelization ####
no_cores <- ifelse(exists("no_cores"), no_cores, 1)
doMC::registerDoMC(no_cores) # for plyr::llply
future::plan(future::multisession) # for furrr::future_map

## folder_names ####
data_dir <- paste0(root, "/data")
results_dir <- paste0(root, "/results")
scores_dir <- paste0(root, "/scores")

xn_dir <- paste0(data_dir,"/nD/x")
yn_dir <- paste0(data_dir,"/nD/y")
plotn_dir <- paste0(data_dir,"/nD/scatterplots")
x2_dir <- paste0(data_dir,"/2D/x")
y2_dir <- paste0(data_dir,"/2D/y")
thres_dir <- paste0(data_dir,"/2D/thresholds")
filt_dir <- paste0(data_dir,"/2D/filters")
