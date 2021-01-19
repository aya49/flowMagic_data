# date created: 2020-01-04
# author: alice yue
# input: score matrices
# output: plots


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 10#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
source("helpers_2D.R")
libr(c(
  "furrr" #"rslurm",
))


## input ####
main_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
sc2_dir <- paste0(main_dir,"/scores/2D/deepCyTOF_labels")
scn_dir <- paste0(main_dir,"/scores/nD/deepCyTOF_labels")





start <- Sys.time()

scoredf_ <- furrr::future_map_dfr(dc2_files, function(dn_files) {
})
time_output(start)






