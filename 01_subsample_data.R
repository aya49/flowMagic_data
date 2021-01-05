# date created: 2020-01-03
# author: alice yue
# input: 2D and nD csv/clr
# output: subsampled versions (because they're taking too much space...); each files should have max m rows; I'll experiment with m


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 5#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
libr(c(
  "furrr", #"rslurm",
  "data.table"
))



