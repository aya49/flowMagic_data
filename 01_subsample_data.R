# date created: 2020-01-03
# author: alice yue
# input: 2D and nD csv/clr
# output: subsampled versions (because they're taking too much space...); each files should have max m rows; I'll experiment with m


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))



