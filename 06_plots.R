# date created: 2020-01-04
# author: alice yue
# input: score matrices
# output: plots


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
sc2_dir <- paste0(root,"/scores/2D")
scn_dir <- paste0(root,"/scores/nD")


## load scores ####
start <- Sys.time()
sc2_files <- list.files(sc2_dir, full.names=TRUE, recursive=TRUE, pattern=".csv.gz")
score2D <- plyr::ldply(sc2_files, data.table::fread, data.table=FALSE)
write.csv(score2D, file=gzfile(paste0(sc2_dir,".csv.gz")))
time_output(start)

start <- Sys.time()
scn_files <- list.files(scn_dir, full.names=TRUE, recursive=TRUE, pattern=".csv.gz")
scorenD <- plyr::ldply(scn_files, data.table::fread, data.table=FALSE)
write.csv(scorenD, file=gzfile(paste0(scn_dir,".csv.gz")))
time_output(start)





