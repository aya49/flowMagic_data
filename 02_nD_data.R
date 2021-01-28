# date created: 2020-01-04
# author: alice yue
# input: nD csv
# output: Rtsne 2D representations


## set directory, load packages, set parallel ####
no_cores <- 10#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/home/ayue/projects/flowMagic_data"
# root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## output ####
xn_dirs <- list_leaf_dirs(xn_dir)
gs_xr <- function(x) gsub("data/nD/x/","results/nD/Rtsne/",x)
plyr::l_ply(gs_xr(xn_dirs), dir.create, recursive=TRUE, showWarnings=FALSE)


## load inputs ####
xn_files <- list.files(xn_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")


## START ####
start <- Sys.time()

res <- plyr::llply(xn_files, function(xn_file) {
  xn <- data.table::fread(xn_file, data.table=FALSE)
  tx <- Rtsne::Rtsne(xn[!duplicated(xn),,drop=FALSE])$Y
  tx_file <- gs_xr(xn_file)
  tx_dir <- stringr::str_split(tx_file,"/")[[1]]
  tx_dir <- paste0(tx_dir[-length(tx_dir)],collapse="/")
  colnames(tx) <- c("tsne 1", "tsne 2")
  write.table(tx, file=gzfile(tx_file), sep=",", row.names=FALSE, col.names=TRUE)
}, .parallel=TRUE)
time_output(start)




