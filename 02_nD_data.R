# date created: 2020-01-04
# author: alice yue
# input: nD csv
# output: Rtsne 2D representations
set.seed(1)

## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
# root <- "/home/ayue/projects/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## output ####
xn_dirs <- list_leaf_dirs(xn_dir)
gs_xr <- function(x) gsub("data/nD/x/","results/nD/Rtsne/",x)
plyr::l_ply(gs_xr(xn_dirs), dir.create, recursive=TRUE, showWarnings=FALSE)


## load inputs ####
xn_files <- list.files(xn_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")

xn_files <- xn_files[sapply(xn_files, function(x) {
  if (file.exists(gs_xr(x))) if (file.size(gs_xr(x))>0) return(FALSE)
  return(TRUE)
})]


## START ####
start <- Sys.time()

# loop_ind <- loop_ind_f(sample(fe), no_cores)
# plyr::l_ply(loop_ind, function(ii) { plyr::l_ply(ii, function(i) {
fe <- seq_len(length(xn_files))
for (i in fe) {
  start1 <- Sys.time()
  xn_f <- xn_files[i]
  
  xn <- data.table::fread(xn_f, data.table=FALSE)
  tx <- Rtsne::Rtsne(xn[!duplicated(xn),,drop=FALSE])$Y
  colnames(tx) <- c("tsne 1", "tsne 2")
  write.table(tx, file=gzfile(gs_xr(xn_f)), sep=",", row.names=FALSE, col.names=TRUE)
  
  time_output(start1, xn_f)
}
# }) }, .parallel=TRUE)
time_output(start)






