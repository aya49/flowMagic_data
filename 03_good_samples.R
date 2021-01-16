# date created: 2020-01-04
# author: alice yue
# input: 2D density (unnormalized) + scatterplot with no density + scatterplot with density
# output: distance matrices between samples of the same 2D gate


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 10#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
libr(c(
  "furrr" #"rslurm",
))


## input ####
out_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"


## output ####
x2_dir <- paste0(out_dir,"/results/2D/x_2Ddensity"); 
for (dl in list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)) {
  dir.create(gsub("/results/2D/x/","/results/2D/cluster/",dl), recursive=TRUE, showWarnings=FALSE)
}


## load inputs ####
x2_dirs <- list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)
x2_dirs <- x2_dirs[sapply(x2_dirs, function(x) sum(grepl(x,x2_dirs))==1)]
x2_files <- plyr::llply(x2_dirs, function(x) list.files(x, full.names=TRUE, pattern="csv"))

ks <- 1:10


## START ####
start <- Sys.time()


# load csv
loop_ind <- loop_ind_f(sample(seq_len(length(x2_files))), no_cores)
res <- furrr::future_map(loop_ind, function(ii) { purrr::map(ii, function(i) {
  
  dist_dir <- gsub(paste0("/",file_name(x2_files[[i]][1])),".Rdata",x2_files[[i]][1])
  dist_dir <- gsub("x_2Ddensity","xtrain_2Ddensity_euclidean",dist_dir)
  dir.create(gsub(paste0("/",file_name(dist_dir)),"",dist_dir), recursive=TRUE, showWarnings=FALSE)
  
  # don't overwrite
  if (file.exists(dist_dir))
    if (file.size(dist_dir)>0) return(NULL)
  
  # load files into matrix
  start1 <- Sys.time()
  x2_fs <- x2_files[[i]]
  all2D <- purrr::map(x2_fs, function(x2_f) {
    a <- data.table::fread(x2_f, data.table=FALSE)
    a <- a/max(a)
    # a <- data.table::fread(gsub("/x_2Ddensity/","/x_2Dscatter/",x2_f), data.table=FALSE)
    a <- as.vector(as.matrix(a))
    return(a)
  })
  names(all2D) <- sapply(x2_fs, file_name)
  all2D <- as.matrix(dplyr::bind_cols(all2D))
  if (nrow(all2D)!=length(x2_fs)) all2D <- t(all2D)
  time_output(start1)
  
  
  # make distance object from matrix
  start1 <- Sys.time()
  d <- dist(all2D, method="euclidean") # "manhattan", "canberra", "binary" or "minkowski"
  
  save(d, file=dist_dir)
  rm(d)
  
  time_output(start1)
}) })
time_output(start)


dist_files <- list.files(gsub("x_2Ddensity","xtrain_2Ddensity_euclidean",x2_dir),
                         recursive=TRUE, pattern=".Rdata", full.names=TRUE)

start1 <- Sys.time()

loop_ind <- loop_ind_f(sample(seq_len(length(dist_files))), no_cores)
res <- furrr::future_map(loop_ind, function(ii) { purrr::map(ii, function(i) {
  # k medoid with k=ks
  start1 <- Sys.time()
  dist_file <- dist_files[i]
  d <- get(load(dist_file))
  ddim <- nrow(as.matrix(d))
  pk <- purrr::map(ks[ks<ddim], function(x) cluster::pam(d,k=x))
  time_output(start1)
  
  # save k medoid results
  cl_dir <- gsub(".Rdata","",gsub("euclidean/","euclidean_pam/",dist_file))
  
  pkm <- purrr::map(pk, function(x) {
    kmm <- purrr::map(seq_len(length(x$medoids)), function(medi) {
      if (length(x$medoids)==0) return(NULL)
      
      filesi <- names(x$clustering[x$clustering==medi])
      filesi <- filesi[!grepl(x$medoids[medi], filesi)]
      
      cli_dir <- paste0(cl_dir,length(x$medoids))
      dir.create(cli_dir, recursive=TRUE, showWarnings=FALSE)
      if (length(filesi)>0) {
        write.table(filesi, file=gzfile(paste0(cli_dir,"/",x$medoids[medi])), 
                    row.names=FALSE, col.names=FALSE)
      }
    })
  })
}) })

time_output(start)

# res <- furrr::future_map(loop_ind, function(ii) { purrr::map(ii, function(i) {
#   # load csv
#   file.remove(paste0(gsub("/x/","/x_scatter/",x2_files[i]),".png"))
#   file.remove(paste0(gsub("/x/","/x_colscat/",x2_files[i]),".png"))
# }) })


