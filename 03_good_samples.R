# date created: 2020-01-04
# author: alice yue
# input: 2D density (unnormalized) + scatterplot with no density + scatterplot with density
# output: distance matrices between samples of the same 2D gate


## set directory, load packages, set parallel ####
no_cores <- 20#parallel::detectCores() - 5
root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/home/ayue/projects/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
ingrid <- "x_2Ddensity"
x2_dir <- paste0(root,"/results/2D/",ingrid); 


## output ####
distn <- "euclidean"
clustn <- "pam"
for (dl in list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)) {
  clus_dir <- gsub(paste0("2D/",ingrid),paste0("2D/",ingrid,"/",distn,"/",clustn),dl)
  dist_dir <- gsub(paste0("2D/",ingrid),paste0("2D/",ingrid,"/",distn),dl)
  dir.create(c(clus_dir, dist_dir), recursive=TRUE, showWarnings=FALSE)
}


## load inputs ####
x2_dirs <- list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)
x2_dirs <- x2_dirs[sapply(x2_dirs, function(x) sum(grepl(x,x2_dirs))==1)]
x2_files <- plyr::llply(x2_dirs, list.files, full.names=TRUE, pattern="csv")

ks <- 1:10


## START ####
start <- Sys.time()


# load csv
# loop_ind <- loop_ind_f(sample(seq_len(length(x2_files))), no_cores)
plyr::l_ply(x2_files, function(x2_fs) {
  dist_dir <- paste0(gsub(paste0("2D/",ingrid),paste0("2D/",ingrid,"/",distn), folder_name(x2_fs[1])),".Rdata")

  # load files into matrix
  start1 <- Sys.time()
  
  all2D <- purrr::map(x2_fs, function(x2_f) {
    a <- data.table::fread(x2_f, data.table=FALSE)
    as.vector(as.matrix(a/max(a)))
  })
  fnames <- names(all2D) <- gsub(".csv.gz","",sapply(x2_fs, file_name))
  all2D <- as.matrix(dplyr::bind_cols(all2D))
  if (nrow(all2D)!=length(x2_fs)) all2D <- t(all2D)
  time_output(start1)
  
  
  # make distance object from matrix
  start1 <- Sys.time()
  d <- dist(all2D, method="euclidean") # "manhattan", "canberra", "binary" or "minkowski"
  save(d, file=dist_dir)
  
  
  # clust
  ddim <- nrow(as.matrix(d))
  pk <- plyr:llply(ks[ks<ddim], function(x) cluster::pam(d,k=x))
  time_output(start1)
  
  # save k medoid results
  cl_dir <- gsub(".Rdata","",gsub(distn,paste0(dist,"_",clustn),dist_dir))
  
  plyr::l_ply(pk, function(x) {
    plyr::l_ply(seq_len(length(x$medoids)), function(medi) {
      if (length(x$medoids)==0) return(NULL)
      
      filesi <- names(x$clustering[x$clustering==medi])
      filesi <- filesi[!grepl(x$medoids[medi], filesi)]
      
      cli_dir <- paste0(cl_dir,length(x$medoids))
      dir.create(cli_dir, showWarnings=FALSE)
      if (length(filesi)>0) {
        write.table(filesi, file=gzfile(paste0(cli_dir,"/",x$medoids[medi])), 
                    row.names=FALSE, col.names=FALSE)
      }
    })
  })
  time_output(start1)
}, .parallel=TRUE)
time_output(start)

