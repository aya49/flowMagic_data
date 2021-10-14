# date created: 2020-01-04
# author: alice yue
# input: method res vector (for all support sizes 1-5,10,15,20; data, scatterplot, sample)
# output: score matrix (f1)

# conda: R 3.6; conda install -c r r; conda install -c r r-essentials 


## set directory, load packages, set parallel ####
no_cores <- 32#parallel::detectCores() - 5
# root <- "/home/ayue/projects/flowMagic_data"
# install.packages("devtools")
root <- "/home/aya43/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))
future::plan(future::multisession, workers=no_cores) # for furrr

## input ####
m2_dir <- paste0(results_dir,"/2D/method"); 


## load inputs ####
dc2_dirs <- list_leaf_dirs(m2_dir)
# dc2_files <- list.files(dc2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")


## output ####
gs_xr_ <- function(x,y) gs_xr(x,y,"scores") 
# plyr::l_ply(append(list_leaf_dirs(dc2_dir),list_leaf_dirs(dcn_dir)), function(x)
#   dir.create(gs_xr_(x,"deepCyTOF_labels"), recursive=TRUE, showWarnings=FALSE) )


# # get train samples
# trs <- gs_xr(dc2_dir,"x_2Ddensity_euclidean_rankkmed")


## START ####
start <- Sys.time()

for (dc_dir_ in dc2_dirs) {
  start1 <- Sys.time()
  print(dc_dir_)
  
  dc_files <- list.files(dc_dir_, full.names=TRUE)

  li <- loop_ind_f(dc_files, no_cores)
  bests <- furrr::map(li, function(dc_fs) plyr::ldply(dc_fs, function(dc_f) {
    x2predicted <- read.csv(dc_f, header=FALSE)
    
    x2discrete_file <- gsub("results","data",gsub("method[/]setr[/]1","x_2Ddiscrete", dc_f))
    x2discrete <- read.csv(x2discrete_file, header=FALSE)
    
    ypred <- apply(x2discrete, 1, function(xy) x2predicted[xy[1], xy[2]])
    
    # yactual_file <- gsub("results","data",gsub("method[/]setr[/]1","y_vector_", dc_f))
    # yactual <- read.csv(yactual_file, header=FALSE)[,1]
    
    yactualfull_file <- gsub("results","raw",gsub("method[/]setr[/]1","y", dc_f))
    actual <- read.csv(yactualfull_file, check.names=FALSE)
    cpops <- colnames(actual)
    
    # get meta data
    # pus <- sort(unique(predicted)) # unique labels
    path_stuff <- stringr::str_split(dc_f,"/")[[1]]
    di <- which(path_stuff=="method")
    method <- path_stuff[di+1]
    shots <- as.numeric(path_stuff[di+2])
    dset <- path_stuff[di+3]
    scat <- path_stuff[di+4]
    fname <- gsub(".csv.gz","",path_stuff[length(path_stuff)])
    
    ## score each cpop ####
    plyr::ldply(seq_len(ncol(actual)), function(cpopi) { 
      cbind(data.frame(
        method=method,
        dataset=dset, scatterplot=scat, cpop=cpops[cpopi], 
        train_no=shots, fcs=fname, 
        train=NA
      ), f1score(actual[,cpopi]==1, ypred==cpopi)) ### +1 ?????
    })
  }))
  bests <- Reduce(rbind, bests)
  bests <- bests[,colnames(bests)!=".id"]
  score_file <- paste0(gs_xr_(dc_dir_,"method"),".csv.gz") ### ????
  dir.create(folder_name(score_file), recursive=TRUE, showWarnings=FALSE)
  write.table(bests, file=gzfile(score_file), sep=",", row.names=FALSE, col.names=TRUE)
  time_output(start1)
}
time_output(start)



