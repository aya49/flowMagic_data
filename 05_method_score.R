# date created: 2020-01-04
# author: alice yue
# input: method res vector (for all support sizes 1-5,10,15,20; data, scatterplot, sample)
# output: score matrix (f1)

# conda: R 3.6; conda install -c r r; conda install -c r r-essentials 


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
# root <- "/home/ayue/projects/flowMagic_data"
# install.packages("devtools")
root <- "/home/aya43/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
m2_dir <- paste0(results_dir,"/method"); 


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
  bests <- plyr::ldply(li, function(dc_fs) plyr::ldply(dc_fs, function(dc_f) {
    predicted <- read.csv(dc_f, header=FALSE)[,1]
    actual_f <- stringr::str_split(predicted,"/")[[1]]
    suppressWarnings({
        actual_f <- paste0(actual_f[is.na(sapply(actual_f, as.numeric))], collapse="/")
    })
    actual_f <- paste0(gs_xr(dc_f,"y","raw"),".gz")
    actual <- data.table::fread(actual_f, data.table=FALSE)
    cpops <- colnames(actual)
    
    # get meta data
    # pus <- sort(unique(predicted)) # unique labels
    path_stuff <- stringr::str_split(dc_f,"/")[[1]]
    di <- which(path_stuff=="method")
    method <- path_stuff[di+1]
    shots <- path_stuff[di+2]
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
      ), f1score(actual[,cpopi]==1, predicted==cpopi)) ### +1 ?????
    })
  }), .parallel=TRUE)
  bests <- bests[,colnames(bests)!=".id"]
  score_file <- paste0(gs_xr_(dc_dir_,method),".csv.gz") ### ????
  dir.create(folder_name(score_file), recursive=TRUE, showWarnings=FALSE)
  write.table(bests, file=gzfile(score_file), sep=",", row.names=FALSE, col.names=TRUE)
  time_output(start1)
}
time_output(start)



