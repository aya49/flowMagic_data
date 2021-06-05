# date created: 2020-01-04
# author: alice yue
# input: deepCyTOF
# output: score matrix (f1)


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
# root <- "/home/ayue/projects/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
dc2_dir <- paste0(results_dir,"/2D/deepCyTOF_labels"); 
dcn_dir <- paste0(results_dir,"/nD/deepCyTOF_labels"); 


## load inputs ####
dc2_dirs <- list_leaf_dirs(dc2_dir)
dcn_dirs <- list_leaf_dirs(dcn_dir)
dc2_files <- list.files(dc2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")
dcn_files <- list.files(dcn_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")


## output ####
gs_xr_ <- function(x,y) gs_xr(x,y,"scores") 
# plyr::l_ply(append(list_leaf_dirs(dc2_dir),list_leaf_dirs(dcn_dir)), function(x)
#   dir.create(gs_xr_(x,"deepCyTOF_labels"), recursive=TRUE, showWarnings=FALSE) )


# get train samples
trs <- gs_xr(dc2_dir,"x_2Ddensity_euclidean_rankkmed")


## START ####
start <- Sys.time()

for (dc_dir_ in dc2_dirs) {
  start1 <- Sys.time()
  print(dc_dir_)
  
  dc_files <- list.files(dc_dir_, full.names=TRUE)
  
  actual_f <- paste0(gs_xr(dc_files[1],"y","raw"),".gz")
  actual_f <- stringr::str_split(actual_f,"/")
  k <- as.numeric(actual_f[[1]][!is.na(sapply(actual_f[[1]], as.numeric))])
  
  nD <- grepl("/nD/",dc_dir_)
  if (!nD) {
    tr_fnames <- gsub(".csv.gz","",list.files(paste0(gs_xr(dc_dir_,"x_2Ddensity_euclidean_rankkmed"),"/", k)))
  } else {
    tr_fnames <- gsub(".csv","",sapply(dc_files[round(seq(from=2, to=length(dc_files), length=10))], file_name))
  }

  li <- loop_ind_f(dc_files, no_cores)
  bests <- plyr::ldply(li, function(dc_fs) plyr::ldply(dc_fs, function(dc_f) {
    predicted <- read.csv(dc_f, header=FALSE)[,1]
    actual_f <- paste0(gs_xr(dc_f,"y","raw"),".gz")
    actual_f <- stringr::str_split(actual_f,"/")
    suppressWarnings({
      actual_f <- paste0(actual_f[[1]][is.na(sapply(actual_f[[1]], as.numeric))], collapse="/")
    })
    actual <- data.table::fread(actual_f, data.table=FALSE)
    cpops <- colnames(actual)
    
    # get meta data
    # pus <- sort(unique(predicted)) # unique labels
    nD <- grepl("/nD/", dc_f)
    path_stuff <- stringr::str_split(dc_f,"/")[[1]]
    di <- which(path_stuff=="deepCyTOF_labels")
    dset <- path_stuff[di+2]
    scat <- ifelse(!nD,path_stuff[di+3],NA)
    fname <- gsub(".csv","",path_stuff[length(path_stuff)])
    
    ## score each cpop ####
    plyr::ldply(seq_len(ncol(actual)), function(cpopi) { 
      cbind(data.frame(
        method="deepCyTOF",
        dataset=dset, scatterplot=scat, cpop=cpops[cpopi], 
        train_no=k, fcs=fname, 
        train=fname%in%tr_fnames
      ), f1score(actual[,cpopi]==1, predicted==cpopi))
    })
  }), .parallel=TRUE)
  bests <- bests[,colnames(bests)!=".id"]
  score_file <- paste0(gs_xr_(dc_dir_,"deepCyTOF"),".csv.gz")
  dir.create(folder_name(score_file), recursive=TRUE, showWarnings=FALSE)
  write.table(bests, file=gzfile(score_file), 
              sep=",", row.names=FALSE, col.names=TRUE)
  time_output(start1)
}
time_output(start)



