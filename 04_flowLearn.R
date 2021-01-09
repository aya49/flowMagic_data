# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 200x200


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 15#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
libr(c(
  "furrr", #"rslurm",
  "flowLearn" # devtools::install_github("mlux86/flowLearn")
))


## input ####
out_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/data"
x2_dir <- paste0(out_dir,"/2D/x"); 
dir.create(x2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr
y2_dir <- paste0(out_dir,"/2D/y"); 
dir.create(y2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr


## output ####
for (dl in list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)) {
  dir.create(gsub("/x/","/results/flowLearn/",dl), recursive=TRUE, showWarnings=FALSE)
}


## load inputs ####
x2_files <- list.files(x2_dir, recursive=TRUE, full.names=TRUE)


## START ####
start <- Sys.time()

loop_ind <- loop_ind_f(seq_len(length(x2_files)), no_cores)
res <- furrr::future_map(loop_ind, function(ii) { purrr::map(ii, function(i) {
  # load csv
  x2discrete <- x2 <- data.table::fread(x2_files[i], data.table=FALSE)
  y2 <- data.table::fread(gsub("/x/","/y/",x2_files[i]), data.table=FALSE)
  

}) })
time_output(start)



