# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 400x400


## parallelization ####
# future::plan(future::multiprocess)
# no_cores <- 15#parallel::detectCores() - 5
# doMC::registerDoMC(no_cores)

## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
setwd(root)


## packages ####
source("src/helpers.R")
libr(c(
  "furrr", #"rslurm",
  "data.table"
))


## input ####
xn_dir <- paste0(root,"/data/nD/x"); 


## output ####



## load inputs ####
xn_files <- list.files(xn_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")


## START ####
start <- Sys.time()

res <- purrr::map(xn_files, function(xn_file) {
  xn <- data.table::fread(xn_file, data.table=FALSE)
  tx <- Rtsne::Rtsne(xn[!duplicated(xn),,drop=FALSE])$Y
  tx_file <- gsub("/x/","/Rtsne/",gsub("data/","results/",xn_file))
  tx_dir <- stringr::str_split(tx_file,"/")[[1]]
  tx_dir <- paste0(tx_dir[-length(tx_dir)],collapse="/")
  dir.create(tx_dir, recursive=TRUE, showWarnings=FALSE)
  write.table(tx, file=gzfile(tx_file), sep=",", row.names=FALSE, col.names=TRUE)
})
time_output(start)




