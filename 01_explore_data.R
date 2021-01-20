# date created: 2020-01-03
# author: alice yue
# input: 2D and nD csv/clr
# output: explore data


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 5#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
libr(c(
  "furrr", #"rslurm",
  "data.table"
))


## input ####
out_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/data"
x2_dir <- paste0(out_dir,"/2D/x")
y2_dir <- paste0(out_dir,"/2D/y")
xn_dir <- paste0(out_dir,"/nD/x")
yn_dir <- paste0(out_dir,"/nD/y")

x2_dirs <- list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)
xn_dirs <- list.dirs(xn_dir, recursive=TRUE, full.names=TRUE)
x_dirs <- append(x2_dirs, xn_dirs)


for (x_dir in x_dirs) {
  x_files <- sample(list.files(x_dir, full.names=TRUE, ".csv.gz"),30)
  if (length(x_files)==0) next
  cat("\n\n", x_dir)
  
  x_lengths <- sapply(x_files, function(x) length(count.fields(x,skip=1)))
  ys <- sapply(x_files, function(x) data.table::fread(gsub("/x/","/y/",x)))

  if (nrow(ys[[1]])!=x_lengths[1])
    cat("\n- diff num cells")
  if (any(rowSums(ys[[1]])==0))
    cat("\n- not all labelled")
  cat("\n- files:",length(x_files))
  cat("\n- mean number of cells:", mean(x_lengths))
  cat("\n- number of cell pops:", ncol(ys[[1]]), "(",paste0(colnames(ys[[1]]), collapse=", "),")")
}



