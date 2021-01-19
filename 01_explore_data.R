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
  x_files <- list.files(x_dir, full.names=TRUE, ".csv.gz")
  if (length(x_files)==0) next
  cat("\n\n", x_dir)
  y1 <- data.table::fread(gsub("/x/","/y/",x_files[1]))
  cat("\n- files:",length(x_files))
  cat("\n- mean number of cells:", mean(sapply(sample(x_files,30), function(x) length(count.fields(x,skip=1)))))
  cat("\n- number of cell pops:", ncol(y1), "(",paste0(colnames(y1), collapse=", "),")")
}



