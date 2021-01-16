# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 200x200


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 10#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
source("helpers_2D.R")
libr(c(
  "furrr", #"rslurm",
  "flowLearn" # devtools::install_github("mlux86/flowLearn")
))


## input ####
out_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/data"
x2_dir <- paste0(out_dir,"/2D/x"); 
y2_dir <- paste0(out_dir,"/2D/y"); 
t2_dir <- paste0(out_dir,"/2D/thresholds"); 


## output ####
score_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/scores"
dir.create(score_dir, recursive=TRUE, showWarnings=FALSE)
fl2_dir <- paste0(gsub("/data","/results",out_dir),"/2D/flowLearn")


## load inputs ####
t2_dirs <- list.dirs(t2_dir)[-1] # includes t2_dir so remove it
x2_folders <- list.dirs(x2_dir, recursive=TRUE)
x2_folders <- x2_folders[sapply(x2_folders, function(x) sum(grepl(x,x2_folders))==1)]



dfn <- 512 # number of density features for flowlearn
ks <- 1:10 # number of training samples



## START ####
start <- Sys.time()

for (t2_dir_ in t2_dirs) {
  print(paste0("data set:",t2_dir_))
  idataset <- file_name(t2_dir_)
  t2_files <- list.files(t2_dir_, full.names=TRUE)
  
  # get file names, all thresholds, and scatterplot names
  fnames_ <- purrr::map(t2_files, function(t2_f)  gsub(".Rdata","",file_name(t2_f)))
  t2sp <- purrr::map(t2_files, function(t2_file) get(load(t2_file)))
  names(t2sp) <- fnames_
  scats <- names(t2sp[[1]])
  
  for (scat in scats) { 
    try ({
    start1 <- Sys.time()
    cat(scat)
    
    # get the scatterplot name
    x2_diri <- paste0(x2_dir,"/",idataset,"/",scat)
    
    # get x, y, threshold files ^
    x2_files <- list.files(x2_diri, full.names=TRUE)
    x2_files_fn <- gsub(".csv.gz","",sapply(x2_files, file_name))
    x2s <- purrr::map(fnames_, function(fname) {
      xind <- which(x2_files_fn==fname); if (length(xind)==0) return(NULL)
      data.table::fread(x2_files[xind], data.table=FALSE) })
    names(x2s) <- fnames_
    x2s <- plyr::compact(x2s)
    fnames <- names(x2s)
    
    t2s <- purrr::map(t2sp[fnames], function(t2si) t2si[[scat]])
    names(t2s) <- fnames
    
    # prepare to save predicted thresholds
    fl_dir <- gsub("/data/2D/x","/results/2D/flowLearn",x2_diri)

    time_output(start1, "loaded files")
    
    # for each threshold
    for (t2name in names(t2s[[1]])) {
      # make a density data object for flowLearn
      densdat <- new('DensityData')
      for (fname in fnames) {
        # some files don't have certain thresholds because they are un-needed
        if (is.na(t2s[[fname]][t2name])) next 
        # there might be error files, I will remove these!
        # if (!t2name%in%colnames(x2s[[fname]])) {
        #   # file.remove(x2_files[grep(fname, x2_files)])
        #   cat("\nactual markers:", paste0(x2s, collapse=", "), "not ", t2name)
        #   next
        # }
        df <- flowLearn::flEstimateDensity(x2s[[fname]][,t2name], dfn)
        densdat <- flowLearn::flAdd(
          densdat, fname, "cpop", 1, df$x, df$y, NaN, t2s[[fname]][t2name])
      }
    #   dd_dir <- gsub("flowLearn","flowLearn_densdat",fl_dir)
    #   dir.create(dd_dir, recursive=TRUE, showWarnings=FALSE)
    #   save(densdat, file=paste0(dd_dir,"/",t2name,".Rdata"))
    # }
    # 
    # for (t2name in names(t2s[[1]])) {
    #   dd_dir <- gsub("flowLearn","flowLearn_densdat",fl_dir)
    #   densdat <- get(load(paste0(dd_dir,"/",t2name,".Rdata")))
      
      # for each number of training data, predict threshold
      for (k in ks) {
        protoIdx <- flowLearn::flSelectPrototypes(densdat, k)
        ddp <- flowLearn::flPredictThresholds(densdat, protoIdx)
        
        ddpt <- ddp@data$gate.high
        ddpt[rownames(densdat@data)%in%protoIdx] <- NA
        names(ddpt) <- densdat@data$fcs
        
        dt_dir <- paste0(gsub("flowLearn","flowLearn_thresholds",fl_dir),"/",t2name)
        dir.create(dt_dir, recursive=TRUE, showWarnings=FALSE)
        save(ddpt, file=paste0(dt_dir,"/",k,".Rdata"))
      }
    }
    })
    time_output(start1)
  }
}
time_output(start)




