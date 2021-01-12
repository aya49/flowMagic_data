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
  fnames <- purrr::map(t2_files, function(t2_f)  gsub(".Rdata","",file_name(t2_f)))
  t2sp <- purrr::map(t2_files, function(t2_file) get(load(t2_file)))
  names(t2sp) <- fnames
  scats <- names(t2sp[[1]])
  
  for (scat in scats) { 
    start1 <- Sys.time()
    cat(scat)
    
    # get the scatterplot name
    x2_diri <- paste0(x2_dir,"/",idataset,"/",scat)
    
    # get x, y, threshold files ^
    x2_files <- list.files(x2_diri, full.names=TRUE)
    x2s <- purrr::map(fnames, function(fname) 
      data.table::fread(x2_files[grep(fname,x2_files)], data.table=FALSE))
    
    t2s <- purrr::map(t2sp, function(t2si) t2si[[scat]])
    names(x2s) <- names(t2s) <- fnames
    
    # prepare to save predicted thresholds
    fl_dir <- gsub("/data/2D/x","/results/2D/flowLearn",x2_diri)

    time_output(start1, "loaded files")
    
    # for each threshold
    for (t2name in names(t2s[[1]])) {
      # make a density data object for flowLearn
      densdat <- new('DensityData')
      for (fname in fnames) {
        if (is.na(t2s[[fname]][t2name])) next # some files don't have certain thresholds because they are un-needed
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
        names(ddpt) <- fnames
        
        dt_dir <- paste0(gsub("flowLearn","flowLearn_thresholds",fl_dir),"/",t2name)
        dir.create(dt_dir, recursive=TRUE, showWarnings=FALSE)
        save(ddpt, file=paste0(dt_dir,"/",k,".Rdata"))
      }
    }
    time_output(start1)
  }
}
time_output(start)


## load predicted thresholds and get f1 ####
f1score <- function(tfactual, tfpred) {
  ap <- sum(tfactual & tfpred)
  ap0 <- ap==0
  
  pred <- sum(tfpred)
  actu <- sum(tfactual)
  
  if (ap0) {
    precision <- recall <- f1 <- 0
  } else {
    precision <- ap/pred
    recall <- ap/actu
    
    f1 <- ifelse(precision+recall == 0, 0, 
                 2 * precision * recall/(precision + recall)) 
  }
  data.frame(
    precision=precision, recall=recall, f1=f1, 
    true_proportion=actu/length(tfactual), 
    predicted_proportion=pred/length(tfpred),
    true_size=actu,
    predicted_size=pred)
}

flPlot <- function(x2, marknames, ft, fto, main) {
  layout(matrix(c(2,2,1,4,4,3,4,4,3), 3,3, byrow = TRUE))
  par(mar=c(1,3,3,3)) # bottom, left, top, right
  plot.new()
  par(mar=c(1,4,5,1))
  plot(density(x2[,1]), xlab=NA, xaxt="n", 
       main=main)
  if (colnames(x2)[1]%in%marknames) {
    abline(v=fto[colnames(x2)[1]], col="red")
    if (!is.na(ft[colnames(x2)[1]])) 
      abline(v=ft[colnames(x2)[1]], col="blue")
  }
  dens <- density(x2[,2])
  densymax <- max(dens$y)
  par(mar=c(4,1,1,1))
  plot(dens$y, dens$x, ylab=NA, yaxt="n",
       type="l", xlab="Density", 
       main=ifelse(is.na(ft[colnames(x2)[2]])," (train sample)",""))
  if (colnames(x2)[1]%in%marknames) {
    abline(h=fto[colnames(x2)[2]], col="red")
    if (!is.na(ft[colnames(x2)[2]])) 
      abline(h=ft[colnames(x2)[2]], col="blue")
  }
  par(mar=c(4,4,1,1))
  plot_dens(
    x2, main=NA, xlab=paste0(
      colnames(x2)[1], ifelse(is.na(ft[colnames(x2)[1]])," (train sample)","")),
    ylab=paste0(paste0(
      colnames(x2)[2], ifelse(is.na(ft[colnames(x2)[2]])," (train sample)",""))))
  if (colnames(x2)[1]%in%marknames) {
    abline(v=fto[colnames(x2)[1]], col="red")
    if (!is.na(ft[colnames(x2)[1]])) 
      abline(v=ft[colnames(x2)[1]], col="blue")
  }
  if (colnames(x2)[1]%in%marknames) {
    abline(h=fto[colnames(x2)[2]], col="red")
    if (!is.na(ft[colnames(x2)[2]])) 
      abline(h=ft[colnames(x2)[2]], col="blue")
  }
  
}

start <- Sys.time()

scoredf_ <- purrr::map(x2_folders, function(x2_fold) {
  ft_dir <- gsub("data/2D/x","results/2D/flowLearn_thresholds",x2_fold)
  if (!dir.exists(ft_dir)) next
  
  marknames <- list.dirs(ft_dir, full.names=FALSE)[-1]
  scat <- file_name(x2_fold)
  
  # get x, y, threshold files ^
  
  fto_dir <- gsub(paste0("/",scat),"",gsub("/x/","/thresholds/",x2_fold))
  dataset <- file_name(fto_dir)
  
  print(paste0("data/scatterplot: ",dataset,"/",scat," (",paste0(marknames,collapse=", "),")"))

  fto_files <- list.files(fto_dir, full.names=TRUE)
  ftos <- purrr::map(fto_files, function(fto_file) get(load(fto_file))[[scat]])
  
  fnames <- gsub(".Rdata","",sapply(fto_files, file_name))
  
  x2_files <- list.files(x2_fold, full.names=TRUE)
  x2s <- purrr::map(fnames, function(fname) 
    data.table::fread(x2_files[grep(fname,x2_files)], data.table=FALSE))
  
  y2_files <- list.files(gsub("/x/","/y/",x2_fold), full.names=TRUE)
  y2s <- purrr::map(fnames, function(fname) 
    data.table::fread(y2_files[grep(fname,y2_files)], data.table=FALSE))
  
  names(x2s) <- names(y2s) <- names(ftos) <- fnames
  

  # for each number of train data used
  # (not very efficient, actual values are the same)
  ks_ <- as.numeric(gsub(".Rdata","",list.files(paste0(ft_dir,"/",marknames[1]), full.names=FALSE)))
  
  tempmark <- rep(NA,length(marknames))
  names(tempmark) <- marknames
  scoredf_k_ <- furrr::future_map_dfr(ks_, function(k) {
    fts <- lapply(fnames, function(x) tempmark)
    names(fts) <- fnames
    for (mn in marknames) {
      markthresi <- get(load(paste0(ft_dir,"/",mn,"/",k,".Rdata")))
      for (fname in fnames) 
        fts[[fname]][mn] <- markthresi[fname]
    }
    
    # plot
    markthresi_fol <- paste0(gsub("thresholds","plots",ft_dir),"/",k)
    dir.create(markthresi_fol, recursive=TRUE, showWarnings=FALSE)
    scoredf_fname_ <- purrr::map(fnames, function(fname) {
      x2 <- x2s[[fname]]
      fto <- ftos[[fname]]
      ft <- fts[[fname]]
      
      png(paste0(markthresi_fol,"/",fname,
                 ifelse(any(is.na(ft)),"_train",""),".png"), 
          width=500, height=500)
      flPlot(x2, marknames, ft, fto, 
             main=paste0("data: ", dataset, "\nscatterplot: ",scat,
                         "\n(blue=predicted, red=actual)"))
      graphics.off()
    })
    
    # for each cell population
    scoredf_cpop_ <- purrr::map_dfr(colnames(y2s[[1]]), function(cpop) {
      cpop_ <- flowLearn::flNormalizePopulationName(cpop)
      
      cpopx <- x2s[[1]][as.logical(y2s[[1]][,cpop]),,drop=FALSE]
      fto <- ftos[[1]]
      ft <- fts[[1]]
      
      mn_tf <- tempmark
      for (mn in marknames) {
        if (fto[mn] <= min(cpopx[,mn])) {
          mn_tf[mn] <- 1 # lower threshold
        } else if (fto[mn] >= max(cpopx[,mn])) {
          mn_tf[mn] <- 0 # higher threshold
        }
      }
      
      # for each fcs file
      scoredf_fname_ <- purrr::map_dfr(fnames, function(fname) {
        x2 <- x2s[[fname]]
        fto <- ftos[[fname]]
        ft <- fts[[fname]]
        ftna <- is.na(ft)
        ft[ftna] <- fto[ftna]

        tfactual <- tfpred <- rep(TRUE, nrow(x2))
        for (mn in marknames) {
          if (is.na(mn_tf[mn])) next
          if (mn_tf[mn]==1) {
            if (!is.na(fto[mn])) tfactual <- tfactual & fto[mn]<=x2[,mn]
            if (!is.na(ft[mn])) tfpred <- tfpred & ft[mn]<=x2[,mn]
          } else if (mn_tf[mn]==0) {
            if (!is.na(fto[mn])) tfactual <- tfactual & fto[mn]>x2[,mn]
            if (!is.na(ft[mn])) tfpred <- tfpred & ft[mn]>x2[,mn]
          }
        }
        
        cbind(data.frame(
          data_set=dataset, scatter_plot=scat, cell_population=cpop, 
          train_samples=k, total_samples=length(fnames), fcs=fname, 
          train=any(is.na(fts[[fname]]))
        ), f1score(tfactual, tfpred))
      })
      # return(dplyr::bind_rows(scoredf_fname_))
    })
    # return(dplyr::bind_rows(scoredf_cpop_))
  })
  # return(dplyr:::bind_rows(scoredf_k_))
})
scoredf <- dplyr::bind_rows(scoredf_)
write.table(scoredf, file=gzfile(paste0(score_dir,"/flowLearn.csv.gz")), sep=",")
time_output(start)



