# date created: 2020-01-04
# author: alice yue
# input: 2D csv/clr + thresholds
# output: flowLearn predicted thresholds, plots, scores


## set directory, load packages, set parallel ####
no_cores <- 14#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## output ####
thres_dirs <- list_leaf_dirs(thres_dir)
plyr::l_ply(c(
  gsub("data/2D/thresholds","results/2D/flowLearn_plots",thres_dirs), 
  gsub("data/2D/thresholds","results/2D/flowLearn_thresholds",thres_dirs)), 
  dir.create, recursive=TRUE, showWarnings=FALSE)


## parameters ####
dfn <- 512 # number of density features for flowlearn
ks <- 1:10 # number of training samples


## plot function: scatterplot + densities ####
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


## START ####
start <- Sys.time()

res <- plyr::llply(thres_dirs, function(thres_dir_) { try ({
  thres_dir_s <- stringr::str_split(thres_dir_,"/")[[1]]
  scat <- thres_dir_s[length(thres_dir_s)]
  dset <- thres_dir_s[length(thres_dir_s)-1]
  
  start1 <- Sys.time()
  print(paste0(dset," > ", scat))
  
  t2_files <- list.files(thres_dir_, full.names=TRUE)
  fnames <- plyr::laply(t2_files, function(t2_f)  gsub(".Rdata","",file_name(t2_f)))

  # get x, thresholds
  x2_diri <- paste0(x2_dir,"/",dset,"/",scat)
  ftos <- purrr::map(t2_files, function(t2_file) get(load(t2_file)))
  x2s <- purrr::map(paste0(x2_diri,"/",fnames,".csv.gz"), 
                    data.table::fread, data.table=FALSE)
  y2s <- purrr::map(paste0(gsub("/x/","/y/",x2_diri),"/",fnames,".csv.gz"), 
                    data.table::fread, data.table=FALSE)
  names(y2s) <- names(x2s) <- names(ftos) <- fnames
  marknames <- names(ftos[[1]])
  
  time_output(start1, "loaded files")
  
  
  ## predict ####
  
  # for each threshold, predict for all files
  fl_dir <- gsub("/data/2D/x","/results/2D/flowLearn_thresholds",x2_diri)
  protoIdxs <- NULL
  for (markname in marknames) {
    # make a density data object for flowLearn
    densdat <- new('DensityData')
    for (fname in fnames) {
      if (is.na(ftos[[fname]][markname])) next 
      df <- flowLearn::flEstimateDensity(x2s[[fname]][,markname], dfn)
      densdat <- flowLearn::flAdd(
        densdat, fname, "cpop", 1, df$x, df$y, NaN, ftos[[fname]][markname])
    }

    # for each number of training data, predict threshold
    if (is.null(protoIdxs)) {
      protoIdxs <- plyr::llply(ks, function(k) 
        flowLearn::flSelectPrototypes(densdat, k))
      names(protoIdxs) <- as.character(ks)
    }
    for (k in ks) {
      if (k>=nrow(densdat@data)) break
      protoIdx <- protoIdxs[[as.character(k)]]
      ddp <- flowLearn::flPredictThresholds(densdat, protoIdx)
      
      ddpt <- ddp@data$gate.high
      ddpt[rownames(densdat@data)%in%protoIdx] <- NA
      names(ddpt) <- densdat@data$fcs
      
      dt_dir <- paste0(fl_dir,"/",markname)
      dir.create(dt_dir, showWarnings=FALSE)
      save(ddpt, file=paste0(dt_dir,"/",k,".Rdata"))
    }
  }
  time_output(start1, "predicted thresholds")
  
  
  ## plot ####
  ks_ <- as.numeric(gsub(".Rdata","",list.files(paste0(fl_dir,"/",names(ftos[[1]])[1]))))
  
  # load prediced thresholds
  fts <- lapply(names(ftos[[1]]), function(x) lapply(ks_, function(k) 
    get(load(paste0(fl_dir,"/",x,"/",k,".Rdata"))) ))
  names(fts) <- names(ftos[[1]])

  pl_dir <- gsub("thresholds","plots",fl_dir)
  dir.create(pl_dir, recursive=TRUE, showWarnings=FALSE)
  for (fname in fnames) 
    for (ki in ks_) {
      x2 <- x2s[[fname]]
      fto <- ftos[[fname]]
      ft <- sapply(fts, function(x) x[[ki]][fname]); names(ft) <- marknames
      png(paste0(pl_dir,"/",fname,"_",ks_[ki],".png"), width=400, height=400)
      flPlot(x2, marknames, ft, fto, main=paste0("data: ", dset, "\nscatterplot: ",scat, "\n(blue=predicted, red=actual)"))
      graphics.off()
    }
  time_output(start1, "plotted")
  
  
  ## score ####
  cpops <- colnames(y2s[[1]])
  tempmark <- rep(NA, length(marknames))
  names(tempmark) <- marknames
  
  # get how thresholds work with cpops
  mn_tfs <- plyr::llply(cpops, function(cpop) {
    cpopx <- x2s[[1]][y2s[[1]][,cpop]==1,,drop=FALSE]
    cxn <- nrow(cpopx)
    cpopx <- cpopx[sample(seq_len(cxn), min(cxn, 300)),,drop=FALSE]
    
    mn_tf <- tempmark
    for (mn in marknames)
      if (ftos[[1]][mn] <= min(cpopx[,mn])) {
        mn_tf[mn] <- 1 # lower threshold
      } else if (ftos[[1]][mn] >= max(cpopx[,mn])) {
        mn_tf[mn] <- 0 # higher threshold
      }
    return(mn_tf)
  })
  names(mn_tfs) <- cpops
  
  # get actual cell populations not from y but from thresholds, just to be safe
  cpop_fname_actual <- plyr::llply(cpops, function(cpop) {
    # for each fcs file
    fname_actual <- plyr::llply(fnames, function(fname) {
      x2 <- x2s[[fname]]
      fto <- ftos[[fname]]
      mn_tf <- mn_tfs[[cpop]]
      
      tfactual <- rep(TRUE, nrow(x2))
      for (mn in marknames) {
        if (is.na(mn_tf[mn]) | is.na(fto[mn])) next
        if (mn_tf[mn]==1) {
          tfactual <- tfactual & fto[mn]<=x2[,mn]
        } else if (mn_tf[mn]==0) {
          tfactual <- tfactual & fto[mn]>x2[,mn]
        }
      }
      return(tfactual)
    })
    names(fname_actual) <- fnames
    return(fname_actual)
  })
  names(cpop_fname_actual) <- cpops
  
  # for each number of train data used; cell population
  fl <- length(fnames)
  for (ki in seq_len(length(ks_))) for (cpop in cpops) {
    k <- ks_[ki]
    mn_tf <- mn_tfs[[cpop]]
    # for each fcs file
    scoredf_fname_ <- purrr::map_dfr(fnames, function(fname) {
      print(fname)
      x2 <- x2s[[fname]]
      fto <- ftos[[fname]]
      ft <- sapply(fts, function(x) x[[ki]][fname]); names(ft) <- marknames
      ftna <- is.na(ft)
      ft[ftna] <- fto[ftna]

      tfpred <- rep(TRUE, nrow(x2))
      for (mn in marknames) {
        if (is.na(mn_tf[mn])) next
        if (mn_tf[mn]==1) {
          if (!is.na(ft[mn])) tfpred <- tfpred & ft[mn]<=x2[,mn]
        } else if (mn_tf[mn]==0) {
          if (!is.na(ft[mn])) tfpred <- tfpred & ft[mn]>x2[,mn]
        }
      }
      tfactual <- cpop_fname_actual[[cpop]][[fname]]
      
      cbind(data.frame(
        method="flowLearn",
        data_set=dset, scatter_plot=scat, cell_population=cpop, 
        train_samples=k, fcs=fname, 
        train=any(is.na(fts[[fname]]))
      ), f1score(tfactual, tfpred))
    })
    scores_dir <- paste0(gsub("results/2D/flowLearn_thresholds", "scores/2D/flowLearn",fl_dir),"/",cpop)
    dir.create(scores_dir, recursive=TRUE, showWarnings=FALSE)
    write.table(scoredf_fname_, file=gzfile(paste0(scores_dir,"/",k,".csv.gz")), sep=",", row.names=FALSE, col.names=TRUE)
  }
  time_output(start1, "scored")
}) }, .parallel=TRUE)
time_output(start)




