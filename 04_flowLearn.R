# date created: 2020-01-04
# author: alice yue
# input: 2D csv/clr + thresholds
# output: flowLearn predicted thresholds, plots, scores


## set directory, load packages, set parallel ####
no_cores <- 1#parallel::detectCores() - 5
# root <- "/home/ayue/projects/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## output ####
thres_dirs <- list_leaf_dirs(thres_dir)
gs_xr_ <- function(x,y) gs_xr(x,y,"results") 
plyr::l_ply(append(gs_xr_(thres_dirs,"flowLearn_plots"),
                   gs_xr_(thres_dirs,"flowLearn_thresholds")), 
            dir.create, recursive=TRUE, showWarnings=FALSE)


## parameters ####
dfn <- 512 # number of density features for flowlearn
ks <- c(1:5,10,15,20) # number of training samples

## plot function: scatterplot + densities ####
flPlot <- function(x2, marknames, ft, fto, filt, main) {
  # layout, margins, new plot
  layout(matrix(c(2,2,1,4,4,3,4,4,3), 3,3, byrow = TRUE))
  par(mar=c(1,3,3,3)) # bottom, left, top, right
  plot.new()
  
  # top density plot
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
  traintitle <- ifelse((is.na(ft[colnames(x2)[colnames(x2)[1]]]) & 
                         colnames(x2)[1]%in%names(ft)) &
    (is.na(ft[colnames(x2)[colnames(x2)[2]]]) & 
       colnames(x2)[2]%in%names(ft)),
  " (train sample)","")
  plot(dens$y, dens$x, ylab=NA, yaxt="n",
       type="l", xlab="Density", 
       main=traintitle)
  if (colnames(x2)[1]%in%marknames) {
    abline(h=fto[colnames(x2)[2]], col="red")
    if (!is.na(ft[colnames(x2)[2]])) 
      abline(h=ft[colnames(x2)[2]], col="blue")
  }
  par(mar=c(4,4,1,1))
  plot_dens(
    x2, main=NA, xlab=colnames(x2)[1],
    ylab=colnames(x2)[2])
  if (colnames(x2)[1]%in%marknames) {
    abline(v=fto[colnames(x2)[1]], col="red")
    if (!is.na(ft[colnames(x2)[1]])) 
      abline(v=ft[colnames(x2)[1]], col="blue")
  }
  cpops <- names(filt)
  colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")[seq_len(length(cpops))]
  names(colours) <- cpops
  for (cpop in cpops)
    lines(filt[[cpop]], lwd=2, col=colours[cpop])
  legend("topright", legend=cpops, col=colours, 
         lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
  
  if (colnames(x2)[1]%in%marknames) {
    abline(h=fto[colnames(x2)[2]], col="red")
    if (!is.na(ft[colnames(x2)[2]])) 
      abline(h=ft[colnames(x2)[2]], col="blue")
  }
  
}


## START ####
start <- Sys.time()
par_scat <- FALSE
overwrite_thresholds <- TRUE
overwrite_plot <- TRUE


# res <- plyr::llply(thres_dirs, function(thres_dir_) { try({
  # CD66CD14_Livecells_
  for (thres_dir_ in thres_dirs[-c(1:3)]) { #try({
  thres_dir_s <- stringr::str_split(thres_dir_,"/")[[1]]
  scat <- thres_dir_s[length(thres_dir_s)]
  dset <- thres_dir_s[length(thres_dir_s)-1]
  
  x2_diri <- paste0(x2_dir,"/",dset,"/",scat)
  fl_dir <- gs_xr_(x2_diri, "flowLearn_thresholds")
  scores_dir_ <- paste0(gs_xr(fl_dir, "flowLearn","scores"),".csv.gz")
  # if (file.exists(scores_dir_)) return(NULL)
  
  start1 <- Sys.time()
  print(paste0(dset," > ", scat))
  
  t2_files <- list.files(thres_dir_, full.names=TRUE)
  fnames <- sapply(t2_files, function(t2_f)  gsub(".Rdata","",file_name(t2_f)))
  
  # get x, thresholds
  ftos <- purrr::map(t2_files, function(t2_file) get(load(t2_file)))
  x2s <- purrr::map(paste0(x2_diri,"/",fnames,".csv.gz"), 
                    data.table::fread, data.table=FALSE)
  y2_diri <- gsub("/x/","/y/",x2_diri)
  y2s <- purrr::map(paste0(y2_diri,"/",fnames,".csv.gz"), 
                    data.table::fread, data.table=FALSE)
  filt2_diri <- gsub("/x/","/filters/",x2_diri)
  filt2s <- purrr::map(paste0(filt2_diri,"/",fnames,".Rdata"),
                       function(x) get(load(x)))
  names(y2s) <- names(x2s) <- names(ftos) <- names(filt2s) <- fnames
  marknames <- names(ftos[[1]])
  
  time_output(start1, "loaded files")
  
  
  ## predict ####
  
  # for each threshold, predict for all files
  ## TEMP
  if (overwrite_thresholds) {
    protoIdxs <- NULL
    for (markname in marknames) {
      dt_dir <- paste0(fl_dir,"/",markname)
      if (!overwrite_thresholds & 
          all(sapply(paste0(dt_dir,"/",ks,".Rdata"), file.exists))) next
      
      # make a density data object for flowLearn
      densdat <- new('DensityData')
      for (fname in fnames) {
        if (is.na(ftos[[fname]][markname])) next 
        df <- flowLearn::flEstimateDensity(x2s[[fname]][,markname], dfn)
        densdat <- flowLearn::flAdd(
          densdat, fname, "cpop", 1, df$x, df$y, NaN, ftos[[fname]][markname])
      }
      
      # for each number of training data, predict threshold
      if (is.null(protoIdxs))
        protoIdxs <- plyr::llply(ks, function(k) 
          which(rownames(densdat@data)%in%flowLearn::flSelectPrototypes(densdat, k)))
      
      dir.create(dt_dir, showWarnings=FALSE, recursive=TRUE)
      for (ki in seq_len(length(ks))) {
        k <- ks[ki]
        if (k>=nrow(densdat@data)) break
        protoIdx <- protoIdxs[[ki]]
        
        ddp <- flowLearn::flPredictThresholds(densdat, protoIdx)
        
        ddpt <- ddp@data$gate.high
        ddpt[rownames(densdat@data)%in%protoIdx] <- NA
        names(ddpt) <- densdat@data$fcs
        
        save(ddpt, file=paste0(dt_dir,"/",k,".Rdata"))
      }
    }
    time_output(start1, "predicted thresholds")
  }
  
  
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
      if (round(ftos[[1]][mn], digits=6) <= 
          round(min(cpopx[,mn]), digits=6)) {
        mn_tf[mn] <- 1 # lower threshold
      } else if (round(ftos[[1]][mn], digits=6) >= 
                 round(min(cpopx[,mn]), digits=6)) {
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
      na_other <- cpop=="other"
      for (mn in marknames) {
        if (is.na(mn_tf[mn]) | is.na(fto[mn])) next
        if (mn_tf[mn]==1) {
          tfactual <- tfactual & fto[mn]<=x2[,mn]
        } else if (mn_tf[mn]==0) {
          tfactual <- tfactual & fto[mn]>x2[,mn]
        }
        na_other <- FALSE
      }
      if (na_other) return(NA)
      return(tfactual)
    })
    names(fname_actual) <- fnames
    return(fname_actual)
  })
  names(cpop_fname_actual) <- cpops
  if (cpops[length(cpops)]=="other")
    for (fi in seq_len(length(cpop_fname_actual$other)))
      if (is.na(cpop_fname_actual$other[[fi]])) 
        cpop_fname_actual$other[[fi]] <- 
    !Reduce("|",lapply(cpops[!cpops%in%"other"], function(ci)
      cpop_fname_actual[[ci]][[fi]]))
  
  
  # for each number of train data used; cell population
  ks_ <- as.numeric(gsub(".Rdata","",list.files(
    paste0(fl_dir,"/",names(ftos[[1]])[1]))))
  
  # load prediced thresholds
  fl_dir <- gs_xr_(x2_diri,"flowLearn_thresholds")
  fts <- lapply(names(ftos[[1]]), function(x) {
    a <- lapply(ks_, function(k) get(load(paste0(fl_dir,"/",x,"/",k,".Rdata"))))
    names(a) <- as.character(ks_)
    a
  })
  names(fts) <- names(ftos[[1]])
  
  # get predicted cell population T/F's
  testpars <- expand.grid(cpops, ks_, fnames, stringsAsFactors=FALSE)
  colnames(testpars) <- c("cpop", "train_no", "fcs")
  
  # for each set of parameters @ index i, get threshold
  fts_ <- plyr::llply(seq_len(nrow(testpars)), function(i) {
    k <- testpars[i,"train_no"]
    fname <- testpars[i,"fcs"]
    ft <- sapply(fts, function(x) x[[as.character(k)]][fname])
    names(ft) <- marknames
    ft
  })
  tfpreds <- plyr::llply(seq_len(nrow(testpars)), function(i) {
    cpop <- testpars[i,"cpop"]
    fname <- testpars[i,"fcs"]
    
    ft <- fts_[[i]]
    x2 <- x2s[[fname]]
    fto <- ftos[[fname]]
    ftna <- is.na(ft)
    ft[ftna] <- fto[ftna]
    
    tfpred <- rep(TRUE, nrow(x2))
    for (mn in marknames) {
      if (is.na(mn_tfs[[cpop]][mn])) next
      if (mn_tfs[[cpop]][mn]==1) {
        if (!is.na(ft[mn])) tfpred <- tfpred & ft[mn]<=x2[,mn]
      } else if (mn_tfs[[cpop]][mn]==0) {
        if (!is.na(ft[mn])) tfpred <- tfpred & ft[mn]>x2[,mn]
      }
    }
    return(tfpred)
  })
  for (k in ks_) {
    fly_dir <- paste0(gsub("thresholds","y",fl_dir),"/",k)
    dir.create(fly_dir, recursive=TRUE, showWarnings=FALSE)
    for (fname in fnames) {
      tpi <- testpars[,"fcs"]==fname & testpars[,"train_no"]==k
      tpm <- Reduce(cbind,tfpreds[tpi])
      if (is.na(dim(tpm))) tpm <- matrix(tpm, ncol=1)
      colnames(tpm) <- testpars[tpi,"cpop"]
      write.table(tpm, file=gzfile(paste0(fly_dir,"/",fname,".csv.gz")), 
                  sep=",", row.names=FALSE, col.names=TRUE)
    }
  }
  time_output(start1, "saved clrs")
  
  scoredf_k_cpop <- plyr::ldply(seq_len(nrow(testpars)), function(i) {
    ft <- fts_[[i]]
    tfactual <- cpop_fname_actual[[testpars[i,"cpop"]]][[testpars[i,"fcs"]]]
    cbind(data.frame(train=is.na(ft)), f1score(tfactual, tfpreds[[i]]))
  })
  scoredf_k_cpop <- cbind("flowLearn", dset, scat, testpars, scoredf_k_cpop)
  colnames(scoredf_k_cpop)[c(1:3)] <- c("method", "dataset", "scatterplot")
  # save(scoredf_k_cpop, file=scores_dir_)
  dir.create(folder_name(scores_dir_), recursive=TRUE, showWarnings=FALSE)
  write.table(scoredf_k_cpop, file=gzfile(scores_dir_), 
              sep=",", row.names=FALSE, col.names=TRUE)
  time_output(start1, "scored")
  
  
  ## plot ####
  pl_dir <- gsub("thresholds","plots",fl_dir)
  a <- sapply(paste0(pl_dir,"/",ks_), dir.create, showWarnings=FALSE, recursive=TRUE)
  for (fname in fnames) {
    for (k in ks_) {
      png_name <- paste0(pl_dir,"/",k,"/",fname,".png")
      if (!overwrite_plot & file.exists(png_name)) next
      x2 <- x2s[[fname]]
      fto <- ftos[[fname]]
      ft <- sapply(fts, function(x) x[[as.character(k)]][fname]); names(ft) <- marknames
      filt <- filt2s[[fname]]
      png(png_name, width=400, height=400)
      flPlot(x2, marknames, ft, fto, filt, main=paste0("data set: ", dset, "\nscatterplot: ",scat, "\n(blue=predicted, red=actual)"))
      graphics.off()
    }
  }
  time_output(start1, "plotted")
}#)}#, .parallel=par_scat)
time_output(start)


