# date created: 2020-01-04
# author: alice yue
# input: flowLearn
# output: score matrix


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
  "flowCore",
  "flowLearn" # devtools::install_github("mlux86/flowLearn")
))


## input ####
main_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
out_dir <- paste0(main_dir,"/data")
res_dir <- paste0(main_dir,"/results")
x2_dir <- paste0(out_dir,"/2D/x"); 
y2_dir <- paste0(out_dir,"/2D/y"); 
t2_dir <- paste0(out_dir,"/2D/thresholds"); 


## output ####
score_dir <- paste0(main_dir,"/scores")
dir.create(score_dir, recursive=TRUE, showWarnings=FALSE)


## load inputs ####
# actual thresholds
t2_dirs <- list.dirs(t2_dir)[-1] # includes t2_dir so remove it

# predicted thresholds
tp2_dirs <- list.dirs(paste0(res_dir,"/2D/flowLearn_thresholds"), recursive=TRUE)
tp2_dirs <- unique(tp2_dirs[sapply(tp2_dirs, function(x) sum(grepl(x,tp2_dirs))>1)])
tp2_dirs <- tp2_dirs[sapply(tp2_dirs, function(x) sum(grepl(x,tp2_dirs))==1)]

# # x data
# x2_folders <- list.dirs(x2_dir, recursive=TRUE)
# x2_folders <- x2_folders[sapply(x2_folders, function(x) sum(grepl(x,x2_folders))==1)]


## load predicted thresholds and get f1 ####

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

scoredf_ <- purrr::map(tp2_dirs, function(ft_dir) {
  x2_fold <- gsub("results/2D/flowLearn_thresholds","data/2D/x",ft_dir)
  marknames <- list.dirs(ft_dir, full.names=FALSE)[-1]
  scat <- file_name(ft_dir)
  fto_dir <- gsub(paste0("/",scat),"",gsub("/x/","/thresholds/",x2_fold))
  dataset <- file_name(fto_dir)
  
  # get x, y, threshold files ^
  print(paste0("data/scatterplot: ",dataset,"/",scat," (",paste0(marknames,collapse=", "),")"))
  
  fto_files <- list.files(fto_dir, full.names=TRUE)
  ftos <- purrr::map(fto_files, function(fto_file) get(load(fto_file))[[scat]])
  
  fnames <- gsub(".Rdata","",sapply(fto_files, file_name))
  
  x2_files <- list.files(x2_fold, full.names=TRUE)
  x2_files_fn <- gsub(".csv.gz","",sapply(x2_files, file_name))
  x2s <- purrr::map(fnames, function(fname) {
    x2_id <- which(x2_files_fn==fname); if (length(x2_id)==0) return(NULL)
    data.table::fread(x2_files[x2_id], data.table=FALSE) })
  names(x2s) <- fnames
  x2s <- plyr::compact(x2s)
  fnames <- names(x2s)
  
  y2_files <- list.files(gsub("/x/","/y/",x2_fold), full.names=TRUE)
  y2_files_fn <- gsub(".csv.gz","",sapply(y2_files, file_name))
  y2s <- purrr::map(fnames, function(fname) 
    data.table::fread(y2_files[y2_files_fn==fname], data.table=FALSE))
  
  names(x2s) <- names(y2s) <- names(ftos) <- fnames
  
  tempmark <- rep(NA,length(marknames))
  names(tempmark) <- marknames
  cpops <- colnames(y2s[[1]])

  # get how thresholds work with cpops
  mn_tfs <- purrr::map(cpops, function(cpop) {
    cpopx <- x2s[[1]][y2s[[1]][,cpop]==1,,drop=FALSE]

    mn_tf <- tempmark
    for (mn in marknames) {
      if (ftos[[1]][mn] <= min(cpopx[,mn])) {
        mn_tf[mn] <- 1 # lower threshold
      } else if (ftos[[1]][mn] >= max(cpopx[,mn])) {
        mn_tf[mn] <- 0 # higher threshold
      }
    }
    return(mn_tf)
  })
  names(mn_tfs) <- cpops
  
  # get actual cell populations not from y but from thresholds just to be safe
  cpop_fname_actual <- purrr::map(cpops, function(cpop) {
    # for each fcs file
    fname_actual <- purrr::map(fnames, function(fname) {
      x2 <- x2s[[fname]]
      fto <- ftos[[fname]]
      mn_tf <- mn_tfs[[cpop]]
      
      tfactual <- rep(TRUE, nrow(x2))
      for (mn in marknames) {
        if (is.na(mn_tf[mn])) next
        if (mn_tf[mn]==1) {
          if (!is.na(fto[mn])) tfactual <- tfactual & fto[mn]<=x2[,mn]
        } else if (mn_tf[mn]==0) {
          if (!is.na(fto[mn])) tfactual <- tfactual & fto[mn]>x2[,mn]
        }
      }
      return(tfactual)
    })
    names(fname_actual) <- fnames
    return(fname_actual)
  })
  names(cpop_fname_actual) <- cpops
  
  # for each number of train data used
  # (not very efficient, but whatever lol)
  ks_ <- as.numeric(gsub(".Rdata","",list.files(paste0(ft_dir,"/",marknames[1]), full.names=FALSE)))
  
  tempmark <- rep(NA,length(marknames))
  names(tempmark) <- marknames
  scoredf_k_ <- purrr::map_dfr(ks_, function(k) {
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
      mn_tf <- mn_tfs[[cpop]]

      # for each fcs file
      scoredf_fname_ <- purrr::map_dfr(fnames, function(fname) {
        x2 <- x2s[[fname]]
        fto <- ftos[[fname]]
        ft <- fts[[fname]]
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
          data_set=dataset, scatter_plot=scat, cell_population=cpop, 
          train_samples=k, total_samples=length(fnames), fcs=fname, 
          train=any(is.na(fts[[fname]]))
          # true_size=actu,
          # predicted_size=pred,
          # true_proportion=actu/length(tfactual), 
          # predicted_proportion=pred/length(tfpred)
        ), f1score(tfactual, tfpred))
      })
      # return(dplyr::bind_rows(scoredf_fname_))
    })
    # return(dplyr::bind_rows(scoredf_cpop_))
  })
  # return(dplyr:::bind_rows(scoredf_k_))
})
scoredf <- dplyr::bind_rows(scoredf_)
dir.create(paste0(score_dir,"/2D"), showWarnings=FALSE)
write.table(scoredf, file=gzfile(paste0(score_dir,"/2D/flowLearn.csv.gz")), sep=",")
time_output(start)

