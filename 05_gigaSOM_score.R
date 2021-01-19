# date created: 2020-01-04
# author: alice yue
# input: GigaSOM output
# output: score matrix (f1)


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 13#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
setwd(root)


## packages ####
source("src/helpers.R")
source("src/helpers_2D.R")
libr(c(
  "combinat", 
  "flowCore",
  "furrr" #"rslurm",
))


## input ####
y2_dir <- paste0(root,"/data/2D/y"); 
gs2_dir <- paste0(root,"/results/2D/GigaSOM_clusters"); 
gsn_dir <- paste0(root,"/results/nD/GigaSOM_clusters"); 


## output ####
gss_dir <- paste0(root,"/2D/GigaSOM")
dir.create(gss_dir, recursive=TRUE, showWarnings=FALSE)


## load inputs ####
gs2_files <- list.files(gs2_dir, recursive=TRUE, pattern=".csv.gz", full.names=TRUE)
gsn_files <- list.files(gsn_dir, recursive=TRUE, pattern=".csv.gz", full.names=TRUE)



gs2_files <- sample(gs2_files)

start <- Sys.time()

# cpop combos
clustn <- 4
cpop_combos_all <- lapply(2:4, function(cpopn) {
  # number of clusters in each cpop to try
  maxcl <- clustn - cpopn + 1
  maxcll <- lapply(seq_len(cpopn), function(x) seq_len(maxcl))
  cpop_size <- expand.grid(maxcll) 
  cpop_size <- cpop_size[rowSums(cpop_size)==clustn,,drop=FALSE]
  
  # all permutations of clusters
  cpop_order <- combinat::permn(seq_len(clustn))
  
  # get all number and perm of clusters
  cpop_combo <- lapply(seq_len(nrow(cpop_size)), function(csi) {
    cs <- cpop_size[csi,]
    lapply(cpop_order, function(co) {
      a <- list()
      clusi <- 1
      for (i in seq_len(length(cs))) {
        csi <- cs[i]
        clusj <- as.numeric(clusi+csi)
        a[[i]] <-  sort(co[clusi:(clusj-1)])
        clusi <- clusj
      }
      return(a)
    })
  })
  cpop_combo <- unlist(cpop_combo, recursive=FALSE)
  cpop_combo <- cpop_combo[!duplicated(cpop_combo)]
  return(cpop_combo)
})
names(cpop_combos_all) <- as.character(2:4)


loop_ind <- loop_ind_f(gs2_files, no_cores)
scoredf2_ <- furrr::future_map(loop_ind, function(gs_files) { purrr::map(gs_files, function(gs_file) { try({
  score_dir <- gsub("results/2D/GigaSOM_clusters", "scores/2D/GigaSOM", gs_file)
  if (file.exists(score_dir))
    if (file.size(score_dir)>0)
      return(NULL)
  
  predicted <- as.data.frame(data.table::fread(gs_file, data.table=FALSE))
  y <- as.data.frame(data.table::fread(gsub("results/2D/GigaSOM_clusters","data/2D/y",gs_file), data.table=FALSE))
  x <- as.data.frame(data.table::fread(gsub("results/2D/GigaSOM_clusters","data/2D/x",gs_file), data.table=FALSE))
  
  # # get rid of other cell population, these aren't important
  # good_ind <- rep(TRUE, nrow(y_))
  # if (any(colnames(y_)%in%"other"))
  #   if (ncol(y_)>2)
  #     good_ind <- y_[,"other"]==0
  # y <- y_[good_ind,!colnames(y_)%in%"other",drop=FALSE]
  # x <- x_[good_ind,,drop=FALSE]
  # predicted <- predicted_[good_ind,,drop=FALSE]
  
  # get meta data
  gs_f <- stringr::str_split(gs_file,"/")[[1]]
  fname <- gsub(".csv.gz","",gs_f[length(gs_f)])
  scat <- gs_f[length(gs_f)-1]
  dataset <- gs_f[length(gs_f)-2]
  
  cpops <- colnames(y)
  cpopsn <- length(cpops)
  clustn <- max(predicted)
  
  
  # for each cpop combo
  cpop_combos <- cpop_combos_all[[as.character(cpopsn)]]
  scoredf_cpop_combos <- purrr::map(cpop_combos, function(cpop_combo) {
    
    # for each cell population
    scoredf_cpop_ <- purrr::map_dfr(seq_len(cpopsn), function(cpopi) {
      cpop <- cpops[cpopi]
      tfactual <- y[,cpop]==1
      tfpred <- predicted$index%in%cpop_combo[[cpopi]]
      
      cbind(data.frame(
        method="gigaSOM",
        data_set=dataset, scatter_plot=scat, cell_population=cpop, 
        train_samples=NA, total_samples=NA, fcs=fname, 
        train=NA
      ), f1score(tfactual, tfpred))
    })
  })
  if ("other"%in%cpops & cpopsn>2) {
    cpop_combos <- cpop_combos_all[[as.character(cpopsn-1)]]
    scoredf_cpop_combos <- append(
      scoredf_cpop_combos, purrr::map(cpop_combos, function(cpop_combo) {
      
      # for each cell population
      scoredf_cpop_ <- purrr::map_dfr(seq_len(cpopsn-1), function(cpopi) {
        cpop <- cpops[cpopi]
        tfactual <- y[,cpop]==1
        tfpred <- predicted$index%in%cpop_combo[[cpopi]]
        
        cbind(data.frame(
          method="gigaSOM",
          data_set=dataset, scatter_plot=scat, cell_population=cpop, 
          train_samples=NA, total_samples=NA, fcs=fname, 
          train=NA
        ), f1score(tfactual, tfpred))
      })
      scoredf_cpop_ <- rbind(scoredf_cpop_, rep(0,ncol(scoredf_cpop_)))
      scoredf_cpop_$cell_population[nrow(scoredf_cpop_)] <- "other"
      return(scoredf_cpop_)
    }))
  }
  
  # get best cpop combo
  meanf1s <- sapply(scoredf_cpop_combos, function(x) mean(as.numeric(x[,"f1"])))
  best_combo <- which.max(meanf1s)
  best <- scoredf_cpop_combos[[best_combo]]
  best <- best[best$cell_population!="other",,drop=FALSE]
  
  score_directory <- stringr::str_split(score_dir,"/")[[1]]
  score_directory <- paste0(score_directory[-length(score_directory)], collapse="/")
  dir.create(gsub(paste0(root,"/"), "", score_directory), recursive=TRUE, showWarnings=FALSE)
  write.table(best, file=gzfile(gsub(paste0(root,"/"), "", score_dir)), sep=',', row.names=FALSE, col.names=TRUE)
  
  ## plot!
  try ({
    tfs <- purrr::map(cpops, function(cpop) 
      predicted$index%in%cpop_combos[[best_combo]][cpop] )
    names(tfs) <- cpops
    
    gs_file_ <- gsub("clusters","plots",gs_file)
    gs_file_ <- stringr::str_split(gs_file_, "/")[[1]]
    gs_file_ <- paste0(gs_file_[-length(gs_file_)], collapse="/")
    dir.create(gs_file_, recursive=TRUE, showWarnings=FALSE)
    png(paste(sub(".csv.gz",".png",gsub("clusters","plots",gs_file),".png")), 
        width=800, height=450)
    par(mfrow=c(1,2))
    plot_dens(x, xlab="tsne 1", ylab="tsne 2", main="actual")
    colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")
    colours <- colours[seq_len(length(cpops))]
    names(colours) <- cpops
    for (cpop in cpops) {
      xcpop <- as.matrix(x[y[,cpop]==1,,drop=FALSE])
      if (nrow(xcpop)==0) next
      ch <- chull(xcpop)#$rescoords
      ch <- xcpop[c(ch,ch[1]),]
      lines(ch, lwd=2, col=colours[cpop])
    }
    legend("topright", legend=cpops, col=colours, 
           lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
    clust_vec <- rep(NA, nrow(predicted))
    for (cpopi in seq_len(length(cpops)))
      for (clust in cpop_combos[[best_combo]][cpopi]) 
        clust_vec[predicted[,1]%in%clust] <- cpops[cpopi]
    cols <- colours[clust_vec]
    # cols[is.na(cols) | !good_ind] <- "black"
    plot(x, cex=.1, col=cols, xlab="tsne 1", ylab="tsne 2", main="clustered")
    legend("topright", legend=cpops, col=colours, 
           lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
    graphics.off()
  })
  # return(dplyr::bind_rows(scoredf_cpop_))
  # return(best)
  
}) }) })
# scoredf2_ <- unlist(scoredf2_)
# names(scoredf2_) <- gs2_files
# scoredf2 <- plyr::compact(scoredf2_)
# gs2_files_error <- gs2_files[!gs2_files%in%names(scoredf2)]
# 
# scoredf2tb <- Reduce(rbind, scoredf2)
# dir.create(paste0(score_dir,"/2D"), showWarnings=FALSE)
# write.table(scoredf2tb, file=gzfile(paste0(score_dir,"/2D/gigaSOM.csv.gz")), sep=",", row.names=FALSE)
time_output(start)


start <- Sys.time()

scoredf <- furrr::future_map(gsn_files, function(gs_file) {
  # get_inds <- sample(seq_len(nrow(predicted)), 400)
  predicted_ <- as.data.frame(data.table::fread(gs_file, data.table=FALSE))
  y_ <- as.data.frame(data.table::fread(gsub("results/nD/GigaSOM_clusters","data/nD/y",gs_file), data.table=FALSE))
  x_ <- as.data.frame(data.table::fread(gsub("results/nD/GigaSOM_clusters","data/nD/x",gs_file), data.table=FALSE))
  
  # get rid of other cell population, these aren't important
  good_ind <- rep(TRUE, nrow(y_))
  if (any(colnames(y_)%in%"other"))
    if (ncol(y_)>2)
      good_ind <- y_[,"other"]==0
  y <- y_[good_ind,!colnames(y_)%in%"other",drop=FALSE]
  x <- x_[good_ind,,drop=FALSE]
  predicted <- predicted_[good_ind,,drop=FALSE]
  
  # get meta data
  gs_f <- stringr::str_split(gs_file,"/")[[1]]
  fname <- gsub(".csv.gz","",gs_f[length(gs_f)])
  scat <- gs_f[length(gs_f)-1]
  dataset <- gs_f[length(gs_f)-2]
  
  cpops <- colnames(y)
  
  cpopsn <- length(cpops)
  clustn <- max(predicted[,1])
  
  # number of clusters in each cpop to try
  maxcl <- clustn - cpopsn + 1
  maxcll <- lapply(seq_len(cpopsn), function(x) seq_len(maxcl))
  cpop_size <- expand.grid(maxcll) 
  cpop_size <- cpop_size[rowSums(cpop_size)==clustn,,drop=FALSE]
  
  # all permutations of clusters
  cpop_order <- combinat::permn(seq_len(clustn))
  
  # get all number and perm of clusters
  cpop_combos <- purrr::map(seq_len(nrow(cpop_size)), function(csi) {
    cs <- cpop_size[csi,]
    lapply(cpop_order, function(co) {
      a <- list()
      clusi <- 1
      for (i in seq_len(length(cs))) {
        csi <- cs[i]
        clusj <- as.numeric(clusi+csi)
        a[[cpops[i]]] <-  sort(co[clusi:(clusj-1)])
        clusi <- clusj
      }
      return(a)
    })
  })
  cpop_combos <- unlist(cpop_combos, recursive=FALSE)
  cpop_combos <- cpop_combos[!duplicated(cpop_combos)]
  
  # for each cpop combo
  scoredf_cpop_combos <- furrr::future_map(cpop_combos, function(cpop_combo) {
    
    # for each cell population
    scoredf_cpop_ <- purrr::map_dfr(cpops, function(cpop) {
      tfactual <- y[,cpop]==1
      tfpred <- predicted$index%in%cpop_combo[cpop]
      
      cbind(data.frame(
        method="gigaSOM",
        data_set=dataset, scatter_plot=NA, cell_population=cpop, 
        train_samples=NA, total_samples=NA, fcs=fname, 
        train=NA
      ), f1score(tfactual, tfpred))
    })
  })
  
  # get best cpop combo
  meanf1s <- sapply(scoredf_cpop_combos, function(x) mean(as.numeric(x$f1)))
  best_combo <- which.max(meanf1s)
  best <- scoredf_cpop_combos[[best_combo]]
  
  score_dir <- gsub("results/nD/GigaSOM_clusters", "scores/nD/GigaSOM", gs_file)
  score_directory <- stringr::str_split(score_dir,"/")[[1]]
  score_directory <- paste0(score_directory[-length(score_directory)], collapse="/")
  dir.create(score_directory, recursive=TRUE, showWarnings=FALSE)
  write.table(best, file=gzfile(score_dir), sep=",", row.names=FALSE, col.names=TRUE)
  
  try({
    tx <- as.matrix(data.table::fread(gsub("GigaSOM_clusters","Rtsne",gs_file), data.table=FALSE))
    
    png_file <- paste(sub(".csv.gz",".png",gsub("clusters","plots",gs_file),".png"))
    dir.create(gsub(file_name(png_file),"",png_file), recursive=TRUE, showWarnings=FALSE)
    png(png_file, 
        width=900, height=400)
    par(mfrow=c(1,2))
    colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")
    colours <- colours[seq_len(length(cpops))]
    names(colours) <- cpops
    clust_vec <- rep(NA, nrow(y_))
    for (cpop in cpops) 
      clust_vec[y_[,cpop]==1] <- cpop
    cols <- colours[clust_vec]
    cols[is.na(cols) | !good_ind] <- "black"
    plot(tx, cex=.1, xlab="tsne 1", ylab="tsne 2", main="actual", col=cols)
    legend("topright", legend=cpops, col=colours, 
           lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
    clust_vec <- rep(NA, nrow(predicted_))
    for (cpop in cpops) 
      for (clust in cpop_combos[[best_combo]][cpop]) 
        clust_vec[predicted_[,1]%in%clust] <- cpop
    cols <- colours[clust_vec]
    cols[is.na(cols) | !good_ind] <- "black"
    plot(tx, cex=.1, xlab="tsne 1", ylab="tsne 2", main="predicted", col=cols)
    legend("topright", legend=cpops, col=colours, 
           lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
    graphics.off()
  })
  
  # return(best)
  # return(dplyr::bind_rows(scoredf_cpop_))
})
# dir.create(paste0(score_dir,"/nD"), showWarnings=FALSE)
# write.table(scoredf, file=gzfile(paste0(score_dir,"/nD/gigaSOM_nD.csv.gz")), sep=",", row.names=FALSE)
time_output(start)




