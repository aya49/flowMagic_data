# date created: 2020-01-04
# author: alice yue
# input: GigaSOM output
# output: score matrix (f1)


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
gs2_dir <- paste0(root,"/results/2D/GigaSOM_clusters"); 
gsn_dir <- paste0(root,"/results/nD/GigaSOM_clusters"); 


## load inputs ####
gs2_dirs <- list_leaf_dirs(gs2_dir)
gsn_dirs <- list_leaf_dirs(gsn_dir)


## output ####
plyr::l_ply(append(gs2_dirs, gsn_dirs), function(x) {
  dir.create(gsub("_clusters","",gsub("results","scores",x)), 
             recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("_clusters","_plots",x), recursive=TRUE, showWarnings=FALSE)
})


## scoring function ####
clust_score <- function(cpop_combos, clusttf, actualtf, cpops, dset, scat, fname) {
  purrr::map(cpop_combos, function(cpop_combo) {
    # for each cell population
    purrr::map_dfr(seq_len(length(cpops)), function(cpopi) {
      cpop <- cpops[cpopi]
      tfactual <- actualtf[[cpopi]]
      tfpred <- Reduce('|',clusttf[cpop_combo[[cpopi]]])
      
      cbind(data.frame(
        method="gigaSOM",
        data_set=dset, scatter_plot=scat, cell_population=cpop, 
        train_samples=NA, fcs=fname, 
        train=NA
      ), f1score(tfactual, tfpred))
    })
  })
}


## cpop combo function ####
get_cpop_combos <- function(clustn, cpopn) {
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
}


## START ####
start <- Sys.time()

clustn <- 4
cpop_combos_all_2D <- lapply(2:4, function(cpopn) 
  get_cpop_combos(clustn, cpopn) )
names(cpop_combos_all_2D) <- as.character(2:4)

# loop_ind <- loop_ind_f(sample(append(gs2_files, gsn_files)), no_cores)
plyr::l_ply(append(gsn_dirs, gs2_dirs), function(gs_dir_) { 
  start1 <- Sys.time()
  gs_files <- list.files(gs_dir_, full.names=TRUE, pattern=".csv.gz")
  
  # get meta data
  gs_f_ <- stringr::str_split(gs_dir_,"/")[[1]]
  scat <- gs_f_[length(gs_f_)]
  dset <- gs_f_[length(gs_f_)-1]
  if (dset=="GigaSOM_clusters") {
    scat <- NA
    dset <- scat
    nD <- TRUE
    
    predicted <- as.data.frame(data.table::fread(gs_files[1], data.table=FALSE))
    y_f <- gsub("results", "data", gsub("GigaSOM_clusters","x",gs_files[1]))
    y <- as.data.frame(data.table::fread(y_f, data.table=FALSE))
    cpopsn <- ifelse("other"%in%cpops, ncol(y)-1, ncol(y))
    clustn <- max(predicted[,1])
    
    # get cpop combos
    cpop_combos <- get_cpop_combos(clustn, cpopn)
  }

  plyr::l_ply(gs_files, function(gs_f) { try({
    score_file <- gsub("results", "scores", gsub("_clusters","",gs_f))
    
    predicted <- as.data.frame(data.table::fread(gs_f, data.table=FALSE))
    y_f <- gsub("results", "data", gsub("GigaSOM_clusters","x",gs_f))
    y <- as.data.frame(data.table::fread(y_f, data.table=FALSE))
    if (nD & "other"%in%colnames(y)) y <- y[,colnames(y)!="other",drop=FALSE]
    x_f <- gsub("results", "data", gsub("GigaSOM_clusters","y",gs_f))
    x <- as.data.frame(data.table::fread(x_f, data.table=FALSE))
    
    fname <- gsub(".csv.gz","",gs_f_[length(gs_f_)])
    cpops <- colnames(y)
    cpopsn <- length(cpops)
    clustn <- max(predicted[,1])
    
    clusttf <- plyr::llply(seq_len(clustn), function(x) predicted[,1]==x)
    actualtf <- plyr::llply(seq_len(cpopsn), function(x) y[,x]==1)
    names(actualtf) <- cpops
    
    # get cluster combinations
    if (!nD) cpop_combos <- cpop_combos_all_2D[[as.character(cpopsn)]]

    # for each cpop combo
    scoredf_cpop_combos <- clust_score(
      cpop_combos, clusttf, actualtf, cpops, dset, scat, fname)
    if (!nD & "other"%in%cpops & cpopsn>2) 
      scoredf_cpop_combos <- append(scoredf_cpop_combos,clust_score(
        cpop_combos_all_2D[[as.character(cpopsn-1)]], 
        clusttf, actualtf, cpops, dset, scat, fname))
    
    # get best cpop combo
    meanf1s <- sapply(scoredf_cpop_combos, function(x) mean(as.numeric(x[,"f1"])))
    best_combo <- which.max(meanf1s)
    best <- scoredf_cpop_combos[[best_combo]]
    best <- best[best$cell_population!="other",,drop=FALSE]
    
    write.table(best, file=gzfile(score_file), sep=',', row.names=FALSE, col.names=TRUE)
    
    
    ## plot! ####
    tx <- x
    ft_ <- NULL
    if (nD) {
      tx <- as.matrix(data.table::fread(gsub("GigaSOM_clusters","Rtsne",gs_file), data.table=FALSE))
      ft_ <- get(load(gsub("results/2D/GigaSOM_clusters","data/2D/filters",gs_file)))
    }
    tfs <- plyr::llply(cpops, function(cpop) Reduce('|',clusttf[cpop_combo[[cpop]]]) )
    names(tfs) <- cpops
    
    colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")[seq_len(length(cpops))]
    names(colours) <- cpops
    
    png(gsub(".csv.gz",".png",gsub("clusters","plots",gs_file)),width=800,height=450)
    par(mfrow=c(1,2))
    
    plot_dens(tx, main="actual")
    for (cpop in cpops) {
      if (!nD & cpop%in%names(ft_)) lines(ft_[[cpop]], lwd=2, col=colours[cpop])
      if (nD) points(tx[actualtf[[cpop]],,drop=FALSE], cex=.1, col=colours[cpop])
    }
    legend("topright", legend=cpops, col=colours, 
           lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
    
    clust_vec <- rep(NA, nrow(predicted))
    for (cpop in cpops) 
      clust_vec[tfs[[cpop]]] <- cpop
    cols <- colours[clust_vec]
    plot(tx, xlab=colnames(tx)[1], ylab=colnames(tx)[2], 
         cex=.1, main="clustered", col=cols)
    legend("topright", legend=cpops, col=cols, 
           lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
    
    graphics.off()
  }) })
}, .parallel=TRUE)
time_output(start)




