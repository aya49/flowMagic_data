# date created: 2020-01-04
# author: alice yue
# input: GigaSOM output
# output: score matrix (f1)


## set directory, load packages, set parallel ####
no_cores <- 14#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
gs2_dir <- paste0(root,"/results/2D/GigaSOM_clusters"); 
gsn_dir <- paste0(root,"/results/nD/GigaSOM_clusters"); 
gl2_dir <- paste0(root,"/results/2D/GigaSOM_labels"); 
gln_dir <- paste0(root,"/results/nD/GigaSOM_labels"); 


## load inputs ####
gs2_dirs <- list_leaf_dirs(gs2_dir)
gsn_dirs <- list_leaf_dirs(gsn_dir)


## output ####
plyr::l_ply(append(gs2_dirs, gsn_dirs), function(x) {
  dir.create(gsub("_clusters","",gsub("results","scores",x)), 
             recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("_clusters","_plots",x), recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("_clusters","_labels",x), recursive=TRUE, showWarnings=FALSE)
})


## scoring function ####
clust_score <- function(cpop_combos, clusttf, actualtf) {
  purrr::map(cpop_combos, function(cpop_combo) {
    # for each cell population
    purrr::map_dfr(seq_len(length(cpops)), function(cpopi) {
      if (cpopi>length(cpop_combo)) return(NULL)
      cpop <- cpops[cpopi]
      tfactual <- actualtf[[cpopi]]
      tfpred <- Reduce('|',clusttf[cpop_combo[[cpopi]]])
      
      f1score(tfactual, tfpred)
    })
  })
}


## cpop combo function ####
# all combination of clusters vs cpops (i.e. which clusters are in which cpops)
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
clustn <- 4
cpop_combos_all_2D <- lapply(2:4, function(cpopn) 
  get_cpop_combos(clustn, cpopn) )
names(cpop_combos_all_2D) <- as.character(2:4)

start <- Sys.time()

overwrite <- TRUE
par_scat <- TRUE# parallelize by scatterplot, not by file

# loop_ind <- loop_ind_f(sample(append(gs2_dirs, gsn_dirs)), no_cores)
res <- furrr::future_map(append(gsn_dirs, gs2_dirs), function(gs_dir_) {
# for (gs_dir_ in gs2_dirs) {
  start1 <- Sys.time()
  gs_files <- list.files(gs_dir_, full.names=TRUE, pattern=".csv.gz")
  
  # get meta data
  gs_f_ <- stringr::str_split(gs_dir_,"/")[[1]]
  scat <- gs_f_[length(gs_f_)]
  dset <- gs_f_[length(gs_f_)-1]
  
  nD <- FALSE
  if (dset=="GigaSOM_clusters") {
    scat <- NA
    dset <- scat
    nD <- TRUE
    
    # make cpop_combos
    cc_file <- paste0(gsub("_clusters","_combos",gs_dir_),".Rdata")
    if (overwrite | !file.exists(cc_file)) {
      predicted <- as.data.frame(data.table::fread(gs_files[1], data.table=FALSE))
      y_f <- gsub("results", "data", gsub("GigaSOM_clusters","y",gs_files[1]))
      y <- as.data.frame(data.table::fread(y_f, data.table=FALSE))
      cpops <- colnames(y)
      cpopn <- ifelse("other"%in%cpops, ncol(y)-1, ncol(y))
      clustn <- max(predicted)
      
      # get cpop combos
      cpop_combos <- get_cpop_combos(clustn, cpopn)
      dir.create(folder_name(cc_file), recursive=TRUE, showWarnings=FALSE)
      save(cpop_combos, file=cc_file)
    }
    cpop_combos <- get(load(cc_file))
  }
  cat("\n",dset,">",scat)

  loop_ind <- loop_ind_f(gs_files, no_cores)
  bests <- plyr::ldply(loop_ind, function(x) plyr::ldply(x, function(gs_f) {
    # for (gs_f in gs_files) {
    # if (file.exists(gsub("_clusters","_labels",gs_f)))
    #   if (file.size(gsub("_clusters","_labels",gs_f))>0)
    #     return(NULL)
    
    score_file <- gsub("results", "scores", gsub("_clusters","",gs_f))
    
    predicted <- as.data.frame(data.table::fread(gs_f, data.table=FALSE))[,1]
    y_f <- gsub("results", "data", gsub("GigaSOM_clusters","y",gs_f))
    y <- as.data.frame(data.table::fread(y_f, data.table=FALSE))
    if (nD & "other"%in%colnames(y)) y <- y[,colnames(y)!="other",drop=FALSE]
    x_f <- gsub("results", "data", gsub("GigaSOM_clusters","x",gs_f))
    x <- as.data.frame(data.table::fread(x_f, data.table=FALSE))
    
    fname <- gsub(".csv.gz","",gs_f_[length(gs_f_)])
    cpops <- colnames(y)
    cpopn <- length(cpops)
    clustn <- max(predicted)
    
    clusttf <- plyr::llply(seq_len(clustn), function(x) predicted==x)
    actualtf <- plyr::llply(seq_len(cpopn), function(x) y[,x]==1)
    names(actualtf) <- cpops
    
    # get cluster combinations
    if (!nD) cpop_combos <- cpop_combos_all_2D[[as.character(cpopn)]]
    if (!nD & "other"%in%cpops & cpopn>2) 
      cpop_combos <- append(cpop_combos, cpop_combos_all_2D[[as.character(cpopn-1)]])

    # for each cpop combo
    if (!nD & clustn>4) cpop_combos <- get_cpop_combos(clustn, cpopn)
    scoredf_cpop_combos <- clust_score(cpop_combos, clusttf, actualtf)

    # get best cpop combo
    meanf1s <- sapply(scoredf_cpop_combos, function(x) mean(as.numeric(x[,"f1"])))
    best_combo <- which.max(meanf1s)
    best <- scoredf_cpop_combos[[best_combo]]
    best <- cbind(cpops[seq_len(nrow(best))], best)
    best <- best[best[,1]!="other",,drop=FALSE]
    
    # write.table(best, file=gzfile(score_file), sep=',', row.names=FALSE, col.names=TRUE)
    
    cpops_ <- cpops
    if (length(cpop_combos[[best_combo]])<length(cpops))
      cpops <- cpops[!cpops%in%"other"]
    tfs <- plyr::llply(seq_len(length(cpops)), function(cpopi) 
      Reduce('|',clusttf[cpop_combos[[best_combo]][[cpopi]]]) )
    if (length(tfs)==1) {
      tfs <- matrix(tfs, ncol=1)
    } else {
      tfs <- Reduce(cbind, tfs)
    }
    colnames(tfs) <- cpops
    if (length(cpops_)>length(cpops)) {
      tfs <- cbind(tfs, rep(FALSE,nrow(tfs)))
      colnames(tfs)[ncol(tfs)] <- "other"
    }
    write.table(tfs, file=gzfile(gsub("_clusters","_labels",gs_f)),
                sep=',', row.names=FALSE, col.names=TRUE)
    # save(tfs, file=gsub(".csv.gz",".Rdata",gsub("results/_clusters","",gs_f)))
    
    return(best)
  }), .parallel=!par_scat)
  bests <- bests[,colnames(bests)!=".id"]
  bests <- cbind("gigaSOM", dset, scat, bests[,1], 0, 
                 gsub(".csv.gz","",sapply(gs_files, file_name)), FALSE, bests[,-1])
  colnames(bests)[c(1:7)] <- c("method","dataset","scatterplot","cpop","train_no","fcs","train")
  # save(bests, file=paste0(gsub("results","scores",gsub("_clusters","",gs_dir_)),".Rdata"))
  write.table(bests, file=gzfile(gsub("results","scores",gsub("_clusters","",gs_dir_)),".csv.gz"), 
              sep=",", row.names=FALSE, col.names=TRUE)
  time_output(start1)
})
time_output(start)


# ## plot! ####
# start <- Sys.time()
# 
# gl2_files <- list.files(gl2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")
# gln_files <- list.files(gln_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")
# 
# loop_ind <- loop_ind_f(sample(append(gl2_files, gln_files)), no_cores)
# plyr::l_ply(loop_ind, function(gl_fs) { plyr::l_ply(gl_fs, function(gl_f) { try({
#   nD <- grepl("/nD/", gl_f)
# 
#   tfs <- as.data.frame(data.table::fread(gl_f, data.table=FALSE))
#   y_f <- gsub("results", "data", gsub("GigaSOM_labels","y",gs_f))
#   y <- as.data.frame(data.table::fread(y_f, data.table=FALSE))
#   cpops <- colnames(y)
# 
#   ft_ <- NULL
#   if (nD) {
#     tx <- as.matrix(data.table::fread(gsub("GigaSOM_labels","Rtsne",gs_f), data.table=FALSE))
#   } else {
#     x_f <- gsub("results", "data", gsub("GigaSOM_labels","x",gl_f))
#     tx <- as.data.frame(data.table::fread(x_f, data.table=FALSE))
#     ft_ <- get(load(gsub(".csv.gz",".Rdata",gsub("results/2D/GigaSOM_labels","data/2D/filters",gs_f))))
#   }
# 
#   colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")[seq_len(length(cpops))]
#   names(colours) <- cpops
# 
#   png(gsub(".csv.gz",".png",gsub("labels","plots",gl_f)),width=800,height=450)
#   par(mfrow=c(1,2))
# 
#   plot_dens(tx, main="actual")
#   for (cpop in cpops) {
#     if (!nD & cpop%in%names(ft_)) lines(ft_[[cpop]], lwd=2, col=colours[cpop])
#     if (nD) points(tx[y[,cpop]==1,,drop=FALSE], cex=.1, col=colours[cpop])
#   }
#   legend("topright", legend=cpops, col=colours,
#          lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
# 
#   clust_vec <- rep(NA, nrow(x))
#   for (cpop in cpops)
#     clust_vec[tfs[,cpop]] <- cpop
#   cols <- colours[clust_vec]
#   plot(tx, xlab=colnames(tx)[1], ylab=colnames(tx)[2],
#        cex=.1, main="clustered", col=cols)
#   legend("topright", legend=cpops, col=colours,
#          lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
# 
#   graphics.off()
# }) }) }, .parallel=TRUE)
# time_output(start)




