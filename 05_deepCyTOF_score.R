# date created: 2020-01-04
# author: alice yue
# input: deepCyTOF
# output: score matrix (f1)


## set directory, load packages, set parallel ####
no_cores <- 20#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
dc2_dir <- paste0(root,"/results/2D/deepCyTOF_labels"); 
dcn_dir <- paste0(root,"/results/nD/deepCyTOF_labels"); 


## load inputs ####
dc2_files <- list.files(dc2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")
dcn_files <- list.files(dcn_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")


## output ####
plyr::l_ply(append(dc2_files, dcn_files), function(x) {
  dir.create(gsub("_labels","",gsub("results","scores",x)), 
             recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("_labels","_plots",x), recursive=TRUE, showWarnings=FALSE)
})


# get train samples
trs <- gsub(".csv.gz","",list.files(paste0(results_dir, "/x_2Ddensity_euclidean_pam/10")))


## START ####
start <- Sys.time()

# loop_ind <- loop_ind_f(dc2_files, no_cores)
l_ply(append(dcn_files, dc2_files), function(dc_file) { try({
  score_file <- gsub("results/", "scores/", gsub("_labels","",dc_file))
  
  predicted <- read.csv(dc_file)[,1]
  actual <- data.table::fread(paste0(gsub("/deepCyTOF_labels","/y",gsub(
    "/results/","/data/",dc_file)),".gz"), data.table=FALSE)
  cpops <- colnames(actual)
  
  # get meta data
  # pus <- sort(unique(predicted)) # unique labels
  nD <- grepl("/nD/", dc_file)
  path_stuff <- stringr::str_split(dc_file,"/")[[1]]
  di <- which(path_stuff=="deepCyTOF_labels")
  dset <- path_stuff[di+1]
  scat <- ifelse(nD,path_stuff[di+2],NA)
  fname <- gsub(".csv","",path_stuff[length(path_stuff)])
  
  
  ## score each cpop ####
  best <- plyr::ldply(seq_len(ncol(actual)), function(cpopi) { 
    cbind(data.frame(
      method="deepCyTOF",
      data_set=dset, scatter_plot=scat, cell_population=cpops[cpopi], 
      train_samples=10, fcs=fname, 
      train=fname%in%trs
    ), f1score(actual[,cpopi]==1, predicted==cpopi))
  })
  write.table(best, file=gzfile(paste0(score_file,".gz")), 
              sep=",", row.names=FALSE, col.names=TRUE)

  
  ## plot ####
  ft_ <- NULL
  if (nD) {
    tx <- data.table::fread(paste0(gsub(
      "/deepCyTOF_labels","/Rtsne",dc_file),".gz"), data.table=FALSE)
  } else {
    tx <- data.table::fread(paste0(gsub("results","data",gsub(
      "/deepCyTOF_labels","/x",dc_file)),".gz"), data.table=FALSE)
    ft_ <- get(load(gsub("results/2D/GigaSOM_clusters","data/2D/filters",gs_file)))
  }
  colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")[seq_len(length(cpops))]
  names(colours) <- cpops
  
  png(gsub(".csv",".png",gsub("labels","plots",dc_file)),width=800,height=450)
  par(mfrow=c(1,2))
  
  plot_dens(tx, main="actual")
  for (cpopi in seq_len(length(cpops))) {
    cpop <- cpops[cpopi]
    if (!nD & cpop%in%names(ft_)) lines(ft_[[cpop]], lwd=2, col=colours[cpop])
    if (nD & cpopi%in%predicted) 
      points(tx[predicted==cpopi,,drop=FALSE], cex=.1, col=colours[cpop])
  }
  legend("topright", legend=cpops, col=colours, 
         lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
  
  clust_vec <- rep(NA, length(predicted))
  for (cpopi in seq_len(length(cpops)))
    clust_vec[predicted==i] <- cpops[cpopi]
  cols <- colours[clust_vec]
  plot(tx, xlab=colanmes(tx)[1], ylab=colnames(tx)[2],
       cex=.1, col=cols, main="classified")
  legend("topright", legend=cpops, col=colours, 
         lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
  
  graphics.off()
}) }, .parallel=TRUE)
time_output(start)
