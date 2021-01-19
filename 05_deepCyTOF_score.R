# date created: 2020-01-04
# author: alice yue
# input: deepCyTOF
# output: score matrix


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 14#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
setwd(root)


## packages ####
source("src/helpers.R")
source("src/helpers_2D.R")
libr(c(
  "furrr" #"rslurm",
))


## input ####
main_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
dc2_dir <- paste0(main_dir,"/results/2D/deepCyTOF_labels"); 
dcn_dir <- paste0(main_dir,"/results/nD/deepCyTOF_labels"); 


## output ####
score_dir <- paste0(main_dir,"/scores")
dir.create(score_dir, recursive=TRUE, showWarnings=FALSE)


## load inputs ####
# actual thresholds
dc2_files <- list.files(dc2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")
dcn_files <- list.files(dcn_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")



start <- Sys.time()

loop_ind <- loop_ind_f(dc2_files, no_cores)
scoredf_ <- furrr::future_map_dfr(loop_ind, function(dn_files) {
  purrr::map_dfr(dn_files, function(dn_file) {
    score_dir <- gsub("results/", "scores/", gsub("deepCyTOF_labels","deepCyTOF",dn_file))
    if (file.exists(score_dir))
      if (file.size(score_dir)>0)
        next
  actual <- data.table::fread(paste0(gsub("/deepCyTOF_labels","/y",gsub("/results/","/data/",dn_file)),".gz"), data.table=FALSE)
  predicted <- read.csv(dn_file)[,1]
  pus <- sort(unique(predicted))
  
  path_stuff <- stringr::str_split(dn_file,"/")[[1]]
  di <- which(path_stuff=="deepCyTOF_labels")
  dataset <- path_stuff[di+1]
  scat <- path_stuff[di+2]
  fname <- gsub(".csv","",path_stuff[length(path_stuff)])
  
  best <- purrr::map_dfr(seq_len(ncol(actual)), function(cpopi) {
    cbind(data.frame(
      data_set=dataset, scatter_plot=scat,
      cell_population=colnames(actual)[cpopi], 
      train_samples=10, 
      # total_samples=length(fnames), 
      fcs=fname
      # train=any(is.na(fts[[fname]]))
    ), f1score(actual[,cpopi]==1, predicted==cpopi))
  })
  score_dir_ <- stringr::str_split(score_dir,"/")[[1]]
  score_dir_ <- paste0(score_dir_[-length(score_dir_)], collapse="/")
  dir.create(score_dir_, recursive=TRUE, showWarnings=FALSE)
  write.table(best, file=gzfile(paste0(score_dir,".gz")), sep=",", row.names=FALSE, col.names=TRUE)
  
  })
})
# scoredf_ <- unlist(scoredf_)
# scoredf <- dplyr::bind_rows(scoredf_)
# dir.create(paste0(score_dir,"/2D"), showWarnings=FALSE)
# write.table(scoredf_, file=gzfile(paste0(score_dir,"/2D/deepCyTOF.csv.gz")), sep=",", row.names=FALSE)
time_output(start)

start <- Sys.time()

scoredf_ <- furrr::future_map_dfr(dcn_files, function(dc_file) {
  score_dir <- gsub("results/", "scores/", gsub("deepCyTOF_labels","deepCyTOF",dc_file))
  if (file.exists(score_dir))
    if (file.size(score_dir)>0)
      next
  
  actual <- data.table::fread(paste0(gsub("/deepCyTOF_labels","/y",gsub("/results/","/data/",dc_file)),".gz"), data.table=FALSE)
  predicted <- read.csv(dc_file)[,1]
  pus <- sort(unique(predicted))
  
  path_stuff <- stringr::str_split(dc_file,"/")[[1]]
  di <- which(path_stuff=="deepCyTOF_labels")
  dataset <- path_stuff[di+1]
  fname <- gsub(".csv","",path_stuff[length(path_stuff)])
  
  best <- purrr::map_dfr(seq_len(ncol(actual)), function(cpopi) {
    cbind(data.frame(
      data_set=dataset, #scatter_plot=scat,
      cell_population=colnames(actual)[cpopi], 
      train_samples=10, 
      # total_samples=length(fnames), 
      fcs=fname
      # train=any(is.na(fts[[fname]]))
    ), f1score(actual[,cpopi]==1, predicted==cpopi))
  })
  score_dir_ <- stringr::str_split(score_dir,"/")[[1]]
  score_dir_ <- paste0(score_dir_[-length(score_dir_)], collapse="/")
  dir.create(score_dir_, recursive=TRUE, showWarnings=FALSE)
  write.table(best, file=gzfile(paste0(score_dir,".gz")), sep=",", row.names=FALSE, col.names=TRUE)
})
# scoredf_ <- unlist(scoredf_)
# scoredf <- dplyr::bind_rows(scoredf_)
# dir.create(paste0(score_dir,"/nD"), showWarnings=FALSE)
# write.table(scoredf_, file=gzfile(paste0(score_dir,"/nD/deepCyTOF.csv.gz")), sep=",", row.names=FALSE)
time_output(start)


## plotting ####
loop_ind <- loop_ind_f(append(dc2_files,dcn_files), no_cores)
scoredf_ <- furrr::future_map(loop_ind, function(dc_files) {
  scoredf_ <- purrr::map(dc_files, function(dc_file) {
    actual <- data.table::fread(paste0(gsub("/deepCyTOF_labels","/y",gsub("/results/","/data/",dc_file)),".gz"), data.table=FALSE)
    predicted <- read.csv(dc_file)[,1]
    
    cpops <- colnames(actual)
    try ({
      if (grepl("2D",dc_file)) {
        x <- data.table::fread(paste0(gsub("/deepCyTOF_labels","/x",dc_file),".gz"), data.table=FALSE)
      } else {
        x <- data.table::fread(paste0(gsub("/deepCyTOF_labels","/Rtsne",dc_file),".gz"), data.table=FALSE)
        colnames(x) <- c("tsne 1", "tsne 2")
      }
      
      gs_file_ <- gsub("labels","plots",dc_file)
      gs_file_ <- stringr::str_split(gs_file_, "/")[[1]]
      gs_file_ <- paste0(gs_file_[-length(gs_file_)], collapse="/")
      dir.create(gs_file_, recursive=TRUE, showWarnings=FALSE)
      png(gsub(".csv.gz",".png",gsub("labels","plots",dc_file)), width=800, height=450)
      par(mfrow=c(1,2))
      plot_dens(x, main="actual")#, xlab="tsne 1", ylab="tsne 2")
      colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")
      colours <- colours[seq_len(length(cpops))]
      names(colours) <- cpops
      for (cpop in cpops) {
        xcpop <- as.matrix(x[actual[,cpop]==1,,drop=FALSE])
        if (nrow(xcpop)==0) next
        ch <- chull(xcpop)#$rescoords
        ch <- xcpop[c(ch,ch[1]),]
        lines(ch, lwd=2, col=colours[cpop])
      }
      legend("topright", legend=cpops, col=colours, 
             lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
      clust_vec <- rep(NA, length(predicted))
      for (cpopi in seq_len(length(cpops)))
        clust_vec[predicted==i] <- cpops[cpopi]
      cols <- colours[clust_vec]
      # cols[is.na(cols) | !good_ind] <- "black"
      plot(x, cex=.1, col=cols, main="classified")#, xlab="tsne 1", ylab="tsne 2")
      legend("topright", legend=cpops, col=colours, 
             lty=rep(1,length(cpop)), lwd=rep(2,length(cpop)))
      graphics.off()
    })
  })
})
time_output(start)
