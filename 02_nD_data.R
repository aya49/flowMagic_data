# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 400x400 ++++ plot! this is to make sure everything's ok


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
# root <- "/home/ayue/projects/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## output ####
x2_folds <- list_leaf_dirs(x2_dir)
x2_folds_ <- gsub("/data/","/results/",x2_folds)
plyr::l_ply(gsub("/x/","/scatterplots/",x2_folds), dir.create, recursive=TRUE, showWarnings=FALSE)
plyr::l_ply(
  paste0("/",c("x_2Ddensity","x_2Dscatter","x_2Ddiscrete","y_2Dncells","y_vector","y_2D"),"/"), 
  function(x) plyr::l_ply(gsub("/x/",x,x2_folds_), dir.create, recursive=TRUE, showWarnings=FALSE))


## load inputs ####
x2_files <- list.files(x2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")


## parameters
dimsize <- c(400,400) # dimension of 2D outputs
overwrite <- FALSE

wi <- 10; hi <- 6 # number of plots per png


## START ####
# fe <- which(!unlist(plyr::llply(gsub("/data/2D/x","/results/2D/y_2D",x2_files), file.exists)))
fe <- 1:length(x2_files)

start <- Sys.time()

cat("out of",length(fe),"\n")
loop_ind <- loop_ind_f(sample(fe), no_cores)
plyr::l_ply(loop_ind, function(ii) { purrr::map(ii, function(i) { try({
  # res <- plyr::llply(loop_ind, function(ii) { plyr::l_ply(ii, function(i) { try({
  cat(i," ")
  x2_file <- x2_files[i]
  
  # load csv
  x2discrete <- x2 <- data.table::fread(x2_file, data.table=FALSE)
  y2 <- data.table::fread(gsub("/x/","/y/",x2_file), data.table=FALSE)
  
  # save a vector version of y
  y2i_file <- gsub("/data/2D/x","/results/2D/y_vector",x2_file)
  y2i <- apply(y2, 1, function(x) which(x==1)[1])
  # if ("other"%in%colnames(y2))
  #   y2i[y2i==which(colnames(y2)=="other")] <- 0
  # write.table(y2i, file=gzfile(y2i_file), col.names=FALSE, row.names=FALSE, sep=",")
  
  
  # discretize x
  xr <- range(x2[,1])
  yr <- range(x2[,2])
  
  x2discrete[,1] <- ceiling((dimsize[1]-1)*(x2[,1]-xr[1])/(xr[2]-xr[1]))+1
  x2discrete[,2] <- ceiling((dimsize[2]-1)*(x2[,2]-yr[1])/(yr[2]-yr[1]))+1
  x2discrete[x2discrete>400] <- 400
  x2discrete[x2discrete<1] <- 1
  x2discrete_ <- as.matrix(x2discrete[!duplicated(x2discrete),,drop=FALSE])
  write.table(x2discrete, file=gzfile(gsub("/data/2D/x","/results/2D/x_2Ddiscrete",x2_file)), 
              col.names=FALSE, row.names=FALSE, sep=",")
  
  # grid, reverse dimensions when plotting please!
  plotsc <- matrix(0, nrow=dimsize[1], ncol=dimsize[2])
  
  # scatterplot
  plotsc[x2discrete_] <- 1
  write.table(plotsc, file=gzfile(gsub("/data/2D/x","/results/2D/x_2Dscatter",x2_file)),
              col.names=FALSE, row.names=FALSE, sep=",")
  # gplots::heatmap.2(plotsc, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
  
  # # density
  # dens2 = KernSmooth::bkde2D(
  #   x2, gridsize=c(dimsize[2],dimsize[1]),
  #   bandwidth=c(max(x2[,1])-min(x2[,1]), max(x2[,2])-min(x2[,2]))/30)$fhat # bandwidth in each coordinate
  # write.table(dens2, file=gzfile(gsub("/data/2D/x","/results/2D/x_2Ddensity",x2_file)), col.names=FALSE, row.names=FALSE, sep=",")
  # # gplots::heatmap.2(dens2, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
  
  # answer: number of cells in each pixel
  plotgn <- table(x2discrete)
  write.table(plotgn, file=gzfile(gsub("/data/2D/x","/results/2D/y_2Dncells",x2_file)),
              col.names=FALSE, row.names=FALSE, sep=",")
  
  # answer: label of each pixel
  plotgs <- plotsc
  x2discrete_y <- apply(x2discrete_, 1, function(x) 
    getmode(y2i[x2discrete[,1]==x[1] & x2discrete[,2]==x[2]]) )
  for (yi in unique(x2discrete_y))
    plotgs[x2discrete_[x2discrete_y==yi,,drop=FALSE]] <- yi
  write.table(plotgs, file=gzfile(gsub("/data/2D/x","/results/2D/y_2D",x2_file)), 
              col.names=FALSE, row.names=FALSE, sep=",")
  # gplots::heatmap.2(plotgs, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
  
}) }) }, .parallel=TRUE)
time_output(start)


## plotting ####
start <- Sys.time()
plyr::l_ply(x2_folds, function(x2_fold) {
  x2_files <- list.files(x2_fold, full.names=TRUE, pattern=".csv.gz")
  plot_no <- ceiling(length(x2_files)/(wi*hi))
  
  # dir.create(gsub("2D/x","2D/scatterplots",x2_fold), recursive=TRUE, showWarnings=FALSE)
  lpf <- paste0(gsub("2D/x","2D/scatterplots",x2_fold),"/",plot_no,".png")
  if (file.exists(lpf)) if (file.size(lpf)>0) return(NULL)
  
  xi <- 1
  for (i in seq_len(plot_no)) {
    png(file=paste0(gsub("2D/x","2D/scatterplots",x2_fold),"/",i,".png"),
        width=wi*350, height=hi*350)
    par(mfcol=c(hi,wi))
    
    for (pi in seq_len(wi*hi)) {
      x2_file <- x2_files[xi]
      if (!file.exists(x2_file)) break
      
      x2 <- data.table::fread(x2_file)
      f2 <- get(load(gsub("/x/","/filters/",gsub(".csv.gz",".Rdata",x2_file))))
      
      f2l <- length(f2)
      plot_dens(x2, main=file_name(x2_file))
      colours <- RColorBrewer::brewer.pal(f2l, "Dark2")[seq_len(f2l)]
      for (f2i in seq_len(f2l))
        lines(f2[[f2i]], lwd=2, col=colours[f2i])
      legend("topright", legend=names(f2), col=colours, 
             lty=rep(1,f2l), lwd=rep(2,f2l))
      
      xi <- xi+1
    }
    graphics.off()
  }
}, .parallel=TRUE)
time_output(start)