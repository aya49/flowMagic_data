# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 400x400


## parallelization ####
# future::plan(future::multiprocess)
no_cores <- 15#parallel::detectCores() - 5
doMC::registerDoMC(no_cores)

## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
libr(c(
  "furrr", #"rslurm",
  "data.table"
))


## input ####
out_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/data"
x2_dir <- paste0(out_dir,"/2D/x"); 
y2_dir <- paste0(out_dir,"/2D/y"); 


## output ####
dls <- list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)
for (dl in dls[!grepl("x$",dls)]) {
  dir.create(gsub("/data/","/results/",gsub("/x/","/x_2Ddensity/",dl)), recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("/data/","/results/",gsub("/x/","/x_2Dscatter/",dl)), recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("/data/","/results/",gsub("/x/","/x_2Dncells/",dl)), recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("/data/","/results/",gsub("/x/","/y_vector/",dl)), recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("/data/","/results/",gsub("/x/","/y_2D/",dl)), recursive=TRUE, showWarnings=FALSE)
}


## load inputs ####
x2_files <- list.files(x2_dir, recursive=TRUE, full.names=TRUE)

dimsize <- c(400,400)
overwrite <- FALSE


## START ####
start <- Sys.time()

loop_ind <- loop_ind_f(sample(seq_len(length(x2_files))), no_cores)
# res <- furrr::future_map(loop_ind, function(ii) { purrr::map(ii, function(i) {
res <- plyr::llply(loop_ind, function(ii) { purrr::map(ii, function(i) {
    # load csv
  x2discrete <- x2 <- data.table::fread(x2_files[i], data.table=FALSE)
  y2 <- data.table::fread(gsub("/x/","/y/",x2_files[i]), data.table=FALSE)
  
  # save a vector version of y
  y2i_file <- gsub("/data/","/results/",gsub("/x/","/y_vector/",x2_files[i]))
  y2i <- apply(y2, 1, function(x) which(x==1)[1])
  if ("other"%in%colnames(y2)) 
    y2i[y2i==which(colnames(y2)=="other")] <- 0
  write.table(y2i, file=gzfile(y2i_file), col.names=FALSE, row.names=FALSE, sep=",")
  
  
  # discretize x
  xr <- range(x2[,1])
  yr <- range(x2[,2])
  
  x2discrete[,1] <- ceiling(dimsize[1]*(x2[,1]-xr[1])/(xr[2]-xr[1]))
  x2discrete[,2] <- ceiling(dimsize[2]*(x2[,2]-yr[1])/(yr[2]-yr[1]))
  x2discrete[x2discrete>400] <- 400
  x2discrete_ <- as.matrix(x2discrete[!duplicated(x2discrete),,drop=FALSE])

  # grid, reverse dimensions when plotting please!
  plotsc <- matrix(0, nrow=dimsize[1], ncol=dimsize[2])
  
  # scatterplot
  plotsc[x2discrete_] <- 1
  write.table(plotsc, file=gzfile(gsub("/data/","/results/",gsub("/x/","/x_2Dscatter/",x2_files[i]))), col.names=FALSE, row.names=FALSE, sep=",")
  # gplots::heatmap.2(plotsc, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')

  # density
  dens2 = KernSmooth::bkde2D(
    x2, gridsize=c(dimsize[2],dimsize[1]),
    bandwidth=c(max(x2[,1])-min(x2[,1]), max(x2[,2])-min(x2[,2]))/30)$fhat # bandwidth in each coordinate 
  write.table(dens2, file=gzfile(gsub("/data/","/results/",gsub("/x/","/x_2Ddensity/",x2_files[i]))), col.names=FALSE, row.names=FALSE, sep=",")
  # gplots::heatmap.2(dens2, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
  
  # answer: number of cells in each pixel
  plotgn <- table(x2discrete)
  write.table(plotgn, file=gzfile(gsub("/data/","/results/",gsub("/x/","/x_2Dncells/",x2_files[i]))), col.names=FALSE, row.names=FALSE, sep=",")
  
  # answer: label of each pixel
  plotgs <- plotsc
  x2discrete_y <- apply(x2discrete_, 1, function(x) 
    getmode(y2i[x2discrete[,1]==x[1] & x2discrete[,2]==x[2]]) )
  for (yi in unique(x2discrete_y))
    plotgs[x2discrete_[x2discrete_y==yi,,drop=FALSE]] <- yi
  write.table(plotgs, file=gzfile(gsub("/data/","/results/",gsub("/x/","/y_2D/",x2_files[i]))), col.names=FALSE, row.names=FALSE, sep=",")
  # gplots::heatmap.2(plotgs, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
  
  
  
  # # contour
  # dens2_ = KernSmooth::bkde2D(
  #   x2, gridsize=dimsize,
  #   bandwidth=c(max(x2[,1])-min(x2[,1]), max(x2[,2])-min(x2[,2]))/30) # bandwidth in each coordinate 
  # dens2 <- round(1000*dens2_$fhat/max(dens2_$fhat), digits=5)
  # write.table(t(dens2), file=gzfile(gsub("/x/","/x_contour/",x2_files[i])), row.names=FALSE, col.names=FALSE, sep=",")
  # 
  # # scatter
  # png(paste0(gsub("/x/","/x_scatter/",x2_files[i]),".png"), height=dimsize[2], width=dimsize[1])
  # grid::grid.points(
  #   x2[,1],x2[,2], default.units="native", 
  #   vp=grid::dataViewport(xData=x2[,1], yData=x2[,2]),
  #   gp=grid::gpar(col="black", cex=.1), pch=6)
  # graphics.off()
  # scat2 = OpenImageR::rgb_2gray(OpenImageR::readImage(paste0(gsub("/x/","/x_scatter/",x2_files[i]),".png"))) * 1000
  # write.table(scat2, file=gzfile(gsub("/x/","/x_scatter/",x2_files[i])), row.names=FALSE, col.names=FALSE, sep=",")
  # 
  # 
  # # coloured scatter
  # dcol <- densCols(x2, map=dens2_, nbin=128, colramp=colorRampPalette(c("white","black")))
  # png(paste0(gsub("/x/","/x_colscat/",x2_files[i]),".png"), height=dimsize[1], width=dimsize[2])
  # grid::grid.points(
  #   x2[,1],x2[,2], default.units="native", 
  #   vp=grid::dataViewport(xData=x2[,1], yData=x2[,2]),
  #   gp=grid::gpar(col=dcol, cex=.1), pch=6)
  # graphics.off()
  # scat2 = OpenImageR::rgb_2gray(OpenImageR::readImage(paste0(gsub("/x/","/x_colscat/",x2_files[i]),".png"))) * 1000
  # write.table(scat2, file=gzfile(gsub("/x/","/x_colscat/",x2_files[i])), row.names=FALSE, col.names=FALSE, sep=",")

}) }, .parallel=TRUE)
time_output(start)

# res <- furrr::future_map(loop_ind, function(ii) { purrr::map(ii, function(i) {
#   # load csv
#   file.remove(paste0(gsub("/x/","/x_scatter/",x2_files[i]),".png"))
#   file.remove(paste0(gsub("/x/","/x_colscat/",x2_files[i]),".png"))
# }) })


