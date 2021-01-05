# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 200x200


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 5#parallel::detectCores() - 5


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
dir.create(x2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr
y2_dir <- paste0(out_dir,"/2D/y"); 
dir.create(y2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr


## output ####
for (dl in list.dirs(x2_dir, recursive=TRUE, full.names=TRUE)) {
  dir.create(gsub("/x/","/x_contour/",dl), recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("/x/","/x_scatter/",dl), recursive=TRUE, showWarnings=FALSE)
  dir.create(gsub("/x/","/x_colscat/",dl), recursive=TRUE, showWarnings=FALSE)
}


## load inputs ####
x2_files <- list.files(x2_dir, recursive=TRUE, full.names=TRUE)

dimsize <- c(200,200)


## START ####
start <- Sys.time()

loop_ind <- loop_ind_f(seq_len(length(x2_files)), no_cores)
res <- furrr::future_map(loop_ind, function(ii) { purrr::map(ii, function(i) {
  # load csv
  x2 <- data.table::fread(x2_files[i], data.table=FALSE)
  
  # contour
  dens2_ = KernSmooth::bkde2D(
    x2, gridsize=dimsize,
    bandwidth=c(max(x2[,1])-min(x2[,1]), max(x2[,2])-min(x2[,2]))/30) # bandwidth in each coordinate 
  dens2 <- round(1000*dens2_$fhat/max(dens2_$fhat), digits=5)
  write.table(dens2, file=gzfile(gsub("/x/","/x_contour/",x2_files[i])), row.names=FALSE, col.names=FALSE, sep=",")
  
  # scatter
  png(paste0(out_dir,"/temp.png"), height=dimsize[1], width=dimsize[2])
  grid::grid.points(
    x2[,1],x2[,2], default.units="native", 
    vp=grid::dataViewport(xData=x2[,1], yData=x2[,2]),
    gp=grid::gpar(col="black", cex=.1), pch=6)
  graphics.off()
  scat2 = OpenImageR::rgb_2gray(OpenImageR::readImage(paste0(out_dir,"/temp.png"))) * 1000
  file.remove(paste0(out_dir,"/temp.png"))
  write.table(scat2, file=gzfile(gsub("/x/","/x_scatter/",x2_files[i])), row.names=FALSE, col.names=FALSE, sep=",")
  
  
  # coloured scatter
  dcol <- densCols(x2, map=dens2_, nbin=128, colramp=colorRampPalette(c("white","black")))
  png(paste0(out_dir,"/temp.png"), height=dimsize[1], width=dimsize[2])
  grid::grid.points(
    x2[,1],x2[,2], default.units="native", 
    vp=grid::dataViewport(xData=x2[,1], yData=x2[,2]),
    gp=grid::gpar(col=dcol, cex=.1), pch=6)
  graphics.off()
  scat2 = OpenImageR::rgb_2gray(OpenImageR::readImage(paste0(out_dir,"/temp.png"))) * 1000
  file.remove(paste0(out_dir,"/temp.png"))
  write.table(scat2, file=gzfile(gsub("/x/","/x_colscat/",x2_files[i])), row.names=FALSE, col.names=FALSE, sep=",")

}) })
time_output(start)



