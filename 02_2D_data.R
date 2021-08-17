# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 400x400 ++++ plot! this is to make sure everything's ok


## set directory, load packages, set parallel ####
no_cores <- 14#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
# root <- "/home/ayue/projects/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## output ####
x2_folds <- list_leaf_dirs(x2_dir)
x2_folds_ <- gsub("/raw/","/data/",x2_folds)
folds <- c("x_2Ddensity","x_2Dscatter","x_2Ddiscrete","x_2Dcontour","x_2Ddenscat","y_2Dncells","y_vector","y_2D")
plyr::l_ply(folds, function(y) plyr::l_ply(
    gs_xr(x2_folds,y), dir.create, recursive=TRUE, showWarnings=FALSE))


## load inputs ####
x2_files <- list.files(x2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")


## parameters
dimsize <- 200 # dimension of 2D outputs
overwrite <- FALSE


# function to turn 2D density matrix into n x 3 matrix
kd2Dto3D <- function(x) {
    as.data.frame(x) %>% #convert the matrix to data frame
        tibble::rownames_to_column() %>% #get row coordinates
        tidyr::gather(key, value, -rowname) %>% #convert to long format
        dplyr::mutate(key = as.numeric(gsub("V", "", key)), #convert the column names to numbers
                      rowname = as.numeric(rowname))
}


## START ####
# fe <- which(!unlist(plyr::llply(gsub("/data/2D/x","/results/2D/y_2D",x2_files), file.exists)))
fe <- 1:length(x2_files)
# fe <- grep("pregnancy",x2_files)
# fe <- which(Reduce("|", lapply(paste0("L0000555",c(80:95)), grepl, x2_files)))

start <- Sys.time()

f <- flowCore::read.FCS("/mnt/FCS_local3/backup/Brinkman group/current/Alice/G69019FF_SEB_CD4.fcs")

cat("out of",length(fe),"\n")
loop_ind <- loop_ind_f(sample(fe), no_cores)
plyr::l_ply(loop_ind, function(ii) { purrr::map(ii, function(i) {
    # res <- plyr::llply(loop_ind, function(ii) { plyr::l_ply(ii, function(i) { try({
    x2_file <- x2_files[i]
    # if (!overwrite & file.exists(gs_xr(x2_file,"y_2D"))) 
    # if (!overwrite & file.exists(gs_xr(x2_file,"x_2Dcontour")))
    #     if (file.info(gs_xr(x2_file,"x_2Dcontour"))$size>0)
    #         return()
    cat(i," ")
    
    # load csv
    x2discrete <- x2 <- data.table::fread(x2_file, data.table=FALSE)
    y2 <- data.table::fread(gsub("/x/","/y/",x2_file), data.table=FALSE)

    # save a vector version of y
    y2i <- apply(y2, 1, function(x) which(x==1)[1])
    if ("other"%in%colnames(y2))
      y2i[y2i==which(colnames(y2)=="other")] <- 0
    write.table(y2i, file=gzfile(gs_xr(x2_file,"400/y_vector")),
                col.names=FALSE, row.names=FALSE, sep=",")


    # discretize x
    xr <- range(x2[,1])
    yr <- range(x2[,2])

    x2discrete[,1] <- ceiling((dimsize-1)*(x2[,1]-xr[1])/(xr[2]-xr[1]))+1
    x2discrete[,2] <- ceiling((dimsize-1)*(x2[,2]-yr[1])/(yr[2]-yr[1]))+1
    x2discrete[x2discrete>dimsize] <- dimsize
    x2discrete[x2discrete<1] <- 1
    x2discrete_ <- as.matrix(x2discrete[!duplicated(x2discrete),,drop=FALSE])
    write.table(x2discrete, file=gzfile(gs_xr(x2_file,"x_2Ddiscrete")),
                col.names=FALSE, row.names=FALSE, sep=",")


    # scatterplot
    # grid, reverse dimensions when plotting please!
    plotsc <- matrix(0, nrow=dimsize, ncol=dimsize)
    plotsc[x2discrete_] <- 1
    write.table(plotsc, file=gzfile(gs_xr(x2_file,"x_2Dscatter")),
                col.names=FALSE, row.names=FALSE, sep=",")
    gplots::heatmap.2(plotsc, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
    plotsc <- as.matrix(data.table::fread(gs_xr(x2_file,"x_2Dscatter")))
    
    # density
    dens2scat <- dens2 <-
        # KernSmooth::bkde2D(
        #   x2, gridsize=rep(dimsize,2),
        #   bandwidth=c(max(x2[,1])-min(x2[,1]), max(x2[,2])-min(x2[,2]))/30)$fhat # bandwidth in each coordinate
        MASS::kde2d(x2[,1],x2[,2], n=200)$z
    write.table(dens2, file=gzfile(gs_xr(x2_file,"x_2Ddensity")),
                col.names=FALSE, row.names=FALSE, sep=",")
    gplots::heatmap.2(dens2, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
    
    min0 <- min(dens2[dens2>0])
    dens2scat[plotsc==0] <- 0
    dens2scat[plotsc==1 & dens2scat==0] <- min0
    dens2scat <- 100 * dens2scat/max(dens2scat)
    write.table(dens2scat, file=gzfile(gs_xr(x2_file,"x_2Ddenscat")),
                col.names=FALSE, row.names=FALSE, sep=",")
    # gplots::heatmap.2(dens2scat, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
    
    # density contour
    f@exprs <- as.matrix(x2)
    hor <- unique(c(
        flowDensity::deGate(f, 2, all.cuts=TRUE),
        flowDensity::deGate(f, 2, use.upper=TRUE, upper=TRUE), 
        flowDensity::deGate(f, 2, use.upper=TRUE, upper=FALSE) ))
    hor <- ((dimsize-1)*(hor-yr[1])/(yr[2]-yr[1]))+1
    ver <- unique(c(
        flowDensity::deGate(f, 1, all.cuts=TRUE),
        flowDensity::deGate(f, 1, use.upper=TRUE, upper=TRUE), 
        flowDensity::deGate(f, 1, use.upper=TRUE, upper=FALSE) ))
    ver <- ((dimsize-1)*(ver-xr[1])/(xr[2]-xr[1]))+1
    
    gp <- ggplot2::ggplot(kd2Dto3D(dens2)) +
        ggplot2::geom_contour(aes(x=key, y=rowname, z=value), colour="black", size=.25) + # switched key and rowname for easier conversion
        ggplot2::scale_x_continuous(limits=c(0, 200), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits=c(0, 200), expand = c(0, 0)) +
        ggplot2::theme(
            panel.grid.major=ggplot2::element_blank(),
            panel.grid.minor=ggplot2::element_blank(),
            panel.background=ggplot2::element_blank(),
            axis.line=ggplot2::element_blank(), axis.text=ggplot2::element_blank(),
            axis.ticks=ggplot2::element_blank(), axis.title=ggplot2::element_blank(),
            plot.margin=grid::unit(c(0,0,0,0), "mm"))
    
    for (veri in ver)
        gp = gp + ggplot2::geom_hline(yintercept=veri, colour="black", size=.25)
    for (hori in hor)
        gp = gp + ggplot2::geom_vline(xintercept=hori, colour="black", size=.25)
    
    ggplot2::ggsave(paste0(gs_xr(x2_file,"x_2Dcontour"),".png"), gp, dpi=dimsize, height=1, width=1)
    
    # load and greyscale
    contpng <- png::readPNG(paste0(gs_xr(x2_file,"x_2Dcontour"),".png"))
    file.remove(paste0(gs_xr(x2_file,"x_2Dcontour"),".png"))
    contpng <- contpng[,,1] + contpng[,,2] + contpng[,,3]
    contpng <- contpng/max(contpng)
    contpng <- 1 - contpng[nrow(contpng):1,,drop=FALSE] # flip values 0 - 1 & y axis
    contpng <- 100 * contpng
    write.table(contpng, file=gzfile(gs_xr(x2_file,"x_2Dcontour")),
                col.names=FALSE, row.names=FALSE, sep=",")
    # gplots::heatmap.2(contpng, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')


    # contour(dens2, drawlabels=FALSE, xlim=c(0,1), ylim=c(0,1))
    # abline(v=1); abline(v=0); abline(h=0); abline(h=1)

    # answer: number of cells in each pixel
    plotgn <- table(x2discrete)
    write.table(plotgn, file=gzfile(gs_xr(x2_file,"y_2Dncells")),
                col.names=FALSE, row.names=FALSE, sep=",")

    # answer: label of each pixel
    plotgs <- plotsc
    x2discrete_y <- apply(x2discrete_, 1, function(x)
        getmode(y2i[x2discrete[,1]==x[1] & x2discrete[,2]==x[2]]) )
    for (yi in unique(x2discrete_y))
        plotgs[x2discrete_[x2discrete_y==yi,,drop=FALSE]] <- yi
    write.table(plotgs, file=gzfile(gs_xr(x2_file,"y_2D")),
                col.names=FALSE, row.names=FALSE, sep=",")
    # gplots::heatmap.2(plotgs, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
    
}) }, .parallel=TRUE)
time_output(start)

# plyr::l_ply(loop_ind, function(ii) { purrr::map(ii, function(i) { try({
#   x2_file <- x2_files[i]
#   
#   dens2scat <- as.matrix(data.table::fread(gs_xr(x2_file,"x_2Ddenscat")))
#   if (max(dens2scat != 100)) {
#     dens2scat <- 100 * dens2scat/max(dens2scat)
#     write.table(dens2scat, file=gzfile(gs_xr(x2_file,"x_2Ddenscat")),
#                 col.names=FALSE, row.names=FALSE, sep=",")
#   }
# 
#   contpng <- as.matrix(data.table::fread(gs_xr(x2_file,"x_2Dcontour")))
#   if (max(contpng) == 1) {
#     contpng <- 100 * contpng
#     write.table(contpng, file=gzfile(gs_xr(x2_file,"x_2Dcontour")),
#                 col.names=FALSE, row.names=FALSE, sep=",")
#   }
# }) }) }, .parallel=TRUE)

