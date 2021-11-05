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
future::plan(future::multisession, workers=no_cores) # for furrr


## output ####
x2_folds <- list_leaf_dirs(x2_dir)
x2_folds_ <- gsub("/raw/","/data/",x2_folds)
folds <- c(
    "x_2Dcontour","x_2Dcontour_plot_",
    "x_2Ddenscat","x_2Ddenscat_plot_", # "x_2Dscatter", "x_2Ddensity",
    "x_2Ddiscrete", 
    "y_vector_","y_2D", #, "y_2Dncells"
    "temp_score"
)
plyr::l_ply(folds, function(y) plyr::l_ply(
    gs_xr(x2_folds,y), dir.create, recursive=TRUE, showWarnings=FALSE))


## load inputs ####
x2_files <- list.files(x2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")


## parameters
dimsize <- 256 # dimension of 2D outputs
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
fe <- fe[sapply(gs_xr(x2_files[fe],"y_2D"), function(x) {
    if (!file.exists(x)) return(TRUE)
    if (file.info(x)$size==0) return(TRUE)
    return(FALSE)
})]
# fe <- fe[sapply(x2_files, function(x2_file) !file.exists(paste0(gs_xr(x2_file,"temp_score"),".Rdata")))]
# fe <- fe[sapply(x2_files[fe], function(x2_file) !file.exists(gs_xr(x2_file,"x_2Ddenscat")))]
# fe <- fe[!grepl("HIPCmyeloid[/]FSCASSCA_Singlets", x2_files[fe]) & !grepl("HIPCmyeloid[/]viabilitydyeSSCA_Allcells", x2_files[fe])]
loop_ind <- loop_ind_f(sample(fe, length(fe), replace=FALSE), no_cores)
# fe <- fe[plyr::llply(loop_ind, function(ii) 
#     sapply(ii, function(i) {
#         a <- data.table::fread(gs_xr(x2_files[i],"x_2Ddenscat"))
#         return(is.na(a[1]))    
#     }), .parallel=TRUE)]
# fe <- grep("pregnancy",x2_files)
# fe <- which(Reduce("|", lapply(paste0("L0000555",c(80:95)), grepl, x2_files)))

start <- Sys.time()

f <- flowCore::read.FCS("/mnt/FCS_local3/backup/Brinkman group/current/Alice/G69019FF_SEB_CD4.fcs")

#7231

cat("out of",length(fe),"\n")
# loop_ind <- loop_ind_f(sample(fe), no_cores)
# plyr::l_ply(loop_ind, function(ii) {
a <- furrr::future_map(loop_ind, function(ii) {
# plyr::l_ply(loop_ind, function(ii) {# tryCatch({
    for (i in ii) { tryCatch({
        # res <- plyr::llply(loop_ind, function(ii) { plyr::l_ply(ii, function(i) { try({
        x2_file <- x2_files[i]
        # if (file.exists(gs_xr(x2_file,"x_2Ddenscat")))
        #     if (!is.na(data.table::fread(gs_xr(x2_file,"x_2Ddenscat"))[1,1]))
        #         next
        # if (!overwrite & file.exists(paste0(gs_xr(x2_file,"temp_score"),".Rdata")))
        #     # if (!overwrite & file.exists(gs_xr(x2_file,"y_2D")))
        #     # if (!overwrite & file.exists(gs_xr(x2_file,"x_2Dcontour")))
        #     #     if (file.info(gs_xr(x2_file,"x_2Dcontour"))$size>0)
        #     return()
        cat(i," ")
        
        # load csv
        # x2discrete <- x2 <- data.table::fread(x2_file, data.table=FALSE)
        # y2 <- data.table::fread(gsub("/x/","/y/",x2_file), data.table=FALSE)
        # 
        # save a vector version of y: a vector of cpops for each cell 0,1,2,3,...
        # y2i <- apply(y2, 1, function(x) which(x==1)[1])
        # if ("other"%in%colnames(y2))
        #     y2i[y2i==which(colnames(y2)=="other")] <- 0
        # # write.table(y2i, file=gzfile(gs_xr(x2_file,"y_vector_")),
        # #             col.names=FALSE, row.names=FALSE, sep=",")

        y2i <- read.csv(gs_xr(x2_file,"y_vector_"), header=FALSE)[,1]
        # # discretize x: cell x 2 markers i.e. which pixel each cell is at
        # xr <- range(x2[,1])
        # yr <- range(x2[,2])
        # 
        # x2discrete[,1] <- ceiling((dimsize-1)*(x2[,1]-xr[1])/(xr[2]-xr[1]))+1
        # x2discrete[,2] <- ceiling((dimsize-1)*(x2[,2]-yr[1])/(yr[2]-yr[1]))+1
        # x2discrete[x2discrete>dimsize] <- dimsize
        # x2discrete[x2discrete<1] <- 1
        # write.table(x2discrete, file=gzfile(gs_xr(x2_file,"x_2Ddiscrete")),
        #             col.names=FALSE, row.names=FALSE, sep=",")
        #
        x2discrete <- read.csv(gs_xr(x2_file,"x_2Ddiscrete"), header=FALSE)
        x2discrete_ <- as.matrix(x2discrete[!duplicated(x2discrete),,drop=FALSE])
        
        # scatterplot
        # grid, reverse y dimension when plotting please!
        plotsc <- matrix(0, nrow=dimsize, ncol=dimsize)
        plotsc[x2discrete_] <- 1
        # write.table(plotsc, file=gzfile(gs_xr(x2_file,"x_2Dscatter")),
        #             col.names=FALSE, row.names=FALSE, sep=",")
        # gplots::heatmap.2(plotsc, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
        # plotsc <- as.matrix(data.table::fread(gs_xr(x2_file,"x_2Dscatter")))
        
        # # density
        # dens2scat <- dens2 <-
        #     # KernSmooth::bkde2D(
        #     #   x2, gridsize=rep(dimsize,2),
        #     #   bandwidth=c(max(x2[,1])-min(x2[,1]), max(x2[,2])-min(x2[,2]))/30)$fhat # bandwidth in each coordinate
        #     MASS::kde2d(x2[,1],x2[,2], n=dimsize)$z
        # # write.table(dens2, file=gzfile(gs_xr(x2_file,"x_2Ddensity")),
        # #             col.names=FALSE, row.names=FALSE, sep=",")
        # # gplots::heatmap.2(dens2, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
        # 
        # min0 <- min(dens2[dens2>0])
        # dens2scat[plotsc==0] <- 0
        # dens2scat[plotsc==1 & dens2scat==0] <- min0
        # dens2scat <- 100 * dens2scat/max(dens2scat)
        # write.table(dens2scat, file=gzfile(gs_xr(x2_file,"x_2Ddenscat")),
        #             col.names=FALSE, row.names=FALSE, sep=",")
        # # gplots::heatmap.2(dens2scat, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
        # # # png(paste0(gs_xr(x2_file,"x_2Ddenscat_plot_"),".png"), width=dimsize, height=dimsize)
        # # # plot_dens(x2)
        # # # graphics.off()
        # #
        # 
        # # density contour
        # f@exprs <- as.matrix(x2)
        # hor <- unique(c(
        #     flowDensity::deGate(f, 2, all.cuts=TRUE),
        #     flowDensity::deGate(f, 2, use.upper=TRUE, upper=TRUE),
        #     flowDensity::deGate(f, 2, use.upper=TRUE, upper=FALSE) ))
        # hor <- ((dimsize-1)*(hor-yr[1])/(yr[2]-yr[1]))+1
        # ver <- unique(c(
        #     flowDensity::deGate(f, 1, all.cuts=TRUE),
        #     flowDensity::deGate(f, 1, use.upper=TRUE, upper=TRUE),
        #     flowDensity::deGate(f, 1, use.upper=TRUE, upper=FALSE) ))
        # ver <- ((dimsize-1)*(ver-xr[1])/(xr[2]-xr[1]))+1
        # 
        # gp <- ggplot2::ggplot(kd2Dto3D(dens2)) +
        #     ggplot2::geom_contour(ggplot2::aes(x=key, y=rowname, z=value), colour="black", size=.15) + # switched key and rowname for easier conversion
        #     ggplot2::scale_x_continuous(limits=c(0, dimsize), expand = c(0, 0)) +
        #     ggplot2::scale_y_continuous(limits=c(0, dimsize), expand = c(0, 0)) +
        #     ggplot2::theme(
        #         panel.grid.major=ggplot2::element_blank(),
        #         panel.grid.minor=ggplot2::element_blank(),
        #         panel.background=ggplot2::element_blank(),
        #         axis.line=ggplot2::element_blank(), axis.text=ggplot2::element_blank(),
        #         axis.ticks=ggplot2::element_blank(), axis.title=ggplot2::element_blank(),
        #         plot.margin=grid::unit(c(0,0,0,0), "mm"))
        # 
        # for (veri in ver)
        #     gp = gp + ggplot2::geom_hline(yintercept=veri, colour="black", size=.15)
        # for (hori in hor)
        #     gp = gp + ggplot2::geom_vline(xintercept=hori, colour="black", size=.15)
        # 
        # ggplot2::ggsave(paste0(gs_xr(x2_file,"x_2Dcontour_plot_"),".png"), gp, dpi=dimsize, height=1, width=1)
        # 
        # # load and greyscale
        # contpng <- png::readPNG(paste0(gs_xr(x2_file,"x_2Dcontour_plot_"),".png"))
        # # file.remove(paste0(gs_xr(x2_file,"x_2Dcontour"),".png"))
        # contpng <- contpng[,,1] + contpng[,,2] + contpng[,,3]
        # contpng <- contpng/max(contpng)
        # contpng <- 1 - contpng[nrow(contpng):1,,drop=FALSE] # flip values 0 - 1 & y axis
        # contpng <- 100 * contpng
        # write.table(contpng, file=gzfile(gs_xr(x2_file,"x_2Dcontour")),
        #             col.names=FALSE, row.names=FALSE, sep=",")
        # # gplots::heatmap.2(contpng, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
        
        
        # # contour(dens2, drawlabels=FALSE, xlim=c(0,1), ylim=c(0,1))
        # # abline(v=1); abline(v=0); abline(h=0); abline(h=1)
        # 
        # 
        # # # answer: number of cells in each pixel
        # # plotgn <- table(x2discrete)
        # # write.table(plotgn, file=gzfile(gs_xr(x2_file,"y_2Dncells")),
        # #             col.names=FALSE, row.names=FALSE, sep=",")
        # 
        # answer: label of each pixel
        plotgs <- plotsc
        x2discrete_y <- apply(x2discrete_, 1, function(x)
            getmode(y2i[x2discrete[,1]==x[1] & x2discrete[,2]==x[2]]) )
        for (yi in unique(x2discrete_y))
            plotgs[x2discrete_[x2discrete_y==yi,,drop=FALSE]] <- yi
        write.table(plotgs, file=gzfile(gs_xr(x2_file,"y_2D")),
                    col.names=FALSE, row.names=FALSE, sep=",")
        # gplots::heatmap.2(plotgs, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')

# 
#         # answer: label of each cell according to its pixel
#         # i.e. convert 2D label to vector label
#         y2dp <- apply(x2discrete, 1, function(xy) plotgs[xy[1], xy[2]])
# 
#         # baseline accuracy in this dimension
#         file_split <- stringr::str_split(x2_file, "/")[[1]]
#         fl <- length(file_split)
#         cpops <- colnames(y2)[colnames(y2)!="other"]
#         a <- cbind(
#             data.table(method="2Dbaseline", dataset=file_split[fl-2], scatpop=file_split[fl-1], cpop=cpops, train_no=0, fcs=file_split[fl], train=FALSE),
#             f1_score(y2i, y2dp, silent=TRUE))
# 
#         save(a, file=paste0(gs_xr(x2_file,"temp_score"),".Rdata"))
        # return()
        
        # }, error = function(e) return(i) ) 
        # })
    }, error = function(e) next ) }
})#, .parallel=TRUE)
# })
time_output(start)

blscore <- Reduce(rbind, plyr::llply(x2_files[fe], function(x) {
    tryCatch({
        get(load(paste0(gs_xr(x,"temp_score"),".Rdata")))
    }, error = function(e) return())
    
}))

write.table(blscore, file=gzfile(paste0(scores_dir,"/2D/pixels_baseline.csv.gz")),
            row.names=FALSE, sep=",")

# f1 scores; top=data, left=method, y=f1, x=scat|cpop; colour=mean_true_prop
xlab <- "(data set) scatterplot > cell population"
gtitle <- "scatterplot || cell population VS F1 (data set vs method)\n(colour=prop, fill=count)"
gtheme <- ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
    ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))
g2f1 <- ggplot2::ggplot(blscore, ggplot2::aes(
    x=reorder(scatpop, f1), y=f1)) + 
    ggplot2::geom_boxplot(outlier.shape=NA) + 
    ggplot2::facet_grid(method~dataset, scales="free_x") + 
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/f1_2D_baseline_",dimsize,".png"), 
                plot=g2f1, dpi=600, units="in", width=18, height=8)

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

