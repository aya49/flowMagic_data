# date created: 2020-01-04
# author: alice yue
# input: score matrices
# output: plots


# increase the encoding layer of deepCyTOF
# synthetic data

## set directory, load packages, set parallel ####
no_cores <- 4 #parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
sc2_dir <- paste0(root,"/scores/2D")
scn_dir <- paste0(root,"/scores/nD")

x2_folds <- list_leaf_dirs(x2_dir)

dir.create(paste0(root,"/plots/2D/scores"), recursive=TRUE, showWarnings=FALSE)


## load scores ####
start <- Sys.time()
sc2_files <- list.files(sc2_dir, full.names=TRUE, recursive=TRUE, pattern=".csv.gz")
loop_ind <- loop_ind_f(sc2_files, no_cores)
score2D <- plyr::ldply(loop_ind, function(sc2_fs)
    plyr::ldply(sc2_fs, data.table::fread, data.table=FALSE) )
write.csv(score2D, file=gzfile(paste0(sc2_dir,".csv.gz")))
time_output(start)
# score2D <- data.table::fread(paste0(sc2_dir,".csv.gz"))

# start <- Sys.time()
# scn_files <- list.files(scn_dir, full.names=TRUE, recursive=TRUE, pattern=".csv.gz")
# loop_ind <- loop_ind_f(scn_files, no_cores)
# scorenD <- plyr::ldply(loop_ind, function(scn_fs) {
#   plyr::ldply(scn_fs, data.table::fread, data.table=FALSE)
# }, .parallel=TRUE)
# write.csv(scorenD, file=gzfile(paste0(scn_dir,".csv.gz")))
# time_output(start)


# remove "other" cell population
score2D <- score2D[score2D$cpop!="other",]
score2D$train_no <- as.factor(score2D$train_no)

# add scat|cpop and data|scat|cpop columns
score2D$scatpop <- paste0(score2D$cpop," || ",score2D$scatterplot)
score2D$dscatpop <- paste0(score2D$dataset," || ",score2D$scatpop)

# add mean count/prop/f1 for each scat|cpop
score2D <- score2D %>% dplyr::group_by(dataset, scatterplot, cpop) %>%
    dplyr::mutate(mean_true_size=mean(true_size))
score2D <- score2D %>% dplyr::group_by(dataset, scatterplot, cpop) %>%
    dplyr::mutate(mean_true_proportion=mean(true_proportion))

# filter
score2D_ <- score2D %>% dplyr::filter(
    !((method=="flowLearn" & train_no!=10) | 
          (method=="deepCyTOF" & train_no!=10)) )


# colour palette for continuous values
colour_palette <- c('blue','cyan','yellow','red')
gtheme <- ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
    ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))


## f1 score plots ####
xlab <- "(data set) scatterplot > cell population"
gtitle <- "scatterplot || cell population VS F1 (data set vs method)\n(colour=prop, fill=count)"

# f1 scores; top=data, left=method, y=f1, x=scat|cpop; colour=mean_true_prop
g2f1 <- ggplot2::ggplot(score2D_, ggplot2::aes(
    x=reorder(scatpop, f1), y=f1, # y
    colour=mean_true_proportion, fill=mean_true_size)) + 
    ggplot2::geom_boxplot(outlier.shape=NA) + 
    ggplot2::facet_grid(method~dataset, scales="free_x") + 
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/f1.png"), 
                plot=g2f1, dpi=600, units="in", width=18, height=8)

# x ordered by scat|cpop
g2f1_ <- ggplot2::ggplot(score2D_, ggplot2::aes(
    x=reorder(scatpop, dscatpop), y=f1, # y
    colour=mean_true_proportion, fill=mean_true_size)) +
    ggplot2::geom_boxplot(outlier.shape=NA) + 
    ggplot2::facet_grid(method~dataset, scales="free_x") + 
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/f1_reordered.png"), 
                plot=g2f1_, dpi=600, units="in", width=18, height=8)


## actual vs pred proportions; top=data, left=method, y=f1, x=scat|cpop
# duplicate number of rows for actual vs predicted
score2D__ <- score2D___ <- score2D_
score2D__$proportion <- score2D_$true_proportion
score2D___$proportion <- score2D_$predicted_proportion
score2D__$size <- score2D_$true_size
score2D___$size <- score2D_$predicted_size
score2D__$actualvspred <- "actual"
score2D___$actualvspred <- "predicted"
score2D_ap <- rbind(score2D__,score2D___)

xlab <- "(data set) scatterplot > cell population"
gtitle <- "scatterplot || cell population VS prop"

g2propap <- ggplot2::ggplot(score2D_ap, ggplot2::aes(
    x=reorder(scatpop, f1), y=proportion, fill=actualvspred)) + 
    ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
    ggplot2::geom_boxplot(outlier.shape=NA) +
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/prop.png"), plot=g2propap, dpi=600, units="in", width=18, height=8)

# x ordered by scat|cpop
g2propap_ <- ggplot2::ggplot(score2D_ap, ggplot2::aes(
    x=reorder(scatpop, dscatpop), y=proportion, fill=actualvspred)) + 
    ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
    ggplot2::geom_boxplot(outlier.shape=NA) +
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/prop_reordered.png"), plot=g2propap_, dpi=600, units="in", width=18, height=8)


## f1 vs prop/size ####
score2D_ <- score2D %>% dplyr::group_by(dataset, scatterplot, cpop, method) %>%
    dplyr::mutate(mean_f1=mean(f1))

# f1 vs proportion
g2_f1prop <- ggplot2::ggplot(score2D_, ggplot2::aes(
    x=mean_true_proportion, y=mean_f1, colour=mean_true_size)) +
    ggplot2::ggtitle("prop vs F1") + ggplot2::geom_point() +
    ggplot2::facet_grid(method~dataset, scales="free_x") + gtheme +
    ggplot2::theme(legend.position="none")

# f1 vs size
g2_f1size <- ggplot2::ggplot(score2D_, ggplot2::aes(
    x=mean_true_size, y=mean_f1, fill=mean_true_proportion)) +
    ggplot2::ggtitle("count vs F1") + ggplot2::geom_point() +
    ggplot2::facet_grid(method~dataset, scales="free_x") + gtheme +
    ggplot2::theme(legend.position="none")

g2_f1ps <- ggpubr::ggarrange(g2_f1prop, g2_f1size, ncol=2)

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/propsizevsf1.png"), plot=g2_f1ps, dpi=600, units="in", width=8, height=4)


## f1 vs number of train samples ####
score2D_ <- score2D %>% 
    dplyr::group_by(method, train_no, dataset, scatterplot, cpop) %>%
    dplyr::summarize(mean_f1=mean(f1))

g2f1_trainno <- ggplot2::ggplot(score2D_, ggplot2::aes(
    x=train_no, y=mean_f1)) +
    ggplot2::geom_boxplot(outlier.shape=NA) + 
    ggplot2::stat_summary(fun=mean, geom="line", ggplot2::aes(group=1))  + 
    # ggplot2::scale_colour_gradientn(colours=colour_palette) +
    ggplot2::facet_grid(method~dataset, scales="free_x") + 
    ggplot2::xlab("number of train samples") + ggplot2::ggtitle("train samples VS F1")

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/f1trainno.png"), 
                plot=g2f1_trainno, dpi=600, units="in", width=18, height=8)

# flowLearn
score2D_ <- score2D %>% dplyr::filter(method=="flowLearn")

g2f1_fl <- ggplot2::ggplot(score2D_, ggplot2::aes(
    x=reorder(scatpop, f1), y=f1, # y
    colour=mean_true_proportion, fill=mean_true_size)) + 
    ggplot2::geom_boxplot(outlier.shape=NA) + 
    ggplot2::facet_grid(train_no~dataset, scales="free_x") + 
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/f1flowLearn.png"), 
                plot=g2f1_fl, dpi=600, units="in", width=16, height=18)

# deepCytOF
score2D_ <- score2D %>% dplyr::filter(method=="deepCyTOF")

g2f1_dc <- ggplot2::ggplot(score2D_, ggplot2::aes(
    x=reorder(scatpop, f1), y=f1, # y
    colour=mean_true_proportion, fill=mean_true_size)) + 
    ggplot2::geom_boxplot(outlier.shape=NA) + 
    ggplot2::facet_grid(train_no~dataset, scales="free_x") + 
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/f1deepCyTOF.png"), 
                plot=g2f1_dc, dpi=600, units="in", width=16, height=8)



## all sample grid plots ####
gs_xr_ <- function(x,y) gs_xr(x,y,"plots") 

wi <- hi <- 10
plotfunc <- function(
    lpf_fold, x2_files, lfiles, x2, f2,
    cpops, fnames, colours_f=NULL, plotf=NULL, 
    sample_cells=30000, wi=10, hi=10, size=400) {
    dir.create(lpf_fold, recursive=TRUE, showWarnings=FALSE)
    xi <- 1
    for (i in seq_len(plot_no)) {
        png(file=paste0(lpf_fold,"/",i,".png"), width=wi*size, height=hi*size)
        par(mfcol=c(hi,wi))
        
        for (pi in seq_len(wi*hi)) {
            if (xi>xl) break
            
            x2 <- x2s[[xi]]
            f2 <- f2s[[xi]]
            subsetx2 <- sample_cells<nrow(x2)
            if (subsetx2) {
                x2ind <- sample(seq_len(nrow(x2)), sample_cells)
                x2 <- x2[x2ind,,drop=FALSE]
            }
            
            lfile <- lfiles[xi]
            if (!file.exists(lfile)) {
                f2l <- length(f2)
                plot_dens(x2, main=fnames[xi])
                for (f2i in names(f2))
                    lines(f2[[f2i]], lwd=2, col=colours[f2i])
                legend("topright", legend=names(f2), col=colours, 
                       lty=rep(1,f2l), lwd=rep(2,f2l))
            } else {
                if (is.null(plotf)) {
                    colours_ <- colours_f(lfile, x2l=nrow(x2), cpops)
                    if (subsetx2) 
                        colours_ <- colours_[x2ind]
                    
                    f2l <- length(f2)
                    plot(x2, main=fnames[xi], col=colours[colours_], cex=.1)
                    for (f2i in names(f2)) {
                        lines(f2[[f2i]], lwd=2, col=colours[f2i])
                        lines(f2[[f2i]], lwd=2, lty=3, col="black")
                    }
                    legend("topright", legend=names(f2), col=colours, 
                           lty=rep(1,f2l), lwd=rep(2,f2l))
                } else {
                    plotf(lfile)
                }
            }
            
            xi <- xi+1
        }
        graphics.off()
    }
}


start <- Sys.time()
for (x2_fold in x2_folds) { try({
    start1 <- Sys.time()
    print(x2_fold)
    
    x2_files <- list.files(x2_fold, full.names=TRUE, pattern=".csv.gz")
    plot_no <- ceiling(length(x2_files)/(wi*hi))
    
    x2s <- plyr::llply(x2_files, data.table::fread, data.table=FALSE)
    f2s <- plyr::llply(gsub("/x/","/filters/",gsub(".csv.gz",".Rdata",x2_files)), 
                       function(x) get(load(x)))
    names(x2s) <- names(f2s) <- fnames <- 
        gsub(".csv.gz", "", sapply(x2_files, file_name))
    xl <- length(x2s)
    time_output(start1,"loaded files")
    
    cpops <- names(f2s[[1]])
    colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")[seq_len(length(cpops))]
    names(colours) <- cpops
    
    ## original ####
    start1 <- Sys.time()
    # plot folder name
    lpf_fold <- gs_xr_(x2_fold,"all/original")
    # label file name
    lfiles <- rep("",length(x2_files))
    # plot
    plotfunc(lpf_fold, x2_files, lfiles, x2s, f2s, cpops, fnames)
    time_output(start1,"plotted original")
    
    
    ## gigaSOM ####
    start1 <- Sys.time()
    # function to get cpop label vector given label file name
    colours_f <- function(lfile, x2l=nrow(x2s[[1]]), cpops) {
        gsl <- data.table::fread(lfile, data.table=FALSE)
        colours_ <- rep(NA,x2l)
        for (gsi in colnames(gsl)) colours_[gsl[,gsi]] <- gsi
        return(colours_)
    }
    # plot folder name
    lpf_fold <- gs_xr_(x2_fold,"all/gigaSOM")
    # label file name
    lfiles <- sapply(x2_files, function(x) gs_xr(x,"gigaSOM_labels","results") )
    # plot
    plotfunc(lpf_fold, x2_files, lfiles, x2s, f2s, cpops, fnames, colours_f)
    time_output(start1,"plotted gigaSOM")
    
    
    ## deepCyTOF ####
    start1 <- Sys.time()
    # function to get cpop label vector given label file name
    colours_f <- function(lfile, x2l=nrow(x2s[[1]]), cpops) {
        colours_ <- rep(NA,x2l)
        gsl <- data.table::fread(lfile, data.table=FALSE)
        for (gsi in unique(gsl[gsl>0])) colours_[gsl==gsi] <- cpops[gsi]
        return(colours_)
    }
    # get the number of training samples used
    a <- stringr::str_split(gs_xr(x2_fold,"deepCyTOF_labels","results"),"/")[[1]]
    ks_ <- list.dirs(paste0(a[-c((length(a)-1):length(a))], collapse="/"), 
                     recursive=FALSE, full.names=FALSE)
    for (k in ks_) {
        # plot folder name
        lpf_fold <- gs_xr_(x2_fold,paste0("all/deepCyTOF/",k))
        # label file name
        lfiles <- sapply(x2_files, function(x)
            gsub("deepCyTOF_labels",paste0("deepCyTOF_labels/",k),
                 gsub(".gz","",gs_xr(x,"deepCyTOF_labels","results"))) )
        # plot
        plotfunc(lpf_fold, x2_files, lfiles, x2s, f2s, cpops, fnames, colours_f)
    }
    time_output(start1,"plotted deepCyTOF")
    
    ## flowLearn ####
    start1 <- Sys.time()
    # function to copy and paste flowLearn plot image
    plotf <- function(lfile) {
        gsl <- png::readPNG(lfile)
        par(mar=rep(0,4), xpd=NA)
        plot(1:10, ty="n", axes=0)
        rasterImage(gsl,1,1,10,10)
    }
    # get the number of training samples used
    ks_ <- list.dirs(gs_xr(x2_fold,"flowLearn_labels","results"), 
                     full.names=FALSE, recursive=FALSE)
    for (k in ks_) {
        # plot folder name
        lpf_fold <- gs_xr_(x2_fold,paste0("all/flowLearn/",k))
        # label file name
        lfiles <- sapply(x2_files, function(x)
            gsub(".csv.gz",paste0("_",k,".png"), gs_xr_(x,"flowLearn_plots")) )
        # plot
        plotfunc(lpf_fold, x2_files, lfiles, x2s, f2s, cpops, fnames, colours_f, plotf)
    }
    time_output(start1,"plotted flowLearn")
    rm(x2s, f2s); gc()
}) }
time_output(start)


