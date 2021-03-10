# date created: 2020-01-04
# author: alice yue
# input: score matrices
# output: plots


## set directory, load packages, set parallel ####
no_cores <- 15 #parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
sc2_dir <- paste0(root,"/scores/2D")
scn_dir <- paste0(root,"/scores/nD")

x2_folds <- list_leaf_dirs(x2_dir)
gs_xr_ <- function(x,y) gs_xr(x,y,"results") 


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


## plot 2D ####

# add scat|cpop and data|scat|cpop columns
score2D_ <- score2D[!(score2D$method=="flowLearn" & score2D$train_no!=10),]
score2D_ <- score2D_[!(score2D$method=="deepCyTOF" & score2D$train_no!=10),]
score2D_$scatpop <- paste0(score2D_$cpop," || ",score2D_$scatterplot)
score2D_$dscatpop <- paste0(score2D_$dataset, " || ", score2D_$cpop," - ",score2D_$scatterplot)

# add mean count/prop/f1 for each scat|cpop
score2D_ <- score2D_ %>% dplyr::group_by(scatpop)%>%
  dplyr::mutate(mean_true_size=mean(true_size))
score2D_ <- score2D_ %>% dplyr::group_by(scatpop)%>%
  dplyr::mutate(mean_true_proportion=mean(true_proportion))
score2D_ <- score2D_ %>% dplyr::group_by(scatpop)%>%
  dplyr::mutate(mean_f1=mean(f1))


# colour palette for continuous values
colour_palette <- c('blue','cyan','red')

## f1 scores; top=data, left=method, y=f1, x=scat|cpop; colour=mean_true_prop
g2f1size <- ggplot2::ggplot(score2D_, ggplot2::aes(
  x=reorder(scatpop, f1), # x
  y=f1, # y
  fill=mean_true_size
  )) +
  ggplot2::scale_colour_gradientn(colours=colour_palette) +
  ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
  ggplot2::geom_boxplot(outlier.shape=NA) + #coord_flip() +
  ggplot2::labs(x="(data set) scatterplot > cell population") +
  # ggplot2::theme_bw() +
  ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
  ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(root,"/scores/2D_f1_size.png"), plot=g2f1size, dpi=600, units="in", width=22, height=10)

# x ordered by scat|cpop
g2f1size_ <- ggplot2::ggplot(score2D_, ggplot2::aes(
  x=reorder(scatpop, scatterplot), # x
  y=f1, # y
  fill=mean_true_size
)) +
  ggplot2::scale_colour_gradientn(colours=colour_palette) +
  ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
  ggplot2::geom_boxplot(outlier.shape=NA) + #coord_flip() +
  ggplot2::labs(x="(data set) scatterplot > cell population") +
  # ggplot2::theme_bw() +
  ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
  ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(root,"/scores/2D_f1_size_.png"), plot=g2f1size_, dpi=600, units="in", width=22, height=10)


## f1 scores; top=data, left=method, y=f1, x=scat|cpop; colour=mean_true_prop
g2f1prop <- ggplot2::ggplot(score2D_, ggplot2::aes(
  x=reorder(scatpop, f1), # x
  y=f1, # y
  fill=mean_true_proportion
)) +
  ggplot2::scale_colour_gradientn(colours=colour_palette) +
  ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
  ggplot2::geom_boxplot(outlier.shape=NA) + #coord_flip() +
  ggplot2::labs(x="(data set) scatterplot > cell population") +
  # ggplot2::theme_bw() +
  ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
  ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(root,"/scores/2D_f1_prop.png"), plot=g2f1prop, dpi=600, units="in", width=22, height=10)

# x ordered by scat|cpop
g2f1prop_ <- ggplot2::ggplot(score2D_, ggplot2::aes(
  x=reorder(scatpop, scatterplot), # x
  y=f1, # y
  fill=mean_true_proportion
)) +
  ggplot2::scale_colour_gradientn(colours=colour_palette) +
  ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
  ggplot2::geom_boxplot(outlier.shape=NA) + #coord_flip() +
  ggplot2::labs(x="(data set) scatterplot > cell population") +
  # ggplot2::theme_bw() +
  ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
  ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(root,"/scores/2D_f1_prop_.png"), plot=g2f1prop_, dpi=600, units="in", width=22, height=10)


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

g2propap <- ggplot2::ggplot(score2D_ap, ggplot2::aes(
  x=reorder(scatpop, f1), # x
  y=proportion, # y
  fill=actualvspred
)) +
  ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
  ggplot2::geom_boxplot(outlier.shape=NA) + #coord_flip() +
  ggplot2::labs(x="(data set) scatterplot > cell population") +
  # ggplot2::theme_bw() +
  ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
  ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(root,"/scores/2D_prop.png"), plot=g2propap, dpi=600, units="in", width=22, height=10)

# x ordered by scat|cpop
g2propap_ <- ggplot2::ggplot(score2D_ap, ggplot2::aes(
  x=reorder(scatpop, scatterplot), # x
  y=proportion, # y
  fill=actualvspred
)) +
  ggplot2::facet_grid(method~dataset, scales="free_x") + # biggest box
  ggplot2::geom_boxplot(outlier.shape=NA) + #coord_flip() +
  ggplot2::labs(x="(data set) scatterplot > cell population") +
  # ggplot2::theme_bw() +
  ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
  ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))

ggplot2::ggsave(filename=paste0(root,"/scores/2D_prop_.png"), plot=g2propap_, dpi=600, units="in", width=22, height=10)


## f1 vs proportion
g2_f1prop <- ggplot2::ggplot(score2D_, ggplot2::aes(
  x=mean_true_proportion, y=mean_f1)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(method~dataset, scales="free_x")
  
ggplot2::ggsave(filename=paste0(root,"/scores/2D_propvsf1.png"), plot=g2_f1prop, dpi=600, units="in", width=10, height=10)

## f1 vs size
g2_f1size <- ggplot2::ggplot(score2D_, ggplot2::aes(
  x=mean_true_size, y=mean_f1)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(method~dataset, scales="free_x")

ggplot2::ggsave(filename=paste0(root,"/scores/2D_sizevsf1.png"), plot=g2_f1size, dpi=600, units="in", width=10, height=10)





## plotting ####
wi <- 10; hi <- 10 # number of plots per png
size <- 400

# HIPCbcell/FSCASSCA_Singlets
start <- Sys.time()
for (x2_fold in x2_folds) {
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

  # original
  lpf_fold <- gs_xr_(x2_fold,"allfilesplots_original")
  dir.create(lpf_fold, recursive=TRUE, showWarnings=FALSE)
  lpf <- paste0(lpf_fold,"/",plot_no,".png")
  xi <- 1
  for (i in seq_len(plot_no)) {
    png(file=paste0(lpf_fold,"/",i,".png"), width=wi*size, height=hi*size)
    par(mfcol=c(hi,wi))
    
    for (pi in seq_len(wi*hi)) {
      if (xi>xl) break
      
      x2 <- x2s[[xi]]
      f2 <- f2s[[xi]]
      
      f2l <- length(f2)
      plot_dens(x2, main=fnames[xi])
      for (f2i in names(f2))
        lines(f2[[f2i]], lwd=2, col=colours[f2i])
      legend("topright", legend=names(f2), col=colours, 
             lty=rep(1,f2l), lwd=rep(2,f2l))
      
      xi <- xi+1
    }
    graphics.off()
  }
  time_output(start1,"plotted original")
  
  
  # gigaSOM
  lpf_fold <- gs_xr_(x2_fold,"allfilesplots_gigaSOM")
  dir.create(lpf_fold, recursive=TRUE, showWarnings=FALSE)
  lpf <- paste0(lpf_fold,"/",plot_no,".png")
  xi <- 1
  for (i in seq_len(plot_no)) {
    png(file=paste0(lpf_fold,"/",i,".png"), width=wi*size, height=hi*size)
    par(mfcol=c(hi,wi))
    
    for (pi in seq_len(wi*hi)) {
      if (xi>xl) break
      
      x2 <- x2s[[xi]]
      f2 <- f2s[[xi]]
      
      lfile <- gs_xr_(x2_files[xi],"gigaSOM_lbls")
      if (!file.exists(lfile)) {
        f2l <- length(f2)
        plot_dens(x2, main=fnames[xi])
        for (f2i in names(f2))
          lines(f2[[f2i]], lwd=2, col=colours[f2i])
        legend("topright", legend=names(f2), col=colours, 
               lty=rep(1,f2l), lwd=rep(2,f2l))
      } else {
        gsl <- data.table::fread(lfile, data.table=FALSE)
        colours_ <- rep(NA,nrow(x2))
        for (gsi in colnames(gsl))
          colours_[gsl[,gsi]] <- gsi
        
        f2l <- length(f2)
        plot(x2, main=fnames[xi], col=colours[colours_], cex=.1)
        for (f2i in names(f2)) {
          lines(f2[[f2i]], lwd=2, col=colours[f2i])
          lines(f2[[f2i]], lwd=2, lty=3, col="black")
        }
        legend("topright", legend=names(f2), col=colours, 
               lty=rep(1,f2l), lwd=rep(2,f2l))
      }
      
      xi <- xi+1
    }
    graphics.off()
  }
  time_output(start1,"plotted gigaSOM")
  
  
  # deepCyTOF
  lpf_fold <- gs_xr_(x2_fold,"allfilesplots_deepCyTOF")
  dir.create(lpf_fold, recursive=TRUE, showWarnings=FALSE)
  lpf <- paste0(lpf_fold,"/",plot_no,".png")
  xi <- 1
  for (i in seq_len(plot_no)) {
    png(file=paste0(lpf_fold,"/",i,".png"), width=wi*size, height=hi*size)
    par(mfcol=c(hi,wi))
    
    for (pi in seq_len(wi*hi)) {
      if (xi>xl) break
      
      x2 <- x2s[[xi]]
      f2 <- f2s[[xi]]
      
      lfile <- gsub(".gz","",gs_xr_(x2_files[xi],"deepCyTOF_labels"))
      if (!file.exists(lfile)) {
        f2l <- length(f2)
        plot_dens(x2, main=fnames[xi])
        for (f2i in names(f2))
          lines(f2[[f2i]], lwd=2, col=colours[f2i])
        legend("topright", legend=names(f2), col=colours, 
               lty=rep(1,f2l), lwd=rep(2,f2l))
      } else {
        gsl <- data.table::fread(lfile, data.table=FALSE)
        colours_ <- rep(NA,nrow(x2))
        for (gsi in unique(gsl[gsl>0]))
          colours_[gsl==gsi] <- cpops[gsi]
        
        f2l <- length(f2)
        plot(x2, main=fnames[xi], col=colours[colours_], cex=.1)
        for (f2i in names(f2)) {
          lines(f2[[f2i]], lwd=2, col=colours[f2i])
          lines(f2[[f2i]], lwd=2, lty=3, col="black")
        }
        legend("topright", legend=names(f2), col=colours, 
               lty=rep(1,f2l), lwd=rep(2,f2l))
      }
      
      xi <- xi+1
    }
    graphics.off()
  }
  time_output(start1,"plotted deepCyTOF")
  
  # flowLearn
  lpf_fold <- gs_xr_(x2_fold,"allfilesplots_flowLearn")
  if (!dir.exists(lpf_fold)) next
  for (k in list.dirs(gsub("allfilesplots_flowLearn","flowLearn_y",lpf_fold), 
                      full.names=FALSE)[-1]) {
    lpf_fold <- gs_xr_(x2_fold,paste0("allfilesplots_flowLearn/",k))
    dir.create(lpf_fold, recursive=TRUE, showWarnings=FALSE)
    lpf <- paste0(lpf_fold,"/",plot_no,".png")
    xi <- 1
    for (i in seq_len(plot_no)) {
      png(file=paste0(lpf_fold,"/",i,".png"), width=wi*size, height=hi*size)
      par(mfcol=c(hi,wi))
      
      for (pi in seq_len(wi*hi)) {
        if (xi>xl) break
        
        x2 <- x2s[[xi]]
        f2 <- f2s[[xi]]
        
        lfile <- gsub(".csv.gz",paste0("_",k,".png"),
                      gs_xr_(x2_files[xi],"flowLearn_plots"))
        if (!file.exists(lfile)) {
          f2l <- length(f2)
          plot_dens(x2, main=fnames[xi])
          for (f2i in names(f2))
            lines(f2[[f2i]], lwd=2, col=colours[f2i])
          legend("topright", legend=names(f2), col=colours, 
                 lty=rep(1,f2l), lwd=rep(2,f2l))
        } else {
          gsl <- png::readPNG(lfile)
          par(mar=rep(0,4), xpd=NA)
          plot(1:10, ty="n", axes=0)
          rasterImage(gsl,1,1,10,10)
        }
        
        xi <- xi+1
      }
      graphics.off()
    }
  }
  time_output(start1,"plotted flowLearn")
  
}
time_output(start)


