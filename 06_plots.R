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

start <- Sys.time()
scn_files <- list.files(scn_dir, full.names=TRUE, recursive=TRUE, pattern=".csv.gz")
loop_ind <- loop_ind_f(scn_files, no_cores)
scorenD <- plyr::ldply(loop_ind, function(scn_fs) {
  plyr::ldply(scn_fs, data.table::fread, data.table=FALSE)
}, .parallel=TRUE)
write.csv(scorenD, file=gzfile(paste0(scn_dir,".csv.gz")))
time_output(start)


# ## plot ####
# score2D <- data.table::fread(paste0(sc2_dir,".csv.gz"))
# dsets <- unique(score2D$data_set)
# indices <- lapply(dsets, function(dset) { # dset > scat > cpop > method
#   dsi <- which(score2D$data_set==dset)
#   sd <- score2D[dsi,,drop=FALSE]
#   scats <- unique(sd$scatter_plot)
#   scis <- lapply(scats, function(scat) {
#     sci <- dsi[sd$scatter_plot==scat]
#     sc <- score2D[sci,,drop=FALSE]
#     cpops <- unique(sc$cell_population)
#     cpis <- lapply(cpops, function(cpop) {
#       cpi <- sci[sc$cell_population==cpop]
#       sp <- score2D[cpi,,drop=FALSE]
#       mthds <- unique(sp$method)
#       mtis <- lapply(mthds, function(mthd) cpi[sp$method==mthd] )
#       names(mtis) <- mthds
#       return(mtis)
#     })
#     names(cpis) <- cpops
#     return(cpops)
#   })
#   names(scis) <- scats
#   return(scis)
# })
# names(indices) <- dsets

score2D$scatpop <- paste0(score2D$scatter_plot," > ", score2D$cell_population)
score2D$dsetscat <- paste0(score2D$data_set," > ", score2D$scatter_plot)

for (dset in unique(score2D$data_set)) {
  g2 <- ggplot2::ggplot(score2D[score2D$data_set==dset,], ggplot2::aes(
    x=reorder(cell_population, f1), # x
    y=f1, # y
    fill=method,
    dodge=method)) +
    ggplot2::facet_grid(~scatter_plot) + # biggest box
    ggplot2::geom_boxplot() + #coord_flip() +
    ggplot2::labs(x="(data set) scatterplot > cell population") +
    # ggplot2::theme_bw() +
    ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
    ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
  ggplot2::ggsave(filename=paste0(root,"/scores/2D_",dset,".png"), plot=g2, dpi=600, width=15, height=18)
}


gg <- ggplot2::ggplot(scorenD, ggplot2::aes(
  x=reorder(cell_population, f1), # x
  y=f1, # y
  fill=method,
  dodge=method)) +
  ggplot2::facet_grid(~data_set, scales="free_x") + # biggest box
  ggplot2::geom_boxplot() + #coord_flip() +
  ggplot2::labs(x="(data set) cell population") +
  # ggplot2::theme_bw() +
  ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
  ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
ggplot2::ggsave(filename=paste0(root,"/scores/nD.png"), plot=gg, dpi=600, width=15, height=20)







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
      
      lfile <- gs_xr_(x2_files[xi],"gigaSOM_labels")
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
    lpf_fold <- gs_xr_(x2_fold,paste0("allfilesplots_flowLearn_",k))
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


