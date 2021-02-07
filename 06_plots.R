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
score2D <- plyr::ldply(loop_ind, function(scn_fs) {
  plyr::ldply(scn_fs, data.table::fread, data.table=FALSE)
}, .parallel=TRUE)
write.csv(scorenD, file=gzfile(paste0(scn_dir,".csv.gz")))
time_output(start)


## plot ####
score2D <- data.table::fread(paste0(sc2_dir,".csv.gz"))
dsets <- unique(score2D$data_set)
indices <- lapply(dsets, function(dset) { # dset > scat > cpop > method
  dsi <- which(score2D$data_set==dset)
  sd <- score2D[dsi,,drop=FALSE]
  scats <- unique(sd$scatter_plot)
  scis <- lapply(scats, function(scat) {
    sci <- dsi[sd$scatter_plot==scat]
    sc <- score2D[sci,,drop=FALSE]
    cpops <- unique(sc$cell_population)
    cpis <- lapply(cpops, function(cpop) {
      cpi <- sci[sc$cell_population==cpop]
      sp <- score2D[cpi,,drop=FALSE]
      mthds <- unique(sp$method)
      mtis <- lapply(mthds, function(mthd) cpi[sp$method==mthd] )
      names(mtis) <- mthds
      return(mtis)
    })
    names(cpis) <- cpops
    return(cpops)
  })
  names(scis) <- scats
  return(scis)
})
names(indices) <- dsets

score2D$scatpop <- paste0(score2D$scatter_plot," > ", score2D$cell_population)


for (dset in names(indices)) {
    dat <- score2D[score2D$data_set==dset & 
                     (score2D$scatter_plot=="CD19SSCA_Nongranulocytes" |
                        score2D$scatter_plot=="FSCASSCA_Singlets"),,drop=FALSE]

    gg <- ggplot2::ggplot(dat, ggplot2::aes(
      x=cell_population, # x
      y=f1, # y
      fill=method)) +
      ggplot2::facet_grid(~scatter_plot) + # biggest box
      ggplot2::geom_boxplot() +
      # ggplot2::labs(x="") +
      # ggplot2::theme_bw() +
      ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
      ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) + 
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
    gg
}



