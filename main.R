# 2020-11-20

library(stringr)
library(data.table)
library(KernSmooth)

golden_path = "/mnt/f/Brinkman group/COVID/data/structure_test/golden_samples"
data_path = "/mnt/f/Brinkman group/COVID/data/structure_test/data"
out_path = "/mnt/f/Brinkman group/current/Alice/flowMagic_data/data"

file_endings = list.files(golden_path, pattern="csv$", recursive=TRUE)
files_tf = sapply(paste0(data_path,"/",file_endings), file.exists)
file_endings = file_endings[files_tf]

# CD27CD10- has 1 file; 
# else 17 has 4-9, 
# CD3CD56, CD4CD8, CD14CD16 and CD3+ has 30-40
# HLADR+CD14+ has 399,
# NKCD16-CD56- and CD19CD20- has 1.3k
length(file_endings)
fil_end_split = stringr::str_split(file_endings, "/")
table(sapply(fil_end_split, function(x) x[1])) 
table(sapply(fil_end_split, function(x) x[2]))
table(sapply(fil_end_split, function(x) paste0(x[1],"/",x[2])))

greyscale = colorRampPalette(c("grey", "black"))
nbin = c(250,250)
par(mar=c(0,0,0,0))
for (fe in file_endings) { try({
  print(fe)
  graphics.off()
  
  fe_f = str_extract(fe,"A[A-Za-z0-9]+[/][^/]+")
  dir.create(paste0(out_path,"/dens_points/",fe_f), recursive=TRUE, showWarnings=FALSE)
  dir.create(paste0(out_path,"/class_points/",fe_f), recursive=TRUE, showWarnings=FALSE)
  
  gf = read.csv(paste0(golden_path,"/",fe), header=TRUE)
  df = read.csv(paste0(data_path,"/",fe), header=TRUE)
  
  # density
  dens_list = densCols(df, colramp=greyscale, nbin=c(250,250))

  # class
  class_l = rep(0, nrow(df))
  for (gfc in seq_len(ncol(gf))) class_l[gf[,gfc]] = gfc
  
  png(paste0(out_path,"/dens_points/",gsub("[.]csv",".png",fe)), width=nbin[1], height=nbin[2])
  plot(df, cex=.1, axes=0, col=dens_list$cols)
  graphics.off()
  
  png(paste0(out_path,"/class_points/",gsub("[.]csv",".png",fe)), width=nbin[1], height=nbin[2])
  plot(df, cex=.1, axes=0, col=class_l)
  graphics.off()
})}

## from the document: https://docs.google.com/document/d/1sq1Iahoha-PUhmNFS2P7h8ltF5FJRf1j2qv8zNJArP0/edit
## given 2D matrix, gets 2D density + colours for each point
densCols = function(x, y=NULL, nbin=c(250,250), bandwidth=NULL, colramp=colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))) {
  xy = xy.coords(x, y, setLab=FALSE)
  select = is.finite(xy$x) & is.finite(xy$y)
  X = cbind(xy$x, xy$y)[select, ]
  
  if (is.null(bandwidth)) {
    bandwidth = diff(apply(X, 2, stats::quantile, probs=c(0.05, 0.95), na.rm=TRUE, names=FALSE))/25
    bandwidth[bandwidth == 0] = 1
  }
  if (!is.numeric(bandwidth)) stop("'bandwidth' must be numeric")
  if (any(bandwidth <= 0)) stop("'bandwidth' must be positive")
  
  map = KernSmooth::bkde2D(X, bandwidth=bandwidth, gridsize=nbin)
  mkBreaks = function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin = cut(X[, 1], mkBreaks(map$x1), labels=FALSE)
  ybin = cut(X[, 2], mkBreaks(map$x2), labels=FALSE)
  dens = map$fhat[cbind(xbin, ybin)]
  dens[is.na(dens)] = 0
  colpal = cut(dens, length(dens), labels=FALSE)
  cols = rep(NA_character_, length(select))
  cols[select] = colramp(length(dens))[colpal]
  
  return(list(cols=cols, map=map))
}




