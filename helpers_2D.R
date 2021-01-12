# library("PBSmapping")

## functions for 02_2dgate.R; made for cmpt822 fall 2018
## 2018-11 aya43@sfu.ca



## input: 3 points, and whether to output degree
## output: angle degree between 3 points
anglefun = function(a,b,c,deg=T) {
  if (identical(a,b) | identical(b,c) | identical(a,c)) return(NA)
  ab = a-b
  bc = c-b
  ang = acos(sum(ab*bc) / ( sqrt(sum(ab^2)) * sqrt(sum(bc^2)) ))
  if (deg) ang = ang*180/pi
  return(ang)
}



## input: 3 points, and whether to output degree
## output: angle degree between 3 points
anglefun2 = function(a,b,deg=T) {
  if (identical(a,b)) return(NA)
  diff = a-b
  ang = atan2(diff[1], diff[2])
  if (deg) ang = ang*180/pi
  return(ang)
}

## input: degrees, weights of each degree (int)
## output: weighted avg of input degrees
angleavg = function(x, wts=NULL, deg=T, outdeg=T) {
  if (deg) x = x*pi/180
  if (!is.null(wts)) x = rep(x,round(wts)) # wts need to be integers; temporary work around!
  n = length(x)
  ss = sum(sin(x))
  sc = sum(cos(x))
  angl = atan2(ss, sc)
  if (outdeg) angl = angl*180/pi
  
  # quadrant adjustment: http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=How%20Linear%20Directional%20Mean%20(Spatial%20Statistics)%20works
  # if ((ss>0 | ss<0) & sc<0) angl = 180-angl
  # if (ss<0 & sc>0) angl = 360-angl
  
  cv = 1 - sqrt(ss^2 + sc^2)/n # circular var

  dat_ = list(angle=angl,cv=cv)
  return(dat_)
}

## input: degrees a,b!
## output: smallest degree difference
anglediff = function(a,b) {
  c = abs(a-b) %% 360
  c = ifelse(c>180, 360-c, c)
  return(c)
} 



## input: 2d points
## output: points converted to PBSmapping format
PBSpoints = function(x)
  return( data.frame(PID=rep(1, nrow(x)), POS=1:nrow(x), 
                     X=x[,1], Y=x[,2]) )

## input: 2d points
## output: area
pts2area = function(x)
  return( calcArea(PBSpoints(x))$area )

## input: vector areas
## output: mean 1/2 sidelength
area2meanside = function(x, sidelen_perc=.5)
  return( mean( floor(sqrt(floor(x))*sidelen_perc) ) )

## input: 2d points
## output: angles between every point
pts2angle = function(x, deg=F) {
  angle = matrix(NA, nrow=nrow(x), ncol=nrow(x))
  for (i in 1:(nrow(x)-1)) {
    for (j in (i+1):nrow(x)) {
      angle[i,j] = anglefun2(x[i,], x[j,], deg=T)
      angle[j,i] = anglefun2(x[j,], x[i,], deg=T)
    }
  }
  return(angle)
}


## input: x is a vector from 0-1 (normalized)!
## output: binned x
binvector = function(x, bins=20, norm=F) {
  minx = min(x)
  maxx = max(x)
  n = length(x)
  if (minx<maxx) {
    # x = (x-minx) / (maxx-minx)
    x = round(x*bins)
    x = round(x/min(x))
    if (norm) x = x/(maxx-minx) + minx
  } else {
    x = round(rep(1,n))
    if (norm) x = rep(1/n,n)
  }
  return(x)
}


ptsonedge = function(xx, xylim) {
  x = matrix(F, nrow=nrow(xx), ncol=nrow(xx))
  x[x<3] = T
  x[x[,1]>(xylim[1]-3), 1] = T
  x[x[,2]>(xylim[2]-3), 2] = T
  return(x)
}



## input: area, xlim, ylim of original image, xlim, ylim of new image
## output: area of new image
areaconv = function(ar, xylim, xymax_new) {
  art = (xylim$xlim[2]-xylim$xlim[1]) * 
    (xylim$ylim[2]-xylim$ylim[1])
  return( xymax_new[1]*xymax_new[2] * ar/art )
}

## input: pts, xlim, ylim of original image, xlim, ylim of new image
## output: pts of new image
ptsconv = function(pts, xylim, xymax_new) {
  xt = (xylim$xlim[2]-xylim$xlim[1])
  yt = (xylim$ylim[2]-xylim$ylim[1])
  pts[,1] = xymax_new[1]*pts[,1]/xt
  pts[,2] = xymax_new[2]*pts[,2]/yt
  return(pts)
}


## convert 2d grayscale image to rgb
grey2rgb = function(x) {
  if (length(dim(x))==2)
    x = lapply(1:3, function(i) array(x, dim=c(dim(x),1)))
  return(Reduce(abind,x))
}


## input: n x p matrix
## output: normalized cross correlation distance matrix between each n
dist_ncc = function(m, m2=NULL) {
  m = as.matrix(m)
  class(m) = "numeric"
  m_centre_ = llply(1:nrow(m), function(xi) m[xi,]-mean(m[xi,]))
  m_centre = as.matrix(do.call(rbind,m_centre_))
  # m_centre = t(apply(m, 1, function(x) x-mean(x)))
  
  if (is.null(m2)) {
    distimc = matrix(0, nrow=nrow(m), ncol=nrow(m))
    for (i in 1:(nrow(m)-1)) {
      for (j in (i+1):nrow(m)) {
        distimc[i,j] = distimc[j,i] = 
          sum(m_centre[i,] * m_centre[j,]) / 
          sqrt(sum(m_centre[i,]^2) * sum(m_centre[j,]^2))
      }
    }
  } else {
    m2 = as.matrix(m2)
    class(m2) = "numeric"
    # m_centre2 = t(apply(m2, 1, function(x) x-mean(x)))
    m_centre2_ = ldply(1:nrow(m2), function(xi) m2[xi,]-mean(m2[xi,]))
    m_centre2 = as.matrix(do.call(rbind,m_centre2_))
    distimc = matrix(0, nrow=nrow(m), ncol=nrow(m2))
    for (i in 1:nrow(m)) {
      for (j in 1:nrow(m2)) {
        distimc[i,j] = sum(m_centre[i,] * m_centre2[j,]) / 
          sqrt(sum(m_centre[i,]^2) * sum(m_centre2[j,]^2))
      }
    }
  }
  return(distimc)
}



## input: flowframe or data frame; 2 columns
## ouput: density intensity plot
# references plotDens from flowDensity
plot_int = function(dat, col, main, pch = ".", ...) {
  if (missing(col)) {
    colPalette = colorRampPalette(blues9[-(1:3)])#colorRampPalette(c("grey", "black"))
    col = densCols_(dat, colramp=colPalette)
  }
  if (nrow(dat) < 2) {
    graphics::plot(1, type = "n", axes = F, ...)
  } else {
    graphics::plot(dat, col = col, pch = pch, ...)
  }
}


densCols_ <- function(x, y=NULL, nbin=128, bandwidth=NULL, map=NULL, 
                     colramp=colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))) {
  xy <- xy.coords(x, y, setLab = FALSE)
  select <- is.finite(xy$x) & is.finite(xy$y)
  x <- cbind(xy$x, xy$y)[select, ]
  
  if (is.null(map)) {
    if (is.null(bandwidth)) {
      bandwidth <- diff(apply(x, 2, stats::quantile, probs = c(0.05, 0.95), na.rm=TRUE, names=FALSE))/25
      bandwidth[bandwidth == 0] <- 1
    } else {
      if (!is.numeric(bandwidth)) 
        stop("'bandwidth' must be numeric")
      if (any(bandwidth <= 0)) 
        stop("'bandwidth' must be positive")
    }
    map <- KernSmooth::bkde2D(x, nbin, bandwidth)
  }
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
  ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
  dens <- map$fhat[cbind(xbin, ybin)]
  dens[is.na(dens)] <- 0
  colpal <- cut(dens, length(dens), labels = FALSE)
  cols <- rep(NA_character_, length(select))
  cols[select] <- colramp(length(dens))[colpal]
  cols
}

plot_dens <- function(df, xlab=colnames(df)[1], ylab=colnames(df)[2], ...) {
  # libr(c("flowCore", "flowDensity"))
  f <- new("flowFrame")
  f@exprs <- as.matrix(df)
  flowDensity::plotDens(f, colnames(df), xlab=xlab, ylab=ylab, ...)
}

# # http://knowledge-forlife.com/r-color-scatterplot-points-density/
# plot_int2 = function(
#   df, ylim=c(min(df[,2]), max(df[,2])),  
#   xlim=c(min(df[,1]),max(df[,1])),
#   xlab=colnames(df)[1], ylab=colnames(df)[2], main="",
#   cols=NULL) {
#   
#   if (is.null(cols)) colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
#   
#   dc <- densCols(df[,1], df[,2], colramp=colorRampPalette(c("black", "white")))
#   df$dens <- col2rgb(dc)[1,] + 1L
#   df$col <- cols[df$dens]
#   plot(df[,2]~df[,1], data=df[order(df$dens),], 
#        ylim=ylim,xlim=xlim,pch=20,col=col,
#        cex=2,xlab=xlab,ylab=ylab,
#        main=main)
# }




# https://math.stackexchange.com/questions/1743995/determine-whether-a-polygon-is-convex-based-on-its-vertices
isconvex = function(vertexlist) {
  vertexlist = vertexlist[!duplicated(vertexlist),] # remove duplicated points
  N = nrow(vertexlist)
  if (N<3) return(F)
  wSign = 0        # First nonzero orientation (positive or negative)
  
  xSign = 0
  xFirstSign = 0   # Sign of first nonzero edge vector x
  xFlips = 0       # Number of sign changes in x
  
  ySign = 0
  yFirstSign = 0   # Sign of first nonzero edge vector y
  yFlips = 0       # Number of sign changes in y
  
  curr = vertexlist[N-1,]   # Second-to-last vertex
  nextv = vertexlist[N,]     # Last vertex
  
  for (vi in 1:N) {       # Each vertex, in order
    v = vertexlist[vi,]
    prev = curr          # Previous vertex
    curr = nextv          # Current vertex
    nextv = v             # Next vertex
    
    # Previous edge vector ("before"):
    bx = curr[1] - prev[1]
    by = curr[2] - prev[2]
    
    # Next edge vector ("after"):
    ax = nextv[1] - curr[1]
    ay = nextv[2] - curr[2]
    
    # Calculate sign flips using the next edge vector ("after"),
    # recording the first sign.
    if (ax > 0) {
      if (xSign == 0) {
        xFirstSign = +1
      } else if (xSign < 0) {
        xFlips = xFlips + 1
      }
      xSign = +1
    } else if (ax < 0) {
      if (xSign == 0) {
        xFirstSign = -1
      } else if (xSign > 0) {
        xFlips = xFlips + 1
      }
      xSign = -1
    }
    
    if (xFlips > 2) {
      return(F)
    }
    
    if (ay > 0) {
      if (ySign == 0) {
        yFirstSign = +1
      } else if (ySign < 0) {
        yFlips = yFlips + 1
      }
      ySign = +1
    } else if (ay < 0) {
      if (ySign == 0) {
        yFirstSign = -1
      } else if (ySign > 0) {
        yFlips = yFlips + 1
      }
      ySign = -1
    }
    
    if (yFlips > 2) {
      return(F)
    }
    
    # Find out the orientation of this pair of edges,
    # and ensure it does not differ from previous ones.
    w = bx*ay - ax*by
    if ((wSign == 0) & (w != 0)) {
      wSign = w
    } else if ((wSign > 0) & (w < 0)) {
      return(F)
    } else if ((wSign < 0) & (w > 0)) {
      return(F)
    }
  }
  
  # Final/wraparound sign flips:
  if ((xSign != 0) & (xFirstSign != 0) & (xSign != xFirstSign)) {
    xFlips = xFlips + 1
  }
  if ((ySign != 0) & (yFirstSign != 0) & (ySign != yFirstSign)) {
    yFlips = yFlips + 1
  }
  
  # Concave polygons have two sign flips along each axis.
  if ((xFlips != 2) | (yFlips != 2)) {
    return(F)
  }
  
  # This is a convex polygon.
  return(T)
}




## input: original points
## output: refined points
coord_refine = function(kptsl0_train, kptsl_train, xylim, gt_adjust=T, ptdim_train_pm, stride, range_cand) {
  feat_pts_x = seq(ptdim_train_pm+1, xylim[1]-ptdim_train_pm, stride)
  feat_pts_y = seq(ptdim_train_pm+1, xylim[2]-ptdim_train_pm, stride)
  
  if (is.null(kptsl_train)) {
    imis = names(kptsl0_train)
    kptsl_train = list()
  } else {
    imis = setdiff(names(kptsl0_train), names(kptsl_train))
  }

  for (imi in imis) {
    kptsii = NULL
    for (kpti in 1:nrow(kptsl0_train[[imi]])) { # !! assume all training files have same number of key points
      # this section was to account for the way i made the ground thruth data set, e.g. non-integer key points, out of bound points etc...
      kpt = kpt_ = round(kptsl0_train[[imi]][kpti,])
      if (gt_adjust) {
        kpt = sapply(kpt_, function(x) max(x,1));
        kpt[1] = min(kpt[1],xylim[1]);
        kpt[2] = min(kpt[2],xylim[2]);
      }
      names(kpt) = c("x","y")
      kptsii = rbind(kptsii, kpt)
    }
    # are keypoints on edge of polygon? (and image)
    types = NULL
    x = kptsii
    for (kpti in 1:nrow(x)) {
      type = rep(0,8); 
      names(type) = c("l", "r", "b", "t", "l0", "r0", "b0", "t0")
      if (x[kpti,"y"]>=xylim[2]-3) type[8] = 1
      if (x[kpti,"y"]==max(x[,"y"])) type[4] = 1
      
      if (x[kpti,"y"]<=3) type[7] = 1
      if (x[kpti,"y"]==min(x[,"y"])) type[3] = 1
      
      if (x[kpti,"x"]>=xylim[1]-3) type[6] = 1
      if (x[kpti,"x"]==max(x[,"x"])) type[2] = 1
      
      if (x[kpti,"x"]<=3) type[5] = 1
      if (x[kpti,"x"]==min(x[,"x"])) type[1] = 1
      types = rbind(types, type)
    }
    kptsii = cbind(kptsii, types)
    
    kptsii_cand = NULL
    for (kpi in 1:nrow(kptsii)) {
      # add this to get actual voxels, minus this to go back to key points
      xplus = 0
      if (kptsii[kpi,"l"]==1 & kptsii[kpi,"r"]==0) xplus = ptdim_train_pm
      if (kptsii[kpi,"l"]==0 & kptsii[kpi,"r"]==1) xplus = -ptdim_train_pm
      yplus = 0
      if (kptsii[kpi,"b"]==1 & kptsii[kpi,"t"]==0) yplus = ptdim_train_pm
      if (kptsii[kpi,"b"]==0 & kptsii[kpi,"t"]==1) yplus = -ptdim_train_pm
      
      # get candidate range from feat_pts for every point
      # if (kptsii[kpi,"r0"]==1) {        xis = rep(feat_pts_x[length(feat_pts_x)], 2)
      # } else if (kptsii[kpi,"l0"]==1) { xis = rep(feat_pts_x[1], 2)
      # } else {
        xis = c(max(feat_pts_x[which.min(abs(feat_pts_x-(kptsii[kpi,"x"]-(range_cand*xylim[1])+xplus)))], 
                    feat_pts_x[1]), 
                min(feat_pts_x[which.min(abs(feat_pts_x-(kptsii[kpi,"x"]+(range_cand*xylim[1])+xplus)))], 
                    feat_pts_x[length(feat_pts_x)]))
      # }
      # if (kptsii[kpi,"t0"]==1) {        yis = rep(feat_pts_y[length(feat_pts_y)], 2)
      # } else if (kptsii[kpi,"b0"]==1) { yis = rep(feat_pts_y[1], 2)
      # } else {
        yis = c(max(feat_pts_y[which.min(abs(feat_pts_y-(kptsii[kpi,"y"]-(range_cand*xylim[2])+yplus)))], 
                      feat_pts_y[1]), 
                  min(feat_pts_y[which.min(abs(feat_pts_y-(kptsii[kpi,"y"]+(range_cand*xylim[2])+yplus)))], 
                      feat_pts_y[length(feat_pts_y)]))
      # }
      
      kptsii_cand = rbind(kptsii_cand, c(xis,yis,xplus,yplus))
    }
    colnames(kptsii_cand) = c("xmin_cand","xmax_cand","ymin_cand","ymax_cand", "xplus","yplus")
    kptsii = cbind(kptsii, kptsii_cand)
    
    kptsl_train[[imi]] = kptsii
  }
  
  return(kptsl_train)
  
  # kptr = c(kpt[1]-ptdim_train_pm, kpt[1]+ptdim_train_pm,
  #              kpt[2]-ptdim_train_pm, kpt[2]+ptdim_train_pm, 
  #              kpt, rep(0,4))
  #     names(kptr) = c("xmin", "xmax", "ymin", "ymax", "x","y", "l", "r", "t", "b")
  #     if (kptr[1]<=3) { 
  #       kptr["l"]=1; 
  #       kptr[1:2]=kptr[1:2]+(1-kptr[1]) }
  #     if (kptr[3]<=3) { 
  #       kptr["b"]=1; 
  #       kptr[3:4]=kptr[3:4]+(1-kptr[3]) }
  #     if (kptr[2]>=xylim[1]-3) { 
  #       kptr["r"]=1; 
  #       kptr[1:2]=kptr[1:2]+(xylim[1]-kptr[2]) }
  #     if (kptr[4]>=xylim[2]-3) { 
  #       kptr["t"]=1; 
  #       kptr[3:4]=kptr[3:4]+(xylim[2]-kptr[4]) }
  #     kptrs = rbind(kptrs, kptr)
  #     filts = rbind(filts, filtl[[imi]][kptr[1]:kptr[2], kptr[3]:kptr[4],])
  #   }
  #   kptrstrl_train[[imi]] = kptrs
  #   filtstrl_train[[imi]] = filts
  # }
}




## input: flow frame
## output: # x max, x min, y max, y min thresholds + filter points to tighten up flow frame values for 2 channels
tighten_flow = function(f, channels, qt=.05, alpha=.05) {
  f_ = f
  dat = na.omit(f_@exprs[,channels])
  qts = bounds = NULL
  bounds[1] = deGate(f_, channels[1], upper=T, use.upper=T, alpha=alpha)
  qts[1] = sum(dat[,1]>bounds[1]) / nrow(dat)
  bounds[2] = deGate(f_, channels[1], upper=F, use.upper=T, alpha=alpha)
  qts[2] = sum(dat[,1]<bounds[2]) / nrow(dat)
  bounds[3] = deGate(f_, channels[2], upper=T, use.upper=T, alpha=alpha)
  qts[3] = sum(dat[,2]>bounds[3]) / nrow(dat)
  bounds[4] = deGate(f_, channels[2], upper=F, use.upper=T, alpha=alpha)
  qts[4] = sum(dat[,2]<bounds[4]) / nrow(dat)
  qbb = qts < qt
  maxdat1 = max(dat[,1])
  mindat1 = min(dat[,1])
  maxdat2 = max(dat[,2])
  mindat2 = min(dat[,2])
  bounds = c(ifelse(qbb[1], min(maxdat1, bounds[1]), maxdat1),
             ifelse(qbb[2], max(mindat1, bounds[2]), mindat1),
             ifelse(qbb[3], min(maxdat2, bounds[3]), maxdat2),
             ifelse(qbb[4], max(mindat2, bounds[4]), mindat2))
  bounds_pt = matrix(c(bounds[2], bounds[2], bounds[1], bounds[1],
                       bounds[3], bounds[4], bounds[4], bounds[3]), nrow=4)
  colnames(bounds_pt) = colnames(dat)
  
  return(list(bounds=bounds, bounds_pt=bounds_pt)) # x max, x min, y max, y min
}



## input: flow frame
## output: 2 column data frame + xylim + save png image
fcs_img = function(f, channels, imgs_dir, save_imgs=T, img_size=c(300,300), tight_all=T, qt_all=.02, alpha=.01, contour=T, contour_alpha=40) {
  fifilename = paste0(imgs_dir, "/", gsub('.fcs','', identifier(f)))
  f@exprs = na.omit(f@exprs)
  
  dat = f@exprs[,channels]
  inds = rep(T,nrow(dat))
  if (tight_all) {
    # if (!is.null(rotate_angle)) {
    #   f = rotate.data(f,channels,-rotate_angle)$data
    #   dat = rotate.data(dat,channels,-rotate_angle)$data
    # }
    bds = tighten_flow(f, channels, qt=qt_all, alpha=alpha)
    bounds_all = bds$bounds
    
    # bounds_pts_all = bds$bounds_pt
    
    # # if have filter
    # inds = dat[,1] < ifelse(maxx < bounds_all[1], bounds_all[1], Inf) & 
    #   dat[,1] > ifelse(bounds_all[2] < minx, bounds_all[2], -Inf) & 
    #   dat[,2] < ifelse(maxy < bounds_all[3], bounds_all[3], Inf) & 
    #   dat[,2] > ifelse(bounds_all[4] < miny, bounds_all[4], -Inf)
    
    # if (!is.null(rotate_angle)) {
    #   f = rotate.data(f,channels,rotate_angle)$data
    #   dat = rotate.data(dat,channels,rotate_angle)$data
    # }
    inds = dat[,1] < bounds_all[1] & 
      dat[,1] > bounds_all[2] & 
      dat[,2] < bounds_all[3] & 
      dat[,2] > bounds_all[4]
  }
  
  dat = f@exprs[inds, channels]
    
  xlim = c(min(dat[,1]), max(dat[,1]))
  ylim = c(min(dat[,2]), max(dat[,2]))

  # dat[,1] = (dat[,1]-xlim[1]) / (xlim[2]-xlim[1])
  # dat[,2] = (dat[,2]-ylim[1]) / (ylim[2]-ylim[1])
  
  png(paste0(fifilename,".png"), width=img_size[2], height=img_size[1])
  par(mar=rep(0,4))
  plot_int(dat, axes=0)
  # points(dat_flt, col='yellow', pch='.')
  # lines(bounds_pts_flt)
  # plot(flt_tight)
  # lines(cflt)
  graphics.off()

  img = rgb_2gray(readImage(paste0(fifilename,".png"))) * 255 
  # img = img[xborders,yborders]
  # img = resize(img,xylim[1])
  img = img[,ncol(img):1]
  
  if (!save_imgs) file.remove(paste0(fifilename,".png"))
  
  # contour
  if (contour) {
    dat_ = dat
    dat_[,1] = img_size[1] * (dat[,1]-xlim[1]) / (xlim[2]-xlim[1])
    dat_[,2] = img_size[2] * (dat[,2]-ylim[1]) / (ylim[2]-ylim[1])
    
    if (save_imgs) {
    png(paste0(fifilename,'_contour.png'))
    est = plotSmoothContour(dat_, alpha=contour_alpha)
    graphics.off()
  }
    
    est_ = array(0,dim=c(dim(est),1,1)); est_[,,1,1] = est
    estr = imager::resize(est_, img_size[1], img_size[2])[,,1,1] # imager
    e_max = max(estr)
    e_min = min(estr)
    estrn = 255 * (estr-e_min) / (e_max-e_min)
    
    return(list(xylim=list(xlim=xlim, ylim=ylim), dat=dat, img=img, contour=estrn))
  }
  return(list(xylim=list(xlim=xlim, ylim=ylim), dat=dat, img=img))
}



## input: flowframe, event indices or filter, filename with no file type extension, quantile (for image, for filtered region), alpha (how tight to make image, i.e. how much white space to delete)
## input: f_flt filter is required if polygon!=4 >2
## output: nothing, just saves the two things
getfilter_single = function(f, f_flt_, channels, polygon=4, envelope=T, tight=T, qt_flt=.05, alpha=.03, anglelimit=179, seed=1) { # boundingalpha=0.075) { # how tight to bound a cell population; don't make too large, long run tme (3 seconds per)
  if (!envelope & polygon<=2) 
    print("polygon variable needs to be >2; it's the number of corners in your polygon!; note if polygon=4 default, rotate your points such that the rectangle is angled upright")
  
  set.seed(seed)
  f@exprs = na.omit(f@exprs)
  temps = flowDensity(f, channels, filter=f_flt_, position=c(T,T), gates=c(-0.,0))
  f_flt_ind = temps@index
  
  f_flt = f_flt_
  
  dat = f@exprs[,channels]
  
  # # if f_flt_ind is not the actual data points, but are the indices
  # if (is.null(dim(f_flt_ind))) {
  f_flt_inds = rep(F,nrow(dat))
  f_flt_inds[f_flt_ind] = T
  dat_flt = f@exprs[f_flt_inds, channels]
  # } else {
  #   f_flt_ind = na.omit(f_flt_ind)
  #   if (ncol(f_flt_ind)>2) f_flt_ind = f_flt_ind[,channels]
  #   dat_flt = f_flt_ind
  #   if (tight) {
  #     f_flt_ind_ = f_flt_ind
  #     if (!is.null(rotate_angle)) {
  #       f_flt_ind_ = rotate.data(f_flt_ind,channels,-rotate_angle)$data
  #     }
  #     flt_inds = f_flt_ind_[,1]<=bounds_all[1] & f_flt_ind_[,1]>=bounds_all[2] & 
  #       f_flt_ind_[,2]<=bounds_all[3] & f_flt_ind_[,2]>=bounds_all[4]
  #     dat_flt = f_flt_ind[flt_inds,]
  #   }
  # }
  if (tight) {
    f__ = f
    f__@exprs = f__@exprs[1:nrow(dat_flt),]
    f__@exprs[,channels] = dat_flt
    bds = tighten_flow(f__, channels, qt=qt_flt, alpha=alpha)
    bounds_flt = bds$bounds
    bounds_pts_flt = bds$bounds_pt
    
    flt_indss = dat_flt[,1]<=bounds_flt[1] & dat_flt[,1]>=bounds_flt[2] &
      dat_flt[,2]<=bounds_flt[3] & dat_flt[,2]>=bounds_flt[4]
    dat_flt = dat_flt[flt_indss,]
  } else {
    dat_flt = na.omit(dat_flt)
    bounds_flt = c(max(dat_flt[,1]), min(dat_flt[,1]),
                   max(dat_flt[,2]), min(dat_flt[,2]))
    bounds_pts_flt = matrix(c(bounds_flt[2], bounds_flt[2], bounds_flt[1], bounds_flt[1],
                              bounds_flt[3], bounds_flt[4], bounds_flt[4], bounds_flt[3]),
                            nrow=4)
  }
  
  maxdat1 = max(dat[,1])
  mindat1 = min(dat[,1])
  maxdat2 = max(dat[,2])
  mindat2 = min(dat[,2])
  bounds_pts_all = matrix(c(mindat1, mindat1, maxdat1, maxdat1,
                            mindat2, maxdat2, maxdat2, mindat2),
                          nrow=4)
  
  
  # tightenn original filter if needed
  # maxdat1 = max(dat_flt[,1])
  # mindat1 = min(dat_flt[,1])
  # maxdat2 = max(dat_flt[,2])
  # mindat2 = min(dat_flt[,2])
  # f_flt[f_flt[,1]>maxdat1, 1] = maxdat1
  # f_flt[f_flt[,1]<mindat1, 1] = mindat1
  # f_flt[f_flt[,2]>maxdat2, 2] = maxdat2
  # f_flt[f_flt[,2]<mindat2, 2] = mindat2
  
  p1_flt = PBSpoints(f_flt)
  p2_flt_box = PBSpoints(bounds_pts_flt)
  p3_all_box = PBSpoints(bounds_pts_all)
  bounds_pts_flt_ = joinPolys(p1_flt, p2_flt_box)
  bounds_pts_flt_temp = bounds_pts_flt_ = 
    joinPolys(bounds_pts_flt_, p3_all_box)
  
  
  if (!envelope) {
    area_orig = calcArea(bounds_pts_flt_)$area
    # plot_int(dat)
    # colind = 1
    while (nrow(bounds_pts_flt_temp) > polygon) {
      area_temps = sapply(1:nrow(bounds_pts_flt_temp), function(pi) 
        abs(area_orig - calcArea(bounds_pts_flt_temp[-pi,])$area) )
      bounds_pts_flt_temp = bounds_pts_flt_temp[-which.min(area_temps),]
      # lines(bounds_pts_flt_temp[,3:4], col=colind)
      # colind = colind+1
    }
    bounds_pts_flt = bounds_pts_flt_temp[,c(3,4)]
    
    
    # angles = rep(0,nrow(f_flt))
    # angles[c(1,nrow(f_flt))] = anglefun(f_flt[nrow(f_flt)-1,], f_flt[1,], f_flt[2,], deg=T)
    # for (ai in 2:(nrow(f_flt)-1)) 
    #   angles[ai] = anglefun(f_flt[ai-1,], f_flt[ai,], f_flt[ai+1,], deg=T)
    # 
    # f_flt = f_flt[angles<anglelimit & !is.na(angles),]
    # angles = angles[angles<anglelimit & !is.na(angles)]
    # 
    # f_flt_km = cutree(hclust(dist(f_flt)), polygon)
    # f_flt_km_ = kmeans(f_flt,min(polygon*2,nrow(f_flt)))
    # f_flt_km = kmeans(as.matrix(f_flt_km_$centres),polygon)
    #
    # f_flti = sapply(unique(f_flt_km), function(ci) {
    #   ci_ind = which(f_flt_km==ci)
    #   ci_ind[which.min(angles[ci_ind])]
    # })
    # f_flt__ = f_flt
    # f_flt = f_flt[f_flti,]
    
    # # density contour
    # datz = kde2d(dat_flt[,1], dat_flt[,2])#, n=ceil(c(xlim[2]-xlim[1], ylim[2]-ylim[1]) / boundingalpha)) #MASS
    # cflt = contourLines(datz) #, nlevels=1, levels=max(datz$datz[, which(datz$y > ssca.bound) ]))
    # cfltmax = which.max(sapply(cflt, function(x){length(x$x)}))
    # cflt = cbind(cflt[[cfltmax]]$x, y=cflt[[cfltmax]]$y)
    
    
    # # envelope
    # flt_tight = structure(as.data.frame(dat_flt[chull(dat_flt[,1],dat_flt[,2]), ]))
    # coordinates(flt_tight) = cbind(flt_tight[,1],flt_tight[,2])
    # proj4string(flt_tight) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    # gEnvelope(flt_tight)
    
    # # rgeos
    # gConvexHull(flt)
    # x = readWKT(paste("POLYGON((0 40,10 50,0 60,40 60,40 100,50 90,60 100,60", "60,100 60,90 50,100 40,60 40,60 0,50 10,40 0,40 40,0 40))"))
    # env = gEnvelope(x)
    # plot(x,col='blue', border='green')
    # plot(env,add=TRUE)
    
    # gPolygonize(linelist,getCutEdges=TRUE) # no cut edges
    
  }
  
  dat__ = list(flt=bounds_pts_flt_[,c(3:4)], kpts=bounds_pts_flt)
  return(dat__)
}



# same as getfilter_single; but with more filters per frame; combined, 
# output filter contain combined key points for all filters
# ## CHANGE: SAVE THE DESNITY MAP INSTEAD, then go and revers the y in the 02_dataprep.R script!!
getfilter_multiple = function(f, f_flt_, channels, 
                              polygon=NULL, envelope=NULL, polygon_total=NULL,
                              tight_all=T, tight=NULL,
                              imgsize, fifilename, 
                              qt_all=.02, qt_flt=.05, alpha=.03, 
                              anglelimit=179, seed=1) {
  require("flowDensity")
  
  if (is.null(dim(f_flt_))) {
    # f_flt_rbind = Reduce('rbind',f_flt_)
    if (is.null(polygon)) polygon = rep(4, length(f_flt_))
    if (is.null(envelope)) envelope = rep(T, length(f_flt_))
    if (is.null(tight)) tight = rep(T, length(f_flt_))
    if (is.null(polygon_total)) polygon_total = sum(polygon)
  } else {
    # f_flt_rbind = f_flt_
    if (is.null(polygon)) polygon = 4
    if (is.null(polygon_total)) polygon_total = polygon
    if (is.null(envelope)) envelope = T
    if (is.null(tight)) tight = T
  }
  
  f_ = f
  f@exprs = na.omit(f@exprs)
  
  dat = f@exprs[,channels]
  inds = rep(T,nrow(dat))
  if (tight_all) {
    # if (!is.null(rotate_angle)) {
    #   f = rotate.data(f,channels,-rotate_angle)$data
    #   dat = rotate.data(dat,channels,-rotate_angle)$data
    # }
    bds = tighten_flow(f, channels, qt=qt_all, alpha=alpha)
    bounds_all = bds$bounds
    
    if (is.null(dim(f_flt_))) {
      maxx = max(sapply(f_flt_, function(x) min(x[,1]))) # _| :( ||_| :) _|
      minx = min(sapply(f_flt_, function(x) max(x[,1]))) # |_ :) |_|| :( |_
      maxy = max(sapply(f_flt_, function(x) min(x[,2])))
      miny = min(sapply(f_flt_, function(x) max(x[,2])))
    } else {
      maxx = min(f_flt_[,1])
      minx = max(f_flt_[,1])
      maxy = min(f_flt_[,2])
      miny = max(f_flt_[,2])
    }
    # bounds_pts_all = bds$bounds_pt
    inds = dat[,1] < ifelse(maxx < bounds_all[1], bounds_all[1], Inf) & 
      dat[,1] > ifelse(bounds_all[2] < minx, bounds_all[2], -Inf) & 
      dat[,2] < ifelse(maxy < bounds_all[3], bounds_all[3], Inf) & 
      dat[,2] > ifelse(bounds_all[4] < miny, bounds_all[4], -Inf)
    # if (!is.null(rotate_angle)) {
    #   f = rotate.data(f,channels,rotate_angle)$data
    #   dat = rotate.data(dat,channels,rotate_angle)$data
    # }
  }
  f@exprs = f@exprs[inds,]
  dat = f@exprs[, channels]
  xlim = c(min(dat[,1]), max(dat[,1]))
  ylim = c(min(dat[,2]), max(dat[,2]))
  # dat[,1] = (dat[,1]-xlim[1]) / (xlim[2]-xlim[1])
  # dat[,2] = (dat[,2]-ylim[1]) / (ylim[2]-ylim[1])
  # flt[,1] = (flt[,1]-xlim[1]) / (xlim[2]-xlim[1])
  # flt[,2] = (flt[,2]-ylim[1]) / (ylim[2]-ylim[1])
  
  # multiple filters?
  if (is.null(dim(f_flt_))) {
    f_flt_ = lapply(f_flt_, function(x) { colnames(x) = channels; x })
    if (length(f_flt_)!=length(polygon) | 
        length(f_flt_)!=length(envelope) | 
        length(f_flt_)!=length(tight)) {
      print("must be same length: f_flt_ filter, polygon, envelope, tight_all")
      return(NULL)
    }
    f_flt_combo = NULL
    f_flt_combol = list()
    f_flt_l = f_flt_
    for (i in 1:length(f_flt_)) {
      single_flt = getfilter_single(f, f_flt_[[i]], channels, 
                                      polygon=polygon[i], envelope=envelope[i], 
                                      tight=tight[i], qt_flt=qt_flt, alpha=alpha)
      f_flt_combo = rbind(f_flt_combo, single_flt$kpts)
      f_flt_combol = append(f_flt_combol, single_flt$kpts)
      f_flt_l[[i]] = single_flt$flt
    }
    # find points that are repeated for the different filters, these are just one point to find
    f_flt_combo = f_flt_combo[!duplicated(f_flt_combo),]
    
    if (nrow(f_flt_combo)>polygon_total) {
      f_flt_combo_clust = cutree(hclust(dist(f_flt_combo)), polygon_total)
      f_flt_combo_ = f_flt_combo
      f_flt_combo = NULL
      for (ci in unique(f_flt_combo_clust))
        f_flt_combo = rbind(f_flt_combo, colMeans(f_flt_combo_[f_flt_combo_clust%in%ci,]))
    }
  } else {
    colnames(f_flt_) = channels
    single_flt = getfilter_single(f, f_flt_, channels, polygon=polygon, envelope=envelope, tight=T, qt_flt=qt_flt, alpha=alpha)
    f_flt_combo = f_flt_combol = single_flt$kpts
    f_flt_l = single_flt$flt
  }
  
  
  ## PNG
  png(paste0(fifilename,".png"), width=imgsize[2], height=imgsize[1])
  par(mar=rep(0,4))
  plot_int(dat, axes=0)
  # points(dat_flt, col='yellow', pch='.')
  # lines(bounds_pts_flt)
  # plot(flt_tight)
  # lines(cflt)
  graphics.off()
  ## PNG
  
  png(paste0(fifilename,"_filter.png"), width=imgsize[2], height=imgsize[1])
  par(mar=rep(0,4))
  plot_int(dat, axes=0)
  if (is.null(dim(f_flt_))) {
    for (fi in 1:length(f_flt_))
      lines(f_flt_[[fi]], lty=2, col=fi)
  } else {
    lines(f_flt_, lty=2)
  }
  points(f_flt_combo)
  graphics.off()
  
  # # area of smallest filter
  # if (is.null(dim(f_flt_))) {
  #   fia = min( unlist(lapply(f_flt_, function(x) {
  #     x = bounds_pts_flt(x)
  #     calcArea(x)$area
  #   })) )
  # } else {
  #   x = bounds_pts_flt(f_flt_)
  #   fia = calcArea(x)$area
  # }
  
  dat_ = list(flt=f_flt_combo, fltl=f_flt_combol, flt_orig=f_flt_l, xylim=list(xlim=xlim, ylim=ylim), dat=dat)
  save(dat_, file=paste0(fifilename,".Rdata"))
}


## DON"T USE NOT FIXED!!
## input: same as getfilter_single, but gates are 2d thresholds
## output: keypoints for 4 corners, and centre point
getquadfilter_2d = function(f, gates, channels, tight=T, tight_all=T, imgsize, fifilename, qt_all=.05, qt_flt=.05, alpha=.03) {
  f_ = f
  f@exprs = na.omit(f@exprs)
  
  dat = f@exprs[,channels]
  inds = rep(T,nrow(dat))
  if (tight_all) {
    # if (!is.null(rotate_angle)) {
    #   f = rotate.data(f,channels,-rotate_angle)$data
    #   dat = rotate.data(dat,channels,-rotate_angle)$data
    # }
    bds = tighten_flow(f, channels, qt=qt_all, alpha=alpha)
    bounds_all = bds$bounds
    # bounds_pts_all = bds$bounds_pt
    inds = dat[,1]<bounds_all[1] & dat[,1]>bounds_all[2] & 
      dat[,2]<bounds_all[3] & dat[,2]>bounds_all[4]
    # if (!is.null(rotate_angle)) {
    #   f = rotate.data(f,channels,rotate_angle)$data
    #   dat = rotate.data(dat,channels,rotate_angle)$data
    # }
  }
  f@exprs = f@exprs[inds,]
  dat = f@exprs[, channels]
  xlim = c(min(dat[,1]), max(dat[,1]))
  ylim = c(min(dat[,2]), max(dat[,2]))
  # dat[,1] = (dat[,1]-xlim[1]) / (xlim[2]-xlim[1])
  # dat[,2] = (dat[,2]-ylim[1]) / (ylim[2]-ylim[1])
  # flt[,1] = (flt[,1]-xlim[1]) / (xlim[2]-xlim[1])
  # flt[,2] = (flt[,2]-ylim[1]) / (ylim[2]-ylim[1])
  
  bounds_flt_pts = gates
  for (tfi in 1:4) {
    # x max min, y max min
    if (tfi==1) dat_flt = dat[dat[,1]>=gates[1] & dat[,2]>=gates[2],] # top right
    if (tfi==2) dat_flt = dat[dat[,1]>=gates[1] & dat[,2]< gates[2],] # bot right
    if (tfi==3) dat_flt = dat[dat[,1]< gates[1] & dat[,2]< gates[2],] # bot left
    if (tfi==4) dat_flt = dat[dat[,1]< gates[1] & dat[,2]>=gates[2],] # top left
    if (tight) {
      f_ = f
      f_@exprs = f_@exprs[1:nrow(dat_flt),]
      f_@exprs[,channels] = dat_flt
      bds = tighten_flow(f_, channels, qt=qt_flt, alpha=alpha)
      bounds = bds$bounds
    } else {
      dat_flt = na.omit(dat_flt)
      bounds = c(max(dat_flt[,1]), min(dat_flt[,1]),
                 max(dat_flt[,2]), min(dat_flt[,2]))
    }
    if (tfi==1) bounds_flt_pts = rbind(bounds_flt_pts, c(bounds[1], bounds[3])) # top right
    if (tfi==2) bounds_flt_pts = rbind(bounds_flt_pts, c(bounds[1],bounds[4])) # bot right
    if (tfi==3) bounds_flt_pts = rbind(bounds_flt_pts, c(bounds[2],bounds[4])) # bot left
    if (tfi==4) bounds_flt_pts = rbind(bounds_flt_pts, c(bounds[2],bounds[3])) # top left
  }
  # bounds_flt_pts = bounds_flt_pts[!duplicated(bounds_flt_pts)]
  
  png(paste0(fifilename,".png"), width=imgsize[2], height=imgsize[1])
  par(mar=rep(0,4))
  plot_int(dat, axes=0)
  dev.off()
  
  png(paste0(fifilename,"_filter.png"), width=imgsize[2], height=imgsize[1])
  par(mar=rep(0,4))
  plot_int(dat, axes=0)
  points(bounds_flt_pts, cex=2, col='red')
  dev.off()
  
  dat_ = list(flt=bounds_flt_pts, flt_orig=gates, xylim=list(xlim=xlim, ylim=ylim), dat=dat)
  save(dat_, file=paste0(fifilename,".Rdata"))
}



# quadtree
library(imager)
library(purrr)

#Divide along x, then y
qsplit <- function(im) {
  imsplit(im,"x",2) %>% map(~ imsplit(.,"y",2)) %>% flatten 
}

qunsplit <- function(l) {
  list(l[1:2],l[3:4]) %>% map(~ imappend(.,"y")) %>% imappend("x")
}

#Max std. dev across channels
imsd <- function(im) {
  imsplit(im,"c") %>% map_dbl(sd) %>% max
}

rebuild <- function(l,borders=FALSE) {
  map(l[-5], ~ ifelse (is.cimg(.)), meanim(.,borders=borders), rebuild(.,borders=borders)) %>% qunsplit
}

refine <- function(l) {
  if (is.cimg(l)) { #We have a leaf
    qs <- qsplit(l) #Split
    if (any(dim(l)[1:2] <= 4)) { #Quadrants are very small
      qs$sds <- rep(0,4) #Prevent further refinement
    } else {
      qs$sds <- map_dbl(qs,imsd) #Store std.dev of children
    }
    qs
  }
  else { #Not a leaf, explore further
    indm <- which.max(l$sds) #Find child with max. std. dev
    l[[indm]] <- refine(l[[indm]]) #Refine
    l$sds[indm] <- max(l[[indm]]$sds) #Update std. dev
    l
  }
}

#Produce an image that's just the average of image im
#Optionally, add borders
meanim <- function(im, borders=F) {
  im <- imsplit(im,"c") %>% map(~ 0*. + mean(.)) %>% imappend("c")
  if (borders) {
    im[px.borders(im)] <- 0
  }
  im
}

iter.refine <- function(im,nIter) {
  for (i in seq_len(nIter)) im <- refine(im)
  im
}

qtree = function(im) {
  # https://dahtah.github.io/imager/quadtrees.html
  
  # qsplit(im) %>% as.imlist %>% plot
  # qsplit(im) %>% qunsplit %>% plot
  
  #The first four iterations of the process
  map_il(1:4,~ iter.refine(im,.) %>% rebuild) %>% plot
  
  #After 200 iterations
  iter.refine(im,200) %>% rebuild(borders=T) %>% plot
}

## input: 2d matrix
## output: smoothed contour + filled
plotSmoothContour = function(dat, alpha=30) {
  # dat = get(load("/home/ayue/projects/2dgate/result/bcell/images/subset1/07.0_igdpos_cd27neg/2017-10-04_1b_Specimen_001_02-525-M_30Jul2009_013.Rdata"))$dat
  dat = na.omit(dat)
  require("MASS")
  require("KernSmooth")
  est = KernSmooth::bkde2D(dat, bandwidth=c(max(dat[,1])-min(dat[,1]), max(dat[,2])-min(dat[,2]))/alpha)$fhat # bandwidth in each coordinate
  filled.contour(est, color.palette=colorRampPalette(c("white", "black")))
  # est = kde2d(dat[,1], dat[,2], n=50)$z
  # filled.contour(est, color.palette=colorRampPalette(c("white", "black")), axes=F)
  return(est)
}






# sybil plot functions -------------------------------------------------------

plotDensContour = function(fframe, channels, main, numlevels = 10, pch ='.',  ...){
  plotDens(fframe, channels = channels, main = main, pch = pch, ...) 
  data.new = na.omit(exprs(fframe)[, channels])
  z = kde2d(data.new[, 1], data.new[, 2], n = 50)
  contour(z, drawlabels = FALSE, add = TRUE, nlevels = numlevels, lty = 2)
}  

plotDens2 = function(ff, channels, main, numlevels = 10,  ...){
  # Automatically increases the point size when there are < 1000 events to plot
  if(nrow(ff@exprs) > 1000){
    plotDens(ff, channels = channels, main = main, ...)
  } else {
    plotDens(ff, channels = channels, main = main, pch = 20, cex = 0.2, ...)
  }
  if (nrow(ff@exprs) > 500){
    data.new = na.omit(exprs(ff)[, channels])
    z = kde2d(data.new[, 1], data.new[, 2], n = 50)
    contour(z, drawlabels = FALSE, add = TRUE, nlevels = numlevels, lty = 2)
  }
}

