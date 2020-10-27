rm(list=ls())
library(alphahull)
library(igraph)
library(plyr)
library(doMC)
registerDoMC(8)
ahull2polygon <- function(n10){
  # The return structure of ashape doesn't make it easy to get the polygon out
  # We can re-order the points using the igraph package. 
  # The ashape function returns an object that has the indexes of the start and end point of each segment. We can feed these in as an edge list for a graph object:
  # Input : n10 is an alphahull
  # Output : 
  n10g = graph.edgelist(cbind(as.character(n10$edges[, 'ind1']), as.character(n10$edges[, 'ind2'])), directed = FALSE)
  if(!is.connected(n10g)){
    return(1)
    # stop("Graph not connected")
  }
  # if(any(degree(n10g) != 2)){
  #   stop("Graph not circular")
  # }
  if(igraph::clusters(n10g)$no > 1){
    stop("Graph composed of more than one circle")
  }
  
  
  cutg = n10g - E(n10g)[1]
  # find chain end points
  ends = names(which(degree(cutg) == 1))
  path = get.shortest.paths(cutg, ends[1], ends[2])[[1]][[1]]
  # this is an index into the points
  pathX = as.numeric(V(n10g)[path]$name)
  # join the ends
  pathX = c(pathX, pathX[1])
  # now show the alpha shape plot with our poly on top
  #plot(n10, lwd = 10, col = "gray")
  # get the points from the ashape object
  #lines(n10$x[pathX, ], lwd = 2)
  
  return(poly = n10$x[pathX, ])
}

generate.gates.usr <- function(dat, clrs, weight=c("No_weight","ARI_weight")[1], 
                               no=c("L100")[1])
{
  # if(nrow(dat)>60000){
  #   idx <- sample(1:nrow(dat), 60000)
  #   dat <- dat[idx,]
  # }
  which.col <- grep(colnames(clrs),pattern = paste0(weight,"_",no))
  clrs.new <- clrs[,which.col]
  len <- sort(setdiff(unique(clrs.new),c(0)),decreasing = F)
  1 <- lapply(len,function(i1){
    X <- dat[which(clrs.new==i1),]
    if(nrow(X)>6000){
      x.idx <- sample(1:nrow(X),6000)
      X <- dat[x.idx,]
    }
    x2.ah <- ashape(X, alpha = 0.5)
    filter <- try(ahull2polygon(n10=x2.ah))

    return(filter)
  })
}

generate.gates.usr.chull <- function(dat, clrs, weight=c("No_weight","ARI_weight")[1], 
                               no=c("L100")[1])
{
  which.col <- grep(colnames(clrs),pattern = paste0(weight,"_",no))
  clrs.new <- clrs[,which.col]
  len <- sort(setdiff(unique(clrs.new),c(0)),decreasing = F)
  gates <- lapply(len,function(i1){
    X <- dat[which(clrs.new==i1),]
    filter <- chull(X)
    filter <- X[c(filter,filter[1]),]
    return(filter)
  })
}

file.info <- as.character(read.table("/code/support_files/files_to_analyze.csv",
                        sep=",",header=T,check.names=F)[,1])
file.info <- list.files("/data/results/")
file.info <- file.info[grep("FR-FCM-Z2N5", file.info)]
file.info <- gsub("\u20A0", "/", file.info)
file.info <- file.info[-grep("_pixel", file.info)]
file.info <- gsub("_results.csv","", file.info)
file.info <- file.info[-1]
res.dir <- "/code/test_users/"
data.dir <- "/data/data/"
user.res.dir <- "/data/results/"
gold.dir <- "/data/golden_samples/"
only.col<-F
which.one <- ifelse(only.col==T,yes = "color",no="gate")
grbg <- llply(file.info[1:length(file.info)], function(f1){
  tryCatch({
      print(f1)
       data <- read.table(paste0(data.dir,f1,".csv"),header = T,sep=",",check.names = F)
       # gold.clrs <- tryCatch(read.table(paste0(gold.dir,f1,".csv"),header = T,sep=","
                                        # ,check.names = F),error=function(ex) {return(NA)})
       # if(is.na(gold.clrs))
       #   gold.clrs <- read.table(paste0("~/project/COVID/data/structure_test/clr/",f1,".csv"),header = T,sep=","
       #                                    ,check.names = F)
      usr.clrs <- read.table(paste0(user.res.dir,
                                    gsub(f1,pattern = "/",replacement = '\u20A0'),"_results.csv")
                             ,header = T,sep=",",check.names = F)
      colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
      col <- densCols(data, colramp = colPalette)
      png(paste0(res.dir,gsub(f1,pattern = "/",replacement = '\u20A0'),which.one,"_.png"), width=1200, height=600)
      par(mfrow=c(1,2))
      # plot(data,col=col,pch=".",main="Gold")
      # g.gates <- generate.gates.gold(data,clrs = gold.clrs)
      # tmp <- lapply(g.gates, function(g1) lines(g1,type="l",lwd=2))
      # plot(data,col=col,pch=".",main="100 usrs, no weight, w clusters")
      weight=c("No_weight","ARI_weight")[1]
      no=c("L111")[1]
      which.col <- grep(colnames(usr.clrs),pattern = paste0(weight,"_",no))
      clrs.new <- usr.clrs[,which.col]
      len <- sort(setdiff(unique(clrs.new),c(0)),decreasing = F)
      # u1.gates <- generate.gates.usr(data,clrs = usr.clrs)
      # 
      # if (only.col==T)
      # {
      #   tmp <- lapply(len, function(l1) points(data[clrs.new==l1,],pch=".",col=l1,cex=2))
      # }else{
      # tmp <- lapply(u1.gates, function(g1) lines(g1,type="l",lwd=2))
      # }
      # plot(data,col=col,pch=".",main="100 usrs, ARI weight w clusters")
      # u2.gates <- generate.gates.usr(data,clrs = usr.clrs,weight ="ARI_weight" )
      # if (only.col==T)
      # {
      #   weight=c("No_weight","ARI_weight")[2]
      #   no=c("L100")[1]
      #   which.col <- grep(colnames(usr.clrs),pattern = paste0(weight,"_",no))
      #   clrs.new <- usr.clrs[,which.col]
      #   tmp <- lapply(len, function(l1) points(data[clrs.new==l1,],pch=".",col=l1,cex=2))
      # }else{
      # tmp <- lapply(u2.gates, function(g1) lines(g1,type="l",lwd=2))
      # }
      plot(data,col=col,pch=".", main="Event assignment")
      which.col <- grep(colnames(usr.clrs),pattern = paste0(weight,"_",no))
      clrs.new <- usr.clrs[,which.col]
      tmp <- lapply(len, function(l1) points(data[clrs.new==l1,],pch=".",col=l1,cex=2))
      ### sample for gates
      u3.gates <- tryCatch({return(generate.gates.usr(dat=data,clrs=usr.clrs,weight ="ARI_weight", no="L111" ))}, error = function(x){return(NA)})
      main.title <- "Gates"
      if(is.na(u3.gates)){
        u3.gates <- generate.gates.usr.chull(data,clrs = usr.clrs,weight ="ARI_weight", no="L111" )
        main.title <- "Convex gates (ahull fail)"
      }
      weight=c("No_weight","ARI_weight")[2]
      no=c("L111")[1]

      plot(data,col=col,pch='.', main=main.title)
      tmp <- lapply(u3.gates, function(g1) lines(g1,type="l",lwd=2))

      
      dev.off()
      return(1)
  },error = function(ex) { 
      t<- dev.list()
      if(!is.null(t)){
        graphics.off()
      }
      print(ex)
      return(ex) 
  })
},.parallel = F)
