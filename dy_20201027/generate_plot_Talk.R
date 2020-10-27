rm(list=ls())
library(alphahull)
library(igraph)
library(plyr)
library(doMC)
registerDoMC(8)

generate.gates.gold <- function(dat, clrs)
{
  gates <- lapply(1:ncol(clrs),function(i1){
    X <- dat[which(clrs[,i1]==TRUE),]
    filter <- chull(X)
    filter <- X[c(filter,filter[1]),]
    return(filter)
  })
}

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

#generate.gates.usr <- function(dat, clrs, weight=c("No_weight","ARI_weight")[1], 
#                               no=c("L111")[1])
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

#generate.gates.usr.chull <- function(dat, clrs, weight=c("No_weight","ARI_weight")[1], 
#                                     no=c("L111")[1])

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

#file.info <- as.character(read.table("~/data/TALK/Razzi_files_to_analyze.csv",
#                                     sep=",",header=T,check.names=F)[,1])

#file.info <- as.character(read.table("~/data/TALK/files_to_analyze_test.csv",
#                                     sep=",",header=T,check.names=F)[,1])

file.info <- as.character(read.table("~/data/TALK/examples_ARI/gs_ARI_poor_examples.csv",
                                    sep=",",header=T,check.names=F)[,1])
res.dir <- "~/data/TALK/plots/poor/"
dir.create(res.dir)
data.dir <- "~/data/data/structure_test/data/"
user.res.dir <- "~/data/data/structure_test/results/"
gold.dir <- "~/data/data/structure_test/golden_samples/"

only.col<-F
which.one <- ifelse(only.col==T,yes = "color",no="gate")
grbg <- llply(file.info[1:length(file.info)], function(f1){
  tryCatch({
    print(f1)
    data <- read.table(paste0(data.dir,f1,".csv"),header = T,sep=",",check.names = F)
   
    usr.clrs <- read.table(paste0(user.res.dir,
                                  gsub(f1,pattern = "/",replacement = '\u20A0'),"_results.csv")
                           ,header = T,sep=",",check.names = F)
    
    gold.clrs <- tryCatch(read.table(paste0(gold.dir,f1,".csv"),header = T,sep=","
                                     ,check.names = F),error=function(ex) {return(NA)})
    colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
    col <- densCols(data, colramp = colPalette)
    no=c("all")[1]
    
    png(paste0(res.dir,gsub(f1,pattern = "/",replacement = '\u20A0'),
               which.one,"_",no,"_.png"),
        width=1800, height=600)
    par(mfrow=c(1,3))
    
    weight=c("No_weight","ARI_weight")[1]
    which.col <- grep(colnames(usr.clrs),pattern = paste0(weight,"_",no))
    clrs.new <- usr.clrs[,which.col]
    len <- sort(setdiff(unique(clrs.new),c(0)),decreasing = F)
    #Original plot
    plot(data,col=col,pch='.', main="Original")
    
    #Gold standard
    plot(data,col=col,pch='.', main="Gold standard")
    g.gates <- generate.gates.gold(data,clrs = gold.clrs)
    tmp <- lapply(g.gates, function(g1) lines(g1,type="l",lwd=2))
    
    #Average user
    plot(data,col=col,pch=".", main="Average users result")
    which.col <- grep(colnames(usr.clrs),pattern = paste0(weight,"_",no))
    clrs.new <- usr.clrs[,which.col]
    tmp <- lapply(len, function(l1) points(data[clrs.new==l1,],pch=".",col=l1,cex=2))
    
    

    
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
