
library(flowCore)
library(flowDensity)
library(FlowSOM)
library(flowPeaks)
library(flowMeans)
Make.FCS<- function(markers, data, f.guid=NULL)
{
  pd <- c()  # 'params' phenoData
  des <- list()  # 'description' list
  
  des[["$DATATYPE"]] <- "F"
  min.instrument <- -111
  max.instrument <- 262143
  for (c in 1:ncol(data)) {
    c_name <- colnames(data)[c]
    c_marker<-markers[c]
    c_min <- floor(min(c(data[,c], min.instrument),na.rm = T))
    c_max <- ceiling(max(c(data[,c],max.instrument),na.rm = T))
    c_rng <- c_max - c_min + 1
    
    pl <- matrix(c(c_name, c_marker, c_rng, c_min, c_max),nrow=1)
    colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
    rownames(pl) <- paste("$P",c,sep="") 
    pd <- rbind(pd, pl)
    
    des[[paste("$P",c,"B",sep="")]] <- "32";      # Number of bits
    des[[paste("$P",c,"R",sep="")]] <- toString(c_rng); # Range
    des[[paste("$P",c,"E",sep="")]] <- "0,0";      # Exponent
    des[[paste("$P",c,"N",sep="")]] <- c_name;	    # Name
    des[[paste("$P",c,"S",sep="")]] <- c_marker;	    # Desc	
  }
  frame<-flowFrame(data, as(data.frame(pd), "AnnotatedDataFrame"), description=des)
  if (!is.null(f.guid))
  {
    frame@description$GUID<- f.guid
  }
  return(frame)
}

args <- commandArgs(trailingOnly = TRUE)
ind1 <- as.numeric(args[1])
print(ind1)
in.dir <-"/project/COVID/data/structure_test/data/" 
res.dir <- "/project/COVID/data/clusters2/"
# info <- vector(mode="character", length=3)
# names(info) <- c("filepath","clusters","clusters_excluded")
# write.table(info, file = "~/project/Brinkman group/COVID/HIPC_EPIC/CLR_data_fully_gated_plots/clustersinfo.csv",
#             sep=",",row.names=F)
info <- read.table("/project/COVID/data/Folders_gating_info.csv"
                   ,sep=",",header=T, check.names = F)
ungated <- as.vector(info[which(!info[,2]),1])
g1<-16
# path <- read.csv("/project/COVID/data/code/index.csv")
#dirs <- list.dirs("~/project/COVID/data/structure_test/clr",
#                  full.names = T,recursive = T)
new.name <- list.files(paste0(in.dir,ungated[g1]),full.names = T,
                       pattern = ".csv")[ind1]
true.clusters.n1 <- c("CD10-CD38+","CD14CD66-","CD19+",
                      "CD27-IgD-", "CD34+SSCA-","CD38+CD27+" ,
                      "CD38+CD138+", "IgD-IgM-",
                      "CD16+CD56+", "CD66-CD45+",
                      "CD123-CD11c-","CD123-HLADR-")

true.clusters.n2 <-  c(2,2,2,4,2,2,2,4,4,3,3,2)

print(new.name)
f.name <- basename(new.name)
dat <- as.matrix(read.table(new.name,sep=",",header=T, check.names = F))[,1:2]
if (nrow(dat)<2)
{
  print("data doesn't have enough events, something went wrong.")
}else{
  na.inds <- which(is.na(dat[,1]))
  if (length(na.inds)>0)
    dat <- dat[-na.inds,]
  ff <- Make.FCS(dat, markers = c("M1","M2"),f.guid =f.name )
  
  fmeans <- tryCatch(flowMeans(x = ff,varNames = colnames(ff)),
                     error=function(x) {return(NA)})
  if (is.na(fmeans)){
    fm.cl<- 0
  }else{
    clusters <- unique(fmeans@Label)
    fm.cl <- sum(sapply(clusters, function(c1) 
      ifelse(length(which(fmeans@Label==c1))>50,yes = 1,no = 0)))
  }
  fpeaks <-flowPeaks::flowPeaks(exprs(ff))
  clusters <- unique(fpeaks$peaks.cluster)
  fp.cl <- sum(sapply(clusters, function(c1) 
    ifelse(length(which(fpeaks$peaks.cluster==c1))>50,yes = 1,no = 0)))
  ind2 <- which( true.clusters.n1==ungated[g1])
  if (length(ind2)==1)
  {
    
    write.table(t(c(f.name, fp.cl,fm.cl,true.clusters.n2[ind2])), 
                file = paste0(res.dir,"/",ungated[g1],'\u20A0',f.name),
                col.names = c("file","flowPeaks","flowMeans","truth"
                ),sep=",",row.names = F)
  }else{
    write.table(t(c(f.name, fp.cl,fm.cl)), 
                file = paste0(res.dir,"/",ungated[g1],'\u20A0',f.name),
                col.names = c("file","flowPeaks","flowMeans"),sep=",",row.names = F)
  }
}
