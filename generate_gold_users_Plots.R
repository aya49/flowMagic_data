## Original code by Daniel Yokosawa
## 2020-10-27 gold standard plots

root <- "/mnt/f/Brinkman group/COVID"
setwd(root)

generate.gates.gold <- function(dat, clrs)
{
  gates <- lapply(1:ncol(clrs),function(i1){
    X <- dat[which(clrs[,i1]==TRUE),]
    filter <- chull(X)
    filter <- X[c(filter,filter[1]),]
    return(filter)
  })
}

generate.gates.usr <- function(
  dat, clrs, weight=c("No_weight","ARI_weight")[1], no=c("L100")[1])
{
  which.col <- grep(colnames(clrs),pattern=paste0(weight,"_",no))
  clrs.new <- clrs[,which.col]
  len <- sort(setdiff(unique(clrs.new),c(0)),decreasing=F)
  gates <- lapply(len,function(i1){
    X <- dat[which(clrs.new==i1),]
    filter <- chull(X)
    filter <- X[c(filter,filter[1]),]
    return(filter)
  })
}


file.info <- as.character(
  read.table("data/code/support_files/files_to_analyze_test_group.csv",
             sep=",",header=T,check.names=F)[,1])
res.dir <- "Test_and_result/test_users/"
data.dir <- "data/structure_test/data/"
user.res.dir <- "data/structure_test/results/"
gold.dir <- "data/structure_test/golden_samples/"
only.col<-F
which.one <- ifelse(only.col==T,yes="color",no="gate")
colPalette <- colorRampPalette(c("white", "black"))
for (f1 in file.info) { try({
  print(f1)
  print(file.exists(paste0(data.dir,f1,".csv")))
  next
  # missing / "data/structure_test/data/"
  data <- read.table(paste0(data.dir,f1,".csv"),header=T,sep=",",check.names=F)
  
  # "data/structure_test/golden_samples/"
  gold.clrs <- tryCatch(read.table(paste0(gold.dir,f1,".csv"),header=T,sep=",",check.names=F),error=function(ex) {return(NA)})
  if(is.na(gold.clrs))
    gold.clrs <- read.table(paste0("data/structure_test/clr/",f1,".csv"),header=T,sep=",",check.names=F)
  
  # missing / "data/structure_test/results/"
  usr.clrs <- read.table(
    paste0(user.res.dir, gsub(f1,pattern="/",replacement='\u20A0'),"_results.csv"), header=T,sep=",",check.names=F)
  
  # colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
  col <- densCols(data, colramp=colPalette)
  png(paste0(res.dir,gsub(f1,pattern="/",replacement='\u20A0'),which.one,"_.png"), width=2300,height=700,res=150)
  par(mfrow=c(1,4))
  
  # 
  plot(data,col=col,pch=".",main="Gold")
  g.gates <- generate.gates.gold(data,clrs=gold.clrs)
  tmp <- lapply(g.gates, function(g1) lines(g1,type="l",lwd=2))
  plot(data,col=col,pch=".",main="100 usrs, no weight, w clusters")
  weight=c("No_weight","ARI_weight")[1]
  no=c("L100")[1]
  which.col <- grep(colnames(usr.clrs),pattern=paste0(weight,"_",no))
  clrs.new <- usr.clrs[,which.col]
  len <- sort(setdiff(unique(clrs.new),c(0)),decreasing=F)
  u1.gates <- generate.gates.usr(data,clrs=usr.clrs)
  
  if (only.col==T)
  {
    tmp <- lapply(len, function(l1) points(data[clrs.new==l1,],pch=".",col=l1,cex=2))
  }else{
    tmp <- lapply(u1.gates, function(g1) lines(g1,type="l",lwd=2))
  }
  plot(data,col=col,pch=".",main="100 usrs, ARI weight w clusters")
  u2.gates <- generate.gates.usr(data,clrs=usr.clrs,weight ="ARI_weight" )
  if (only.col==T)
  {
    weight=c("No_weight","ARI_weight")[2]
    no=c("L100")[1]
    which.col <- grep(colnames(usr.clrs),pattern=paste0(weight,"_",no))
    clrs.new <- usr.clrs[,which.col]
    tmp <- lapply(len, function(l1) points(data[clrs.new==l1,],pch=".",col=l1,cex=2))
  }else{
    tmp <- lapply(u2.gates, function(g1) lines(g1,type="l",lwd=2))
  }
  plot(data,col=col,pch=".",main="100 usrs, ARI weight no clusters")
  u3.gates <- generate.gates.usr(data,clrs=usr.clrs,weight ="ARI_weight", no="L111" )
  if (only.col==T){
    weight=c("No_weight","ARI_weight")[2]
    no=c("L111")[1]
    which.col <- grep(colnames(usr.clrs),pattern=paste0(weight,"_",no))
    clrs.new <- usr.clrs[,which.col]
    tmp <- lapply(len, function(l1) points(data[clrs.new==l1,],pch=".",col=l1,cex=2))
  }else{
    tmp <- lapply(u3.gates, function(g1) lines(g1,type="l",lwd=2))
  }
  dev.off()
  
}) }

# 2D matrices
data_files <- list.files(data.dir, recursive=TRUE, pattern="csv$", ignore.case=TRUE)
clr_dir <- "data/structure_test/clr/"
clr_files <- list.files(clr_dir, recursive=TRUE, pattern="csv$", ignore.case=TRUE)
