library(digest)
library(flowDensity)
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

all.files2 <- read.csv(file="/data/Brinkman group/COVID/data/code/index.csv", stringsAsFactors = F)
fname <- all.files2[index,1]
fname <- gsub("/home/rstudio", "", fname)

data.file <- read.csv(fname)
bad <- F
motive <- ""
if(ncol(data.file)==2){
  if(colnames(data.file)[1]==colnames(data.file)[2]){
    bad <- T
    motive <- "colnames"
  } else {
    col.one <- data.file[,1] 
    col.two <- data.file[,2]
    equal.vals <- which(col.one %in% col.two)
    if(length(equal.vals)>0.95*length(col.one)){
      bad <- T
      motive <- "equal columns"
    }
  }
  if(nrow(data.file)<1000){
      bad <- T
      motive <- paste0("only: ", nrow(data.file), " events")
  }
  width1 <- max(data.file[,1])-min(data.file[,1])
  upper11 <- deGate(data.file[,1], upper=T, use.upper = T, alpha = .1, tinypeak.removal = .9)
  upper12 <- deGate(data.file[,1], upper=F, use.upper = T, alpha = .1, tinypeak.removal = .9)
  width2 <- max(data.file[,2])-min(data.file[,2])
  upper21 <- deGate(data.file[,2], upper=T, use.upper = T, alpha = .1, tinypeak.removal = .9)
  upper22 <- deGate(data.file[,2], upper=F, use.upper = T, alpha = .1, tinypeak.removal = .9)
  if(abs(upper11-upper12)/width1<0.1){
    bad <- T
    motive <- paste0("thin distribution: ", colnames(data.file)[1], " dist=",abs(upper11-upper12)/width1)
  }
  if(abs(upper21-upper22)/width2<0.1){
    bad <- T
    motive <- paste0("thin distribution: ", colnames(data.file)[2], " dist=",abs(upper21-upper22)/width2)
  }
  
  
}
if(bad){
  sha.name <- sha1(fname)
  report.val <- cbind(fname, motive)
  colnames(report.val) <- c("file","motive") 
  write.csv(report.val, paste0("/data/Brinkman group/COVID/data/data_check_fail/",sha.name, ".csv"))
}
print("success")
  