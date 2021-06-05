rm(list=ls())
##### Routine to get the output from the merged polygons and create more plots:
library(plyr)
library(flowDensity)
source("/code/3.2.flowPrep.R")

load("/code/support_files/l1.files.Rdata")
### Folder structure:
results.path <- "/data/results/"
base.path <- "/data/"
exprs.path <- "/data/exprs/" 
clr.path <- "/data/clr/"
data.path <- "/data/data/"
to_transfer.path <- "/data/to_transfer/"
storage <- "/mount/"

### CHEAT CODE FOR THE HIPC FILES:
hipc.pops.myeloid <- c("CD3[+]", "CD11b[+]CD16[+]", 
                       "CLR-CD16-", "CLR-gd[+]", 
                       "HLADR[+]CD14[+]", "HLADR[+]CD14[+]]]")
hipc.myeloid.folder <-"HIPC_M_Full"

hipc.pops.b <- c("CD27CD10-", "CD19CD20-")
hipc.b.folder <-"HIPC_M_Full"
hipc3.pop <- c("CCRDR")
hipc3.folder <- c("HIPC-3")


### Steps:
### 1. Find files: 1 should be replaced for an argument from slurm ####
# read result file
args <- commandArgs(trailingOnly = TRUE)
idx.to.analyze <- as.numeric(args[1])

results.files <- list.files(results.path, full.names = T, pattern = ".csv", recursive = T)
results.files <- results.files[-grep("pixel", results.files)]
# Change this to the data set want to be analyzed . R.M. - or comment it out to run for everything
results.files <- results.files[grep("Z2N5", results.files)] 

print(idx.to.analyze)

result.file.to.analyze <- results.files[idx.to.analyze]
print(result.file.to.analyze)
new.assingments <- read.csv(result.file.to.analyze, check.names = F)

#temporary, Razzi's code has multiple outputs, selecting the weighted no clusters column
idxs <- grep("ARI_weight",colnames(new.assingments))
idxs <- idxs[grep("L111", colnames(new.assingments)[idxs])]
print(colnames(new.assingments)[idxs])
new.assingments <- new.assingments[,idxs]

#get base name for all files, remve "_result" string and replace ₠ with a /
clean.base.name <- strsplit(result.file.to.analyze, "/")
clean.base.name <- clean.base.name[[1]][length(clean.base.name[[1]])]
clean.base.name <- gsub("\u20A0", "/", clean.base.name)
base.name.idx <- regexpr("_result", clean.base.name)
base.name <- substr(clean.base.name, 1, base.name.idx-1)
new.file.name <- base.name

#find the data file, we will use this to read the channels used. If they are in the filename, we'll need to remove them to find the exprs files.
data.fname <- paste0(base.name, ".csv")
print(data.fname)
data.file <- read.csv(file=paste0(data.path, data.fname), check.names = F, stringsAsFactors = F)
marker.names <- colnames(data.file)
marker.names <- gsub("<","", marker.names)
marker.names <- gsub(">","", marker.names)
marker.names <- gsub(" ", "_",marker.names)
group.markers <- paste0("_",marker.names[1], "_", marker.names[2])
marker.index <- gregexpr(group.markers, data.fname)[[1]][1]
## we'll need to fix this for level 2
exprs.fname <- data.fname
if(marker.index < 0 ){
  group.markers <- paste0("_",marker.names[2], "_", marker.names[1])
  marker.index <- gregexpr(group.markers, data.fname)[[1]][1]
  if(marker.index > 0 ){
    exprs.fname <- gsub(group.markers, "", exprs.fname)
  }
} else {
  exprs.fname <- gsub(group.markers, "", exprs.fname)
}

# if(gregexpr("/", exprs.fname)[[1]][1]>0){
#   exprs.fname <- strsplit(exprs.fname, "/")[[1]][2]  
# }

#### CHEAT CODE: ####
for(k in hipc.pops.myeloid){
  exprs.fname <- gsub(k, hipc.myeloid.folder, exprs.fname)
}
for(k in hipc.pops.b){
  exprs.fname <- gsub(k, hipc.b.folder, exprs.fname)
}
for(k in hipc3.pop){
  exprs.fname <- gsub(k, hipc3.folder, exprs.fname)
}



#sanity check #1: number of rows of data and result ####
if(length(new.assingments)==nrow(data.file)){
  #read clr file ####
  fnames.clr <- data.fname
  print(fnames.clr)
  if(file.exists(paste0(clr.path,fnames.clr))){
    clr.file <- read.csv(paste0(clr.path,fnames.clr), check.names = F)  
  }else{
    ## will need to fix this for level2.
    fnames.clr <- exprs.fname
    clr.file <- read.csv(paste0(clr.path,fnames.clr), check.names = F)  
  }
  
  #read exprs file ####
  print(exprs.fname)
  g.idx <- gregexpr('_g[1-9]', exprs.fname)[[1]][1]-1
  if(g.idx>0){
    exprs.fname <- strtrim(exprs.fname,g.idx)
    exprs.fname <- paste0(exprs.fname, ".csv")
  }
  exprs.file <- read.csv(paste0(exprs.path,exprs.fname), check.names = F)
  # FOR NOW WE WILL KEEP ALL THE MARKERS
  #markers to remove 
  # used.markers <- gsub(exprs.fname, "", base.name)
  # used.markers <- strsplit(used.markers,"_")
  # used.markers <- unlist(used.markers)
  # idx.to.remove <- grep("g[1-9]", used.markers)
  # idx.to.remove <- c(idx.to.remove, which(used.markers==""))
  # if(length(idx.to.remove)>0){
  #   used.markers <- used.markers[-idx.to.remove]
  # }

  
  # sanity check #2: clr and exprs files' number of rows ####
  if(nrow(clr.file)==nrow(exprs.file)){
    ### 2. Apply CLR to get the correct level ####  
    f.level <- exprs.file[which(clr.file[,1]=="TRUE"),]
    
    # sanity check #3: check level exprs and data file's number of rows ####
    if(nrow(f.level)==nrow(data.file)){
      n.new.pops <- max(new.assingments)
      print(n.new.pops)
      ### 3. Apply gate and Generate plot combinations ####
      for(i in 1:n.new.pops){
        markers <- colnames(f.level)
        markers <- gsub('<','', markers)
        markers <- gsub('>','', markers)
        marker.ind.list <- c(1:ncol(f.level))
        names(marker.ind.list) <- markers
        new.clr.mask <- which(clr.file=="TRUE")
        
        new.file.name.pop <- paste0(new.file.name,"_g",i)
        ### 3.1. Create one more CLR ####
        new.clr <- rep(F,nrow(exprs.file))
        new.clr.ind <- which(clr.file[,1]=='TRUE')
        new.clr.ind <- new.clr.ind[which(new.assingments==i)]
        
        ##### SAVE NEW CLR
        new.clr[new.clr.ind] <- TRUE
        write.csv(new.clr, file=paste0(clr.path,new.file.name.pop,".csv"))
        
        while(length(marker.ind.list)>0){
          f.new.pop <- f.level[which(new.assingments==i),]
          main.ind <- marker.ind.list[1]
          marker.ind.list <- marker.ind.list[-1]
          if(length(marker.ind.list)>0){
            # x.gate <- deGate(f.new.pop[,main.ind], percentile = .98, use.percentile = T)
            # x.gate2 <- deGate(f.new.pop[,main.ind], percentile = .02, use.percentile = T)
            # new.clr.ind <- new.clr.ind[which(f.new.pop[,main.ind]<x.gate)]
            # f.new.pop   <-   f.new.pop[which(f.new.pop[,main.ind]<x.gate),]
            # new.clr.ind <- new.clr.ind[which(f.new.pop[,main.ind]>x.gate2)]
            # f.new.pop   <-   f.new.pop[which(f.new.pop[,main.ind]>x.gate2),]
            for(j in 1:length(marker.ind.list)){
              # y.gate <- deGate(f.new.pop[,marker.ind.list[j]], percentile = .98, use.percentile = T)
              # y.gate2 <- deGate(f.new.pop[,marker.ind.list[j]], percentile = .02, use.percentile = T)
              # new.clr.ind <- new.clr.ind[which(f.new.pop[,marker.ind.list[j]]<y.gate)]
              # f.new.pop   <-   f.new.pop[which(f.new.pop[,marker.ind.list[j]]<y.gate),]
              # new.clr.ind <- new.clr.ind[which(f.new.pop[,marker.ind.list[j]]>y.gate2)]
              # f.new.pop   <-   f.new.pop[which(f.new.pop[,marker.ind.list[j]]>y.gate2),]
              csv.data <- f.new.pop[,c(main.ind,marker.ind.list[j])]
              colnames(csv.data) <- c(names(main.ind), names(marker.ind.list[j]))
              n.pops <- get.number.of.populations(csv.data)
              #### remomving n.pops check
              if(sum(n.pops)==2){
                
                densx <- density(csv.data[,1])
                upperxhigh <- deGate(csv.data[,1], use.upper = T, upper = T, alpha = .3)
                upperxlow <- deGate(csv.data[,1], use.upper = T, upper = F, alpha = .3)
                percxhigh <- deGate(csv.data[,1], use.percentile = T, percentile = .99)
                percxlow <- deGate(csv.data[,1], use.percentile = T, percentile = .01)
                densy <- density(csv.data[,2])
                upperyhigh <- deGate(csv.data[,1], use.upper = T, upper = T, alpha = .3)
                upperylow <- deGate(csv.data[,1], use.upper = T, upper = F, alpha = .3)
                percyhigh <- deGate(csv.data[,1], use.percentile = T, percentile = .99)
                percylow <- deGate(csv.data[,1], use.percentile = T, percentile = .01)
                if(!(abs(upperxhigh-percxhigh)<0.2 & abs(upperyhigh-percyhigh)<0.2 & abs(upperxlow-percxlow)<0.2 & abs(upperylow-percylow)<0.2)){
                  n.pops <- c(2,2)
                } else {
                  print("only one apparent population")
                }
                
              }
              if(sum(n.pops)>2){
                if(nrow(csv.data)>800){
                  new.fname <- paste0(new.file.name.pop,"_", names(main.ind), "_", names(marker.ind.list[j]))
                  if(gregexpr('/',new.fname)[[1]][1]>0){
                    folder.name <- strsplit(new.fname,'/')[[1]][1]
                    dir.create(paste0(clr.path, folder.name))
                    dir.create(paste0(to_transfer.path, folder.name))
                  }
                  #### SAVE NEW CLR FILE ####
                  # new.clr[new.clr.ind] <- TRUE
                  # write.csv(new.clr, file=paste0(clr.path, new.fname,".csv"), row.names = F)
                  #### SAVE NEW DATA FILE ####
                  write.csv(csv.data, file=paste0(to_transfer.path, new.fname, ".csv"), row.names = F)
                } else {
                  new.fname <- paste0(new.file.name.pop,"_", names(main.ind), "_", names(marker.ind.list[j]))
                  print(new.fname)
                  print("not enough events")  
                }
              } else {
                new.fname <- paste0(new.file.name.pop,"_", names(main.ind), "_", names(marker.ind.list[j]))
                print(new.fname)
                print("not enough peaks")
              }
            }
          }
        }
      }
      ### 4. Move files from succesfull files ####
      foldername <- strsplit(data.fname, "/")[[1]][1]
      dir.create(paste0(storage,'/clr/',folder.name), recursive = T)
      dir.create(paste0(storage,'/data/',folder.name), recursive = T)
      dir.create(paste0(storage,'/results/',folder.name), recursive = T)
      
      #clrs 
      t <- file.copy(paste0(clr.path,fnames.clr),paste0(storage,"/clr/",fnames.clr),recursive=T)
      if(t){
        file.remove(paste0(clr.path,fnames.clr))
      }
      #data
      t <- file.copy(paste0(data.path, data.fname),paste0(storage,"/data/",data.fname),recursive = T)
      if(t){
        file.remove(paste0(data.path, data.fname))
      }
      #results
      storage.result.fname <- gsub("/data", "/mount", result.file.to.analyze)
      t <- file.copy(result.file.to.analyze,storage.result.fname,recursive = T)
      if(t){
        file.remove(result.file.to.analyze)
      }
    }
    else{
      print(nrow(f.level))
      print(nrow(data.file))
    }
  }
}



