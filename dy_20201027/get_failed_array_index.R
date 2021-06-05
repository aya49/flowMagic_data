args = commandArgs(trailingOnly=TRUE)
jobID <-  args[1]
slurm.report <- read.table("/mnt/f/Brinkman group/COVID/data/code/failed_jobs.csv", header = T, stringsAsFactors = F)
slurm.report <- slurm.report[-1,]
slurm.report <- slurm.report[-grep(".batch", slurm.report[,1]),]

jobs.to.re.run <- slurm.report[grep(jobID, slurm.report[,1]),]
jobs.to.re.run <- jobs.to.re.run[grep("127:0", jobs.to.re.run[,3]),]
array.index <- jobs.to.re.run[,1]
array.index <- lapply(array.index, FUN = function(x){
  return(strsplit(x,"_")[[1]][2])
})
array.index <- unlist(array.index)
string.to.print <- array.index[1]
for(i in 2:length(array.index)){
  string.to.print <- paste0(string.to.print, ",",array.index[i])
}
print(string.to.print)
