rm(list=ls())
library(rjson)
print("list files")
result <- fromJSON(file="/mnt/f/Brinkman group/COVID/data/code/support_files/filepaths_count.json")
files_to_analyze <- c()
if(length(result$over)>0){
  files_to_analyze <- c(files_to_analyze, result$over)  
}
if(length(result$exact)>0){
  files_to_analyze <- c(files_to_analyze, result$exact)  
}
cop.idx <- grep("COP[1-9]/", files_to_analyze)
files_to_analyze[cop.idx] <- gsub("COP[1-9]/", "COP/",files_to_analyze[cop.idx])

analyzed.files <- list.files("/mnt/f/Brinkman group/COVID/data/structure_test/results/")
analyzed.files <- c(analyzed.files, list.files("/mnt/FCS_local3/COVID/results/"))
analyzed.files <- gsub("\u20A0", "/",analyzed.files)
analyzed.files <- gsub("_results", "",analyzed.files)
idx <- which(files_to_analyze %in% analyzed.files)

if(length(idx)>0){
  files_to_analyze <- files_to_analyze[-idx]
}
files_to_analyze <- gsub(".csv", "", files_to_analyze)

write.csv(files_to_analyze, file="/mnt/f/Brinkman group/COVID/data/code/support_files/files_to_analyze.csv", row.names = F)
print("done ")
