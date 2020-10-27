library(rjson)
result <- fromJSON(file="/mnt/f/Brinkman group/COVID/data/code/filepaths_count.json")
files_to_analyze <- c()
if(length(result$over)>0){
  files_to_analyze <- c(files_to_analyze, result$over)  
}
if(length(result$exact)>0){
  files_to_analyze <- c(files_to_analyze, result$exact)  
}
analyzed.files <- list.files("/mnt/f/Brinkman group/COVID/data/structure_test/results/")
analyzed.files <- gsub("\u20A0", "/",analyzed.files)
idx <- which(files_to_analyze %in% analyzed.files)

if(length(idx)>0){
  files_to_analyze <- files_to_analyze[-idx]
}
files_to_analyze <- gsub(".csv", "", files_to_analyze)
write.csv(files_to_analyze, file="/mnt/f/Brinkman group/COVID/data/code/files_to_analyze.csv", row.names = F)
