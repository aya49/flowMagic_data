## Developed by Albina Rahim
## Date: October 29, 2018
## This is the main preProcessing script which calls the function:
## preProcessingFunc.R : for pre-Processing of the CyTOF datasets for the Immnune Clock of Human Pregnancy


remove(list=ls())

setwd("/home/rstudio/code/Projects/Immune_Clock_Human_Pregnancy/Codes/")



###############################################################################################
## NOTE: Run this section functions either separately before you run the rest of the script or you can run the entire script in just one go
## The helperFunc.R prompts user for input: Training data or Validation data

source("helperFunc.R")
if (interactive() ){
  dataType <- readDataFunc()
}

dataInfo <- NULL

if(dataType == "1"){
  dataInfo <- "Training"
}else if(dataType == "2"){
  dataInfo <- "Validation"
}

#############################################################################################
##############################################################################################

library("plyr")
library("doMC")


## Function for the pre-Processing of the datasets
source("preProcessingFunc.R")
## Function for creating a Global Frame
source("globalFrameFunc.R")



## This part of the script was taken from Sibyl for the purpose of parallelizing the execution of this script
## no_cores to determine how many CPUs to use while implenting the parallelization
no_cores <- detectCores() - 2
registerDoMC(no_cores)

start <- Sys.time()

# ## Creating the Results folder for the Training & Validation cohort
# suppressWarnings(dir.create ("/home/rstudio/results/Immune_Clock_Human_Pregnancy/"))

## Paths to the FCS files, metadata spreadsheets, and output folder in the Bioinformatics drive
if(dataInfo == "Training"){
    inPath <- "/home/rstudio/data/Immune_Clock_Human_Pregnancy/FR-FCM-ZY3Q"
    inputPath <- list(inPath)
    #suppressWarnings(dir.create ("/mnt/f/FCS data/Immune_Clock_Human_Pregnancy/Results/FR-FCM-ZY3Q/"))
    outputPath <- "/home/rstudio/results/Immune_Clock_Human_Pregnancy/FR-FCM-ZY3Q/"
    
}else if(dataInfo == "Validation"){
    inPath <- "/home/rstudio/data/Immune_Clock_Human_Pregnancy/FR-FCM-ZY3R"
    inputPath <- list(inPath)
    #suppressWarnings(dir.create ("/mnt/f/FCS data/Immune_Clock_Human_Pregnancy/Results/FR-FCM-ZY3Q/"))
    outputPath <- "/home/rstudio/results/Immune_Clock_Human_Pregnancy/FR-FCM-ZY3R/"
}


## Calling the preProcessing function
preProcessing.Output <- preProcessingFunc(inputPath, dataInfo)
store.allFCS <- preProcessing.Output$store.allFCS
corrupted.FCS <- preProcessing.Output$corrupted.FCS
lessCells.FCS <- preProcessing.Output$lessCells.FCS
duplicate.FCS <- preProcessing.Output$duplicate.FCS
numChannels <- preProcessing.Output$numChannels



sink(file = paste0(outputPath,"preProcessing-Summary.txt"), split = TRUE)
if(dataInfo == "Training"){
  print(paste0("There are in total ", nrow(store.allFCS), " Training files for analysis."))
}else if(dataInfo == "Validation"){
  print(paste0("There are in total ", nrow(store.allFCS), " Validation files for analysis."))
}
print(paste0("Number of Corrupted files: ", nrow(corrupted.FCS)))
print(paste0("Number of files with < 20,000 cells: ", nrow(lessCells.FCS)))
print(paste0("Number of Duplicate FCS files: ", nrow(duplicate.FCS)))
if(length(numChannels) == 1){
  print(paste0("All files have the same number of channels: ", numChannels))
} else{
  #paste(c("The first three notes are: ", notes), collapse=" ")
  print(paste0(c("Files have different number of channels:", numChannels), collapse = " "))
}
sink()

## Saving the results
save (store.allFCS, file =  paste0(outputPath,"store.allFCS.Rdata") )
if(length(corrupted.FCS) != 0){
  save(corrupted.FCS , file = paste0(outputPath, "corrupted.FCS.Rdata"))
}
if(length(lessCells.FCS) != 0){
  save(lessCells.FCS, file = paste0(outputPath,"lessCells.FCS.Rdata"))
}
if(length(duplicate.FCS) != 0){
  save(duplicate.FCS, file = paste0(outputPath,"duplicate.FCS.Rdata"))
}

cat("Total time is: ",TimeOutput(start),sep="")


############################################################################################################
## We don't need this part for CyTOF data analysis.
## For CyTOF we are using Arcsinhtransformation in the mainGate.R script
# ## Calling the Function for creating the Global Frame
# globalFrame.Output <- globalFrameFunc(store.allFCS, outputPath)
# ## Saving the Global Frame
# globalFrame <- globalFrame.Output$globalFrame
# ## Saving the Expression Matrix of the Global Frame with the information of each FCS file
# ## This matrix can later be used if we need to remove any FCS files and its corresponding expression matrix values from the global frame
# globalFrame.Matrix <- globalFrame.Output$globalFrame.Matrix
# ##Saving the computed transform from estimateLogicle()
# lgl <- globalFrame.Output$lgl
# ## Saving the channel index which will be used for gating
# Transform.index = globalFrame.Output$Transform.index
# ## Saving the list of files which had NAs in some of their channels. These files were already removed from the global frame
# files.w.NAs <- globalFrame.Output$files.w.NAs
# save(globalFrame, file = paste0(outputPath,"globalFrame.Rdata"))
# save(globalFrame.Matrix, file = paste0(outputPath, "globalFrame.Matrix.Rdata"))
# save(lgl, file = paste0(outputPath, "lgl.Rdata"))
# save(Transform.index, file = paste0(outputPath, "Transform.index.Rdata"))
# if(length(files.w.NAs) != 0){
#   save(files.w.NAs, file = paste0(outputPath, "files.w.NAs.Rdata"))
# }

