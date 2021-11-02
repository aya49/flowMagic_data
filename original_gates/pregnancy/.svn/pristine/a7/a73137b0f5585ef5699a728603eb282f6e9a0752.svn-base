## Developed by Albina Rahim
## Date: November 02, 2018
## This function does the Pre-Processing of the CyTOF dataset of the Immnune Clock of Human Pregnancy 
## As input it requires the paths of the raw FCS files (Training & Validation cohorts) & information on the type of data
## It will remove FCS files:
## 1. which are Corrupted
## 2. with less than 20,000 cells
## 3. which are Duplicates
## As output it returns a large matrix- store.allFCS, which contains information of the paths, data type (Training or Validation),
## Names of the FCS files, Number of Channels, and Number of Cells.

#############################################################################################################

## preProcessingFunc <- function(inputPath, dataInfo){
##            ....
## }
## inputPath is a list which contains paths of all the datasets sent at various times
## dataInfo states if the data is Training data or Validation data

##############################################################################################################

preProcessingFunc <- function(inputPath, dataInfo){
  library("flowCore")
  library("stringr")
  
  
 
  
  store.allFCS <- NULL
  corrupted.FCS <- NULL
  lessCells.FCS <- NULL
  duplicate.FCS <- NULL

  
  
  numPaths <- length(inputPath) # Length of inputPath to determine the number of paths from where the files need to be retrieved
  
  
  for(i in 1:numPaths){
    # Path to the FCS files
    pathFCS <- unlist(inputPath[i]) 
    # Reads all folders and files in current path folder and makes a list of all of their paths
    allFCS <- dir(pathFCS, full.names=T, recursive=T, pattern = "*.fcs") 
    
    store.allFCS.temp <- sapply(1:length(allFCS), function(x){pathFCS})
    store.allFCS.temp <- cbind(store.allFCS.temp, dataInfo)
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))-1]}))
    
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(allFCS), function(x){unlist(strsplit(allFCS[x], split = "/"))[length(unlist(strsplit(allFCS[x], split = "/")))]}))
  
    
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:nrow(store.allFCS.temp), function(x){unlist(strsplit(store.allFCS.temp[x,4], split = "_"))[3]}))
    
    temp <- sapply(1:nrow(store.allFCS.temp), function(x){unlist(strsplit(store.allFCS.temp[x,4], split = "[.]"))[1]})
    store.allFCS.temp <- cbind(store.allFCS.temp, sapply(1:length(temp), function(x){unlist(strsplit(temp[x], split = "_"))[4]}))
    
    indexChange <- which(store.allFCS.temp[,5] == "Repeat")
    if(length(indexChange) != 0){
      store.allFCS.temp[indexChange,5] <- sapply(1:length(indexChange), function(x){unlist(strsplit(store.allFCS.temp[indexChange[x],4], split = "_"))[4]}) 
      store.allFCS.temp[indexChange,6] <- sapply(1:length(indexChange), function(x){unlist(strsplit(temp[indexChange[x]], split = "_"))[5]}) 
      
    }
 
  
    colnames(store.allFCS.temp) <- c("Path", "Data Info", "Folder", "FCS files", "Time Points", "Cell Stimulation")
    
   ########################################################################################################
    
    ## Code for removing duplicate FCS files based on their FCS file names 
 
      duplicate.index <- which(duplicated(store.allFCS.temp[,c('FCS files')])==TRUE)
      if(length(duplicate.index) != 0){
        duplicate.FCS <- rbind(duplicate.FCS, store.allFCS.temp[duplicate.index,])
        store.allFCS.temp <- store.allFCS.temp[!duplicated(store.allFCS.temp[,c('FCS files')]),]
      }
     
    
    ################################################################################################################
    ## Checking for files which are Corrupted and removing them. We parallelize this part of the code using ddply()
    ## We also record the number of Channels for all the files and the number of Cells.
    
    print("Start finding the Corrupted Files")
    file.names <- data.frame(store.allFCS.temp, stringsAsFactors = F)
    
    corrupted.cols.cells <- ddply(file.names, "FCS.files", function(x){
      index.Corrupted <- matrix(nrow = 1, ncol = 1, data = NA)
      NumberOfCols <- matrix(nrow = 1, ncol = 1, data = NA)
      NumberOfCells <- matrix(nrow = 1, ncol = 1, data = NA)
      
      f <- try(read.FCS(filename = paste0(x$Path, "/", x$FCS.files)), silent = TRUE)
      
      if(class(f)=="try-error"){
        index.Corrupted[1] <- "Corrupted"
      }else{
        NumberOfCols[1] <- ncol(f@exprs)
        NumberOfCells[1] <- nrow(f@exprs)
      }
      data.frame(index.Corrupted, NumberOfCols, NumberOfCells)
    }, .parallel = TRUE) # end ddply
    
    ## Arranging the output from ddply() in order of the file.names
    corrupted.cols.cells <- join(file.names, corrupted.cols.cells)
    
    ## Combinging the Number of Channels and Number of Cells with the temporary storage matrix.
    store.allFCS.temp <- cbind(store.allFCS.temp, corrupted.cols.cells[,c('NumberOfCols')], corrupted.cols.cells[,c('NumberOfCells')])
    colnames(store.allFCS.temp) <- c("Path", "Data Info", "Folder", "FCS files", "Time Points", "Cell Stimulation", "Number of Channels", "Number of Cells")
    
    # Locating the indices of the Corrupted files
    index.Corrupted <- which(!is.na(corrupted.cols.cells[,c('index.Corrupted')]))
    
    
    ## Storing the information for the Corrupted files, so we can send the information to Centre for resending these files through flowRepository
    ## Removing the Corrupted FCS files from store.allFCS.temp
    if(length(index.Corrupted) != 0){
      corrupted.FCS <- rbind(corrupted.FCS, store.allFCS.temp[index.Corrupted,])
      store.allFCS.temp <- store.allFCS.temp[-index.Corrupted,]
    }
    
    print("End of finding the Corrupted Files and removing them from the stored matrix")
    
    
    ##########################################################################################################
    ## Checking for files which has less than 20,000 cells and storing information for such files separately and then removing them from the main storage matrix
    index.lessCells <- 0
    index.lessCells <- which(as.numeric(store.allFCS.temp[,c('Number of Cells')]) < 20000)
    index.lessCells <- index.lessCells[index.lessCells !=0]
    
    if(length(index.lessCells) != 0){
      lessCells.FCS <- rbind(lessCells.FCS, store.allFCS.temp[index.lessCells,])
      store.allFCS.temp <- store.allFCS.temp[-index.lessCells,]
    }
    
    
    ########################################################################################################
    
 
    ## Combining the Paths, Name of the Files, Genotypes, Barcodes, Assay Dates, Gender, Number of Channels, Number of Cells with each FCS file in a large matrix store.allFCS
    store.allFCS <- rbind(store.allFCS, store.allFCS.temp) 
    
  } # end of outer for-loop
  
 
  ########################################################################################################

  ## Checking if all the FCS files have the same number of Channels
  numChannels <- unique(store.allFCS[,c('Number of Channels')])
  
  return(list(store.allFCS = store.allFCS, corrupted.FCS = corrupted.FCS, lessCells.FCS = lessCells.FCS, duplicate.FCS = duplicate.FCS, numChannels = numChannels))
  
}
