## Developed by Albina Rahim
## Date: November 06, 2018
## This function creates GlobalFrame by storing 1000 random cells from each FCS file. 
## As input it requires the store.allFCS matrix and the outputPath (for saving the Results)


globalFrameFunc <- function(store.allFCS, outputPath){
    library("flowCore")
    library("flowBin")
    
    
    print("Start creating the Global Frame")
    #file.names <- data.frame(store.allFCS[200:nrow(store.allFCS),], stringsAsFactors = F)
    file.names <- data.frame(store.allFCS, stringsAsFactors = F)
    gFrame <- ddply(file.names, "FCS.files", function(x){
     
        f <- try(read.FCS(filename = paste0(x$Path, "/", x$FCS.files)), silent = TRUE)
        if(!is.null(f@description$SPILL)){
          f <- compensate(f, f@description$SPILL)
        }
        
        
      # markers1 <- c("BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "CD235ab_CD61", "CD45", "CD66",
      #              "CD7", "CD19", "CD45RA", "CD11b", "CD4", "CD8a", "CD11c", "CD123", "CREB", "STAT5", "p38", 
      #              "TCRgd", "STAT1", "STAT3", "S6", "CXCR3", "CD161", "CD33", "MAPKAPK2", "Tbet", "FoxP3",
      #              "IkB", "CD16", "NFkB", "ERK", "CCR9", "CD25", "CD3", "CCR7", "CD15", "CCR2", "HLADR",
      #              "CD14", "CD56", "DNA1", "DNA2", "beadDist")
      
      # Make colnames readable using information in the parameter data slot
      #markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(f@parameters@data$desc))
        
      
      
      
      #pregating_channels <- c("Bead", "DNA1", "DNA2", "Dead", "Event_length")
      gating_channels <- c("CD57", "CD19", "CD4", "CD8", "IgD", "CD11c", "CD16", "CD3", "CD38", "CD27", "CD14", "CXCR5", "CCR7", "CD45RA", "CD20", "CD127", "CD33", "CD28", "CD161", "TCRgd", "CD123", "CD56", "HLADR", "CD25", 
                           "CD235ab_CD61", "CD66", "CD45", "Tbet", "CD7", "FoxP3", "CD11b")
      #instrument_channels <- c("Time", "Event_length", "Center", "Offset", "Width", "Residual", "sample")
      
      channels.ind <-sort(Find.markers(f, gating_channels))
      
      temp <- f@exprs[sample(1:length(f@exprs[,1]), 1000), channels.ind ]
      
      return(data.frame(temp, check.names = F))
      
    }, .parallel = TRUE) # end ddply
    
    # ## Arranging the output from ddply() in order of the file.names
    gFrame <- join(file.names, gFrame)
   
    ## Remove rows in gFrame with NAs and also saving files which has NA in their flurophore/channels
    fileswNAs <- NULL
    colswNAs <- which(is.na(gFrame), arr.ind = TRUE)
    if(length(colswNAs) > 0){
      rowswNAs <- unique(colswNAs[,'row'])
      colswNAs <- unique(colswNAs[,'col'])
      fileswNAs <- unique(gFrame[rowswNAs,1])
      gFrame <- gFrame[-rowswNAs,]
    }
    
    ## Reading the first FCS file in the storage matrix as a template for creating the global frame
    
    g <- read.FCS(filename = paste0(store.allFCS[1,c('Path')], "/", store.allFCS[1,c('FCS files')]))
    
    Transform.idx <- unlist(sapply(colnames(gFrame), function(x) {grep(x, colnames(g))})) 
    gexprs.temp <- matrix(0, nrow = nrow(gFrame), ncol = ncol(g@exprs))
    gexprs.temp[, Transform.idx] <- as.matrix(gFrame[, 9:ncol(gFrame)])
    g@exprs <- gexprs.temp
    colnames(g@exprs) <- colnames(g) 
    
    print("End of creating the Global Frame")
    
    ## Tranformation
    print("Start computing the transform using the estimateLogicle()")
    lgl <- estimateLogicle(g, channels = colnames(g)[Transform.idx])
    print("End of computing the transform using the estimateLogicle()")
    
    return(list(globalFrame = g, globalFrame.Matrix = gFrame, lgl = lgl, Transform.index = Transform.idx, files.w.NAs = fileswNAs))
}



