## generic functions ################

## input: full path
## output: file name
file_name <- function(full_path) {
  a <- strsplit(full_path,"/")[[1]]
  a[length(a)]
}

## input: full path
## output: folder path
folder_name <- function(full_path) {
  a <- strsplit(full_path,"/")[[1]]
  paste0(a[-length(a)], collapse="/")
}

## input: full path
## output: all leaf directories in path
list_leaf_dirs <- function(full_path) {
  a <- list.dirs(full_path, recursive=TRUE)
  a[grepl(a, function(x) sum(grepl(x,a))==1)]
}


#' @title Formats time into string.
#' @description Formats time into a string HH:MM:SS given time zone.
#' @param time A time variable of class \code{POSIXct}, \code{POSIXt}.
#' @return Time formatted as a string; used in \code{time_output} function.
#' @examples
#'  # NOT EXPORTED
#'  flowGraph:::tstr(Sys.time())
#'
#' @rdname tstr
tstr <- function(time) format(.POSIXct(time), "%H:%M:%S")


#' @title Outputs elapsed time.
#' @description Given a time, prints the time elapsed from that time until now.
#' @param start A time variable of class \code{POSIXct}, \code{POSIXt}.
#' @param msg A string with a message to print out after the elapsed time.
#' @return Prints to console, the time from which process
#'  started \code{start} - ended, and > time elapsed from
#'  \code{start} until now.
#' @examples
#'
#'  start <- Sys.time()
#'  flowGraph:::time_output(start,'start - now > time elapsed')
#'
#' @rdname time_output
time_output <- function(start, msg="") {
  start <- as.POSIXct(start)
  end <- Sys.time()
  time_elapsed <- difftime(end, start, units="secs")
  message(msg, ifelse(msg == "", "", ": "),
          tstr(start), "-", tstr(end), " > ", tstr(time_elapsed))
}


#' @title Prepares parallel loop indices.
#' @description \code{loop_ind_f} is a helper function that splits
#'  a vector of loop indices into a list of multiple loop indices
#'  for use in parallel processes within the flowGraph package.
#' @param x A vector of loop indices.
#' @param n An integer, or the number of vectors to split \code{x} into.
#' @return list of \code{n} vectors with elements from \code{x}.
#' @examples
#'
#'  old_loop_inds <- 1:10
#'  no_cores <- 5
#'
#'  new_loop_inds <- flowGraph:::loop_ind_f(old_loop_inds, no_cores)
#'  # future::plan(future::multiprocess)
#'  # example_indices <- furrr::future_map(new_loop_inds, function(ii) {
#'  #     purrr::map(ii, function(i) i )
#'  # s})
#'
#' @rdname loop_ind_f
loop_ind_f <- function(x, n) {
  if (n == 1) return(base::list(x))
  return(base::split(x, ceiling(seq_along(x)/ceiling(base::length(x)/n))))
}


#' @title Loads libraries
libr = function(pkgs) {
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) 
    install.packages(setdiff(pkgs, rownames(installed.packages())), verbose=F)
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
    if (!requireNamespace("BiocManager", quietly=T)) install.packages("BiocManager")
    BiocManager::install(setdiff(pkgs, rownames(installed.packages())), ask=F)
  }
  sapply(pkgs, require, character.only=T)
}


## input: vector of numbers
## output: mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


## input: two T/F vectors indicating whether each object is in a class
## output: F1 + (see data frame)
f1score <- function(tfactual, tfpred) {
  ap <- sum(tfactual & tfpred)
  ap0 <- ap==0
  
  pred <- sum(tfpred)
  actu <- sum(tfactual)
  
  if (ap0) {
    precision <- recall <- f1 <- 0
  } else {
    precision <- ap/pred
    recall <- ap/actu
    
    f1 <- ifelse(precision+recall == 0, 0, 
                 2 * precision * recall/(precision + recall)) 
  }
  data.frame(
    precision=precision, recall=recall, f1=f1, 
    true_proportion=actu/length(tfactual), 
    predicted_proportion=pred/length(tfpred),
    true_size=actu,
    predicted_size=pred)
}


## gating functions ####################

## input: 2D matrix
## output: flowDensity::plotDens plot
plot_dens <- function(df, xlab=colnames(df)[1], ylab=colnames(df)[2], ...) {
  libr(c("flowCore", "flowDensity"))
  f <- new("flowFrame")
  f@exprs <- as.matrix(df)
  flowDensity::plotDens(f, colnames(df), xlab=xlab, ylab=ylab, ...)
}


#' @title Rotates values in a fcs or matrix
rotate_fcs <- function(data, chans=NULL, theta=NULL) {
  if (class(data)!= "flowFrame" | is.null(chans)) 
    return(list(data=rotate_matrix(data, theta)$data, theta=theta))
  
  return_data <- rotate_matrix(flowCore::exprs(data)[,chans], theta)
  flowCore::exprs(data)[,chans] <- return_data$data
  return(list(data=data, theta=return_data$theta))
}

rotate_matrix <- function(m, theta=NULL) {
  if (is.null(theta)) {
    reg_slope <- atan(stats::lm(n[,1] ~ m[,2])$coefficients[2])
    theta <- pi/2 - reg_slope
  }
  m_new <- m %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2, byrow=T)
  return(list(data=m_new, theta=theta))
}


##Finds markers in the FCS file
Find.markers <- function(frame,marker.list)
{
  #Parameters:
  #*frame: a flowFrame in the flowSet
  #**marker.list: A vector of characters
  #Output:
  #*channels.ind: a vector of channels numbers in the frame  corresponding to marker.list
  channels.ind <- unlist(lapply(marker.list, function(x) {
    ind <- grep(x, frame@parameters@data[,2], ignore.case=T)
    ind_store <- ind
    if(length(ind)==0){
      warning(paste (x, "not found, check markers!"))
      return(NA)
    } else {
      if(length(ind)>1) {
        cnt <- 0
        repeat{
          cnt <- cnt + 1
          fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x," "))[cnt]))
          ind<-match(x,fs.markers)
          if (is.na(ind))
          {
            fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x,"-"))[cnt]))
            ind<-match(x,fs.markers)
            if(!is.na(ind))
              break;
          } else {
            break;
          }
          if(cnt >= 10) {
            
            if (length(ind_store) >= 2){
              ind <- ind_store[1]
              warning(paste (x, "found more than one, choosing first. Check markers!"))
            } else {
              warning(paste (x, "not found, check markers!"))
            }
            break;
          }
        }
      }
    }
    return(ind)
  }))
  names(channels.ind)<-marker.list
  #Removing NAs in the channel vector, as not all centres have all of these markers
  #Note that most centres should have Live/CD4/CD8/CD44/CD62/CD5/CD161/CD25 in their channels
  ind <- which (is.na(channels.ind))
  if (length(ind)!=0)
    channels.ind <- channels.ind[-ind]
  return(channels.ind)
}
