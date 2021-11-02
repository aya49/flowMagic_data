## Prompts user for Data type

readDataFunc <- function()
{
  dataType <- readline(prompt = "Enter 1 for Training data & 2 for Validation data:")
  if(dataType == "1" | dataType == "2")
  {
    return(dataType)
  }else{
    dataType <- readline(prompt = "Invalid input. Please enter again:")
    return(dataType)
  }
}


##################################################################################################################

# # # removes all cells that exist on the edges of FSC-A, FSC-H, SSC-A, and SSC-H. 
# 
# removeMargins<- function(f,chans,sens=1, debris=FALSE,return.ind=F,neg=500, verbose = T)
# {
#   neg <-cbind(chans,neg)[,2]
#   #Size is a vector of size 2, to be passed to mfrow in case of plotting
#   data <- exprs(f)
#   margins <- c()
#   marg.list <-list()
#   if(!debris)
#   {
#     for(chan in chans)
#     {
#       stain.max <-max(data[,chan])
#       margins <- which ( data[, chan] >= stain.max*sens)
#       marg.list <- c(marg.list, list(margins))
#       data <- data[ -margins, ]
#       if(verbose == T){print(paste(length(margins), "margin events in",colnames(f)[chan], "will be removed.",sep =" "))}
#     }
#     
#   }else
#   {
#     for(i in 1:length(chans))
#     {
#       stain.min <-min(data[,chans[i]])
#       margins <- which ( data[, chans[i]] <= stain.min*sens)
#       if (neg[i]<500)
#       {
#         negs <- which ( data[, chans[i]] < neg[i])
#         margins <- negs
#       }
#       marg.list <- c(marg.list, list(margins))
#       if (length(margins)!=0){
#         data <- data[ -margins, ]
#       }
#       if(verbose == T){print(paste(length(margins), "debris events in",colnames(f)[chans[i]], "will be removed.",sep =" "))}
#     }
#   }
#   exprs(f) <- data
#   if (!return.ind)
#     return(f)
#   else
#     return(list(frame=f,ind=marg.list))
#   
# }
# 

############################################################################################

rotate.data <- function(data, chans=NULL, theta=NULL)
{
  if (class(data)== "flowFrame" & !is.null(chans))
  {
    data.new <- exprs(data)[,chans]
    if (is.null(theta))
    {
      reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
      theta <- pi/2 - reg.slope
    }
    data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
    exprs(data)[,chans] <- data.new
  }else{
    data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
  }
  return(list(data=data,theta=theta))
}






#################################################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
TimeOutput(Sys.Date())


#################################################################################
# #SVD reduction for flowType to be used for grouping patients based on reduced matrix
# #Kmeans used for grouping patients, any other method can be used
# 
# #Author: Mehrnoush Malek
# #Date: April 2012
# 
# svd.reduction <- function(cell.prop,kmean.centres=3,kmeans.start=1000)
#   #cell.prop is a matrix of size 3^k by m, where k is number of markers used in flowType and m is number of samples
# {
#   
#   svdr.Result<-svd(t(cell.prop))
#   #Find a threshold where samples get far from others
#   x11();plot(sort(svdr.Result$d))
#   
#   inds<-which(svdr.Result$d > locator()$y)
#   acf<-t(cell.prop) %*% (svdr.Result$v[,inds])
#   kmeans.cluster<-kmeans(acf, centers=kmean.centres, nstart=kmeans.start)
#   x11();plot(data.frame(acf), pch=16,col=kmeans.cluster$cluster)
#   return(kmeans.cluster)
# }
# 

#######################################################################################

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


###################################################################################################


