library(flowCore)
library(flowDensity)
getGlobalFrame <- function(fs, sample.length=NA, all.cells=F){
  set.seed(123)
  if (is(fs, 'flowFrame')){
    return (frame)
  }
  n <- length(fs)
  sample.n <- ifelse(is.na(sample.length),n,sample.length)
  global.frame <- fsApply(fs[sample(n, sample.n)],
                          function(frame){
                            m <- nrow(frame)
                            sample.size <- ifelse (all.cells,yes = m,no =min(m,ceiling(m/sample.n )))
                            
                            exprs(frame )<- exprs(frame)[sample(m, sample.size),]
                            return (frame)
                          })
  global.frame <- as(global.frame, 'flowFrame')
  return (global.frame)
}
singlet.gate <- function(f, channels, angle=-pi/4,alpha=0.05)
{
  rot <- rotate.data(f, channels,theta = angle)$data
  gate.2 <- deGate(rot, channels[1],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .1, upper=T, alpha=alpha,verbose = F)
  singlets<- rotate.fd(flowDensity(rot,channels,position = c(F,NA),gates=c(gate.2,NA),verbose = F),angle = angle)
  return(singlets)
}

#################################################################
#Data rotation
#################################################################
rotate.data <- function(data, chans=NULL, theta=NULL,min.max=F)
{
  if(nrow(data)<3)
  {
    print("Cannot rotate a matrix with less than 3 rows, returning back the original matrix")
    return(list(data=data,theta=theta))
  }else{
    if (class(data)== "flowFrame" & !is.null(chans))
    {
      if(all(is.na(exprs(data)[,1])))
        return("Cannot rotate a flowFrame with all cells NA")
      no.nas <- which(!is.na(exprs(data)[,chans[1]]))
      data.new <- exprs(data)[no.nas,chans]
      if (is.null(theta))
      {
        reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
        slope <-  atan((max(data.new[,2])-min(data.new[,2]))/(max(data.new[,1])-min(data.new[,1])))
        theta <- ifelse(min.max,no = pi/2-reg.slope,yes = slope-pi/2)
      }
      data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
      exprs(data)[no.nas,chans] <- data.new
    }else{
      col.names<- colnames(data)
      data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
      colnames(data)<-col.names
    }
    return(list(data=data,theta=theta))
  }
}

files <- list.files("~/project/COVID/data/not_gated/FR-FCM-Z2XC/",pattern = ".fcs",full.names = T)
fs <- read.flowSet(files)
fsApply()
#-----------------------------------Transformation-----------------------------------
#-------------------------------------------------------------------------------------

#tl <- transformList(colnames(fs[[1]]),arcsinhTransform(transformationId="cytofTransform",a=0,b=(1/5),c=0))
#fs<-fsApply(fs, transform, tl)
lgl<-estimateLogicle(getGlobalFrame(fs),channels=colnames(fs)[2:73])
fs <- fsApply(fs,function(x) return(transform(x,lgl)))
channels <- setdiff(colnames(fs),colnames(fs)[c(1,2,5:11,56:73)])
pairs <- combn(channels,2)
y <- 1
apply(pairs,2,function(x){
  
  png(paste0("~/project/COVID/data/not_gated/FR-FCM-Z2XC/plots/",gsub(identifier(fs[[y]]),pattern = ".fcs",
                                                                      replacement = paste0(x[1],"_",x[2], ".png"))
             ),width = 800,height = 800)
  plotDens(fs[[y]],main="",channels = x)
  dev.off()
})
