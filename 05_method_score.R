# date created: 2020-01-04
# author: alice yue
# input: method res vector (for all support sizes 1-5,10,15,20; data, scatterplot, sample)
# output: score matrix (f1)

# conda: R 3.6; conda install -c r r; conda install -c r r-essentials 
# conda activate r_env

## set directory, load packages, set parallel ####
no_cores <- 32#parallel::detectCores() - 5
# root <- "/home/ayue/projects/flowMagic_data"
# install.packages("devtools")
root <- "/home/aya43/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))
future::plan(future::multisession, workers=no_cores) # for furrr

## input ####
m2_dir <- paste0(results_dir,"/2D/method"); 


## load inputs ####
dc2_dirs <- list_leaf_dirs(m2_dir)
# dc2_files <- list.files(dc2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv")

dc2_dirs <- dc2_dirs[sapply(dc2_dirs, function(x) grepl("pregnancy", x) &
    (grepl("[/]unetPRETRAINmaskDICE[-]HIPCbcell[-]sangerP2[-]HIPCmyeloid[/]", x) |
        grepl("[/]unetBASEmaskDICE[/]", x) ))]

## output ####
gs_xr_ <- function(x,y) gs_xr(x,y,"scores") 
# plyr::l_ply(append(list_leaf_dirs(dc2_dir),list_leaf_dirs(dcn_dir)), function(x)
#   dir.create(gs_xr_(x,"deepCyTOF_labels"), recursive=TRUE, showWarnings=FALSE) )


# # get train samples
# trs <- gs_xr(dc2_dir,"x_2Ddensity_euclidean_rankkmed")
wi <- hi <- 10
size <- 400


dc_files <- lapply(dc2_dirs, function(x) list.files(x, full.names=TRUE))
dc_fls <- lapply(dc_files, function(dc_files_) {
    if (length(dc_files_)<=wi*hi)
        return(list(dc_files))
    li <- list()
    for (ii in 1:ceiling(length(dc_files_)/(wi*hi))) {
        num <- min(wi*hi, length(dc_files_))
        li[[ii]] <- dc_files_[1:num]
        dc_files_ = dc_files_[-c(1:num)]
    }
    return(li)
})
dc_fls <- unlist(dc_fls, recursive=FALSE)

allind <- which(matrix(0, 256,256)==0, arr.ind=TRUE)

## START ####
start <- Sys.time()
bests <- furrr::future_map_dfr(dc_fls, function(dc_fs) {
    best <- NULL
    methl <- stringr::str_split(dc_fs[1],"/")[[1]]
    methi <- which(methl=="method")
    meth <- paste0(methl[methi:(methi+2)], collapse="[/]")
    yactualfull_file <- gsub("results","raw",gsub(meth,"y", dc_fs[1]))
    
    actual <- read.csv(yactualfull_file, check.names=FALSE)
    cpops <- colnames(actual)[colnames(actual)!="other"]
    colours <- RColorBrewer::brewer.pal(max(length(cpops),3), "Dark2")
    
    scores_file <- gsub("raw","scores",
                        gsub("/y/", paste0("/", stringr::str_extract(
                            dc_fs[1], paste0(methl[(methi+1):(methi+2)], collapse="[/]"))), yactualfull_file))
    png_file <- gsub(".csv.gz",".png", gsub("scores","plots",scores_file))
    print(png_file)
    
    dir.create(folder_name(png_file), showWarnings=FALSE, recursive=TRUE)
    # dir.create(folder_name(scores_file), showWarnings=FALSE, recursive=TRUE)
    png(png_file, width=wi*size, height=hi*size)
    par(mfcol=c(hi,wi))
    for (i in seq_len(length(dc_fs))) { tryCatch({
        dc_f <- dc_fs[i]
        cat(i, " ")
        x2predicted <- read.csv(dc_f, header=FALSE)
        
        # get meta data
        # pus <- sort(unique(predicted)) # unique labels
        path_stuff <- stringr::str_split(dc_f,"/")[[1]]
        di <- which(path_stuff=="method")
        method <- path_stuff[di+1]
        shots <- as.numeric(path_stuff[di+2])
        dset <- path_stuff[di+3]
        scat <- path_stuff[di+4]
        fname <- gsub(".csv.gz","",path_stuff[length(path_stuff)])
        
        x2discrete_file <- gsub(
            "results","data", gsub(paste0(methl[methi:(methi+2)], collapse="[/]"), "x_2Ddiscrete", dc_f))
        x2discrete <- read.csv(x2discrete_file, header=FALSE)
        
        y2actual_file <- gsub("x_2Ddiscrete","y_2D",x2discrete_file)
        y2actual <- read.csv(y2actual_file, header=FALSE)
        
        ypred <- apply(x2discrete, 1, function(xy) x2predicted[xy[1], xy[2]])
        yactual <- apply(x2discrete, 1, function(xy) y2actual[xy[1], xy[2]])
        
        yactualfull_file <- gsub("results","raw",gsub(stringr::str_extract(dc_f,paste0(methl[methi:(methi+2)], collapse="[/]")),"y", dc_f))
        
        actual <- read.csv(yactualfull_file, check.names=FALSE)
        cpops <- colnames(actual)[colnames(actual)!="other"]
        cl <- length(cpops)
        
        ## score each cpop ####
        a <- plyr::ldply(seq_len(cl), function(cpopi) {
            cbind(data.frame(
                method=method,
                dataset=dset, scatterplot=scat, cpop=cpops[cpopi],
                train_no=shots, fcs=fname,
                train=NA
            ), f1score(yactual==cpopi, ypred==cpopi)) ### +1 ?????
        })
        best <- rbind(best, a)
        
        # plot
        x2predicted_ <- which(x2predicted>0, arr.ind=TRUE)
        x2predicted_c <- apply(x2predicted_, 1, function(xy) x2predicted[xy[1], xy[2]])
        
        plot(x2predicted_, main=fname, col=colours[x2predicted_c], cex=.55, pch=16)
        for (cpopi in seq_len(cl)) {
            y2D <- which(y2actual==cpopi, arr.ind=TRUE)
            if (is.null(y2D) | nrow(y2D)<3) next
            y2D_chulli <- chull(y2D)
            y2D_chull <- y2D[append(y2D_chulli, y2D_chulli[1]),]
            lines(y2D_chull, lwd=2, col=colours[cpopi])
            lines(y2D_chull, lwd=2, lty=3, col="black")
        }
        legend("topright", legend=cpops, col=colours, 
               lty=rep(1,cl), lwd=rep(2,cl))
        
    }, error = function(e) {
        print(dc_f)
    })}
    graphics.off()
    # write.table(best, file=gzfile(scores_file),
    #             sep=',', row.names=FALSE, col.names=TRUE)
})
bests <- bests[,colnames(bests)!=".id"]
score_file <- paste0(gs_xr_(m2_dir,"method"),"/SCORE_unetBASEandPRETRAINmaskDICE-HIPCbcell-sangerP2-HIPCmyeloid.csv.gz") ### ????
dir.create(folder_name(score_file), recursive=TRUE, showWarnings=FALSE)
write.table(bests, file=gzfile(score_file), sep=",", row.names=FALSE, col.names=TRUE)
time_output(start)



