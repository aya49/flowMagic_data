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


## output ####
gs_xr_ <- function(x,y) gs_xr(x,y,"scores") 
# plyr::l_ply(append(list_leaf_dirs(dc2_dir),list_leaf_dirs(dcn_dir)), function(x)
#   dir.create(gs_xr_(x,"deepCyTOF_labels"), recursive=TRUE, showWarnings=FALSE) )


# # get train samples
# trs <- gs_xr(dc2_dir,"x_2Ddensity_euclidean_rankkmed")
wi <- hi <- 10
size <- 400


## START ####
start <- Sys.time()
dc_files <- lapply(dc2_dirs, function(x) list.files(x, full.names=TRUE))
li <- list()
lii <- 0
for (dc_files_ in dc_files) {
    if (length(dc_files_) < 200) {
        lii <- lii + 1
        li[[lii]] <- dc_files
        next
    }
    il <- floor(length(dc_files_)/(wi*hi))
    for (i in seq_len(il)) {
        lii <- lii + 1
        if (i==il) {
            li[[lii]] <- dc_files_
            next
        }
        li[[lii]] <- dc_files_[1:(wi*hi)]
        dc_files_ <- dc_files_[-c(1:(wi*hi))]
    }
}

bests <- furrr::map(li, function(dc_fs) {
    best <- NULL
    
    yactualfull_file <- gsub("results","raw",gsub(stringr::str_extract(dc_fs[1],"method/[a-z]+[A-Z]*/[0-9]+"),"y", dc_fs[1]))

    actual <- read.csv(yactualfull_file, check.names=FALSE)
    cpops <- colnames(actual)[colnames(actual)!="other"]
    colours <- RColorBrewer::brewer.pal(length(cpops), "Dark2")[seq_len(length(cpops))]
    
    scores_file <- gsub("raw","scores",
                        gsub("/y/", stringr::str_extract(
                            dc_fs[1], "[/]unet[a-z]*[A-Z]*[/][0-9]+[/]"), yactualfull_file))
    png_file <- gsub(".csv.gz",".png", gsub("scores","plots",scores_file))
    
    dir.create(folder_name(png_file), showWarnings=FALSE, recursive=TRUE)
    dir.create(folder_name(scores_file), showWarnings=FALSE, recursive=TRUE)
    png(png_file, width=wi*size, height=hi*size)
    par(mfcol=c(hi,wi))
    for (dc_f in dc_fs) {
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
            "results","data", gsub(stringr::str_extract(dc_f,"method[/][a-zA-Z_]+[/][0-9]+[/]"),"x_2Ddiscrete/", dc_f))
        x2discrete <- read.csv(x2discrete_file, header=FALSE)
        
        ypred <- apply(x2discrete, 1, function(xy) x2predicted[xy[1], xy[2]])
        
        yactualfull_file <- gsub("results","raw",gsub(stringr::str_extract(dc_f,"method/[a-z]+[A-Z]*/[0-9]+"),"y", dc_f))
        
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
            ), f1score(actual[,cpopi]==1, ypred==cpopi)) ### +1 ?????
        })
        best <- rbind(best, a)
        
        # plot
        x2predicted_ <- which(x2predicted>0, arr.ind=TRUE)
        x2predicted_c <- apply(x2predicted_, 1, function(x) x2predicted[x[1],x[2]])
        y2D <- read.csv(gsub("x_2Ddiscrete","y_2D",x2discrete_file), header=FALSE)
        y2D_ <- which(y2D>0, arr.ind=TRUE)
        y2D_c <- apply(y2D_, 1, function(x) y2D[x[1],x[2]])
        
        plot(x2predicted_, main=fname, col=colours[x2predicted_c], cex=.55, pch=16)
        for (cpopi in seq_len(cl)) {
            cpii <- which(y2D_c==cpopi)
            y2D_chulli <- chull(y2D_[cpii,])
            y2D_chull <- y2D_[cpii[append(y2D_chulli, y2D_chulli[1])],]
            lines(y2D_chull, lwd=2, col=colours[cpopi])
            lines(y2D_chull, lwd=2, lty=3, col="black")
        }
        legend("topright", legend=cpops, col=colours, 
               lty=rep(1,cl), lwd=rep(2,cl))
        
    }
    graphics.off()
    write.table(best, file=gzfile(scores_file),
                sep=',', row.names=FALSE, col.names=TRUE)
})
bests <- Reduce(rbind, bests)
bests <- bests[,colnames(bests)!=".id"]
score_file <- paste0(gs_xr_(dc_dir_,"method"),".csv.gz") ### ????
dir.create(folder_name(score_file), recursive=TRUE, showWarnings=FALSE)
write.table(bests, file=gzfile(score_file), sep=",", row.names=FALSE, col.names=TRUE)
time_output(start)



