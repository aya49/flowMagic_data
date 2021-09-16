# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 400x400 ++++ plot! this is to make sure everything's ok

# program arguments: the user should input arguments:
# - root directory (/src/ is where scripts are stored, /data/2D where data is stored)
# - results folder name (in root/data/2D, hopefully)
# - name of the method used for recording purposes
args = commandArgs(trailingOnly=TRUE)
# run command on console: Rscript --vanilla score.R ROOT RESFOLD METHOD


## set directory, load packages, set parallel ####
no_cores <- 14#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
# root <- "/home/ayue/projects/flowMagic_data"
root <- args[1] # "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## load inputs ####
x2_files <- list.files(paste0(data_dir,"/2D/", args[2]), recursive=TRUE, full.names=TRUE, pattern=".csv.gz")


## START ####
fe <- 1:length(x2_files)

start <- Sys.time()

cat("out of",length(fe),"\n")
loop_ind <- loop_ind_f(sample(fe), no_cores)
blscore <- plyr::llply(loop_ind, function(ii) purrr::map(ii, function(i) {
    x2_file <- x2_files[i]
    cat(i," ")
    
    # load csv
    x2r <- data.table::fread(gsub("/x_2Ddiscrete/","/x_result/",x2_file), data.table=FALSE)
    x2discrete <- data.table::fread(x2_file, data.table=FALSE)
    y2i <- data.table::fread(gsub("/x_2Ddiscrete/","/y_vector_/",x2_file), data.table=FALSE)
    
    
    # result: label of each cell according to its pixel
    # i.e. convert 2D label to vector label
    y2dp <- apply(x2discrete, 1, function(xy) x2r[xy[1], xy[2]])
    
    # baseline accuracy in this dimension
    file_split <- stringr::str_split(x2_file, "/")[[1]]
    fl <- length(file_split)
    cpops <- colnames(y2)[colnames(y2)!="other"]
    cbind(
        data.table(method=args[3], dataset=file_split[fl-2], scatterplot=file_split[fl-1], cpop=cpops, train_no=0, fcs=file_split[fl], train=FALSE),
        f1_score(y2i, y2dp, silent=TRUE))
    
}), .parallel=TRUE)
time_output(start)

blscore <- unlist(blscore)
blscore <- Reduce(rbind, blscore)

write.table(blscore, file=gzfile(paste0(scores_dir,"/2D/method_",args[3],".csv.gz")),
            row.names=FALSE, sep=",")

# f1 scores; top=data, left=method, y=f1, x=scat|cpop; colour=mean_true_prop
xlab <- "(data set) scatterplot > cell population"
gtitle <- "scatterplot || cell population VS F1 (data set vs method)\n(colour=prop, fill=count)"
gtheme <- ggplot2::theme(strip.background=ggplot2::element_rect(fill="black")) + 
    ggplot2::theme(strip.text=ggplot2::element_text(color="white", face="bold")) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust = 0.5, hjust=1))
g2f1 <- ggplot2::ggplot(blscore, ggplot2::aes(
    x=reorder(scatpop, f1), y=f1)) + 
    ggplot2::geom_boxplot(outlier.shape=NA) + 
    ggplot2::facet_grid(method~dataset, scales="free_x") + 
    ggplot2::xlab(xlab) + ggplot2::ggtitle(gtitle) + gtheme

ggplot2::ggsave(filename=paste0(root,"/plots/2D/scores/f1_2D_method_",args[3],".png"), 
                plot=g2f1, dpi=600, units="in", width=18, height=8)

time_output(start)


