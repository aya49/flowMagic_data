# date created: 2020-01-04
# author: alice yue
# input: 2D  csv/clr
# output: density (unnormalized) + scatterplot with no density + scatterplot with density | 400x400 ++++ plot! this is to make sure everything's ok


## set directory, load packages, set parallel ####
no_cores <- 32#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
# root <- "/home/ayue/projects/flowMagic_data"
#root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
root <- "/home/aya43/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))
future::plan(future::multisession, workers=no_cores) # for furrr


## output ####
y2_dir <- paste0(root,"/data/2D/y_2D")
y2_folds <- list_leaf_dirs(y2_dir)
folds <- c("y_2D_rough")
a <- sapply(y2_folds, function(y2_fold) dir.create(gsub('y_2D',folds,y2_fold), recursive=TRUE, showWarnings=FALSE))


## load inputs ####
y2_files <- list.files(y2_dir, recursive=TRUE, full.names=TRUE, pattern=".csv.gz")


## parameters
dimsize <- 256 # dimension of 2D outputs
overwrite <- FALSE


rotate.data <- function(data, chans=NULL, theta=NULL, min.max=F) {
    if (nrow(data)<3) {
        print("Cannot rotate a matrix with less than 3 rows, returning back the original matrix")
        return(list(data=data,theta=theta))
    } else {
        if (class(data)== "flowFrame" & !is.null(chans)) {
            if(all(is.na(exprs(data)[,1])))
                return("Cannot rotate a flowFrame with all cells NA")
            no.nas <- which(!is.na(exprs(data)[,chans[1]]))
            data.new <- exprs(data)[no.nas,chans]
            if (is.null(theta)) {
                reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
                slope <-  atan((max(data.new[,2])-min(data.new[,2]))/(max(data.new[,1])-min(data.new[,1])))
                theta <- ifelse(min.max,no = pi/2-reg.slope,yes = slope-pi/2)
            }
            data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
            exprs(data)[no.nas,chans] <- data.new
        } else {
            if (is.null(theta)) {
                reg.slope <- atan(lm(data[,1] ~ data[,2])$coefficients[2])
                slope <-  atan((max(data[,2])-min(data[,2]))/(max(data[,1])-min(data[,1])))
                theta <- ifelse(min.max,no = pi/2-reg.slope,yes = slope-pi/2)
            }

            col.names <- colnames(data)
            data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
            colnames(data)<-col.names
        }
        return(list(data=data,theta=theta))
    }
}

## START ####
fe <- 1:length(y2_files)
fe <- fe[grepl("HIPCbcell[/]CD10CD27[_]CD19[+]CD20[+][_]", y2_files[fe]) |
             grepl("HIPCbcell[/]CD10CD38[_]IGM", y2_files[fe]) |
             grepl("HIPCbcell[/]CD19CD20[_]CD19[+]Bcells[_]", y2_files[fe]) | 
             grepl("HIPCbcell[/]CD19SSCA[_]Nongranulocytes", y2_files[fe])]
loop_ind <- loop_ind_f(fe, no_cores)

y2 <- as.matrix(data.table::fread(y2_files[1], data.table=FALSE))
allind <- which(y2<10, arr.ind=TRUE)
y20 <- matrix(0, dimsize, dimsize)

start <- Sys.time()

cat("out of",length(fe),"\n")
a <- furrr::future_map(loop_ind, function(ii) {
    # plyr::l_ply(ii, function(i) {# tryCatch({
    for (i in ii) { tryCatch({
        # i <- grep("HIPCbcell[/]CD10CD38[_]IGM", y2_files)[1]
        y2_file <- y2_files[i]
        
        # load csv
        y2 <- as.matrix(data.table::fread(y2_file, data.table=FALSE))
        background <- y2==0
        # gplots::heatmap.2(y2, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
        
        if (grepl("HIPCbcell[/]CD10CD27[_]CD19[+]CD20[+][_]", y2_file)) {
            cpop1 <- which(y2==1, arr.ind=TRUE)
            thres1r <- min(cpop1[,1])
            
            y2[thres1r:dimsize,] = 1
            y2[1:(thres1r-1),] = 2
        } else if (grepl("HIPCbcell[/]CD10CD38[_]IGM", y2_file)) {
            cpop1 <- which(y2==1, arr.ind=TRUE)
            thres1r <- min(cpop1[,1])
            thres2c <- min(cpop1[,2])
            
            y2[thres1r:dimsize,thres2c:dimsize] = 1
            y2[y2!=1] = 2
        } else if (grepl("HIPCbcell[/]CD19CD20[_]CD19[+]Bcells[_]", y2_file)) {
            cpop1 <- which(y2==2, arr.ind=TRUE)
            thres1c <- max(cpop1[,2])
            
            y2[,1:thres1c] = 2
            y2[,(thres1c+1):dimsize] = 1
            y2_ <- y2
            y2[background] <- 0
        } else if (grepl("HIPCbcell[/]CD19SSCA[_]Nongranulocytes", y2_file)) {
            cpop1 <- which(y2==1, arr.ind=TRUE)

            theta0 = -.5 # the original angle was found with min.max
            rota <- rotate.data(allind, theta=theta0)$data
            rot <- rotate.data(cpop1, theta=theta0)$data
            thres1r <- min(rot[,1])
            
            y2 <- y20
            y2[allind[rota[,1]>=thres1r,, drop=FALSE]] <- 1
        } else if (grepl("HIPCbcell[/]CD27IgD[_]CD10[-][_]", y2_file)) {
            cpop1 <- which(y2==3, arr.ind=TRUE)
            thres1r <- max(cpop1[,1])
            thres1c <- min(cpop1[,2])
            
            y2 <- matrix(1, dimsize, dimsize)
            y2[thres1r:dimsize,1:thres1c] = 3
            y2[thres1r:dimsize,(thres1c+1):dimsize] = 2
            y2[1:(thres1r-1),1:thres1c] = 4
        } else if (grepl("HIPCbcell[/]CD34SSCA[_]Livecells[_]", y2_file)) {
            
        } else if (grepl("HIPCbcell[/]CD38CD27[_]CD19[+]Bcells[_]", y2_file)) {
            
        } else if (grepl("HIPCbcell[/]CD38CD138[_]CD19[+]Bcells[_]", y2_file)) {
            
        } else if (grepl("HIPCbcell[/]CD66CD14[_]Livecells[_]", y2_file)) {
            
        } else if (grepl("HIPCbcell[/]FSCASSCA[_]Singlets", y2_file)) {
            
        } else if (grepl("HIPCbcell[/]IgDIgM[_]CD19[+]Bcells[_]", y2_file)) {
            
        } else if (grepl("HIPCbcell[/]viabilitydyeSSCA[_]Allcells", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD3gd[_]CD3[+]Tcells[_]", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD3SSCA[_]HLADR[-]CD14[-][_]", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD11bCD16[.]FITCA[_]Granulocytes[_]", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD14CD16[.]FITCA[_]HLADR[+]CD14[+]Monocytes[_]", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD16[.]FITCACD56[_]CD3[-]", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD16[.]FITCACD56[_]CD3[+]Tcells[_]", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD64CD11b[_]CD11b[+]CD16[+]MatureNeutrophils", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD66CD16[_]Livecells", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD123CD11c[_]HLADR[+]CD14[-][_]", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]CD123HLADR[_]CD56[-]CD16[-]cells", y2_file)) {
            
        } else if (grepl("HIPCmyeloid[/]HLADRCD14[_]CD45[+]CD66[-][(]NonGranulocytes[)]", y2_file)) {
            
        } else if (grepl("pregnancy[/]01_CD66CD45[_]leukocyte", y2_file)) {
            
        } else if (grepl("pregnancy[/]02[_]CD3CD19[_]mononuclear[_]", y2_file)) {
            
        } else if (grepl("pregnancy[/]03[_]CD14CD7[_]NKLinNeg[_]", y2_file)) {
            
        } else if (grepl("pregnancy[/]04[_]CD14CD16[_]lin[-][_]", y2_file)) {
            
        } else if (grepl("pregnancy[/]05[_]CD4CD8[_]Tcell[_]", y2_file)) {
            
        } else if (grepl("pregnancy[/]06[_]CD4CD45RA[_]Tcell[_]", y2_file)) {
            
        } else if (grepl("pregnancy[/]07[_]FoxP3CD25[_]CD4Tcell", y2_file)) {
            
        } else if (grepl("pregnancy[/]08[_]TCRgdCD3[_]CD4Tcell", y2_file)) {
            
        } else if (grepl("pregnancy[/]09[_]CD8CD45RA[_]notCD4CD8Tcell[_]", y2_file)) {
            
        } else if (grepl("pregnancy[/]10[_]TbetCD45RA2[_]CD8Tcell[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]01[_]CD5CD11b[_]CD11b[+]lymphocyte[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]02[_]Ly6CCD11b[_]notgranulocyte[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]03[_]CD11bSSCh[_]notmonocyte[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]04[_]CD161CD19[_]noteosinophils[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]05[_]CD5CD11b[_]CD161[+]", y2_file)) {
            
        } else if (grepl("sangerP2[/]06[_]CD11bLy6C[_]NK[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]07[_]CD11bLy6CT[_]NKTcell[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]08[_]mhciiCD5[_]CD161[-]", y2_file)) {
            
        } else if (grepl("sangerP2[/]09[_]CD19CD11c[_]notTcell", y2_file)) {
            
        } else if (grepl("sangerP2[/]10[_]CD11bmhcii[_]cDC[_]", y2_file)) {
            
        } else if (grepl("sangerP2[/]11[_]CD5CD21[_]Bcell", y2_file)) {
            
        } else if (grepl("sangerP2[/]12[_]CD23CD21[_]B2Bcell", y2_file)) {
            
        }
        y2_ <- y2
        y2[background] <- 0
        
        write.table(y2, file=gzfile(y2_file),
                    col.names=FALSE, row.names=FALSE, sep=",")
        write.table(y2_, file=gzfile(gsub("y_2D",folds,y2_file)),
                    col.names=FALSE, row.names=FALSE, sep=",")
    }) }
    return()
})
time_output(start)

