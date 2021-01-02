# date created: 2020-12-31
# author: alice yue
# input: HIPC (sebastiano) gating sets x Bcell + Myeloid panels
# output: cleaned FCS csv files (cell x marker) + their clr (cell x cell pop)


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 4#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
libr(c(
  "furrr", #"rslurm",
  "flowWorkspace", "flowDensity",
  "stringr", 
  "pracma", "quantmod",
  "colorspace"
))


## input ####

# results from gating
gs_dir <- "/mnt/f/FCS data/Tobias Kollmann/shiny_app_data/B cells Gambia"
# gs_dir <- "/mnt/f/FCS data/Tobias Kollmann/shiny_app_data/Myeloid Gambia"


## ouput ####
out_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/data"
x_dir <- paste0(out_dir,"/nD/x/HIPCbcell"); dir.create(x_dir, recursive=TRUE, showWarnings=FALSE)
y_dir <- paste0(out_dir,"/nD/y/HIPCbcell"); dir.create(y_dir, recursive=TRUE, showWarnings=FALSE)
plot_dir <- paste0(out_dir,"/nD/scatterplots/HIPCbcell"); dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
x2_dir <- paste0(out_dir,"/2D/x/HIPCbcell"); dir.create(x2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr
y2_dir <- paste0(out_dir,"/2D/y/HIPCbcell"); dir.create(y2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr


## gating names ###

leaf_cpops <- c(
  "monocyte","granulocyte","eosinophil", 
  "NK immature Ly6c+", "NK mature Ly6c+", "NK immature Ly6c-", "NK mature Ly6c-", 
  "NKT CD11b-Ly6C+", "NKT CD11b+Ly6C+", "NKT CD11b-Ly6C-", "NKT CD11b+Ly6C-", 
  "Tcell",#"Tcell Ly6C+", "Tcell Ly6C-",
  "not Tcell", 
  "cDC CD8+", "cDC CD11b+", 
  "Bcell B1", "Bcell B2 preB", "Bcell B2 MZB", "Bcell B2 folB")

# Live cells > Non granulocytes, Granulocytes*, Blasts*
# Non granulocytes > CD19+ B cells > 
#   plasma cells*, Plasmablasts*, 
#   IGD-IGM- B cells*, IGD+IGM- B cells*, IGD-IGM+ B cells*, IGD+IGM+ B cells*, 
#   CD19+CD20+, CD19+CD20-*,
#   IGM > Immature transition B cells*
#   
rootc <- "Live cells"
root2D <- "Singlets"


## load inputs ####
gs_folders <- list.dirs(gs_dir, full.names=TRUE, recursive=TRUE)
gs_folders <- gs_folders[grepl("GatingSet_files",gs_folders)]

fcs_files <- list.files(gs_dir, full.names=TRUE, recursive=TRUE, pattern="fcs$")

fcs <- flowCore::read.FCS(fcs_files[1])
markers <- fcs@parameters@data$desc
markers[is.na(markers)] <- fcs@parameters@data$name[is.na(markers)]
names(markers) <- fcs@parameters@data$name


## gating ####
saveandrm <- TRUE # set to FALSE to check, TRUE to run through everything
start <- Sys.time()

loop_ind <- loop_ind_f(seq_len(length(gs_folders)), no_cores)
res <- furrr::future_map(loop_ind, function(ii) { 
  gs <- flowWorkspace::load_gs(gs_folders[ii])
  # gs_get_pop_paths(gs, path = "full") # all cell poulations
  # pData(gs) # no metadata
  
  # live cells
  # gsl <- flowWorkspace::gs_pop_get_gs(gs, node="Live cells")
  
  # nodes
  leaves <- flowWorkspace::gs_get_leaf_nodes(gs)
  leaves <- leaves[grepl(rootc,leaves)]
  leaves_short <- sapply(stringr::str_split(leaves, "/"), function(x) x[length(x)])
  
  # parent nodes
  parents <- unique(unlist(sapply(stringr::str_split(leaves, "/"),
                                  function(x) x[-length(x)])))
  parents_short <- parents[which(parents==root2D):length(parents)]
  
  # fcs file names in the gating set; we'll upload the cleaned fcs file to be safe
  fids <- rownames(flowWorkspace::pData(gs))
  
  res1 <- purrr::map(seq_len(length(gs)), function(i) {
    fcs <- flowCore::read.FCS(fcs_files[grepl(fids[i], fcs_files)])
    gh <- gs[[i]]
    
    # can we directly use this clr (T/F cell x cell population matrix)?
    il <- lapply(leaves, function(leaf) flowWorkspace::gh_pop_get_indices(gh, leaf))
    il <- Reduce(cbind, il)
    colnames(il) <- leaves_short
    
    
    ## Gating #####
    png(file=paste0(plot_dir, "/", fids[i], ".png"), width=2200, height=1800)
    par(mfrow=c(4,5),mar=(c(5,5,4,2)+0.1))
    
    # for each nonleaf cell population
    for (prnt in parents_short) {
      prnt_id <- which(flowWorkspace::gh_pop_get_indices(gh, prnt))
      
      # get each of its child cell population
      chlds <- flowWorkspace::gs_pop_get_children(gh, prnt)
      chlds_short <- sapply(stringr::str_split(chlds, "/"), function(x) x[length(x)])
      chlds_gates <- lapply(chlds, function(x) flowWorkspace::gh_pop_get_gate(gh, x))
      names(chlds_gates) <- names(chlds)
      
      # make fcs for plotting
      fcs_temp <- fcs
      fcs_temp@exprs[!seq_len(nrow(fcs_temp@exprs))%in%prnt_id,] <- NA

      # get the marker pairs on which each of the children were gated on
      dims <- matrix(sapply(chlds_gates, function(x) names(x@parameters)), nrow=2)
      dims <- t(dims)
      dims_id <- apply(dims, 1, function(x) paste0(x, collapse="_"))
      
      # for each marker pair, get csv and clr; 
      for (dim_id in unique(dims_id)) {
        plt_id <- which(dims_id==dim_id)
        
        mrk2 <- dims[plt_id[1],]
        gates <- chlds_gates[plt_id] # for plot
        
        csv_prnt <- fcs@exprs[prnt_id,mrk2]
        clr_chld <- sapply(plt_id, function(x) 
          prnt_id %in% which(flowWorkspace::gh_pop_get_indices(gh, chlds[x])) )
        colnames(clr_chld) <- chlds_short[plt_id]
        
        no_cp <- rowSums(clr_chld)==0
        if (sum(no_cp)!=0) {
          cbind(clr_chld, no_cp)
          colnames(clr_chld)[ncol(clr_chld)] <- "other"
        }
        
        gate_id <- paste0(paste0(markers[mrk2],collapse=""),"_",gsub("[ ]","",prnt))
        
        # save
        flowDensity::plotDens(
          fcs_temp, mrk2, 
          main=paste0("*2D* ",gate_id,"\n",paste0(colnames(clr_chld), collapse=", ")), cex.lab=2, cex.axis=2, cex.main=2)
        for (gate in gates) lines(gate@boundaries, lwd=2)
        
        dir.create(paste0(x2_dir,"/",gate_id), showWarnings=FALSE)
        write.csv(csv_lin, file=gzfile(paste0(x2_dir,"/",gate_id,"/",fids[i],".csv.gz")), row.names=FALSE)
        dir.create(paste0(y2_dir,"/",gate_id), showWarnings=FALSE)
        write.csv(clr_lin, file=gzfile(paste0(y2_dir,"/",gate_id,"/",fids[i],".csv.gz")), row.names=FALSE)
      }
    }
    graphics.off()
  })
  time_output(start1)
  
}) 

time_output(start)