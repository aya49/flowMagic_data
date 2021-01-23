# date created: 2020-12-31
# author: alice yue
# input: HIPC (sebastiano) gating sets x Bcell + Myeloid panels
# output: cleaned FCS csv files (cell x marker) + their clr (cell x cell pop) + straight thresholds for some gates (for flowLearn)


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
source("/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src/RUNME.R")


## input ####

# results from gating
gs_dir <- "/mnt/f/FCS data/Tobias Kollmann/shiny_app_data/B cells Gambia"
#gs_dir <- "/mnt/f/FCS data/Tobias Kollmann/shiny_app_data/Myeloid Gambia"

panel <- "myeloid"
if (grepl("B cells",gs_dir))
  panel <- "bcell"


## ouput ####
dset <- paste0("HIPC",panel)
plyr::a_ply(
  paste0(c(xn_dir, yn_dir, plotn_dir, x2_dir, y2_dir, thres_dir, filt_dir),"/",dset),
  dir.create, recursive=TRUE, showWarnings=FALSE)


## gating names ###

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
dupm <- is.na(markers)
if (any(dupm)) 
  markers[dupm] <- paste0(fcs@parameters@data$name[dupm])
dupm <- duplicated(markers)
if (any(dupm)) 
  markers[dupm] <- paste0(markers[dupm],".",fcs@parameters@data$name[dupm])
markers <- gsub("[- ]","",markers)

# # accidentally included duplicated markers, the following is a fix
# hmfs <- list.files("/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/data/nD/x/HIPCmyeloid", full.names=TRUE)
# for (hmf in hmfs) {
#   csv_ <- data.table::fread(hmf, data.table=FALSE)
#   if (sum(duplicated(colnames(csv_)))==0) break
#   colnames(csv_)[duplicated(colnames(csv_))] <- "CD16.FITC-A"
#   write.csv(csv_, file=gzfile(hmf), row.names=FALSE)
# }

markers[is.na(markers)] <- fcs@parameters@data$name[is.na(markers)]
names(markers) <- fcs@parameters@data$name

gs <- flowWorkspace::load_gs(gs_folders[1])
pdf(file=paste0(data_dir,"/HIPC",panel,".pdf"))
plot(gs)
graphics.off()
rm(gs)

error_files <- c()

## gating ####
# saveandrm <- TRUE # set to FALSE to check, TRUE to run through everything
start <- Sys.time()

# res <- furrr::future_map(seq_len(length(gs_folders)), function(ii) { 
for (ii in seq_len(length(gs_folders))) {
  cat("\n\n",gs_folders[ii])
  start1 <- Sys.time()
  suppressMessages( gs <- flowWorkspace::load_gs(gs_folders[ii]) )
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
  
  # nodes <- flowWorkspace::gs_get_leaf_nodes(gs)
  # nodes <- nodes[grepl(paste0(rootc,"/"),nodes)]
  # all_gates <- lapply(nodes, function(x) flowWorkspace::gh_pop_get_gate(gs[[1]], x))
  # mrks <- unique(as.vector(sapply(all_gates, function(x) names(x@parameters))))
  
  
  # fcs file names in the gating set; we'll upload the cleaned fcs file to be safe
  fids <- rownames(flowWorkspace::pData(gs))
  gsl <- lapply(seq_len(length(gs)), function(x) gs[[x]])
  rm(gs)
  
  loop_ind <- loop_ind_f(seq_len(length(gsl)), no_cores)
  res <- plyr::llply(loop_ind, function(j) {
    plyr::l_ply(j, function(i) {
      # for (i in seq_len(length(gsl))) {
      # if (fids[i]%in%error_files) next
      # if (file.exists((paste0(plot_dir, "/", fids[i], ".png"))))
      #   if (file.size(paste0(plot_dir, "/", fids[i], ".png"))>500000) next
      fid <- fids[i]
      cat("\n- ",i, fid)
      # gh <- gsl[[i]]
      # fcs <- flowCore::read.FCS(fcs_files[grepl(fids[i], fcs_files)])
      
      # can we directly use clr (T/F cell x cell population matrix)?
      ir <- flowWorkspace::gh_pop_get_indices(gsl[[i]], rootc)
      il <- lapply(leaves, function(leaf) flowWorkspace::gh_pop_get_indices(gsl[[i]], leaf))
      il <- do.call(cbind, il)
      colnames(il) <- leaves_short
      other <- which(ir)%in%which(rowSums(il)==0)
      ilm <- cbind(il[ir,,drop=FALSE], other)
      
      if (panel=="bcell") {
        ilm_ <- ilm[,!colnames(ilm)%in%c(
          "Plasmablasts","plasma cells","Blasts",
          "IGD-IGM- B cells","IGD-IGM+ B cells","IGD+IGM- B cells","IGD+IGM+ B cells",
          "Immature transition B cells")]
      } else {
        ilm_ <- ilm[,!colnames(ilm)%in%"CD56+CD16+ NKT cells"]
      }
      rowi <- rowSums(ilm_)==1
      ilm_ <- ilm_[rowi,]
      ilm_[which(ilm_)] <- 1
      
      clm <- flowWorkspace::gh_pop_get_data(gsl[[i]], rootc)@exprs[,!markers%in%"Time" & !markers%in%"viability dye"]
      clm_ <- clm[rowi,]
      
      colnames(clm_) <- markers[colnames(clm_)]
      
      write.table(clm_, file=gzfile(paste0(x_dir,"/",dset,"/",fid,".csv.gz")), sep="," ,row.names=FALSE)
      write.table(ilm_, file=gzfile(paste0(y_dir,"/",dset,"/",fid,".csv.gz")), sep=",", row.names=FALSE)
      
      rm(clm,clm_,ilm,ilm_)
      
      ## Gating #####
      gthres <- list()
      
      png(file=paste0(plot_dir, "/", dset,"/",fid, ".png"), width=2200, height=1800)
      par(mfrow=c(4,5),mar=(c(5,5,4,2)+0.1))
      
      # for each nonleaf cell population
      for (prnt in parents_short) {
        prnt_id <- which(flowWorkspace::gh_pop_get_indices(gsl[[i]], prnt))
        if (length(prnt_id)==0) next
        
        # get each of its child cell population
        chlds <- flowWorkspace::gs_pop_get_children(gsl[[i]], prnt)
        chlds_short <- sapply(stringr::str_split(chlds, "/"), function(x) x[length(x)])
        chlds_gates <- lapply(chlds, function(x) flowWorkspace::gh_pop_get_gate(gsl[[i]], x))
        names(chlds_gates) <- names(chlds)
        
        # make fcs for plotting
        fcs_temp <- flowWorkspace::gh_pop_get_data(gsl[[i]], prnt)
        
        # get marker pairs on which each of the children were gated on
        dims <- matrix(sapply(chlds_gates, function(x) names(x@parameters)), nrow=2)
        dims <- t(dims)
        dims_id <- apply(dims, 1, function(x) paste0(x, collapse="_"))
        
        # for each marker pair, get csv and clr; 
        for (dim_id in unique(dims_id)) {
          plt_id <- which(dims_id==dim_id)
          
          mrk2 <- dims[plt_id[1],]
          gates <- chlds_gates[plt_id] # for plot
          
          scat <- paste0(paste0(markers[mrk2],collapse=""),"_",gsub("[ ]","",prnt))
          
          # flowLearn straight gates?
          fltf <- any(sapply(c(
            # myeloid
            "CD14[.A-Za-z]*CD16[.A-Za-z]*[_]HLADR[+]CD14[+]Monocytes",
            "CD123[.A-Za-z]*CD11c[.A-Za-z]*[_]HLADR[+]CD14[-]",
            "CD3[.A-Za-z]*SSCA[.A-Za-z]*[_]HLADR[-]CD14[-]", 
            "CD11b[.A-Za-z]*CD16[.A-Za-z]*[_]Granulocytes",
            # "CD64CD11b_CD11b+CD16+MatureNeutrophils",
            # b cell
            "CD66[.A-Za-z]*CD14[.A-Za-z]*[_]Livecells", 
            "CD19[.A-Za-z]*CD20[.A-Za-z]*[_]CD19+Bcells", 
            "IgD[.A-Za-z]*IgM[.A-Za-z]*[_]CD19[+]Bcells", 
            "CD10[.A-Za-z]*CD27[.A-Za-z]*[_]CD19+CD20+",
            "CD7[.A-Za-z]*IgD[.A-Za-z]*[_]CD10[-]"
          ), function(x) grepl(x, scat, ignore.case=TRUE)))
          if (fltf) scat <- paste0(scat,"_")
          
          # get 2D csv and clr
          csv_prnt <- fcs_temp@exprs[,stringr::str_split(dim_id,"_")[[1]],drop=FALSE]
          colnames(csv_prnt) <- markers[mrk2]
          nacol <- is.na(colnames(csv_prnt))
          if (any(nacol))
            colnames(csv_prnt)[nacol] <- fcs_temp@parameters@data$desc[
              fcs_temp@parameters@data$name==mrk2[nacol]]
          clr_chld <- sapply(plt_id, function(x) 
            prnt_id %in% which(flowWorkspace::gh_pop_get_indices(gsl[[i]], chlds[x])) )
          if (is.na(dim(clr_chld))) clr_chld <- matrix(clr_chld, ncol=1)
          colnames(clr_chld) <- chlds_short[plt_id]
          
          # if there are too little cells in "other", remove them
          no_cp <- rowSums(clr_chld)
          if (sum(no_cp==0)>0) {
            if (sum(no_cp==0)<50) {
              c_inds <- no_cp!=0
              csv_prnt <- csv_prnt[c_inds,,drop=FALSE]
              clr_chld <- clr_chld[c_inds,,drop=FALSE]
            } else {
              clr_chld <- cbind(clr_chld, no_cp==0)
              colnames(clr_chld)[ncol(clr_chld)] <- "other"
            }
          }
          
          # handle some errors: don't do IGM
          if (any(colnames(clr_chld)%in%"IGM"))
            clr_chld <- clr_chld[,!colnames(clr_chld)%in%"IGM",drop=FALSE]
          
          # handle some errors: convex hulled cell populations overlap i.e. "CD45+CD66-_NonGranulocytes"
          if (sum(rowSums(clr_chld)>1)>0) { #scat=="CD45+CD66-_NonGranulocytes") {
            for (j in which(rowSums(clr_chld)>1)) {
              a <- which(clr_chld[j,])[1]
              clr_chld[j,] <- rep(FALSE, ncol(clr_chld))
              clr_chld[j,a] <- TRUE
            }
          }
          
          if (sum(rowSums(clr_chld)>1)>0) print(paste0(prnt,": ",dim_id))
          # while (sum(no_cp>1)>0) {
          #   a <- clr_chld[no_cp>1,]
          #   a <- a[,colSums(a)>0]
          #   suscol <- colnames(a)[colSums(a[a[,1],])>0] # suspicious cols haha
          #   supersuscol <- suscol[which.max(colSums(clr_chld[,suscol]))]
          #   clr_chld[no_cp>1,supersuscol] <- FALSE
          #   no_cp <- rowSums(clr_chld)
          # }
          
          # plot + save filters
          flowDensity::plotDens(
            fcs_temp, mrk2, 
            main=paste0(scat,"\n",paste0(colnames(clr_chld), collapse=", ")), cex.lab=2, cex.axis=2, cex.main=2)
          filt_ <- list()
          for (gate in gates) {
            lines(gate@boundaries, lwd=2)
            filt_[[gate@filterId]] <- gate@boundaries
          }
          dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
          save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))

          clr_chld <- clr_chld[,colSums(clr_chld)>0,drop=FALSE]
          if (ncol(clr_chld)==1) next
          
          # get thresholds
          if (fltf) {
            gate_ <- c()
            for (mrki in seq_len(length(mrk2))) {
              mrkr <- colnames(csv_prnt)[mrki]
              ranges <- apply(clr_chld, 2, function(x) range(csv_prnt[x,mrkr]))
              ranges <- ranges[,order(ranges[1,])]
              gth <- c()
              for (scol in seq_len(ncol(ranges)-1)) {
                if (any(ranges[1,]>ranges[2,scol])) {
                  gth <- append(gth, min(ranges[1,ranges[1,]>ranges[2,scol]]))
                }
              }
              gth <- unique(gth)
              if (length(gth)>0) {
                gthres[[scat]] <- append(gthres[[scat]], gth[1])
                gate_ <- append(gate_, gth[1])
                names(gthres[[scat]])[length(gthres[[scat]])] <- names(gate_)[length(gate_)] <- mrkr
              }
            }
            dir.create(paste0(gate_,"/",dset,"/",scat), showWarnings=FALSE)
            save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
          }
          
          clr_chld[which(clr_chld)] <- 1
          
          # save
          dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
          write.csv(csv_prnt, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
          dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
          write.csv(clr_chld, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
        }
      }
      graphics.off()
      
      # gthres <- gthres[!names(gthres)%in%"CD16CD56_CD3-"]
      # save(gthres, file=paste0(thres_dir,"/",fids[i],".Rdata"))
      
      rm(fcs)
    })
  }, .parallel=TRUE)
  time_output(start1)
  rm(gsl)
}

time_output(start)