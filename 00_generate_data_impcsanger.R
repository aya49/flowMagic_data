# date created: 2020-12-12
# author: alice yue
# input: IMPC Sanger P2 data set's raw FCS files + their gates
# output: cleaned FCS csv files (cell x marker) + their clr (cell x cell pop) + straight thresholds for some gates (for flowLearn)


## parallelization ####
future::plan(future::multiprocess)
no_cores <- 15#parallel::detectCores() - 5


## directory ####
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src"
setwd(root)


## packages ####
source("helpers.R")
libr(c(
  "furrr", #"rslurm",
  "flowDensity", "flowCore",
  "stringr", 
  "pracma", "quantmod",
  "colorspace"
))


## input ####

# results from gating
results_dir <- "/mnt/FCS_local/IMPC/IMPC-Results/3i/Panel_P2-cell/Results_20170121"

# indices
id_dir <- paste0(results_dir, "/indices")
gt_dir <- paste0(results_dir, "/gates")

# fcs metadata; obtained from 01_preProcessing
store.allFCS_dir <- paste0(results_dir, "/store.allFCS.Rdata") 

# matrix used to transform fcs files; obtained from 02_Clean, or create by applying logicleTransform on globalFrame obtained from 01_preProcessing
lgl_dir <- paste0(results_dir, "/lgl.Rdata")

#indices to delete from each file; obtained from 02_Clean
ind.marg.neg.clean.all_dir <- paste0(results_dir, "/ind.marg.neg.clean.all.Rdata") 

# raw fcs files
fcs_dir <- "/mnt/FCS_local/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_Sanger"


## ouput ####
out_dir <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/data"
x_dir <- paste0(out_dir,"/nD/x/sangerP2"); dir.create(x_dir, recursive=TRUE, showWarnings=FALSE)
y_dir <- paste0(out_dir,"/nD/y/sangerP2"); dir.create(y_dir, recursive=TRUE, showWarnings=FALSE)
plot_dir <- paste0(out_dir,"/scatterplots/sangerP2"); dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
x2_dir <- paste0(out_dir,"/2D/x/sangerP2"); dir.create(x2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr
y2_dir <- paste0(out_dir,"/2D/y/sangerP2"); dir.create(y2_dir, recursive=TRUE, showWarnings=FALSE) #csv, clr
thres_dir <- paste0(out_dir,"/2D/thresholds/sangerP2"); dir.create(thres_dir, recursive=TRUE, showWarnings=FALSE)


## gating names ###

markers <- c(
  "FSCa","FSCh","FSCw",
  "SSCa","SSCh","SSCw",
  "mhcii","CD5","Ly6C",
  "Live","CD23","CD11b",
  "CD161","CD11c","CD21",
  "CD19")
fluor <- c(
  "FSC-A", "FSC-H", "FSC-W",
  "SSC-A", "SSC-H", "SSC-W",
  "FITC-A", "APC-A", "Alexa Fluor 700-A",
  "APC-Cy7-A", "BV421-A", "BV510-A",
  "BV650-A", "BV786-A", "PE-A", 
  "PE-Cy7-A"
)

leaf_cpops <- c(
  "monocyte","granulocyte","eosinophil", 
  "NK immature Ly6c+", "NK mature Ly6c+", "NK immature Ly6c-", "NK mature Ly6c-", 
  "NKT CD11b-Ly6C+", "NKT CD11b+Ly6C+", "NKT CD11b-Ly6C-", "NKT CD11b+Ly6C-", 
  "Tcell",#"Tcell Ly6C+", "Tcell Ly6C-",
  "not Tcell", 
  "cDC CD8+", "cDC CD11b+", 
  "Bcell B1", "Bcell B2 preB", "Bcell B2 MZB", "Bcell B2 folB")

CD161slant <- TRUE # I will use slanted gates here because it is moer natural, not used in original


## load inputs ####
store.allFCS <- get(load(store.allFCS_dir))
ind.marg.neg.clean.all <- get(load(ind.marg.neg.clean.all_dir))
lgl <- get(load(lgl_dir))

error_files <- list.files(paste0(results_dir,"/Figures/gating_flagged"),
                          recursive=TRUE, pattern="png$")
error_bc <- stringr::str_extract(error_files,"L[0-9]+")

# store.allFCS <- store.allFCS[!store.allFCS[,"Barcodes"]%in%error_bc,]

fcs_files <- list.files(fcs_dir, recursive=TRUE, pattern="fcs$", full.names=TRUE)
fcs_files <- fcs_files[grepl("P2",fcs_files)]
fcs_files <- fcs_files[unlist(sapply(store.allFCS[,"Barcodes"], grep, fcs_files))]



## gating ####
saveandrm <- TRUE # set to FALSE to check, TRUE to run through everything
start <- Sys.time()

loop_ind <- loop_ind_f(seq_len(nrow(store.allFCS)), no_cores)
res <- furrr::future_map(loop_ind, function(ii) { res1 <- purrr::map(ii, function(i) {
  # res <- rslurm::slurm_apply(
  #   params=seq_len(nrow(store.allFCS)),
  #   jobname=""
  #   f=function(ii) { purrr::map(ii, function(i) {
  if (store.allFCS[i,"Barcodes"]%in%error_bc) return(NULL)
  
  fid <- paste0(stringr::str_pad(i, 4, pad="0"), "_", 
                store.allFCS[i,"Genotype"], "_",
                store.allFCS[i,"Barcodes"])
  
  cat("\n",fid)
  
  # load gates
  # idi <- get(load(file=paste0(id_dir, "/", fid, ".Rdata")))
  gti <- get(load(file=paste0(gt_dir, "/", fid, ".Rdata")))
  
  # load fcs
  cat("loading... ")
  f <- flowCore::read.FCS(fcs_files[i])
  
  cat("Preprocessing... ")
  if (length(ind.marg.neg.clean.all[[i]]$ind.marg.neg)>0) 
    f@exprs <- f@exprs[-ind.marg.neg.clean.all[[i]]$ind.marg.neg,]
  if (length(ind.marg.neg.clean.all[[i]]$ind.clean)>0) 
    f@exprs <- f@exprs[-ind.marg.neg.clean.all[[i]]$ind.clean,]
  if (det(f@description$SPILL)!=1) 
    f <- flowCore::compensate(f, f@description$SPILL) 
  f <- flowCore::transform(f, lgl)
  f@exprs <- f@exprs[,fluor]
  
  cat("gating... ")
  start1 <- Sys.time() # time to gate one file
  
  # errored <- FALSE
  try({
    gthres <- list()
    
    graphics.off()
    png(file=paste0(plot_dir, "/", fid, ".png"), width=2200, height=1800)
    par(mfrow=c(4,5),mar=(c(5,5,4,2)+0.1))
    
    ## Gating All Events. for singlets #####
    temp <- flowDensity::flowDensity(
      f, channels=c("FSC-A", "SSC-W"), 
      position=c(NA, F), gates=c(NA, gti["singlets.gate"]))
    singlets.flowD.l <- flowDensity::flowDensity(
      temp, channels=c("FSC-A", "SSC-W"), 
      position=c(NA, T), gates=c(NA, gti["singlets.gate.l"]))
    
    suppressWarnings({
      flowDensity::plotDens(
        f, c("FSC-A","SSC-W"), main="All events", 
        xlab="01. FSC-A", ylab="06. SSC-W", 
        xlim=c(0, 275000), ylim=c(0, 275000), devn=FALSE)
      lines(singlets.flowD.l@filter,lwd=1)
    })
    
    rm(temp)
    
    
    ## Gating Singlets. Plotting Live/Dead_SSC-A for live ####
    live.flowD <- flowDensity::flowDensity(
      singlets.flowD.l, channels=c(10,4), 
      position=c(F,NA), gates=c(gti["live.gate"],NA))
    
    suppressWarnings({
      flowDensity::plotDens(
        singlets.flowD.l@flow.frame,  c(10,4), main="Singlets", 
        xlab="10. live", ylab="04. SSC-A", 
        xlim=c(0, 4.5), ylim=c(0, 275000), devn=FALSE); 
      abline(v=gti["live.gate"], lwd=2); 
      lines(live.flowD@filter, lwd=1)
    })
    
    if (saveandrm)
      rm(singlets.flowD.l)
    
    
    ## Gating Live. Plotting FSC-A_SSC-A for Lymphocytes ####
    temp <- flowDensity::flowDensity(
      live.flowD, channels=c(1,4), 
      position=c(T, F), gates=c(gti["fsc.a.gate"], gti["ssc.a.gate"]))
    lymph.flowD <- flowDensity::flowDensity(
      temp, channels=c(1,4), 
      position=c(F,NA), gates=c(gti["fsc.a.gate.high"], NA))
    
    suppressWarnings({
      flowDensity::plotDens(
        live.flowD@flow.frame, c(1,4), main="Live", 
        xlab="01. FSC-A", ylab="04. SSC-A", 
        xlim=c(0, 200000), ylim=c(0, 200000), devn=FALSE); 
      abline(v=gti["fsc.a.gate"], lwd=2); 
      abline(v=gti["fsc.a.gate.high"], lwd=2); 
      abline(h=gti["ssc.a.gate"], lwd=2)
      lines(lymph.flowD@filter, lwd=1)
    })
    
    id_lymph <- lymph.flowD@index
    csv_f <- lymph.flowD@flow.frame@exprs
    colnames(csv_f) <- markers
    csv_lymph <- csv_f[id_lymph,,drop=FALSE]
    
    clr_lymph <- matrix(0, nrow=nrow(csv_lymph), ncol=length(leaf_cpops))
    colnames(clr_lymph) <- leaf_cpops
    
    
    if (saveandrm) 
      rm(live.flowD, temp)
    
    
    ## Gating Lymphocytes. Plotting CD5(8)_CD11b(12) for Not/Granulocytes ####
    # Create CD5 gate to left (2 ovals); CD11b gate
    lymph1.flowD <- flowDensity::flowDensity(
      lymph.flowD, channels=c(8,12), 
      position=c(NA,T), gates=c(NA, gti["CD11b.gate"]))
    
    tempp <- flowDensity::flowDensity(
      lymph1.flowD, channels=c(8,9), 
      position=c(F,F), gates=c(gti["CD5.gate.high"], gti["Ly6C.gate.high"]))    
    temp <- flowDensity::flowDensity(
      lymph1.flowD, channels=c(8,9), position=c(NA,T), 
      gates=c(NA,gti["Ly6C.gate1"]))
    gran.flowD <- flowDensity::flowDensity(
      tempp, channels=c(8,9), 
      position=c(T,T), gates=c(gti["CD5.gate"],gti["Ly6C.gate1"]))
    
    temp1 <- flowDensity::flowDensity(
      lymph.flowD, channels=c(8,12), 
      position=c(NA,F), gates=c(NA, gti["CD11b.gate"]))
    notGran.flowD <- flowDensity::notSubFrame(
      lymph1.flowD@flow.frame, channels=c(8,9), 
      position="logical", gates="missing", gran.flowD@filter)
    notGran.flowD@index <- unique(sort(c(notGran.flowD@index, temp1@index)))
    notGran.flowD@flow.frame@exprs[temp1@index,] <- 
      temp1@flow.frame@exprs[temp1@index,]
    inter_ngran <- intersect(notGran.flowD@index, gran.flowD@index)
    if (length(inter_ngran)>0) {
      notGran.flowD@index <- setdiff(notGran.flowD@index, inter_ngran)
      notGran.flowD@flow.frame@exprs[inter_ngran,] <- NA
    }
    
    suppressWarnings({
      flowDensity::plotDens(
        lymph1.flowD@flow.frame, channels=c(8,9), main="*2D* CD5CD11b: \nLymphocytes-CD11b+ [leaf: granulocyte]", 
        xlab="08. CD5", ylab="09. Ly6C", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["CD5.gate"], lwd=2); 
      abline(h=gti["Ly6C.gate1"], lwd=2); 
      abline(v=gti["CD5.gate.high"], lwd=2); 
      abline(h=gti["Ly6C.gate.high"], lwd=2)
      
      flowDensity::plotDens(
        notGran.flowD@flow.frame, channels=c(8,12), main="not Granulocytes", 
        xlab="08. CD5", ylab="12. CD11b", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["CD5.gate"], lty=2); 
      abline(v=gti["CD5.gate.high"], lwd=2); 
      abline(h=gti["CD11b.gate"], lty=2); 
      abline(h=gti["CD11b.gate.high"], lwd=2)
    })
    
    gthres[["CD5Ly6C_lymphocyte"]] <- gti[c("CD5.gate","Ly6C.gate1")]
    names(gthres[["CD5Ly6C_lymphocyte"]]) <- c("CD5","Ly6C")

    clr_lymph[id_lymph%in%gran.flowD@index,"granulocyte"] <- 1
    
    
    csv_lymph1 <- csv_f[lymph1.flowD@index, c(8,9), drop=FALSE]
    clr_lymph1 <- matrix(0, nrow=nrow(csv_lymph1), ncol=2)
    colnames(clr_lymph1) <- c("granulocyte","not granulocyte")
    clr_lymph1[lymph1.flowD@index%in%gran.flowD@index,1] <- 1
    clr_lymph1[clr_lymph1[,1]==0,2] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD5Ly6C_lymphocyte"), showWarnings=FALSE)
      write.csv(csv_lymph1, file=gzfile(paste0(x2_dir,"/CD5Ly6C_lymphocyte/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD5Ly6C_lymphocyte"), showWarnings=FALSE)
      write.csv(clr_lymph1, file=gzfile(paste0(y2_dir,"/CD5Ly6C_lymphocyte/",fid,".csv.gz")), row.names=FALSE)
      
      rm(lymph1.flowD, temp, tempp, csv_lymph1, clr_lymph1)
    }
    
    
    ## Gating Not Granulocytes. Plotting Ly6C(9)_CD11b(12) for Not/Monocytes ####
    temp <- flowDensity::flowDensity(
      notGran.flowD, channels=c(9,12), 
      position=c(T,T), gates=c(gti["Ly6C.gate"], gti["CD11b.gate"]))
    mon.flowD <- flowDensity::flowDensity(
      temp, channels=c(9,12), 
      position=c(NA, F), gates=c(NA, gti["CD11b.gate.high"]))
    
    notMon.flowD <- flowDensity::notSubFrame(
      notGran.flowD@flow.frame, channels=c(9,12), 
      position= "logical", gates="missing", mon.flowD@filter)
    inter_nmon <- intersect(notMon.flowD@index, mon.flowD@index)
    if (length(inter_nmon)>0) {
      notMon.flowD@index <- setdiff(notMon.flowD@index, inter_nmon)
      notMon.flowD@flow.frame@exprs[inter_nmon,] <- NA
    }
    
    suppressWarnings({
      flowDensity::plotDens(
        notGran.flowD@flow.frame, channels=c(9,12), main="*2D* Ly6cCD11b: \nNot Granulocytes [leaf: monocyte]", 
        xlab="09. Ly6c", ylab="12. CD11b", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["Ly6C.gate"], lwd=2); 
      abline(h=gti["CD11b.gate"], lwd=2); 
      abline(h=gti["CD11b.gate.high"], lwd=2)
      lines(mon.flowD@filter, lwd=1) #; points(notGran.flowD@flow.frame@exprs[notMon.flowD@index,c(9,12)], col=2, cex=0.1)
    })
    
    gthres[["Ly6cCD11b_NotGranulocyte"]] <- gti[c("Ly6C.gate","CD11b.gate")]
    names(gthres[["Ly6cCD11b_NotGranulocyte"]]) <- c("Ly6C","CD11b")
    
    
    clr_lymph[id_lymph%in%mon.flowD@index,"monocyte"] <- 1
    
    
    csv_notgran <- csv_f[notGran.flowD@index, c(9,12), drop=FALSE]
    clr_notgran <- matrix(0, nrow=nrow(csv_notgran), ncol=2)
    colnames(clr_notgran) <- c("monocyte","not monocyte")
    clr_notgran[notGran.flowD@index%in%mon.flowD@index,1] <- 1
    clr_notgran[clr_notgran[,1]==0,2] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/Ly6cCD11b_NotGranulocyte"), showWarnings=FALSE)
      write.csv(csv_notgran, file=gzfile(paste0(x2_dir,"/Ly6cCD11b_NotGranulocyte/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/Ly6cCD11b_NotGranulocyte"), showWarnings=FALSE)
      write.csv(clr_notgran, file=gzfile(paste0(y2_dir,"/Ly6cCD11b_NotGranulocyte/",fid,".csv.gz")), row.names=FALSE)
      
      rm(lymph.flowD, gran.flowD, temp, csv_notgran, clr_notgran) 
    }
    
    
    ## Gating Not Monocytes. Plotting CD11b(12)_SSC-H(5) for Not/Eosinophils ####
    temp <- flowDensity::flowDensity(
      notMon.flowD, channels=c(12,5), 
      position=c(T,T), gates=c(gti["CD11b.gate"], gti["ssc.h.gate"]))
    eos.flowD <- flowDensity::flowDensity(
      temp, channels=c(12,5), 
      position=c(F,F), gates=c(gti["CD11b.gate.high"], gti["ssc.h.gate.high"]))
    
    notEos.flowD <- flowDensity::notSubFrame(
      notMon.flowD@flow.frame, channels=c(12,5), 
      position="logical", gates="missing", eos.flowD@filter)
    inter_neos <- intersect(notEos.flowD@index, eos.flowD@index)
    if (length(inter_neos)>0) {
      notEos.flowD@index <- setdiff(notEos.flowD@index, inter_neos)
      notEos.flowD@flow.frame@exprs[inter_neos,] <- NA
    }
    
    suppressWarnings({
      flowDensity::plotDens(
        notMon.flowD@flow.frame, channels=c(12,5), main="*2D* CD11bSSCh: \nNot Monocytes [leaf: eosinophil]", 
        xlab="12. CD11b", ylab="05. SSC-H", 
        xlim=c(0, 4.5), ylim=c(0, 200000), devn=FALSE); 
      abline(v=gti["CD11b.gate"], lwd=2); 
      abline(v=gti["CD11b.gate.high"], lwd=2); 
      abline(h=gti["ssc.h.gate"], lwd=2); 
      abline(h=gti["ssc.h.gate.high"], lwd=2)
      lines(eos.flowD@filter, lwd=1) #; points(notMon.flowD@flow.frame@exprs[notEos.flowD@index,c(12,5)], col=2, cex=0.1)
    })
    
    gthres[["CD11bSSCh_NotMonocyte"]] <- gti[c("CD11b.gate","ssc.h.gate")]
    names(gthres[["CD11bSSCh_NotMonocyte"]]) <- c("CD11b","SSCh")
    
    clr_lymph[id_lymph%in%eos.flowD@index,"eosinophil"] <- 1
    
    csv_notmon <- csv_f[notMon.flowD@index, c(12,5), drop=FALSE]
    clr_notmon <- matrix(0, nrow=nrow(csv_notmon), ncol=2)
    colnames(clr_notmon) <- c("eosinophil","not eosinophil")
    clr_notmon[notMon.flowD@index%in%eos.flowD@index,1] <- 1
    clr_notmon[clr_notmon[,1]==0,2] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD11bSSCh_NotMonocyte"), showWarnings=FALSE)
      write.csv(csv_notmon, file=gzfile(paste0(x2_dir,"/CD11bSSCh_NotMonocyte/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD11bSSCh_NotMonocyte"), showWarnings=FALSE)
      write.csv(clr_notmon, file=gzfile(paste0(y2_dir,"/CD11bSSCh_NotMonocyte/",fid,".csv.gz")), row.names=FALSE)
      
      rm(notGran.flowD, mon.flowD, temp, csv_notmon, clr_notmon)
    }
    
    
    ## Gating Not Eosinophils. Plotting CD161(13)_CD19(16) for CD161+/- ####
    temp <- flowDensity::flowDensity(
      notEos.flowD, channels=c(13,16), 
      position=c(T, NA), gates=c(gti["CD161.gate"], NA))
    CD161.flowD <- flowDensity::flowDensity(
      temp, channels=c(13,16), 
      position=c(F, F), gates=c(gti["CD161.gate.high"], gti["CD19.gate"]))
    
    notCD161.flowD <- flowDensity::notSubFrame(
      notEos.flowD@flow.frame, channels=c(13,16), 
      position= "logical", gates="missing", CD161.flowD@filter)
    inter_nCD161 <- intersect(notCD161.flowD@index, CD161.flowD@index)
    if (length(inter_nCD161)>0) {
      notCD161.flowD@index <- setdiff(notCD161.flowD@index, inter_nCD161)
      notCD161.flowD@flow.frame@exprs[inter_nCD161,] <- NA
    }
    
    suppressWarnings({
      flowDensity::plotDens(
        notEos.flowD@flow.frame, channels=c(13,16), main="*2D* CD161CD19: Not Eosinophils", 
        xlab="13. CD161", ylab="16. CD19", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["CD161.gate"], lwd=2); 
      abline(v=gti["CD161.gate.high"], lwd=2); 
      abline(h=gti["CD19.gate"], lwd=2)
      lines(CD161.flowD@filter, lwd=1) #; points(notEos.flowD@flow.frame@exprs[CD161.flowD@index,c(13,16)], col=2, cex=0.1)
    })
    
    gthres[["CD161CD19_NotEosinophils"]] <- gti[c("CD161.gate","CD19.gate")]
    names(gthres[["CD161CD19_NotEosinophils"]]) <- c("CD161","CD19")
    
    
    csv_noteos <- csv_f[notEos.flowD@index, c(13,16), drop=FALSE]
    clr_noteos <- matrix(0, nrow=nrow(csv_noteos), ncol=2)
    colnames(clr_noteos) <- c("CD161+","CD161-")
    clr_noteos[notEos.flowD@index%in%CD161.flowD@index,1] <- 1
    clr_noteos[clr_noteos[,1]==0,2] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD161CD19_NotEosinophils"), showWarnings=FALSE)
      write.csv(csv_noteos, file=gzfile(paste0(x2_dir,"/CD161CD19_NotEosinophils/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD161CD19_NotEosinophils"), showWarnings=FALSE)
      write.csv(clr_noteos, file=gzfile(paste0(y2_dir,"/CD161CD19_NotEosinophils/",fid,".csv.gz")), row.names=FALSE)
      
      rm(notMon.flowD, eos.flowD, csv_noteos, clr_noteos)
    }
    
    
    ## Gating CD161+. Plotting CD5(8)_CD11b(13) for NK T-/Cells ####
    temp <- flowDensity::flowDensity(
      CD161.flowD,channels=c(8,13), 
      position=c(F,NA), gates=c(3.3,NA))
    temps <- flowDensity::flowDensity(
      temp,channels=c(8,13), 
      position=c(F,NA), gates=c(2.8,NA))
    temps@flow.frame <- rotate_fcs(temp@flow.frame, c(8,13), theta=-pi/4)$data
    
    NK.flowD <- flowDensity::flowDensity(
      temps,channels=c(8,13), 
      position=c(F,NA), gates=c(gti["CD5.gate.slant2"],NA))
    NK.flowD@filter <- rotate_fcs(NK.flowD@filter,theta=pi/4)$data; 
    NK.flowD@flow.frame <- rotate_fcs(NK.flowD@flow.frame,c(8,13),theta=pi/4)$data
    
    temp <- flowDensity::flowDensity(
      temp, channels=c(8,13), 
      position=c(NA,F), gates=c(NA, gti["CD161.gate.high"]))
    NKT.flowD <- flowDensity::notSubFrame(
      temp@flow.frame, channels=c(8,13), 
      position= "logical", gates="missing", NK.flowD@filter)
    inter_NKT <- intersect(NKT.flowD@index, NK.flowD@index)
    if (length(inter_NKT)>0) {
      NKT.flowD@index <- setdiff(NKT.flowD@index, inter_NKT)
      NKT.flowD@flow.frame@exprs[inter_NKT,] <- NA
    }
    
    suppressWarnings({
      flowDensity::plotDens(
        CD161.flowD@flow.frame, channels=c(8,13), main="*2D* CD5CD11b: CD161+", 
        xlab="08. CD5 #2~", ylab="13. CD161", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["CD5.gate"], lty=2); 
      abline(v=gti["CD5.gate2"], lwd=2); 
      abline(v=gti["CD5.gate.high"], lwd=2); 
      abline(h=gti["CD161.gate"], lwd=2); 
      abline(h=gti["CD161.gate.high"], lwd=2)
      lines(NK.flowD@filter, lwd=1)
      lines(NKT.flowD@filter, lwd=1)
    })
    
    csv_cd161 <- csv_f[CD161.flowD@index, c(8,13), drop=FALSE]
    clr_cd161 <- matrix(0, nrow=nrow(csv_cd161), ncol=2)
    colnames(clr_cd161) <- c("NK","NKT")
    clr_cd161[CD161.flowD@index%in%NK.flowD@index,1] <- 1
    clr_cd161[CD161.flowD@index%in%NKT.flowD@index,2] <- 1
    cd161_rows <- rowSums(clr_cd161)>0
    csv_cd161 <- csv_cd161[cd161_rows,,drop=FALSE]
    clr_cd161 <- clr_cd161[cd161_rows,,drop=FALSE]

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD5CD11b_lymphocyte"), showWarnings=FALSE)
      write.csv(csv_cd161, file=gzfile(paste0(x2_dir,"/CD5CD11b_lymphocyte/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD5CD11b_lymphocyte"), showWarnings=FALSE)
      write.csv(clr_cd161, file=gzfile(paste0(y2_dir,"/CD5CD11b_lymphocyte/",fid,".csv.gz")), row.names=FALSE)
      
      rm(notEos.flowD, temp, temps, csv_cd161, clr_cd161)
    }
    
    
    ## Gating NK Cells. Plotting CD11b(12)_Ly6C(9) ####
    NKimmatLy6C.flowD <- flowDensity::flowDensity(
      NK.flowD, channels=c(9,12), 
      position=c(T,F), gates=c(gti["Ly6C.gate2"],gti["CD11b.gate2"]))
    
    NKmatLy6C.flowD <- flowDensity::flowDensity(
      NK.flowD, channels=c(9,12), 
      position=c(T,T), gates=c(gti["Ly6C.gate2"],gti["CD11b.gate2"]))
    
    NKimmatNotLy6C.flowD <- flowDensity::flowDensity(
      NK.flowD, channels=c(9,12), 
      position=c(F,F), gates=c(gti["Ly6C.gate2"],gti["CD11b.gate2"]))
    
    NKmatNotLy6C.flowD <- flowDensity::flowDensity(
      NK.flowD, channels=c(9,12), 
      position=c(F,T), gates=c(gti["Ly6C.gate2"],gti["CD11b.gate2"]))
    
    suppressWarnings({
      flowDensity::plotDens(
        NK.flowD@flow.frame, channels=c(9,12), 
        main="*2D* CD11bLy6c: NK-cell [LEAVES]", 
        ylab="12. CD11b #2", xlab="09. Ly6c #2", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(h=gti["CD11b.gate"], lty=2); 
      abline(h=gti["CD11b.gate2"], lwd=2); 
      abline(v=gti["Ly6C.gate"], lty=2); 
      abline(v=gti["Ly6C.gate2"], lwd=2)
      lines(NKimmatLy6C.flowD@filter, lwd=1)
      lines(NKmatNotLy6C.flowD@filter, lwd=1)
    })
    
    gthres[["CD11bLy6c_NK"]] <- gti[c("Ly6C.gate2","CD11b.gate2")]
    names(gthres[["CD11bLy6c_NK"]]) <- c("Ly6C","CD11b")
    
    
    clr_lymph[id_lymph%in%NKimmatLy6C.flowD@index,"NK immature Ly6c+"] <- 1
    clr_lymph[id_lymph%in%NKmatLy6C.flowD@index,"NK mature Ly6c+"] <- 1
    clr_lymph[id_lymph%in%NKimmatNotLy6C.flowD@index,"NK immature Ly6c-"] <- 1
    clr_lymph[id_lymph%in%NKmatNotLy6C.flowD@index,"NK mature Ly6c-"] <- 1
    
    
    csv_nk <- csv_f[NK.flowD@index, c(9,12), drop=FALSE]
    clr_nk <- matrix(0, nrow=nrow(csv_nk), ncol=4)
    colnames(clr_nk) <- c("NK immature Ly6c+","NK mature Ly6c+","NK immature Ly6c-","NK mature Ly6c-")
    clr_nk[NK.flowD@index%in%NKimmatLy6C.flowD@index,1] <- 1
    clr_nk[NK.flowD@index%in%NKmatLy6C.flowD@index,2] <- 1
    clr_nk[NK.flowD@index%in%NKimmatNotLy6C.flowD@index,3] <- 1
    clr_nk[NK.flowD@index%in%NKmatNotLy6C.flowD@index,4] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD11bLy6c_NK"), showWarnings=FALSE)
      write.csv(csv_nk, file=gzfile(paste0(x2_dir,"/CD11bLy6c_NK/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD11bLy6c_NK"), showWarnings=FALSE)
      write.csv(clr_nk, file=gzfile(paste0(y2_dir,"/CD11bLy6c_NK/",fid,".csv.gz")), row.names=FALSE)
      
      rm(NKimmatNotLy6C.flowD, NKimmatLy6C.flowD, NKmatNotLy6C.flowD, NKmatLy6C.flowD)
    }
    
    
    ## Gating NK T-cells. Plotting CD11b(12)_Ly6C(9) ####
    NKTnotCD11bLy6C.flowD <- flowDensity::flowDensity(
      NKT.flowD, channels=c(9,12), 
      position=c(T,F), gates=c(gti["Ly6C.gate3"],gti["CD11b.gate2"]))
    
    NKTCD11bLy6C.flowD <- flowDensity::flowDensity(
      NKT.flowD, channels=c(9,12), 
      position=c(T,T), gates=c(gti["Ly6C.gate3"],gti["CD11b.gate2"]))
    
    NKTnotCD11bNotLy6C.flowD <- flowDensity::flowDensity(
      NKT.flowD, channels=c(9,12), 
      position=c(F,F), gates=c(gti["Ly6C.gate3"],gti["CD11b.gate2"]))
    
    NKTCD11bNotLy6C.flowD <- flowDensity::flowDensity(
      NKT.flowD, channels=c(9,12), 
      position=c(F,T), gates=c(gti["Ly6C.gate3"],gti["CD11b.gate2"]))
    
    suppressWarnings({
      flowDensity::plotDens(
        NKT.flowD@flow.frame, channels=c(9,12), 
        main="*2D* CD11bLy6cT: NKT-cell [LEAVES]", 
        ylab="12. CD11b #2", xlab="09. Ly6c #3", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(h=gti["CD11b.gate"], lty=2); 
      abline(h=gti["CD11b.gate2"], lwd=2); 
      abline(v=gti["Ly6C.gate"], lty=2); 
      abline(v=gti["Ly6C.gate2"], lty=2); 
      abline(v=gti["Ly6C.gate3"], lwd=2)
      lines(NKTCD11bLy6C.flowD@filter, lwd=1)
      lines(NKTnotCD11bNotLy6C.flowD@filter, lwd=1)
    })
    
    gthres[["CD11bLy6cT_NKTcell"]] <- gti[c("Ly6C.gate3","CD11b.gate2")]
    names(gthres[["CD11bLy6cT_NKTcell"]]) <- c("Ly6C","CD11b")
    
    
    clr_lymph[id_lymph%in%NKTnotCD11bLy6C.flowD@index,"NKT CD11b-Ly6C+"] <- 1
    clr_lymph[id_lymph%in%NKTCD11bLy6C.flowD@index,"NKT CD11b+Ly6C+"] <- 1
    clr_lymph[id_lymph%in%NKTnotCD11bNotLy6C.flowD@index,"NKT CD11b-Ly6C-"] <- 1
    clr_lymph[id_lymph%in%NKTCD11bNotLy6C.flowD@index,"NKT CD11b+Ly6C-"] <- 1
    
    
    csv_nkt <- csv_f[NKT.flowD@index, c(9,12), drop=FALSE]
    clr_nkt <- matrix(0, nrow=nrow(csv_nkt), ncol=4)
    colnames(clr_nkt) <- c("NKT CD11b-Ly6C+","NKT CD11b+Ly6C+","NKT CD11b-Ly6C-","NKT CD11b+Ly6C-")
    clr_nkt[NKT.flowD@index%in%NKTnotCD11bLy6C.flowD@index,1] <- 1
    clr_nkt[NKT.flowD@index%in%NKTCD11bLy6C.flowD@index,2] <- 1
    clr_nkt[NKT.flowD@index%in%NKTnotCD11bNotLy6C.flowD@index,3] <- 1
    clr_nkt[NKT.flowD@index%in%NKTCD11bNotLy6C.flowD@index,4] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD11bLy6cT_NKTcell"), showWarnings=FALSE)
      write.csv(csv_nkt, file=gzfile(paste0(x2_dir,"/CD11bLy6cT_NKTcell/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD11bLy6cT_NKTcell"), showWarnings=FALSE)
      write.csv(clr_nkt, file=gzfile(paste0(y2_dir,"/CD11bLy6cT_NKTcell/",fid,".csv.gz")), row.names=FALSE)
      
      rm(CD161.flowD, NK.flowD, NKT.flowD, NKTCD11bNotLy6C.flowD, 
         NKTnotCD11bNotLy6C.flowD, NKTCD11bLy6C.flowD, NKTnotCD11bLy6C.flowD, csv_nkt, clr_nkt)
    }
    
    
    ## Gating not CD161+. Plotting mhcii(7)_CD5(8) for Not/T-cells ####
    temps <- notCD161.flowD
    temps@flow.frame <- rotate_fcs(notCD161.flowD@flow.frame,c(7,8), theta=-pi/8)$data
    temp <- flowDensity::flowDensity(
      temps, channels=c(7,8), 
      position=c(F,NA), gates=c(gti["mhcii.gate.slant"], NA))
    temp@filter <- rotate_fcs(temp@filter,c(7,8),theta=pi/8)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(7,8),theta=pi/8)$data
    
    Tcell.flowD <- flowDensity::flowDensity(
      temp, channels=c(7,8), 
      position=c(T,T), gates=c(gti["mhcii.gate.low"], gti["CD5.gate2"]))
    
    notTcell.flowD <- flowDensity::notSubFrame(
      notCD161.flowD@flow.frame, channels=c(7,8), 
      position="logical", gates="missing", Tcell.flowD@filter)
    inter_nT <- intersect(notTcell.flowD@index, Tcell.flowD@index)
    if (length(inter_nT)>0) {
      notTcell.flowD@index <- setdiff(notTcell.flowD@index, inter_nT)
      notTcell.flowD@flow.frame@exprs[inter_nT,] <- NA
    }
    
    suppressWarnings({
      flowDensity::plotDens(
        notCD161.flowD@flow.frame, channels=c(7,8), 
        main="*2D* mhciiCD5: Not CD161+",
        xlab="07. MHCII", ylab="08. CD5 #2", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["mhcii.gate.low"], lwd=2); 
      abline(v=gti["mhcii.gate"], lty=2); 
      abline(h=gti["CD5.gate"], lty=2); 
      abline(h=gti["CD5.gate2"], lwd=2)
      lines(Tcell.flowD@filter, lwd=1)
    })
    
    clr_lymph[id_lymph%in%Tcell.flowD@index,"Tcell"] <- 1
    
    csv_ncd161 <- csv_f[notCD161.flowD@index, c(7,8), drop=FALSE]
    clr_ncd161 <- matrix(0, nrow=nrow(csv_ncd161), ncol=2)
    colnames(clr_ncd161) <- c("Tcell","not Tcell")
    clr_ncd161[notCD161.flowD@index%in%Tcell.flowD@index,1] <- 1
    clr_ncd161[clr_ncd161[,1]==0,2] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/mhciiCD5_NotCD161"), showWarnings=FALSE)
      write.csv(csv_ncd161, file=gzfile(paste0(x2_dir,"/mhciiCD5_NotCD161/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/mhciiCD5_NotCD161"), showWarnings=FALSE)
      write.csv(clr_ncd161, file=gzfile(paste0(y2_dir,"/mhciiCD5_NotCD161/",fid,".csv.gz")), row.names=FALSE)
      
      rm(csv_ncd161, clr_ncd161)
    }
    
    
    ## Gating T-cells. Plotting Ly6C(9) for Ly6C+ T-cells ####
    temp <- flowDensity::flowDensity(
      Tcell.flowD@flow.frame, channels=c(9,8), 
      position=c(T,NA), gates=c(gti["Ly6C.gate4"],NA))
    TcellLy6C.flowD <- flowDensity::flowDensity(
      temp, channels=c(9,8), 
      position=c(F,NA), gates=c(gti["Ly6C.gate.high"],NA))
    
    TcellNOTLy6C.flowD <- flowDensity::flowDensity(
      Tcell.flowD, channels=c(9,8), 
      position=c(F,NA), gates=c(gti["Ly6C.gate4"],NA))
    
    d <- density(Tcell.flowD@flow.frame@exprs[,9], na.rm=T)
    plot(d, main = "T-cell", xlab="09. Ly6C #4", xlim=c(0, 4.5), ylab="Count/Density"); abline(v=gti["Ly6C.gate.high"], lwd=2); abline(v=gti["Ly6C.gate"], lty=2); abline(v=gti["Ly6C.gate4"], lwd=2)
    
    # clr_lymph[id_lymph%in%TcellLy6C.flowD@index,"Tcell Ly6C+"] <- 1
    # clr_lymph[id_lymph%in%TcellNOTLy6C.flowD@index,"Tcell Ly6C+"] <- 1
    
    
    if (saveandrm) 
      rm(TcellLy6C.flowD, temp, Tcell.flowD)
    
    
    ## Gating Not T-cells. Plotting CD19(16)_CD11c(14) for cDC/B-cell ####
    temps <- notTcell.flowD
    temps@flow.frame <- rotate_fcs(notTcell.flowD@flow.frame, c(16,14), theta=pi/4)$data
    temp <- flowDensity::flowDensity(
      temps, channels=c(16,14), 
      position=c(NA,T), gates=c(NA, gti["CD11c.gate.slant.high"]))
    temp@filter <- rotate_fcs(temp@filter,c(16,14),theta=-pi/4)$data
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(16,14),theta=-pi/4)$data
    temp <- flowDensity::flowDensity(
      temp, channels=c(16,14), 
      position=c(NA,T), gates=c(NA, gti["CD11c.gate"]))
    cDC.flowD <- flowDensity::flowDensity(
      temp, channels=c(16,14), 
      position=c(F,F), gates=c(gti["CD19.gate"], gti["CD11c.gate.high"]))
    
    temp <- flowDensity::flowDensity(
      temps, channels=c(16,14), 
      position=c(NA,F), gates=c(NA, gti["CD11c.gate.slant.low"]))
    temp@filter <- rotate_fcs(temp@filter,c(16,14),theta=-pi/4)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(16,14),theta=-pi/4)$data
    Bcell.flowD <- flowDensity::flowDensity(
      temp, channels=c(16,14), 
      position=c(T,F), gates=c(gti["CD19.gate2"], gti["CD11c.gate"]))
    
    suppressWarnings({
      flowDensity::plotDens(
        notTcell.flowD@flow.frame, channels=c(16,14), 
        main="*2D* CD19CD11c: Not T-cell", 
        xlab="16. CD19", ylab="14. CD11c", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["CD19.gate"], lwd=2); 
      abline(v=gti["CD19.gate2"], lwd=2); 
      abline(h=gti["CD11c.gate"] , lwd=2); 
      abline(h=gti["CD11c.gate.high"], lwd=2)
      lines(Bcell.flowD@filter, lwd=1)
    })
    
    clr_lymph[id_lymph%in%notTcell.flowD@index & 
                !id_lymph%in%Bcell.flowD@index & 
                !id_lymph%in%cDC.flowD@index,"not Tcell"] <- 1
    
    
    csv_nTcell <- csv_f[notTcell.flowD@index, c(16,14), drop=FALSE]
    clr_nTcell <- matrix(0, nrow=nrow(csv_nTcell), ncol=2)
    colnames(clr_nTcell) <- c("Bcell","not Bcell")
    clr_nTcell[notTcell.flowD@index%in%Bcell.flowD@index,1] <- 1
    clr_nTcell[clr_nTcell[,1]==0,2] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD19CD11c_NotTcell"), showWarnings=FALSE)
      write.csv(csv_nTcell, file=gzfile(paste0(x2_dir,"/CD19CD11c_NotTcell/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD19CD11c_NotTcell"), showWarnings=FALSE)
      write.csv(clr_nTcell, file=gzfile(paste0(y2_dir,"/CD19CD11c_NotTcell/",fid,".csv.gz")), row.names=FALSE)
      
      rm(notCD161.flowD, notTcell.flowD, csv_nTcell, clr_nTcell)
    }
    
    
    ## Gating cDC. Plotting CD11b(12)_mhcii(7) ####
    temp <- flowDensity::flowDensity(
      cDC.flowD, channels=c(12,7), 
      position=c(T,T), gates=c(gti["CD11b.gate.lowlow"], gti["mhcii.gate"]))
    cDCCD8.flowD <- flowDensity::flowDensity(
      temp, channels=c(12,7), 
      position=c(F,NA), gates=c(gti["CD11b.gate3"], NA))
    
    temp <- flowDensity::flowDensity(
      cDC.flowD, channels=c(12,7), 
      position=c(T,T), gates=c(gti["CD11b.gate3"], gti["mhcii.gate"]))
    cDCCD11b.flowD <- flowDensity::flowDensity(
      temp, channels=c(12,7), 
      position=c(F,NA), gates=c(gti["CD11b.gate.high"], NA))
    
    suppressWarnings({
      flowDensity::plotDens(
        cDC.flowD@flow.frame, channels=c(12,7), 
        main="*2D* CD11bmhcii: cDC", 
        xlab="12. CD11b #3", ylab="07. MHCII", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn = F); 
      abline(v=gti["CD11b.gate.lowlow"], lwd=2); 
      abline(v=gti["CD11b.gate2"], lty=2); 
      abline(v=gti["CD11b.gate3"], lwd=2); 
      abline(v=gti["CD11b.gate"], lty=2); 
      abline(v=gti["CD11b.gate.high"], lwd=2); 
      abline(h=gti["mhcii.gate"], lwd=2)
      lines(cDCCD8.flowD@filter, lwd=1)
      lines(cDCCD11b.flowD@filter, lwd=1)
    })
    
    clr_lymph[id_lymph%in%cDCCD8.flowD@index,"cDC CD8+"] <- 1
    clr_lymph[id_lymph%in%cDCCD11b.flowD@index,"cDC CD11b+"] <- 1
    
    
    csv_cdc <- csv_f[cDC.flowD@index, c(12,7), drop=FALSE]
    clr_cdc <- matrix(0, nrow=nrow(csv_cdc), ncol=3)
    colnames(clr_cdc) <- c("cdC CD8+","cdC CD11b+", "cdC")
    clr_cdc[cDC.flowD@index%in%cDCCD8.flowD@index,1] <- 1
    clr_cdc[cDC.flowD@index%in%cDCCD11b.flowD@index,2] <- 1
    clr_cdc[clr_cdc[,1]==0 & clr_cdc[,2]==0,3] <- 1

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD11bmhcii_cDC"), showWarnings=FALSE)
      write.csv(csv_cdc, file=gzfile(paste0(x2_dir,"/CD11bmhcii_cDC/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD11bmhcii_cDC"), showWarnings=FALSE)
      write.csv(clr_cdc, file=gzfile(paste0(y2_dir,"/CD11bmhcii_cDC/",fid,".csv.gz")), row.names=FALSE)
      
      rm(cDCCD8.flowD, cDCCD11b.flowD, csv_cdc, clr_cdc)
    }
    
    
    ## Gating B-cell. Plotting CD5(8)_CD21(15) for B1/2 B-cells ####
    temps <- Bcell.flowD
    temps@flow.frame <- rotate_fcs(Bcell.flowD@flow.frame, c(8,15), theta=-pi/16)$data
    temp <- flowDensity::flowDensity(
      temps, channels=c(8,15), 
      position=c(T,NA), gates=c(gti["CD5.gate.slant"], NA))
    temp@filter <- rotate_fcs(temp@filter,c(8,15),theta=pi/16)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(8,15),theta=pi/16)$data
    BcellB1.flowD <- flowDensity::flowDensity(
      temp, channels=c(8,15), 
      position=c(F,F), gates=c(gti["CD5.gate.high"], gti["CD21.gate.high"]))
    
    BcellB2.flowD <- flowDensity::notSubFrame(
      Bcell.flowD@flow.frame, channels=c(8,15), 
      position= "logical", gates="missing", BcellB1.flowD@filter)
    inter_B1B2 <- intersect(BcellB2.flowD@index, BcellB1.flowD@index)
    if (length(inter_B1B2)>0) {
      BcellB2.flowD@index <- setdiff(BcellB2.flowD@index, inter_B1B2)
      BcellB2.flowD@flow.frame@exprs[inter_B1B2,] <- NA
    }
    
    suppressWarnings({
      flowDensity::plotDens(
        Bcell.flowD@flow.frame, channels=c(8,15), 
        main="*2D* CD5CD21: B-cell", 
        xlab="08. CD5 #2", ylab="15. CD21", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["CD5.gate2"], lty=2); 
      abline(v=gti["CD5.gate"], lty=2); 
      abline(v=gti["CD5.gate.high"], lwd=2); 
      abline(h=gti["CD21.gate.high"], lwd=2)
      lines(BcellB1.flowD@filter, lwd=1)
    })
    
    clr_lymph[id_lymph%in%BcellB1.flowD@index,"Bcell B1"] <- 1
    
    
    csv_bcell <- csv_f[Bcell.flowD@index, c(8,15), drop=FALSE]
    clr_bcell <- matrix(0, nrow=nrow(csv_bcell), ncol=2)
    colnames(clr_bcell) <- c("Bcell B1","Bcell B2")
    clr_bcell[Bcell.flowD@index%in%BcellB1.flowD@index,1] <- 1
    clr_bcell[Bcell.flowD@index%in%BcellB2.flowD@index,2] <- 1
    bcell_rows <- rowSums(clr_bcell)>0
    csv_bcell <- csv_bcell[bcell_rows,,drop=FALSE]
    clr_bcell <- clr_bcell[bcell_rows,,drop=FALSE]

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD5CD21_Bcell"), showWarnings=FALSE)
      write.csv(csv_bcell, file=gzfile(paste0(x2_dir,"/CD5CD21_Bcell/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD5CD21_Bcell"), showWarnings=FALSE)
      write.csv(clr_bcell, file=gzfile(paste0(y2_dir,"/CD5CD21_Bcell/",fid,".csv.gz")), row.names=FALSE)
      
      rm(cDC.flowD, csv_bcell, clr_bcell)
    }
    
    
    ## Gating B2 B-cells. Plotting CD23(11)_CD21(15) ####
    temp <- flowDensity::flowDensity(
      BcellB2.flowD, channels=c(11,15), 
      position=c(T,T), gates=c(gti["CD23.gate.low"], gti["CD21.gate.low"]))
    frame.flowD <- flowDensity::flowDensity(
      temp, channels=c(11,15), 
      position=c(F,NA), gates=c(gti["CD23.gate.high"], NA))
    
    temps <- temps2 <- BcellB2.flowD
    temps@flow.frame <- rotate_fcs(BcellB2.flowD@flow.frame,c(11,15), theta=-pi/8)$data
    temps2@flow.frame <- rotate_fcs(BcellB2.flowD@flow.frame,c(11,15), theta=pi/8)$data
    temp <- flowDensity::flowDensity(
      temps2, channels=c(11,15), 
      position=c(NA,T), gates=c(NA, gti["CD21.gate.slant.high"]))
    temp@filter <- rotate_fcs(temp@filter,c(11,15),theta=-pi/8)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(11,15),theta=-pi/8)$data
    temp <- flowDensity::flowDensity(
      temp, channels=c(11,15), 
      position=c(T,NA), gates=c(gti["CD23.gate.lowlow"], NA))
    MZB.flowD <- flowDensity::flowDensity(
      temp, channels=c(11,15), 
      position=c(F,NA), gates=c(gti["CD23.gate.high"], NA))
    
    temp <- flowDensity::flowDensity(
      temps, channels=c(11,15), 
      position=c(F,F), gates=c(gti["CD23.gate.slant"], gti["CD21.gate.slant.low"]))
    temp@filter <- rotate_fcs(temp@filter,c(11,15),theta=pi/8)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(11,15),theta=pi/8)$data
    temp@filter <- rotate_fcs(temp@filter,c(11,15),theta=pi/8)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(11,15),theta=pi/8)$data
    temp <- flowDensity::flowDensity(
      temp, channels=c(11,15), 
      position=c(NA,F), gates=c(NA, gti["CD21.gate.slant.high"]))
    temp@filter <- rotate_fcs(temp@filter,c(11,15),theta=-pi/8)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(11,15),theta=-pi/8)$data
    preB.flowD <- flowDensity::flowDensity(
      temp, channels=c(11,15), 
      position=c(T,T), gates=c(gti["CD23.gate.low"], gti["CD21.gate.low"]))
    
    temp <- flowDensity::flowDensity(
      temps, channels=c(11,15), 
      position=c(NA,T), gates=c(NA, gti["CD21.gate.slant.low"]))
    temp@filter <- rotate_fcs(temp@filter,c(11,15),theta=(pi/8))$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(11,15),theta=(pi/8))$data
    temp@filter <- rotate_fcs(temp@filter,c(11,15),theta=(pi/8))$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(11,15),theta=(pi/8))$data
    temp <- flowDensity::flowDensity(
      temp, channels=c(11,15), 
      position=c(NA,F), gates=c(NA, gti["CD21.gate.slant.high"]))
    temp@filter <- rotate_fcs(temp@filter,c(11,15),theta=-pi/8)$data; 
    temp@flow.frame <- rotate_fcs(temp@flow.frame,c(11,15),theta=-pi/8)$data
    folB.flowD <- flowDensity::flowDensity(
      temp, channels=c(11,15), 
      position=c(F,NA), gates=c(gti["CD23.gate.high"], NA))
    
    suppressWarnings({
      flowDensity::plotDens(
        BcellB2.flowD@flow.frame, channels=c(11,15), 
        main="*2D* CD23CD21: B-cell B2",
        xlab="11. CD23", ylab="15. CD21", 
        xlim=c(0, 4.5), ylim=c(0, 4.5), devn=FALSE); 
      abline(v=gti["CD23.gate.lowlow"], lwd=2); 
      abline(v=gti["CD23.gate.low"], lwd=2); 
      abline(v=gti["CD23.gate.high"], lwd=2); 
      abline(v=gti["CD23.gate"], lty=2); 
      abline(h=gti["CD21.gate.low"], lwd=2); 
      abline(h=gti["CD21.gate"], lty=2) 
      lines(preB.flowD@filter, lwd=1)
      lines(folB.flowD@filter, lwd=1)
      lines(MZB.flowD@filter, lwd=1)
    })
    
    clr_lymph[id_lymph%in%preB.flowD@index,"Bcell B2 preB"] <- 1
    clr_lymph[id_lymph%in%MZB.flowD@index,"Bcell B2 MZB"] <- 1
    clr_lymph[id_lymph%in%folB.flowD@index,"Bcell B2 folB"] <- 1
    clymph_rows <- rowSums(clr_lymph)>0
    clr_lymph <- clr_lymph[clymph_rows,]
    csv_lymph <- csv_lymph[clymph_rows,]
    
    csv_bcellb2 <- csv_f[BcellB2.flowD@index, c(11,15), drop=FALSE]
    clr_bcellb2 <- matrix(0, nrow=nrow(csv_bcellb2), ncol=4)
    colnames(clr_bcellb2) <- c("preB","MZB","folB", "other")
    clr_bcellb2[BcellB2.flowD@index%in%preB.flowD@index,1] <- 1
    clr_bcellb2[BcellB2.flowD@index%in%MZB.flowD@index,2] <- 1
    clr_bcellb2[BcellB2.flowD@index%in%folB.flowD@index,3] <- 1
    clr_bcellb2[BcellB2.flowD@index%in%frame.flowD@index & 
                  clr_bcellb2[,1]==0 & clr_bcellb2[,2]==0 & 
                  clr_bcellb2[,3]==0, 4] <- 1
    bcellb2_rows <- rowSums(clr_bcellb2)>0
    csv_bcellb2 <- csv_bcellb2[bcellb2_rows,,drop=FALSE]
    clr_bcellb2 <- clr_bcellb2[bcellb2_rows,,drop=FALSE]

    if (saveandrm) {
      dir.create(paste0(x2_dir,"/CD23CD21_B2Bcell"), showWarnings=FALSE)
      write.csv(csv_bcellb2, file=gzfile(paste0(x2_dir,"/CD23CD21_B2Bcell/",fid,".csv.gz")), row.names=FALSE)
      dir.create(paste0(y2_dir,"/CD23CD21_B2Bcell"), showWarnings=FALSE)
      write.csv(clr_bcellb2, file=gzfile(paste0(y2_dir,"/CD23CD21_B2Bcell/",fid,".csv.gz")), row.names=FALSE)
      
      write.csv(csv_lymph, file=gzfile(paste0(x_dir,"/",fid,".csv.gz")), row.names=FALSE)
      write.csv(clr_lymph, file=gzfile(paste0(y_dir,"/",fid,".csv.gz")), row.names=FALSE)
      
      rm(Bcell.flowD, BcellB1.flowD, BcellB2.flowD, MZB.flowD, preB.flowD, folB.flowD, csv_bcellb2, clr_bcellb2)
    }
    
    graphics.off()
    
    save(gthres, file=paste0(thres_dir,"/",fid,".Rdata"))
    
  })
  time_output(start1)
  
}) })

time_output(start)