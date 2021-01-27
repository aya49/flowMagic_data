# date created: 2020-12-12
# author: alice yue
# input: immune pregnancy full csv files + clrs
# output: 2D csv files + clrs + straight thresholds for some gates (for flowLearn)


## set directory, load packages, set parallel ####
no_cores <- 15#parallel::detectCores() - 5
# root <- "/mnt/FCS_local2/Brinkman group/Alice/flowMagic_data"
# root <- "/home/ayue/projects/flowMagic_data"
root <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data"
source(paste0(root,"/src/RUNME.R"))


## input ####
# cleaned csv + clr files
csv_dir1 <- "/mnt/FCS_local3/backup/FCS data/Immune_Clock_Human_Pregnancy/Results/FR-FCM-ZY3Q/LeukocytesRdata"
clr_dir1 <- gsub("LeukocytesRdata","clr",csv_dir1)
csv_dir2 <- gsub("3Q","3R",csv_dir1)
clr_dir2 <- gsub("3Q","3R",clr_dir1)
gate_file1 <- "/mnt/FCS_local3/backup/Brinkman group/current/Alice/gating_projects/pregnancy/gates_original/FR-FCM-ZY3Q.Rdata"
gate_file2 <- gsub("3Q","3R",gate_file1)


## ouput ####
dset <- "pregnancy"
plyr::l_ply(
  paste0(c(xn_dir, yn_dir, plotn_dir, x2_dir, y2_dir, thres_dir, filt_dir),"/",dset),
  dir.create, recursive=TRUE, showWarnings=FALSE)


## prep inputs ####
csv_files <- append(
  list.files(csv_dir1, pattern="csv", full.names=TRUE),
  list.files(csv_dir2, pattern="csv", full.names=TRUE))
csv_files <- csv_files[grepl("Unstim",csv_files)]

gates_tr <- rbind(get(load(gate_file1)), get(load(gate_file2)))
gates_id <- c(1:11,18,19,25,26,27,29:31,35:37)

png(file=paste0(data_dir, "/pregnancy_gates.png"), width=1000, height=600)
par(mfrow=c(5,5),mar=(c(5,5,4,2)+0.1))
for (gates_i in gates_id) 
  plot(density(gates_tr[,gates_i]), main=gates_i)
graphics.off()


## gating markers ####
gating_channels <- c(
  "CD57", "CD19", "CD4", "CD8", "IgD", "CD11c", "CD16", "CD3", "CD38", "CD27", "CD14", "CXCR5", "CCR7", "CD45RA", "CD20", "CD127", "CD33", "CD28", "CD161", "TCRgd", "CD123", "CD56", "HLADR", "CD25", 
  "CD235ab_CD61", "CD66", "CD45", "Tbet", "CD7", "FoxP3", "CD11b")
leaf_cpops <- c(
  "leukocyte",
  "granulocyte",
  "Bcell",
  "NK",
  "ncMC",
  "intMC",
  "notMC",
  "cMC",
  "Tcell other",
  "Tregs naive Tcell", "Tregs memory Tcell",  
  "CD4+ memory Tcell", "CD4+ naive Tcell other", 
  "gamma-delta Tcell", "CD4-CD8- Tcell other",
  # "CD25+CD8+ memory Tcell", "CD25+CD8+ naive Tcell", "CD8+ Tcell other",
  "CD8+ memory Tcell","CD8+ naive Tcell",
  "CD4+ Th1 naive Tcell", "CD4+ Th1 memory Tcell")

# error_files <- c("Gates_PTLG021_1")


## gating ####
saveandrm <- TRUE # set to FALSE to check, TRUE to run through everything
start <- Sys.time()

# files with NA as gates
# nfile <- c("PTLG003_1","PTLG003_2","PTLG003_BL","PTLG012_2","PTLG012_BL","PTLG028_2","PTLG028_3","PTLG030_2","PTLG030_BL","PTLG031_2")

loop_ind <- loop_ind_f(seq_len(length(csv_files)), no_cores)
res <- plyr::llply(loop_ind, function(ii) { plyr::l_ply(ii, function(i) { try({

  fid <- stringr::str_extract(csv_files[i], "Gates[_A-Za-z0-9]+.fcs")
  cat("\n",i,"/",length(csv_files), fid)
  all.gthres <- gates_tr[fid,]
  cat(fid,"\n")
  
  cat("loading... ")
  
  # load and prep fcs
  f <- flowCore::read.FCS(gsub(".csv|LeukocytesRdata/","",gsub("Results","data",csv_files[i])))

  channels.ind <- sort(Find.markers(f, gating_channels))
  
  if(is.na(channels.ind["Time"])){  
    channels.ind["Time"] <- grep('Time', flowWorkspace::pData(flowCore::parameters(f))$name)
  }
  
  pregating.channels <- sort(c(grep(colnames(f),pattern = "Bead|bead*"), grep(f@parameters$desc, pattern = "DNA*"), grep(colnames(f), pattern = "Dead"), grep(colnames(f), pattern = "Event_length")))
  names(pregating.channels) <- colnames(f)[pregating.channels]
  
  channels.to.transform <- c(pregating.channels[c("Ir191Di", "Ir193Di")], channels.ind)
  channels.to.transform <- setdiff(channels.to.transform, channels.to.transform['Time'])
  
  ## Arcsinhtransformation on CyTOF data
  asinhTrans <- flowCore::arcsinhTransform(transformationId="ln-transformation", a=1, b=1, c=1)
  translist <- flowCore::transformList(colnames(f)[channels.to.transform], asinhTrans) 
  
  f <- flowCore::transform(f, translist)
  
  
  cat("gating... ")
  start1 <- Sys.time() # time to gate one file
  
  # errored <- FALSE
  
  png(file=paste0(plotn_dir, "/", dset,"/",fid, ".png"), width=2200, height=1800)
  par(mfrow=c(4,5),mar=(c(5,5,4,2)+0.1))
  
  gthres <- list()
  
  singlets.flowD.temp <- flowDensity::flowDensity(
    f, channels=c('Event_length', 'Ir191Di'), 
    position=c(NA,T),  gates=c(NA, all.gthres[1]))
  singlets.flowD <- flowDensity::flowDensity(
    singlets.flowD.temp, channels=c('Event_length', 'Ir191Di'), 
    position=c(NA,F),  gates=c(NA, all.gthres[2]))
  
  flowDensity::plotDens(
    f, c(pregating.channels["Event_length"], pregating.channels["Ir191Di"]),
    main="All cells", cex.lab=2, cex.axis=2, cex.main=2)
  lines(singlets.flowD@filter, lwd=2)
  
  leukocytes.flowD.temp <- flowDensity::flowDensity(
    singlets.flowD, channels=c('In113Di', 'Ir191Di'), 
    position=c(NA, T), gates=c(NA, all.gthres[3]))
  leukocytes.flowD <- flowDensity::flowDensity(
    leukocytes.flowD.temp, channels=c('In113Di', 'Ir191Di'), 
    position=c(F, F), gates=c(all.gthres[4], all.gthres[2]))
  
  flowDensity::plotDens(
    singlets.flowD, c(channels.ind["CD235ab_CD61"],pregating.channels["Ir191Di"]),
    main="Singlets", cex.lab=2, cex.axis=2, cex.main=2, 
    ylim=c(0, 10), xlim=c(0, 10))
  lines(leukocytes.flowD@filter, lwd=2)
  
  id_leuk <- leukocytes.flowD@index
  csv_f <- leukocytes.flowD@flow.frame@exprs[,channels.ind]
  colnames(csv_f) <- names(channels.ind)
  csv_leuk <- csv_f[id_leuk,,drop=FALSE]
  
  clr_leuk <- matrix(0, nrow=nrow(csv_leuk), ncol=length(leaf_cpops))
  colnames(clr_leuk) <- leaf_cpops
  
  if (saveandrm) 
    rm(leukocytes.flowD.temp, singlets.flowD, singlets.flowD.temp)
  
  
  
  ## gating leukocytes > mononuclear, granulocyte, other ####
  theta=atan(atan(tan(pi/4)))
  R <- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)),2 ,2)
  
  leukocytes.temp <- rotate_fcs(flowDensity::getflowFrame(leukocytes.flowD),c(channels.ind["CD66"],channels.ind["CD45"]),theta=pi/4)$data
  
  mononuclear.flowD.temp <- flowDensity::flowDensity(
    leukocytes.temp, channels=c('La139Di','In115Di'), 
    position=c(F,T), gates=c(all.gthres[6]-0.75, all.gthres[7]+0.5), ellip.gate=T)
  mononuclear.temp <- flowDensity::getflowFrame(mononuclear.flowD.temp)
  mononuclear.temp@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])] <- t(t(R) %*% t(mononuclear.temp@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])]))
  mononuclear.flowD.temp@filter <- rotate_fcs(mononuclear.flowD.temp@filter,c(channels.ind["CD66"],channels.ind["CD45"]),theta=-pi/4)$data
  mononuclear.flowD.temp@flow.frame <- rotate_fcs(flowDensity::getflowFrame(mononuclear.flowD.temp),c(channels.ind["CD66"],channels.ind["CD45"]),theta=-pi/4)$data
  mononuclear.flowD <- flowDensity::flowDensity(
    mononuclear.flowD.temp, channels=c('La139Di','In115Di'), 
    position=c(T,NA), gates=c(all.gthres[5], NA))
  mononuclear.flowD <- flowDensity::flowDensity(
    leukocytes.flowD, channels=c('La139Di','In115Di'), 
    position=TRUE, filter=mononuclear.flowD@filter)
  
  granulocytes.flowD <- flowDensity::flowDensity(
    leukocytes.temp, channels=c('La139Di','In115Di'), 
    position=c(F,F), gates=c(all.gthres[6], all.gthres[7]-0.5), ellip.gate=T)
  granulocytes <- flowDensity::getflowFrame(granulocytes.flowD)
  granulocytes@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])] <- t(t(R) %*% t(granulocytes@exprs[,c(channels.ind["CD66"],channels.ind["CD45"])]))
  granulocytes.flowD@filter <- rotate_fcs(granulocytes.flowD@filter,c(channels.ind["CD66"],channels.ind["CD45"]),theta=-pi/4)$data
  granulocytes.flowD@flow.frame <- rotate_fcs(flowDensity::getflowFrame(granulocytes.flowD),c(channels.ind["CD66"],channels.ind["CD45"]),theta=-pi/4)$data
  granulocytes.flowD <- flowDensity::flowDensity(
    leukocytes.flowD, channels=c('La139Di','In115Di'), 
    position=TRUE, filter=granulocytes.flowD@filter)
  
  scat <- "01_CD66CD45_leukocyte"
  flowDensity::plotDens(
    leukocytes.flowD, c(channels.ind["CD66"],channels.ind["CD45"]), 
    main=paste0(scat, "\n[leaf: granulocyte], mononuclear, other"), 
    cex.lab=2, cex.axis=2, cex.main=2)
  lines(mononuclear.flowD@filter, lwd=2)
  lines(granulocytes.flowD@filter, lwd=2)
  
  clr_leuk[id_leuk%in%granulocytes.flowD@index,"granulocyte"] <- 1
  clr_leuk[!id_leuk%in%granulocytes.flowD@index &
             !id_leuk%in%mononuclear.flowD@index,"leukocyte"] <- 1
  
  csv_ <- csv_f[leukocytes.flowD@index, c("CD66","CD45"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=3)
  colnames(clr_) <- c("granulocyte","mononuclear","other")
  clr_[id_leuk%in%granulocytes.flowD@index,1] <- 1
  clr_[id_leuk%in%mononuclear.flowD@index,2] <- 1
  clr_[rowSums(clr_)==0,3] <- 1
  
  filt_ <- list(mononuclear.flowD@filter, granulocytes.flowD@filter)
  names(filt_) <- colnames(clr_)[1:2]
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(leukocytes.flowD, leukocytes.temp, mononuclear.flowD.temp, mononuclear.temp, granulocytes, granulocytes.flowD) 
  }
  
  
  ## gating mononuclear > Bcell, NKLinNeg, Tcell, other ####
  Bcells.flowD <- flowDensity::flowDensity(
    mononuclear.flowD, channels=c('Er170Di', 'Nd142Di'), 
    position=c(F, T), gates=c(all.gthres[8], all.gthres[9]))
  Tcells.flowD <- flowDensity::flowDensity(
    mononuclear.flowD, channels=c('Er170Di', 'Nd142Di'), 
    position=c(T, F), gates=c(all.gthres[8], all.gthres[9]))
  NK.LinNeg.flowD <- flowDensity::flowDensity(
    mononuclear.flowD, channels=c('Er170Di', 'Nd142Di'), 
    position=c(F,F), gates=c(all.gthres[8], all.gthres[9]))
  
  scat <- "02_CD3CD19_mononuclear_"
  flowDensity::plotDens(
    mononuclear.flowD, c(channels.ind["CD3"],channels.ind["CD19"]), 
    main=paste0(scat,"\n[leaf:Bcell], Tcell, NK lin-, other"), 
    cex.lab=2, cex.axis=2, cex.main=2); 
  abline(h=all.gthres[9], lwd=2); 
  abline(v=all.gthres[8], lwd=2)
  lines(Bcells.flowD@filter, lwd=2)
  lines(Tcells.flowD@filter, lwd=2)
  
  clr_leuk[id_leuk%in%Bcells.flowD@index,"Bcell"] <- 1
  
  csv_ <- csv_f[mononuclear.flowD@index, c("CD3","CD19"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=4)
  colnames(clr_) <- c("Bcell","Tcell","NK lin-","other")
  clr_[mononuclear.flowD@index%in%Bcells.flowD@index,1] <- 1
  clr_[mononuclear.flowD@index%in%Tcells.flowD@index,2] <- 1
  clr_[mononuclear.flowD@index%in%NK.LinNeg.flowD@index,3] <- 1
  clr_[rowSums(clr_)==0,4] <- 1
  
  gthres[[scat]] <- gate_<- all.gthres[8:9]
  names(gthres[[scat]]) <- names(gate_) <- c("CD3","CD19")
  
  filt_ <- list(Bcells.flowD@filter, Tcells.flowD@filter, NK.LinNeg.flowD@filter)
  names(filt_) <- colnames(clr_)[1:3]
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    dir.create(paste0(thres_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(mononuclear.flowD) 
  }
  
  
  ## gating NKLinNeg > NK, lin-, other ####
  NKcells.flowD <- flowDensity::flowDensity(
    NK.LinNeg.flowD, channels=c('Lu175Di', 'Pr141Di'), 
    position=c(T, T), gates=c(all.gthres[10], all.gthres[11]))
  LinNeg.flowD <- flowDensity::flowDensity(
    NK.LinNeg.flowD, channels=c('Lu175Di', 'Pr141Di'), 
    position=c(T, F), gates=c(all.gthres[10], all.gthres[11]))
  
  scat <- "03_CD14CD7_NKLinNeg_"
  flowDensity::plotDens(
    NK.LinNeg.flowD, c(channels.ind["CD14"],channels.ind["CD7"]), 
    main=paste0(scat,"\n[leaf: NK], lin-"), 
    cex.lab=2, cex.axis=2, cex.main=2)
  lines(NKcells.flowD@filter, lwd=2)
  lines(LinNeg.flowD@filter, lwd=2)
  
  clr_leuk[id_leuk%in%NKcells.flowD@index,"NK"] <- 1
  
  csv_ <- csv_f[NK.LinNeg.flowD@index, c("CD14","CD7"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=2)
  colnames(clr_) <- c("lin-","NK")
  clr_[NK.LinNeg.flowD@index%in%LinNeg.flowD@index,1] <- 1
  clr_[NK.LinNeg.flowD@index%in%NKcells.flowD@index,2] <- 1
  nklin_rows <- rowSums(clr_)==1
  csv_ <- csv_[nklin_rows,]
  clr_ <- clr_[nklin_rows,]
  
  filt_ <- list(LinNeg.flowD@filter, NKcells.flowD@filter)
  names(filt_) <- colnames(clr_)[1:2]
  
  gthres[[scat]] <- gate_ <- all.gthres[10]
  names(gthres[[scat]]) <- names(gate_) <- "CD7"
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    dir.create(paste0(thres_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(NK.LinNeg.flowD) 
  }
  
  
  # ## gating NK > CD56-CD16+, CD56+CD16-, other ####
  # cd56low.cd16pos.temp <- flowDensity::flowDensity(
  #   NKcells.flowD, channels=c('Yb176Di', 'Ho165Di'), 
  #   position=c(F,T), gates=c(all.gthres[13], all.gthres[16]))
  # cd56low.cd16pos.flowD <- flowDensity::flowDensity(
  #   cd56low.cd16pos.temp, channels=c('Yb176Di', 'Ho165Di'), 
  #   position=c(T,F), gates=c(all.gthres[12], all.gthres[17]))
  # 
  # cd56pos.cd16neg.flowD.temp <- flowDensity::flowDensity(
  #   NKcells.flowD, channels=c('Yb176Di', 'Ho165Di'), 
  #   position=c(F,F), gates=c(all.gthres[14], all.gthres[16]))
  # cd56pos.cd16neg.flowD <- flowDensity::flowDensity(
  #   cd56pos.cd16neg.flowD.temp, channels=c('Yb176Di', 'Ho165Di'), 
  #   position=c(T,T), gates=c(all.gthres[12], all.gthres[15]))
  # 
  # flowDensity::plotDens(
  #   NKcells.flowD, c(channels.ind["CD56"],channels.ind["CD16"]), 
  #   main="*2D* CD56CD16: NK", cex.lab=2, cex.axis=2, cex.main=2)
  # lines(cd56low.cd16pos.flowD@filter)
  # lines(cd56pos.cd16neg.flowD@filter)
  # 
  # csv_nk <- csv_f[NKcells.flowD@index, c("CD56","CD16"), drop=FALSE]
  # clr_nk <- matrix(0, nrow=nrow(csv_nk), ncol=3)
  # colnames(clr_nk) <- c("NK CD56-CD16+","NK CD56+CD16-","other")
  # clr_nk[NKcells.flowD@index%in%cd56low.cd16pos.flowD@index,1] <- 1
  # clr_nk[NKcells.flowD@index%in%cd56pos.cd16neg.flowD@index,2] <- 1
  # clr_nk[rowSums(clr_nk)==0,3] <- 1
  
  if (saveandrm) {
    # dir.create(paste0(x2_dir,"/CD56CD16"), showWarnings=FALSE)
    # write.csv(csv_nk, file=gzfile(paste0(x2_dir,"/CD56CD16/",fid,".csv.gz")), row.names=FALSE)
    # dir.create(paste0(y2_dir,"/CD56CD16"), showWarnings=FALSE)
    # write.csv(clr_nk, file=gzfile(paste0(y2_dir,"/CD56CD16/",fid,".csv.gz")), row.names=FALSE)
    
    rm(NKcells.flowD, cd56low.cd16pos.temp, cd56pos.cd16neg.flowD.temp, cd56low.cd16pos.flowD, cd56pos.cd16neg.flowD)
  }
  
  
  ## gating lin- > cMC, ncMC, intMC, notMC ####
  cMC.flowD <- flowDensity::flowDensity(
    LinNeg.flowD, channels=c('Lu175Di', 'Ho165Di'), 
    position=c(T, F), gates=c(all.gthres[18], all.gthres[19]))
  ncMC.flowD <- flowDensity::flowDensity(
    LinNeg.flowD, channels=c('Lu175Di', 'Ho165Di'), 
    position=c(F, T), gates=c(all.gthres[18], all.gthres[19]))
  intMC.flowD <- flowDensity::flowDensity(
    LinNeg.flowD, channels=c('Lu175Di', 'Ho165Di'), 
    position=c(T, T), gates=c(all.gthres[18], all.gthres[19]))
  
  NOT.MC.flowD <- flowDensity::flowDensity(
    LinNeg.flowD, channels=c('Lu175Di', 'Ho165Di'), 
    position=c(F, F), gates=c(all.gthres[18], all.gthres[19]))
  
  scat <- "04_CD14CD16_lin-_"
  flowDensity::plotDens(
    LinNeg.flowD, c(channels.ind["CD14"],channels.ind["CD16"]), 
    main=paste0(scat,"\n[leaf: cMC, ncMC, intMC, notMC]"), 
    cex.lab=2, cex.axis=2, cex.main=2)
  lines(cMC.flowD@filter, lwd=2)
  lines(ncMC.flowD@filter, lwd=2)
  lines(intMC.flowD@filter, lwd=2)
  lines(NOT.MC.flowD@filter, lwd=2)
  
  clr_leuk[id_leuk%in%ncMC.flowD@index,"ncMC"] <- 1
  clr_leuk[id_leuk%in%intMC.flowD@index,"intMC"] <- 1
  clr_leuk[id_leuk%in%NOT.MC.flowD@index,"notMC"] <- 1
  clr_leuk[id_leuk%in%cMC.flowD@index,"cMC"] <- 1
  
  csv_ <- csv_f[LinNeg.flowD@index, c("CD14","CD16"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=4)
  colnames(clr_) <- c("cMC","ncMC","intMC","notMC")
  clr_[LinNeg.flowD@index%in%cMC.flowD@index,1] <- 1
  clr_[LinNeg.flowD@index%in%ncMC.flowD@index,2] <- 1
  clr_[LinNeg.flowD@index%in%intMC.flowD@index,3] <- 1
  clr_[LinNeg.flowD@index%in%NOT.MC.flowD@index,4] <- 1
  
  filt_ <- list(cMC.flowD@filter, ncMC.flowD@filter, intMC.flowD@filter, NOT.MC.flowD@filter)
  names(filt_) <- colnames(clr_)[1:4]
  
  gthres[[scat]] <- gate_ <- all.gthres[18:19]
  names(gthres[[scat]]) <- names(gate_) <- c("CD14","CD16")
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    dir.create(paste0(thres_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(cMC.flowD, ncMC.flowD, intMC.flowD) 
  }
  
  
  # ## gating notMC > pDC, something, other ####
  # pDC.flowD <- flowDensity::flowDensity(
  #   ncMC.flowD, channels=c('Nd148Di', 'Yb174Di'), 
  #   position=c(T, T), gates=c(all.gthres[20], all.gthres[21]))
  # 
  # flowDensity::plotDens(
  #   ncMC.flowD, c(channels.ind["CD123"],channels.ind["HLADR"]), 
  #   main="pDCs", cex.lab=2, cex.axis=2, cex.main=2)
  # lines(pDC.flowD@filter, lwd=2)
  # 
  # clr_leuk[id_leuk%in%pDC.flowD@index,"pDC"] <- 1
  # 
  # csv_ncMC <- csv_f[ncMC.flowD@index, c("CD123","HLADR"), drop=FALSE]
  # clr_ncMC <- matrix(0, nrow=nrow(csv_ncMC), ncol=2)
  # colnames(clr_ncMC) <- c("pDC","other")
  # clr_ncMC[ncMC.flowD@index%in%pDC.flowD@index,1] <- 1
  # clr_ncMC[clr_ncMC[,1]==0,2] <- 1
  
  if (saveandrm) {
    #   dir.create(paste0(x2_dir,"/CD123HLADR"), showWarnings=FALSE)
    #   write.csv(csv_ncMC, file=gzfile(paste0(x2_dir,"/CD123HLADR/",fid,".csv.gz")), row.names=FALSE)
    #   dir.create(paste0(y2_dir,"/CD123HLADR"), showWarnings=FALSE)
    #   write.csv(clr_ncMC, file=gzfile(paste0(y2_dir,"/CD123HLADR/",fid,".csv.gz")), row.names=FALSE)
    
    rm(ncMC.flowD, pDC.flowD)
  }
  
  
  # ## gating notMC > mDC, other ####
  # mDCs.flowD <- flowDensity::flowDensity(
  #   NOT.MC.flowD, channels=c('Sm147Di', 'Yb174Di'), 
  #   position=c(T,T), gates=c(all.gthres[22], all.gthres[21]))
  # 
  # flowDensity::plotDens(
  #   NOT.MC.flowD, c(channels.ind["CD11c"],channels.ind["HLADR"]), 
  #   main="mDCs", cex.lab=2, cex.axis=2, cex.main=2)
  # lines(mDCs.flowD@filter, lwd=2)
  # 
  # csv_nMC <- csv_f[NOT.MC.flowD@index, c("CD11c","HLADR"), drop=FALSE]
  # clr_nMC <- matrix(0, nrow=nrow(csv_nMC), ncol=2)
  # colnames(clr_ncMC) <- c("pDC","other")
  # clr_nMC[NOT.MC.flowD@index%in%mDCs.flowD@index,1] <- 1
  # clr_nMC[clr_nMC[,1]==0,2] <- 1
  
  if (saveandrm) {
    # dir.create(paste0(x2_dir,"/CD123HLADR"), showWarnings=FALSE)
    # write.csv(csv_nMC, file=gzfile(paste0(x2_dir,"/CD123HLADR/",fid,".csv.gz")), row.names=FALSE)
    # dir.create(paste0(y2_dir,"/CD123HLADR"), showWarnings=FALSE)
    # write.csv(clr_nMC, file=gzfile(paste0(y2_dir,"/CD123HLADR/",fid,".csv.gz")), row.names=FALSE)
    
    rm(ncMC.flowD, pDC.flowD) 
  }
  
  
  # ## gating cMC > M-MDSC, other, densecluster ####
  # MMDSC.flowD <- flowDensity::flowDensity(
  #   cMC.flowD, channels=c('Nd144Di', 'Yb174Di'), 
  #   position=c(T, F), gates=c(all.gthres[23], all.gthres[24]))
  # 
  # flowDensity::plotDens(
  #   cMC.flowD, c(channels.ind["CD11b"],channels.ind["HLADR"]), 
  #   main="M-MDSCs", cex.lab=2, cex.axis=2, cex.main=2)
  # lines(MMDSC.flowD@filter, lwd=2)
  # 
  # csv_cMC <- csv_f[cMC.flowD@index, c("CD11c","HLADR"), drop=FALSE]
  # clr_cMC <- matrix(0, nrow=nrow(csv_cMC), ncol=2)
  # colnames(clr_ncMC) <- c("MMDSC","other")
  # clr_cMC[cMC.flowD@index%in%MMDSC.flowD@index,1] <- 1
  # clr_cMC[clr_cMC[,1]==0,2] <- 1
  # 
  if (saveandrm) {
    #   dir.create(paste0(x2_dir,"/CD123HLADR2"), showWarnings=FALSE)
    #   write.csv(csv_cMC, file=gzfile(paste0(x2_dir,"/CD123HLADR2/",fid,".csv.gz")), row.names=FALSE)
    #   dir.create(paste0(y2_dir,"/CD123HLADR2"), showWarnings=FALSE)
    #   write.csv(clr_cMC, file=gzfile(paste0(y2_dir,"/CD123HLADR2/",fid,".csv.gz")), row.names=FALSE)
    
    rm(MMDSC.flowD, cMC.flowD) 
  }
  
  
  
  ## gating Tcell > CD8+, CD4+, notCD4CD8Tcells ####
  cd4.Tcells.flowD <- flowDensity::flowDensity(
    Tcells.flowD, channels=c('Nd145Di','Nd146Di'), 
    position=c(T,F), gates=c(all.gthres[25], all.gthres[26]))
  cd8.Tcells.flowD <- flowDensity::flowDensity(
    Tcells.flowD, channels=c('Nd145Di','Nd146Di'), 
    position=c(F,T), gates=c(all.gthres[25], all.gthres[26]))
  NOT.cd4.cd8.Tcells.flowD <- flowDensity::flowDensity(
    Tcells.flowD, channels=c('Nd145Di','Nd146Di'), 
    position=c(F,F), gates=c(all.gthres[25], all.gthres[26]))
  cd4.cd8.Tcells.flowD <- flowDensity::flowDensity(
    Tcells.flowD, channels=c('Nd145Di','Nd146Di'), 
    position=c(T,T), gates=c(all.gthres[25], all.gthres[26]))
  
  scat <- "05_CD4CD8_Tcell_"
  flowDensity::plotDens(
    Tcells.flowD, c(channels.ind["CD4"],channels.ind["CD8"]), 
    main=paste0(scat,"\n[leaf: Tcell other], CD4+, CD8+, CD4-CD8-, CD4+CD8+ Tcell"), 
    cex.lab=2, cex.axis=2, cex.main=2,
    xlim=c(0,10), ylim=c(0,10))
  lines(cd4.Tcells.flowD@filter, lwd=2)
  lines(cd8.Tcells.flowD@filter, lwd=2)
  lines(NOT.cd4.cd8.Tcells.flowD@filter, lwd=2)

  # clr_leuk[id_leuk%in%cd4.Tcells.flowD@index,"CD4+ Tcell"] <- 1
  # clr_leuk[id_leuk%in%cd8.Tcells.flowD@index,"CD8+ Tcell"] <- 1
  # clr_leuk[id_leuk%in%NOT.cd4.cd8.Tcells.flowD@index,"CD4-CD8- Tcell"] <- 1
  clr_leuk[id_leuk%in%cd4.cd8.Tcells.flowD@index,"Tcell other"] <- 1
  
  csv_ <- csv_f[Tcells.flowD@index, c("CD4","CD8"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=4)
  colnames(clr_) <- c("CD4+ Tcell", "CD8+ Tcell", "CD4-CD8- Tcell", "CD4+CD8+ Tcell")
  clr_[Tcells.flowD@index%in%cd4.Tcells.flowD@index,1] <- 1
  clr_[Tcells.flowD@index%in%cd8.Tcells.flowD@index,2] <- 1
  clr_[Tcells.flowD@index%in%NOT.cd4.cd8.Tcells.flowD@index,3] <- 1
  clr_[Tcells.flowD@index%in%cd4.cd8.Tcells.flowD@index,4] <- 1
  
  filt_ <- list(cd4.Tcells.flowD@filter, cd8.Tcells.flowD@filter, NOT.cd4.cd8.Tcells.flowD@filter, cd4.cd8.Tcells.flowD@filter)
  names(filt_) <- colnames(clr_)[1:4]
  
  gthres[[scat]] <- gate_ <- all.gthres[25:26]
  names(gthres[[scat]]) <- names(gate_) <- c("CD4","CD8")
  if (is.na(all.gthres[25])) 
    gthres[["CD4CD8_Tcell"]][1] <- gate_[1] <- 
    min(cd4.Tcells.flowD@flow.frame@exprs[!is.na(cd4.Tcells.flowD@flow.frame@exprs[,"Nd145Di"]),"Nd145Di"])
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    dir.create(paste0(thres_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(Tcells.flowD) 
  }
  
  
  ## gating Tcell > CD8+T, CD4+T, notCD4CD8Tcells ####
  cd4.T.Naive.flowD <- flowDensity::flowDensity(
    cd4.Tcells.flowD, channels=c('Nd145Di', 'Nd143Di'), 
    position=c(NA, T), gates=c(NA, all.gthres[27]))
  cd4.T.Memory.flowD <- flowDensity::flowDensity(
    cd4.Tcells.flowD, channels=c('Nd145Di', 'Nd143Di'), 
    position=c(NA, F), gates=c(NA, all.gthres[27]))
  
  scat <- "06_CD4CD45RA_Tcell_"
  flowDensity::plotDens(
    cd4.Tcells.flowD, c(channels.ind["CD4"],channels.ind["CD45RA"]), 
    main=paste0(scat,"\n CD4+ naive, CD4+ memory Tcell"), 
    cex.lab=2, cex.axis=2, cex.main=2,
    xlim=c(1,10), ylim=c(1,10))
  lines(cd4.T.Naive.flowD@filter, lwd=2)
  lines(cd4.T.Memory.flowD@filter, lwd=2)
  
  clr_leuk[id_leuk%in%cd4.T.Memory.flowD@index,"CD4+ memory Tcell"] <- 1
  
  csv_ <- csv_f[cd4.Tcells.flowD@index, c("CD4","CD45RA"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=2)
  colnames(clr_) <- c("CD4+ naive Tcell", "CD4+ memory Tcell")
  clr_[cd4.Tcells.flowD@index%in%cd4.T.Naive.flowD@index,1] <- 1
  clr_[cd4.Tcells.flowD@index%in%cd4.T.Memory.flowD@index,2] <- 1
  
  filt_ <- list(cd4.T.Naive.flowD@filter, cd4.T.Memory.flowD@filter)
  names(filt_) <- colnames(clr_)[1:2]
  
  gthres[[scat]] <- gate_ <- all.gthres[27]
  names(gthres[[scat]]) <- names(gate_) <- "CD45RA"
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    dir.create(paste0(thres_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(cd4.T.Memory.flowD) 
  }
  
  
  # ## gating CD4+Tcell > CD4+Tcell Th1 naive, CD4+Tcell Th1 memory ####
  # cd4.Th1.Naive.flowD <- flowDensity::flowDensity(
  #   cd4.Tcells.flowD, channels=c('Gd160Di', 'Nd143Di'), 
  #   position=c(T, T), gates=c(all.gthres[28], all.gthres[27]))
  # cd4.Th1.Memory.flowD <- flowDensity::flowDensity(
  #   cd4.Tcells.flowD, channels=c('Gd160Di', 'Nd143Di'), 
  #   position=c(T, F), gates=c(all.gthres[28], all.gthres[27]))
  # 
  # flowDensity::plotDens(
  #   cd4.Tcells.flowD, c(channels.ind["Tbet"],channels.ind["CD45RA"]), 
  #   main="*2D* TbetCD45RA: CD4+ Th1 Tcells", cex.lab=2, cex.axis=2, cex.main=2)
  # lines(cd4.Th1.Naive.flowD@filter, lwd=2)
  # lines(cd4.Th1.Memory.flowD@filter, lwd=2)
  # 
  # # clr_leuk[id_leuk%in%cd4.Th1.Naive.flowD@index,"CD4+ Th1 naive Tcell"] <- 1
  # # clr_leuk[id_leuk%in%cd4.Th1.Memory.flowD@index,"CD4+ Th1 memory Tcell"] <- 1
  # 
  # csv_t4Th <- csv_f[cd4.Tcells.flowD@index, c("Tbet","CD45RA"), drop=FALSE]
  # clr_t4Th <- matrix(0, nrow=nrow(csv_t4Th), ncol=3)
  # colnames(clr_t4Th) <- c("CD4+ naive Tcell", "CD4+ memory Tcell", "other")
  # clr_t4Th[cd4.Tcells.flowD@index%in%cd4.Th1.Naive.flowD@index,1] <- 1
  # clr_t4Th[cd4.Tcells.flowD@index%in%cd4.Th1.Memory.flowD@index,2] <- 1
  # clr_t4Th[rowSums(clr_t4Th)==0,3] <- 1
  
  if (saveandrm) {
    # dir.create(paste0(x2_dir,"/TbetCD45RA"), showWarnings=FALSE)
    # write.csv(csv_t4Th, file=gzfile(paste0(x2_dir,"/TbetCD45RA/",fid,".csv.gz")), row.names=FALSE)
    # dir.create(paste0(y2_dir,"/TbetCD45RA"), showWarnings=FALSE)
    # write.csv(clr_t4Th, file=gzfile(paste0(y2_dir,"/TbetCD45RA/",fid,".csv.gz")), row.names=FALSE)
    
    rm(cd4.Tcells.flowD) 
  }
  
  
  ## gating CD4+Tcell > Tregs naive ####
  cd4.T.Naive.temp <- rotate_fcs(flowDensity::getflowFrame(cd4.T.Naive.flowD),c(channels.ind["FoxP3"],channels.ind["CD25"]),theta=-pi/3)$data
  
  Tregs.Naive.flowD.temp <- flowDensity::flowDensity(
    cd4.T.Naive.temp, channels=c('Dy162Di', 'Tm169Di'), 
    position=c(F,T), gates=c(all.gthres[30], all.gthres[31]))
  Tregs.Naive.flowD <- flowDensity::flowDensity(
    Tregs.Naive.flowD.temp, channels=c('Dy162Di', 'Tm169Di'), 
    position=c(T,NA), gates=c(all.gthres[29], NA))
  
  Tregs.Naive.flowD@filter <- rotate_fcs(Tregs.Naive.flowD@filter,c(channels.ind["FoxP3"],channels.ind["CD25"]),theta=pi/3)$data
  Tregs.Naive.flowD <- flowDensity::flowDensity(
    cd4.T.Naive.flowD, channels=c('Dy162Di', 'Tm169Di'), 
    position=TRUE, filter=Tregs.Naive.flowD@filter)
  
  scat <- "07_FoxP3CD25_CD4Tcell"
  suppressWarnings({
    flowDensity::plotDens(
      cd4.T.Naive.flowD, c(channels.ind["FoxP3"],channels.ind["CD25"]), 
      main=paste0(scat,"\n[leaf: Tregs, other], CD4+ naive Tcell"), 
      cex.lab=2, cex.axis=2, cex.main=2)
    lines(Tregs.Naive.flowD@filter, lwd=2)
    
    clr_leuk[id_leuk%in%Tregs.Naive.flowD@index,"Tregs naive Tcell"] <- 1
    clr_leuk[id_leuk%in%cd4.T.Naive.flowD@index &
               !id_leuk%in%Tregs.Naive.flowD@index,"CD4+ naive Tcell other"] <- 1
  })
  
  csv_ <- csv_f[cd4.T.Naive.flowD@index, c("FoxP3","CD25"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=2)
  colnames(clr_) <- c("Tregs naive Tcell","CD4+ naive Tcell other")
  clr_[cd4.T.Naive.flowD@index%in%Tregs.Naive.flowD@index,1] <- 1
  clr_[clr_[,1]==0,2] <- 1
  
  filt_ <- list(Tregs.Naive.flowD@filter)
  names(filt_) <- colnames(clr_)[1]
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(cd4.T.Naive.flowD, Tregs.Naive.flowD) 
  }
  
  
  ## gating CD4+Tcell naive > Tregs naive, other ####
  gammaDelta.Tcells.flowD <- flowDensity::flowDensity(
    NOT.cd4.cd8.Tcells.flowD, channels=c('Sm152Di', 'Er170Di'), 
    position=c(T,NA), gates=c(all.gthres[35]+0.5,NA), ellip.gate=T)
  
  scat <- "08_TCRgdCD3_CD4Tcell"
  flowDensity::plotDens(
    NOT.cd4.cd8.Tcells.flowD, c(channels.ind["TCRgd"],channels.ind["CD3"]), 
    main=paste0(scat,"\n[leaf: gamma-delta Tcell], CD4-CD8- Tcell other"), 
    cex.lab=2, cex.axis=2, cex.main=2)
  lines(gammaDelta.Tcells.flowD@filter, lwd=2)
  
  clr_leuk[id_leuk%in%gammaDelta.Tcells.flowD@index,"gamma-delta Tcell"] <- 1
  clr_leuk[id_leuk%in%NOT.cd4.cd8.Tcells.flowD@index & 
             !id_leuk%in%gammaDelta.Tcells.flowD@index,"CD4-CD8- Tcell other"] <- 1
  
  csv_ <- csv_f[NOT.cd4.cd8.Tcells.flowD@index, c("TCRgd","CD3"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=2)
  colnames(clr_) <- c("gamma-delta Tcell","CD4-CD8- Tcell other")
  clr_[NOT.cd4.cd8.Tcells.flowD@index%in%gammaDelta.Tcells.flowD@index,1] <- 1
  clr_[clr_[,1]==0,2] <- 1
  
  filt_ <- list(gammaDelta.Tcells.flowD@filter)
  names(filt_) <- colnames(clr_)[1]

  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(NOT.cd4.cd8.Tcells.flowD, gammaDelta.Tcells.flowD) 
  }
  
  
  ## gating notCD4CD8Tcells > gammaDeltaTcells, other ####
  cd8.T.Naive.flowD <- flowDensity::flowDensity(
    cd8.Tcells.flowD, channels=c('Nd146Di', 'Nd143Di'), 
    position=c(NA, T), gates=c(NA, all.gthres[36]))
  cd8.T.Memory.flowD <- flowDensity::flowDensity(
    cd8.Tcells.flowD, channels=c('Nd146Di', 'Nd143Di'), 
    position=c(NA, F), gates=c(NA, all.gthres[36]))
  
  scat <- "09_CD8CD45RA_notCD4CD8Tcell_"
  flowDensity::plotDens(
    cd8.Tcells.flowD, c(channels.ind["CD8"],channels.ind["CD45RA"]), 
    main=paste0(scat,"\n[leaf: naive, memory CD8+ Tcell]"), cex.lab=2, cex.axis=2, cex.main=2)
  lines(cd8.T.Naive.flowD@filter, lwd=2)
  lines(cd8.T.Memory.flowD@filter, lwd=2)

  clr_leuk[id_leuk%in%cd8.T.Naive.flowD@index,"CD8+ naive Tcell"] <- 1
  clr_leuk[id_leuk%in%cd8.T.Memory.flowD@index,"CD8+ memory Tcell"] <- 1
  
  csv_ <- csv_f[cd8.Tcells.flowD@index, c("CD8","CD45RA"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=2)
  colnames(clr_) <- c("CD8+ naive Tcell","CD8+ memory Tcell")
  clr_[cd8.Tcells.flowD@index%in%cd8.T.Naive.flowD@index,1] <- 1
  clr_[cd8.Tcells.flowD@index%in%cd8.T.Memory.flowD@index,2] <- 1
  
  filt_ <- list(cd8.T.Naive.flowD@filter, cd8.T.Memory.flowD@filter)
  names(filt_) <- colnames(clr_)[1:2]
  
  gthres[[scat]] <- gate_ <- all.gthres[36]
  names(gthres[[scat]]) <- names(gate_) <- "CD45RA"
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    dir.create(paste0(thres_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    rm(cd8.T.Naive.flowD, cd8.T.Memory.flowD) 
  }
  
  
  ## gating CD8+Tcell > CD8+Tcell naive, CD8+Tcell memory ####
  cd25.cd8.T.Naive.flowD <- flowDensity::flowDensity(
    cd8.Tcells.flowD, channels=c('Gd160Di', 'Nd143Di'), 
    position=c(T, T), gates=c(all.gthres[37], all.gthres[36]))
  cd25.cd8.T.Memory.flowD <- flowDensity::flowDensity(
    cd8.Tcells.flowD, channels=c('Gd160Di', 'Nd143Di'), 
    position=c(T, F), gates=c(all.gthres[37], all.gthres[36]))
  
  scat <- "10_TbetCD45RA2_CD8Tcell_"
  flowDensity::plotDens(
    cd8.Tcells.flowD, c(channels.ind["Tbet"], channels.ind["CD45RA"]), 
    main=paste0(scat,"\nCD25+CD8+ T Naive & CD25+CD8+ T Memory, other"), 
    cex.lab=2, cex.axis=2, cex.main=2)
  lines(cd25.cd8.T.Naive.flowD@filter, lwd=2)
  lines(cd25.cd8.T.Memory.flowD@filter, lwd=2)
  
  # clr_leuk[id_leuk%in%cd25.cd8.T.Naive.flowD@index,"CD25+CD8+ naive Tcell"] <- 1
  # clr_leuk[id_leuk%in%cd25.cd8.T.Memory.flowD@index,"CD25+CD8+ memory Tcell"] <- 1
  # clr_leuk[id_leuk%in%cd8.Tcells.flowD &
  #            !id_leuk%in%cd25.cd8.T.Memory.flowD@index &
  #            !id_leuk%in%cd25.cd8.T.Naive.flowD@index,"CD8+ Tcell other"] <- 1
  
  csv_ <- csv_f[cd8.Tcells.flowD@index, c("Tbet","CD45RA"), drop=FALSE]
  clr_ <- matrix(0, nrow=nrow(csv_), ncol=3)
  colnames(clr_) <- c("CD25+CD8+ naive Tcell","CD25+CD8+ memory Tcell", "CD25+CD8+ Tcell other")
  clr_[cd8.Tcells.flowD@index%in%cd25.cd8.T.Naive.flowD@index,1] <- 1
  clr_[cd8.Tcells.flowD@index%in%cd25.cd8.T.Memory.flowD@index,2] <- 1
  clr_[rowSums(clr_)==0,3] <- 1
  
  filt_ <- list(cd25.cd8.T.Naive.flowD@filter, cd25.cd8.T.Memory.flowD@filter)
  names(filt_) <- colnames(clr_)[1:2]
  
  gthres[[scat]] <- gate_ <- all.gthres[37:36]
  names(gthres[[scat]]) <- names(gate_) <- c("Tbet","CD45RA")
  
  if (saveandrm) {
    dir.create(paste0(x2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(csv_, file=gzfile(paste0(x2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    dir.create(paste0(y2_dir,"/",dset,"/",scat), showWarnings=FALSE)
    write.csv(clr_, file=gzfile(paste0(y2_dir,"/",dset,"/",scat,"/",fid,".csv.gz")), row.names=FALSE)
    
    dir.create(paste0(filt_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(filt_, file=paste0(filt_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    dir.create(paste0(thres_dir,"/",dset,"/",scat), showWarnings=FALSE)
    save(gate_, file=paste0(thres_dir,"/",dset,"/",scat,"/",fid,".Rdata"))
    
    leuk_rows <- rowSums(clr_leuk)>0
    csv_leuk <- csv_leuk[leuk_rows,]
    clr_leuk <- clr_leuk[leuk_rows,]
    
    write.csv(csv_leuk, file=gzfile(paste0(xn_dir,"/",dset,"/",fid,".csv.gz")), row.names=FALSE)
    write.csv(clr_leuk, file=gzfile(paste0(yn_dir,"/",dset,"/",fid,".csv.gz")), row.names=FALSE)
    
    rm(cd8.Tcells.flowD, cd25.cd8.T.Naive.flowD, cd25.cd8.T.Memory.flowD, csv_leuk, clr_leuk)
  }
  
  # save(gthres, file=paste0(thres_dir,"/",fid,".Rdata"))
  # if (any(is.na(unlist(gthres)))) 
  #   cat("\n",fid,": ",unlist(gthres)[is.na(unlist(gthres))])
  
  graphics.off()
  time_output(start1)
  
}) }) }, .parallel=TRUE)

time_output(start)
