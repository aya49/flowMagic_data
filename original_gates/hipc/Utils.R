
# -------------------- importiamo le funzioni del flowPrep code -------------
import_flowPrep <- function(){
  path_flowPre_R <- "./flowPrep.R" 
  source(path_flowPre_R) 
}


#---------------------- import files on which we work: version 1-------------------
import_files <- function(n_sample,sample_type){
  if(sample_type=="Myeloid_Gambia"){
    files_path <- paste0("/home/rstudio/data/Gambia validation cohort (Aug2018)/Myeloid_Gambia_confcohort/Myeloid_samples_fcsfiles")
  }
  else if(sample_type=="Myeloid_PNG"){
    files_path <- paste0("/home/rstudio/data/PNG pilot/")
  }
  else if(sample_type=="Myeloid_PNG_HIPC_pilot"){
    files_path <- paste0("/home/rstudio/data/HIPC pilots 20160307/Myeloid Panel/HIPC Myeloid Flowjo")
  }
  else if(sample_type=="Bcells_Gambia"){
    files_path <- paste0("/home/rstudio/data/Gambia validation cohort (Aug2018)/Bcell_Gambia_confcohort/Bcells_fcsfiles_samples/")
  } 
  else if(sample_type=="Bcells_PNG"){
    files_path <- paste0("/home/rstudio/data/PNG B cell panel (Extended)/PNG B cells expended panel 01032018/")
  }else{
    stop("please select valid sample_type: Myeloid_Gambia Myeloid_PNG Bcells_Gambia Bcells_PNG")
  }
  
    
  FCSfiles_path_list <- list.files(path = files_path, recursive = F, pattern = ".fcs", full.names = T)
  if(sample_type=="Myeloid_PNG_HIPC_pilot"){
    FCSfiles_path_list <- list.files(path = files_path, recursive = F, pattern = "Specimen", full.names = T)
  }
  # recursive=Should the listing recurse into directories?
  # pattern=".fcs" we take only the fcs fies
  # full_names=T, we take the full path of the files
  if(class(n_sample)=="character"){
    if(n_sample=="All"){
      n_sample<-1:length(FCSfiles_path_list)
    }else{
      stop("please specify a correct value for n_sample in import_files() function")
    }
  }
  fcs_files_list <-lapply(FCSfiles_path_list[n_sample],read.FCS) 
  fs <- as(fcs_files_list,"flowSet")
  return(list(fs=fs,FCSfiles_path_list=FCSfiles_path_list[n_sample])) 
}

# ---------- import files version 2: general importing ---------------
# files_path specify the directory that contains the files of interest
import_files_v2 <- function(n_sample,files_path,n_cores){
  if(is.null(files_path)){
    stop("Directory path that contains the fcs files must be specified")
  }
  FCSfiles_path_list <- list.files(path = files_path, recursive = F, pattern = ".fcs", full.names = T)
  # be sure that no compensation control is taken (especially when we analyze the comp files folder)
  inds<-grep("Stained Control",FCSfiles_path_list)
  if(length(inds)>0){
    FCSfiles_path_list<-FCSfiles_path_list[-inds]
  }
  if(is.null(FCSfiles_path_list)){
    stop("No .fcs file founded in the files_path,please specify a correct path")
  }
  if(class(n_sample)=="character"){
    if(n_sample=="All"){
      n_sample<-1:length(FCSfiles_path_list)
    }else{
      stop("please specify a correct value for n_sample in import_files() function")
    }
  }
  fcs_files_list <-mclapply(FCSfiles_path_list[n_sample],read.FCS,mc.cores = n_cores)
  #print(fcs_files_list)
  fs <- as(fcs_files_list,"flowSet")

  return(list(fs=fs,FCSfiles_path_list=FCSfiles_path_list[n_sample]))
}

# -------------------- Make the directories to save the results ------------------
create_path_output_dir<-function(sample_type){
  if(sample_type=="Myeloid_Gambia"){
    path.output <- "/home/rstudio/results/results_Myeloid_Gambia"
  }
  else if(sample_type=="Myeloid_PNG"){
    path.output <- "/home/rstudio/results/results_Myeloid_PNG"
  }
  else if(sample_type=="Bcells_Gambia"){
    path.output <- "/home/rstudio/results/results_Bcells_Gambia"
  } 
  else if(sample_type=="Bcells_PNG"){
    path.output <- "/home/rstudio/results/results_Bcells_PNG"
  } else{
    stop("please select valid sample_type: Myeloid_Gambia Myeloid_PNG Bcells_Gambia Bcells_PNG")
  }
  path.output_2<-"/home/rstudio/results/flowType_files"
  suppressWarnings ( dir.create ( paste0(path.output,"/Plots/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output,"/Cleaning/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output,"/Stats/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output,"/WSP/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output_2,"/df_thresholds/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output_2,"/fcs_files/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output_2,"/filter_objects/"),recursive = T ))
  return(path.output)
}

# ----------------- function to make the scatter plot with density and counter lines --------------
# the input is  a flowFrame,channels indicates the channels based on which we want to make the counter plot.
# main indicate the plot title
# numlevels indicate the number of "levels", meaning counter lines
#  plotDens = Generate a scatter dot plot with colors based on the distribution of the density 
# of the provided channels. 
plotDensContour <- function(frame, channels, main, numlevels = 10){
  plotDens(frame, channels = channels, main = main, cex.lab = 2, cex.axis = 2, cex.main=2)
  data.new <- na.omit(exprs(frame)[, channels])
  # get the bivariate density  of this cleaned matrix considering as variable the first and second column. n indicates the points la grid matrix
  z <- kde2d(data.new[, 1], data.new[, 2], n = 50)
  # disegnamo le counter lines (sopra al plot di plotDens),fornendo a counter l'oggetto z ottenuto prima,
  # che contiene le coordinate della grid matrix() e il valore della densità bivariata.
  contour(z, drawlabels = FALSE, add = TRUE, nlevels = numlevels, lty = 2)
} 






# basename removes all of the path up to and including the last path separator (if any).
#----- make the GatingSet ------
get_GatingSet<- function(fs,FCSfiles_path_list){
  fsnames <-fname.unique <- basename(FCSfiles_path_list)
  gs <- GatingSet(fs)
  sampleNames(gs)<-fsnames
  return(gs)
}
#------------- import comp_matrix -----------
import_comp_matrix <- function(sample_type){
  if(sample_type=="Myeloid_Gambia"){
    path_comp_matrix <- paste0("/home/rstudio/data/Gambia validation cohort (Aug2018)/Myeloid_Gambia_confcohort")
    comp.mat <- read.csv(list.files(paste0(path_comp_matrix),pattern="Myeloid_comp.csv",full.names = T),check.names = F)
  }else if(sample_type=="Myeloid_PNG"){
    path_comp_matrix <- "/home/rstudio/data/Comp Matrix_PNG_Feb2018-Updated.csv"
    comp.mat <- read.csv(path_comp_matrix,check.names = F)
  }
  else if(sample_type=="Myeloid_PNG_HIPC_pilot"){
    path_comp_matrix <- paste0("/home/rstudio/data/HIPC pilots 20160307/Myeloid Panel/CompMatrix.csv")
    comp.mat <- read.csv(path_comp_matrix,check.names = F)
    comp.mat <- comp.mat[2:length(colnames(comp.mat))]
  }
  else if(sample_type=="Bcells_Gambia"){
    path_comp_matrix <- paste0("/home/rstudio/data/Gambia validation cohort (Aug2018)/Bcell_Gambia_confcohort/")
    comp.mat <- read.csv(paste0(path_comp_matrix,"Bcell_comp.csv"),check.names = F)
    comp.mat <-comp.mat[2:length(comp.mat)]
    ind <- which(colnames(comp.mat)=="PE-Cy5-A :: CD38")
    colnames(comp.mat)[ind]<-"PE-Cy5-A"
    #colnames(comp)<-colnames(fs[[1]])[c(7:19)]
  }
  else if(sample_type=="Bcells_PNG"){
    path_comp_matrix <- paste0("/home/rstudio/data/Sample info/")
    comp.mat <- read.csv(paste0(path_comp_matrix,"Bcell-comp.csv"),check.names = F)
  }
  else{
    stop("please select a valid sample_type for the comp matrix")
  }
  return(comp.mat)
}

# ---------- general version of import comp matrix function ---------------
# path_comp_matrix specify the absolute path for the compensation matrix inside the docker
import_comp_matrix_v2<-function(path_comp_matrix){
  print("Import compensation matrix...")
  if(is.null(path_comp_matrix)){
    stop("Please specify a valid path for the compensation matrix")
  }
  comp.mat<-read.csv(path_comp_matrix,check.names = F)
  classes_v<-as.vector(sapply(comp.mat,class))
  indx<-which(classes_v!="numeric")
  if(length(indx)!=0){
    comp.mat<-comp.mat[,-indx]
  }
  ind <- grep("::",colnames(comp.mat))
  if(length(ind)!=0){
    warning("incorrect format of compensation matrix colnames detected,attempt to correct it...",immediate. = T)
    v<-strsplit(colnames(comp.mat)[ind]," :: ")
    v2<-sapply(1:length(v),function(i){
      ind2<-grep("-A",v[[i]])
      new_name<-v[[i]][ind2]
      return(new_name)
    })
    colnames(comp.mat)[ind]<-v2
    ind <- grep("::",colnames(comp.mat))
    if(length(ind)!=0){
      stop("Attempt to correct the colnames failed. Remove manually :: format")
    }
    return(comp.mat)
   }else{
     return(comp.mat)
   }
  
}
 

#------------------- find markers indices myeloid Gambia --------------------------------
find_markers_indices<-function(fs){
  # We find the markers'indices indicated by the user
  # grep(value = FALSE) returns a vector of the indices of the elements of x that yielded a match in y
  f<-fs[[1]]
  Time.channel <- grep('Time', pData(parameters(f))$name)
  # The which() function will return the position of the elements(i.e., row number/column number/array index) 
  # in a logical vector which are TRUE.
  fluorochrome.chans <- which(pData(parameters(f))$desc != "<NA>")
  names(fluorochrome.chans) <- pData(parameters(f))$desc[fluorochrome.chans]
  names(fluorochrome.chans)[which(tolower(names(fluorochrome.chans))=="viability dye")]<-"Live"
  #----------------------------------------------------------
  # Checking to be sure that specific channels names have always the same structure
  ind<-grep(names(fluorochrome.chans),pattern="HLA-*DR") # with * we consider zero or occurence of the "-" character
  names(fluorochrome.chans)[ind]<-"HLA-DR"
  ind<-which(names(fluorochrome.chans)=="HLA-DR")
  if(length(ind)==0){
    # I suppose that the HLA-DR is associated always to the eF605-A channel
    ind<-which(pData(parameters(f))$name == "eF605-A")
    fluorochrome.chans["HLA-DR"]<-ind
    warning("The channel that refers to the marker name 'HLA-DR' does not exist, the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="gd")
  names(fluorochrome.chans)[ind]<-"gd"
  ind<-which(names(fluorochrome.chans)=="gd")
  if(length(ind)==0){
    # I suppose that the gd is associated always to the PE-A channel
    ind<-which(pData(parameters(f))$name == "PE-A")
    fluorochrome.chans["gd"]<-ind
    warning("The channel that refers to the marker name 'gd' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD11c")
  if(length(ind)==0){
    # I suppose that the CD11c is associated always to the APC-A channel
    ind<-which(pData(parameters(f))$name == "APC-A")
    fluorochrome.chans["CD11c"]<-ind
    warning("The channel that refers to the marker name 'CD11c' does not exist,the entry has been fixed")
    
  } 
  ind<-grep(names(fluorochrome.chans),pattern="CD64")
  if(length(ind)==0){
    # I suppose that the CD64 is associated always to the Alexa fluor channel
    ind<-which(pData(parameters(f))$name == "Alexa Fluor 700-A")
    fluorochrome.chans["CD64"]<-ind
    warning("The channel that refers to the marker name 'CD64' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD14")
  if(length(ind)==0){
    # I suppose that the CD14 is associated always to the V500-A channel
    ind<-which(pData(parameters(f))$name == "V500-A")
    fluorochrome.chans["CD14"]<-ind
    warning("The channel that refers to the marker name 'CD14' does not exist,the entry has been fixed")
    
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD56")
  if(length(ind)==0){
    # I suppose that the CD56 is associated always to the BV650-A channel
    ind<-which(pData(parameters(f))$name == "BV650-A")
    fluorochrome.chans["CD56"]<-ind
    warning("The channel that refers to the marker name 'CD56' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD11b")
  if(length(ind)==0){
    # I suppose that the CD11b is associated always to the BV786-A channel
    ind<-which(pData(parameters(f))$name == "BV786-A")
    fluorochrome.chans["CD11b"]<-ind
    warning("The channel that refers to the marker name 'CD11b' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD3")
  if(length(ind)==0){
    # I suppose that the CD3 is associated always to the PE-CF594-A channel
    ind<-which(pData(parameters(f))$name == "PE-CF594-A")
    fluorochrome.chans["CD3"]<-ind
    warning("The channel that refers to the marker name 'CD3' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD123")
  if(length(ind)==0){
    # I suppose that the CD123 is associated always to the PE-Cy7-A channel
    ind<-which(pData(parameters(f))$name == "PE-Cy7-A")
    fluorochrome.chans["CD123"]<-ind
    warning("The channel that refers to the marker name 'CD123' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD66")
  if(length(ind)==0){
    # I suppose that the CD66 is associated always to the BV711-A channel
    ind<-which(pData(parameters(f))$name == "BV711-A")
    fluorochrome.chans["CD66"]<-ind
    warning("The channel that refers to the marker name 'CD66' does not exist,the entry has been fixed")
  }  
  ind<-grep(names(fluorochrome.chans),pattern="Live")
  if(length(ind)==0){
    # I suppose that the Viability die is associated always to the APC-eF780-A channel
    ind<-which(pData(parameters(f))$name == "APC-eF780-A")
    fluorochrome.chans["Live"]<-ind
    warning("The channel that refers to the marker name 'Viability die' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD45")
  if(length(ind)==0){
    # I suppose that the CD45 is associated always to the V450-A channel
    ind<-which(pData(parameters(f))$name == "V450-A")
    fluorochrome.chans["CD45"]<-ind
    warning("The channel that refers to the marker name 'CD45' does not exist,the entry has been fixed")
  }
  ind<-grep(names(fluorochrome.chans),pattern="CD16")
  if(length(ind)==0){
    # I suppose that the CD16 is associated always to the FITC-A channel
    ind<-which(pData(parameters(f))$name == "FITC-A")
    fluorochrome.chans["CD16"]<-ind
    warning("The channel that refers to the marker name 'CD16' does not exist,the entry has been fixed")
  }
  # --------------------------------------------------------------------------------------
  scat.chans <- c(grep(colnames(f),pattern = "FSC*"), grep(colnames(f),pattern = "SSC*"))
  names(scat.chans) <- colnames(f)[scat.chans]
  # if lenght is more than 1,it means that there is more than 1 index for CD16 marker(and this is no possible,so we need to change the name of one of the two markers)
  if (length(which(names(fluorochrome.chans)=="CD16"))>1)
    names(fluorochrome.chans)[which(names(fluorochrome.chans)=="CD16")[1]]<-"CD45"
  return(list(Time.channel=Time.channel,fluorochrome.chans=fluorochrome.chans,scat.chans=scat.chans))
}

#------------------- find markers indices B cell Gambia --------------------------------

find_markers_indices_bcells<-function(){
  markers<- c("viability dye","CD66","CD14","CD34","HLADR","CD19","CD38","CD138","CD20","IgM","IgD","CD27","CD10")
  channels.ind <- Find.markers(fs[[1]],markers) 
  return(channels.ind)
}

find_markers_indices_bcells_v2<-function(fs){
  f<-fs[[1]]
  Time.channel <- grep('Time', pData(parameters(f))$name)
  # The which() function will return the position of the elements(i.e., row number/column number/array index) 
  # in a logical vector which are TRUE.
  fluorochrome.chans <- which(pData(parameters(f))$desc != "<NA>")
  names(fluorochrome.chans) <- pData(parameters(f))$desc[fluorochrome.chans]
  names(fluorochrome.chans)[which(tolower(names(fluorochrome.chans))=="viability dye")]<-"Live"
  #----------------------------------------------------------
  # Checking to be sure that specific channels names have always the same structure
  ind<-grep(names(fluorochrome.chans),pattern="HLA-*DR") # with * we consider zero or occurence of the "-" character
  names(fluorochrome.chans)[ind]<-"HLADR"
  ind<-which(names(fluorochrome.chans)=="HLADR")
  if(length(ind)==0){
    # I suppose that the HLA-DR is associated always to the BV605-A channel
    ind<-which(pData(parameters(f))$name == "BV605-A")
    fluorochrome.chans["HLA-DR"]<-ind
    warning("The channel that refers to the marker name 'HLA-DR' does not exist, the entry has been fixed")
  }
  #-----------------------------------------------------------
  scat.chans <- c(grep(colnames(f),pattern = "FSC*"), grep(colnames(f),pattern = "SSC*"))
  names(scat.chans) <- colnames(f)[scat.chans]
  return(list(Time.channel=Time.channel,fluorochrome.chans=fluorochrome.chans,scat.chans=scat.chans))
  
}


#------------------------ Check fcs quality ---------------------------
check_fcs_quality<-function(){
  markers_fs_1<-as.vector(fs[[1]]@parameters@data$desc)
  if(all(is.na(markers_fs_1))){
    stop("Description column of your input fcs file contains only NA,please check the files")
  }
}

################################################################################################
################################ function to generate plots and final data ####################
################################################################################################
#------ Data to return --------------------------------------------------------------------
generate_final_data<-function(path.output,name_output,gs){
  freqs<-gs_pop_get_count_fast(gs, format="wide",statistic="freq",path = "auto")
  counts<-gs_pop_get_count_fast(gs,path = "auto")
  write.csv(counts, paste0(path.output,"/","Stats","/",sprintf("pops_counts_%s.csv",name_output)))
  write.csv(freqs, paste0(path.output,"/","Stats","/",sprintf("pops_freqs_%s.csv",name_output)))
  path.output_final_data <- paste0(path.output,"/","WSP","/")
  GatingSet2flowJo(gs, outFile = paste0(path.output_final_data,name_output,"_AutomatedWSP.wsp"))
}


generate_plots_myeloid<-function(gs,n_cores,fs.marg,fs.clean,singlets,fs.sngl,
                                 beads_poly,size_poly,fs.size,live_poly,
                                 fs.live,cd66_45_poly,gran_poly,fs.11p16p,CD64_gate,
                                 fs.gran,CD11bgranulo_gate,fs.66n45p,list_DRn_14n,list_DRp_14n,
                                 list_monocytes,fs.drn.14n,list_CD3_gate,fs.tcell,
                                 list_gdt_poly,list_gdneg_poly,list_all_nkt_related_pops,
                                 fs.nk.56n16n,baso_poly,fs.cd3n,list_nk_56n16n,list_NK_temp,
                                 list_nk_56hi,list_mDC,list_pDC,pdc_poly,mdc_poly,bcell_poly,
                                 fs.mono,CD16mon_gate,path.output,Time.channel,fluorochrome.chans,
                                 scat.chans,clean.inds,CD16granulo_gate,fs){
  
  #--------- generate the plots
  n_plots_processed<-mclapply(1:length(gs),function(i){
    png(file = paste(paste0(path.output,"/Plots/"), gs[[i]]@name, '.png',sep=""),  width = 3000, height = 3000*3.6/4)
    # par regulates the graphic device parameters and make a graphic device.
    #  mfrow = A vector of the form c(nr, nc). Subsequent figures will be drawn in an nr-by-nc array on the device by columns (mfcol), or rows (mfrow), respectively.
    # mar = A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1
    par(mfrow=c(4,4),mar=(c(5, 5, 4, 2) + 0.1))
    # "par()" works only with images produced by plot, you cannot use it with autoplot or ggcyto
    # -------- first row plots
    ### 1
    plotDens(fs.marg[[i]], channels = c(Time.channel, fluorochrome.chans['CD11c']), main = "All Events", cex.lab = 2, cex.axis = 2, cex.main=2)
    points(fs.marg[[i]]@exprs[clean.inds[[i]]$ind, c(Time.channel, fluorochrome.chans['CD11c'])], pch ='.')
    ###2
    plotDens(fs.clean[[i]], channels = c(scat.chans['FSC-A'], scat.chans['FSC-H']), main = "Time(cleaned)", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(singlets[[i]],lwd=2) # ricordati che singlets e' direttamente il filter non il Cell population object
    # add labels
    text(x=max(singlets[[i]][,1])*0.90,y=max(singlets[[i]][,2])*0.90,labels="Singlets",cex = 3.5) # moltiplichiamo per 0.90 per spostare il testo rispetto al confine
    ###3
    plotDens(fs.sngl[[i]], channels = c(scat.chans['FSC-A'], scat.chans['SSC-A']), main = "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(beads_poly[[i]]@boundaries,lwd=2)
    lines(size_poly[[i]]@boundaries,lwd=2) # al posto di usare size.filter usiamo size.poly@boundaries che e' la stessa cosa
    # add labels
    text(x=max(beads_poly[[i]]@boundaries[,1])*0.90,y=max(beads_poly[[i]]@boundaries[,2])*0.95,labels = "Beads",cex=4)
    text(x=max(size_poly[[i]]@boundaries[,1])*0.80,y=max(size_poly[[i]]@boundaries[,2])*0.90,labels = "All cells",cex=4)
    
    ###4
    plotDens(fs.size[[i]], channels = c(fluorochrome.chans['Live'], scat.chans['SSC-A']), main = "All cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(live_poly[[i]]@boundaries, lty=2,lwd=2) # al posto di live.filter usiamo live_poly
    # add labels
    text(x=max(live_poly[[i]]@boundaries[,1])*0.80,y=max(live_poly[[i]]@boundaries[,2])*0.70,labels = "Live cells",cex=4)
    # ------ second row plots
    ###1
    plotDens(fs.live[[i]], channels = c(fluorochrome.chans['CD66'], fluorochrome.chans['CD45']), main = "Live cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd66_45_poly[[i]]@boundaries,lwd=2)
    lines(gran_poly[[i]]@boundaries,lwd=2)
    
    # text(coordinates of start position x and y,other arguments....)
    text(x=median(gran_poly[[i]]@boundaries[, 1]), y=median(gran_poly[[i]]@boundaries[, 2]), labels = 'Granulocytes', cex = 2)
    text(x=median(cd66_45_poly[[i]]@boundaries[,1]),y=median(cd66_45_poly[[i]]@boundaries[,2]),labels="Non Granulocytes",cex=2)

    
    ###2
    # not present in original gating strategy
    dens_x<-density(exprs(fs.11p16p[[i]][,fluorochrome.chans["CD64"]]))$x
    dens_y<-density(exprs(fs.11p16p[[i]][,fluorochrome.chans["CD64"]]))$y
    plot(dens_x,dens_y,type="l",main="CD11b+CD16+ Mature Neutrophils",xlab="CD64",ylab="density",cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=CD64_gate[[i]],lwd=2)
    # add labels
    text(x=CD64_gate[[i]]*1.50,y=0.8,labels="CD64+",cex=4)
    text(x=CD64_gate[[i]]*0.50,y=0.8,labels="CD64-",cex=4)
    ###3
    plotDens(fs.gran[[i]], channels = c(fluorochrome.chans['CD11b'], fluorochrome.chans['CD16']), main = "Granulocytes", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = CD16granulo_gate[[i]],lwd=2)
    abline(v = CD11bgranulo_gate[[i]],lwd=2)
    # add labels
    text(x=CD11bgranulo_gate[[i]]*1.30,y=CD16granulo_gate[[i]]*1.50,labels="Mature Neutrophils",cex=1.7)
    text(x=CD11bgranulo_gate[[i]]*0.80,y=CD16granulo_gate[[i]]*0.50,labels="Immature Neutrophils 1",cex=1.7)
    text(x=CD11bgranulo_gate[[i]]*1.30,y=CD16granulo_gate[[i]]*0.50,labels="CD16-CD11b+ Granulocytes",cex=1.7)
    text(x=CD11bgranulo_gate[[i]]*0.80,y=CD16granulo_gate[[i]]*1.50,labels="Immature Neutrophils 2 ",cex=1.7)
    ###4
    plotDensContour(fs.66n45p[[i]], channels = c(fluorochrome.chans['HLA-DR'], fluorochrome.chans['CD14']), main = "CD45+CD66-")
    lines(list_DRn_14n[[i]]@filter,lwd=2)
    lines(list_DRp_14n[[i]]@filter,lwd=2)
    lines(list_monocytes[[i]]@filter,lwd=2)
    # add labels
    text(x=min(list_DRn_14n[[i]]@filter[,1])*0.80,y=max(list_DRn_14n[[i]]@filter[,2])*0.80,labels="HLADR- CD14-",cex=2) # in questo caso il minimo di x e' negativo quindi per falro avanzare "avanti",cioe' piu' a destr dobbiamo moltiplicare per 0<x<1.
    text(x=min(list_DRp_14n[[i]]@filter[,1])*1.30,y=median(list_DRp_14n[[i]]@filter[,2]),labels="HLADR+ CD14-",cex=2)
    text(x=median(list_monocytes[[i]]@filter[,1]),y=min(list_monocytes[[i]]@filter[,2])*1.30,labels="HLADR+ CD14+",cex=2)
    
    
    # ---- third row plots
    ###1
    plotDens(fs.drn.14n[[i]], channels = c(fluorochrome.chans['CD3'], scat.chans['SSC-A']), main = "HLADR-CD14-", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = list_CD3_gate[[i]],lwd=2)
    # add labels
    text(x=list_CD3_gate[[i]]*0.50,y=50000,labels="CD3-",cex=2)
    text(x=list_CD3_gate[[i]]*1.50,y=50000,labels="CD3+",cex=2)
    ###2
    plotDens(fs.tcell[[i]], channels = c(fluorochrome.chans['CD3'], fluorochrome.chans['gd']),  main = "CD3+ T cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(list_gdt_poly[[i]]@boundaries,type="l",lwd=2)
    lines(list_gdneg_poly[[i]]@boundaries,type="l",lwd=2)
    # add labels
    text(x=median(list_gdt_poly[[i]]@boundaries[,1]),y=median(list_gdt_poly[[i]]@boundaries[,2]),labels="gdTCells",cex=2)
    text(x=median(list_gdneg_poly[[i]]@boundaries[,1]),y=median(list_gdneg_poly[[i]]@boundaries[,2]),labels="gd- T cells",cex=2)
    
    ###3
    plotDens(fs.tcell[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']),  main = "CD3+ T cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = list_all_nkt_related_pops[[i]]$CD16NKT.gate,lwd=2)
    abline(h = list_all_nkt_related_pops[[i]]$CD56NKT.gate,lwd=2)
    # add labels
    text(x=list_all_nkt_related_pops[[i]]$CD16NKT.gate*1.30,y=list_all_nkt_related_pops[[i]]$CD56NKT.gate*1.50,labels="CD56+16+ NKT cells",cex=2)
    text(x=list_all_nkt_related_pops[[i]]$CD16NKT.gate*0.80,y=list_all_nkt_related_pops[[i]]$CD56NKT.gate*0.50,labels="CD56-CD16-",cex=2)
    text(x=list_all_nkt_related_pops[[i]]$CD16NKT.gate*1.30,y=list_all_nkt_related_pops[[i]]$CD56NKT.gate*0.50,labels="CD56-16+",cex=2)
    text(x=list_all_nkt_related_pops[[i]]$CD16NKT.gate*0.80,y=list_all_nkt_related_pops[[i]]$CD56NKT.gate*1.50,labels="CD56+CD16-",cex=2)
    ###4
    plotDens(fs.nk.56n16n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['HLA-DR']),  main = "CD56-CD16- cells", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(baso_poly[[i]]@boundaries,lwd=2) # al posto di basophils@filter usiamo baso_poly@boundaries
    text(x=median(baso_poly[[i]]@boundaries[,1]),y=median(baso_poly[[i]]@boundaries[,2]),labels="Basophils",cex=2)
    
    
    
    # ----- fourth row plots
    ###1
    plotDens(fs.cd3n[[i]], channels = c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56']),  main = "CD3-", cex.lab = 2, cex.axis = 2, cex.main=2)
    data.new <- fs.cd3n[[i]]@exprs[, c(fluorochrome.chans['CD16'], fluorochrome.chans['CD56'])]
    z <- kde2d(data.new[, 1], data.new[, 2])
    contour(z, drawlabels = FALSE, add = TRUE, lty = 2,lwd=2)
    abline(h = list_nk_56n16n[[i]]@gates[2],lwd=2)
    segments(par('usr')[1],  list_nk_56hi[[i]]@gates[2], list_NK_temp[[i]]@gates[1],list_nk_56hi[[i]]@gates[2],lwd=2)
    segments(list_NK_temp[[i]]@gates[1],list_nk_56hi[[i]]@gates[2],list_NK_temp[[i]]@gates[1], par('usr')[4],lwd = 2)
    segments(list_nk_56n16n[[i]]@gates[1],list_nk_56n16n[[i]]@gates[2],list_nk_56n16n[[i]]@gates[1],par('usr')[3],lwd=2)
    segments(list_NK_temp[[i]]@gates[1],list_nk_56n16n[[i]]@gates[2],list_NK_temp[[i]]@gates[1],list_nk_56hi[[i]]@gates[2],lwd=2)
    # add text
    text(par('usr')[1] + 0.35, par('usr')[4] - 0.25, labels = 'CD56 high NK', cex = 2)
    text(par('usr')[2] - 0.65, par('usr')[4] - 0.25, labels = ' CD56dim CD16+ NK', cex = 2)
    text(par('usr')[1] + 0.35, list_nk_56hi[[i]]@gates[2] - 0.25, labels = ' CD56dim CD16- NK', cex = 2)
    text(par('usr')[1] + 0.35, par('usr')[3] + 0.25, labels = 'CD16-CD56- cells', cex = 2)
    text(par('usr')[2] - 0.5,  par('usr')[3] + 0.25, labels = 'CD56-CD16+ NK', cex = 2)
    
    ###2
    # Draw the mother pop HLADR+CD14- with the children pops: mDC,pDC,Bcells.
    # calculate cells proportion 
    #ncells <- round(list_DRp_14n[[i]]@cell.count/beads[[i]]@cell.count, 4)
    plotDens(list_DRp_14n[[i]], channels = c(fluorochrome.chans['CD123'], fluorochrome.chans['CD11c']),
             main = "HLADR+CD14-", cex.lab = 2, cex.axis = 2, cex.main=2,
             xlab =sprintf("%s:CD123", pData(parameters(fs[[1]]))$name[fluorochrome.chans['CD123']]),
             ylab = sprintf("%s:CD11c", pData(parameters(fs[[1]]))$name[fluorochrome.chans['CD11c']]))
             
    abline(h=list_mDC[[i]]@gates[2])
    segments(x0=list_pDC[[i]]@gates[1],y0=list_mDC[[i]]@gates[2],x1=list_pDC[[i]]@gates[1],y1=par('usr')[3]) #  par('usr')[3] indicate the ymin of the y -axis
    segments(x0=list_mDC[[i]]@gates[1],y0=list_mDC[[i]]@gates[2],x1=list_mDC[[i]]@gates[1],y1=par('usr')[4]) #  par('usr')[4] indicate the ymax of the y -axis
    # add text
    text(median(pdc_poly[[i]]@boundaries[,1]), median(pdc_poly[[i]]@boundaries[,2]), labels = 'pDC', cex = 2) # posizione testo poco prima del massimo valore del boundaries di CD123 e poco prima del massimo valore dei boundareis di CD11c
    text(median(mdc_poly[[i]]@boundaries[,1]) - 0.2, median(mdc_poly[[i]]@boundaries[,2]), labels = 'mDC', cex = 2) # stessa cosa per mDC
    # add text
    text(median(bcell_poly[[i]]@boundaries[,1]), median(bcell_poly[[i]]@boundaries[,2]), labels = 'B cells', cex = 2)
    
    ###3
    
    plotDensContour(fs.mono[[i]], channels = c(fluorochrome.chans['CD14'], fluorochrome.chans['CD16']), main = "HLADR+ CD14+ Monocytes") #, cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = CD16mon_gate[[i]],lwd=2)
    # usr = A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of 
    # the plotting region. In other words: to get the range of the x and y axes in the plot
    text(par('usr')[1] + 0.25, par('usr')[3] + 0.45, labels = 'Classical monocytes', cex = 2)
    text(par('usr')[1] + 0.25, CD16mon_gate[[i]] + 0.45, labels = 'Non classical monocytes', cex = 2)
    
    # NOTE!: Remember to put dev.off() at the end or you will obtain an empty file
    dev.off()
    return(i)
  },mc.cores=n_cores)
}



generate_plots_Bcells <- function(gs,n_cores,path.output,fs.marg,clean.inds,fs.clean,singlets,fs.sngl,beads_poly,size_poly,
                                  fs.size,live_poly,fs.live,gran_poly,non_gran_poly,fs.non_gran,Bcells_poly,
                                  fs.bcells,pc_poly,blast_poly,CD20pos_poly,CD20neg_poly,quad1,PB_poly,
                                  fs.CD20pos,CD10neg_poly,CD10pos_poly,fs.igm,trans_poly,fs.cd10neg,quad.3,
                                  fluorochrome.chans,scat.chans){
  n_plots_processed<-mclapply(1:length(gs),function(i){
    #--------- make the matrix of the plots of the figure
    png(paste(paste0(path.output,"/Plots/"), sampleNames(gs)[i],".png",sep=""),width=1800,height=1820,pointsize=18)
    par(mfrow=c(5,5))
    par(mar = c(4,5,2,1)+0.1)
    par(oma=c(0,0,3,0))
    
    plot(1, axes=F, xlab="",ylab="",col="white")
    plot(1, axes=F, xlab="",ylab="",col="white")
    #--------- plot margin with cleaned cells
    plotDens(fs.marg[[i]], c("Time","BV605-A"))
    points(fs.marg[[i]]@exprs[clean.inds[[i]]$ind, c("Time", "BV605-A")], pch ='.')
    # ------- plot of cleaned cells with singlets gate
    plotDens(fs.clean[[i]], c("FSC-A","FSC-H"),main="Cleaned")
    lines(singlets[[i]], type="l",lwd=2)
    text(median(singlets[[i]][,1]), median(singlets[[i]][,2]), labels = 'Singlets', cex = 1)
    plot(1, axes=F, xlab="",ylab="",col="white")
    # ------ plot of singlets with size and beads gates
    plotDens(fs.sngl[[i]],c("FSC-A","SSC-A"),main="Singlet")
    lines(beads_poly[[i]]@boundaries,lwd=2)
    lines(size_poly[[i]]@boundaries,lwd=2) # al posto di usare size.filter usiamo size.poly@boundaries che e' la stessa cosa
    text(median(size_poly[[i]]@boundaries[,1]), median(size_poly[[i]]@boundaries[,2]), labels = 'All cells', cex = 1)
    text(median(beads_poly[[i]]@boundaries[,1]), median(beads_poly[[i]]@boundaries[,2]), labels = 'Beads', cex = 1)
    
    # ------ plot size with live gate 
    plotDens(fs.size[[i]],c("APC-eF780-A","SSC-A"),main="Size")
    lines(live_poly[[i]]@boundaries,type="l",lwd=2)
    text(median(live_poly[[i]]@boundaries[,1]), median(live_poly[[i]]@boundaries[,2]), labels = 'Live cells', cex = 1)
    
    #----- plot live with nonGran gate
    plotDens(fs.live[[i]],c(fluorochrome.chans["CD66"],fluorochrome.chans["CD14"]),main="Live")
    lines(gran_poly[[i]]@boundaries)
    lines(non_gran_poly[[i]]@boundaries)
    text(median(non_gran_poly[[i]]@boundaries[,1]), median(non_gran_poly[[i]]@boundaries[,2]), labels = 'Non granulocytes', cex = 1)
    text(median(gran_poly[[i]]@boundaries[,1])+0.5, median(gran_poly[[i]]@boundaries[,2]), labels = 'Granulocytes', cex = 1)
    
    # ------ plot nonGran with Bcell gate
    plotDens(fs.non_gran[[i]],c(fluorochrome.chans["CD19"],scat.chans["SSC-A"]),main="Not Granulocytes")
    lines(Bcells_poly[[i]]@boundaries,type="l",lwd=2)
    text(median(Bcells_poly[[i]]@boundaries[,1]), median(Bcells_poly[[i]]@boundaries[,2]), labels = 'CD19+ B cells', cex = 1)
    
    # ------ plot Cd19+ with PC gate
    plotDens(fs.bcells[[i]],c(fluorochrome.chans["CD38"],fluorochrome.chans["CD138"]),main="B cells")
    lines(pc_poly[[i]]@boundaries,type="l",lwd=2)
    text(median(pc_poly[[i]]@boundaries[,1]), median(pc_poly[[i]]@boundaries[,2]), labels = 'plasma cells', cex = 1)
    
    # ------ plot live with blast gate
    plotDens(fs.live[[i]],c(fluorochrome.chans["CD34"],scat.chans["SSC-A"]),main="Live")
    lines(blast_poly[[i]]@boundaries)
    text(median(blast_poly[[i]]@boundaries[,1]), median(blast_poly[[i]]@boundaries[,2]), labels = 'Blasts', cex = 1)
    
    # --------- plot Cd19+ with gates CD19+CD20+ and CD19+Ce20-
    plot(1, axes=F, xlab="",ylab="",col="white")
    
    plotDens(fs.bcells[[i]],c(fluorochrome.chans["CD19"],fluorochrome.chans["CD20"]),main="B cells")
    lines(CD20pos_poly[[i]]@boundaries,type="l",lwd=2)
    lines(CD20neg_poly[[i]]@boundaries,type="l",lwd=2)
    text(median(CD20pos_poly[[i]]@boundaries[,1]), median(CD20pos_poly[[i]]@boundaries[,2]), labels = 'CD19+CD20+', cex = 1)
    text(median(CD20neg_poly[[i]]@boundaries[,1]), median(CD20neg_poly[[i]]@boundaries[,2]), labels = 'CD19+CD20-', cex = 1)
    
    # ---------- plot CD19+ with quad1
    plotDens(fs.bcells[[i]],channels=c(fluorochrome.chans["IgD"],fluorochrome.chans["IgM"]),main="B cells")#,nbin=3000
    abline(v=quad1[[i]]@gates[1],h=quad1[[i]]@gates[2],lwd=2)
    text(quad1[[i]]@gates[1]*0.60, quad1[[i]]@gates[2]*1.40, labels = 'IGD-IGM+ B cells', cex = 0.7)
    text(quad1[[i]]@gates[1]*0.60, quad1[[i]]@gates[2]*0.60, labels = 'IGD-IGM- B cells', cex = 0.7)
    text(quad1[[i]]@gates[1]*1.40, quad1[[i]]@gates[2]*0.60, labels = 'IGD+IGM- B cells', cex = 0.7)
    text(quad1[[i]]@gates[1]*1.40, quad1[[i]]@gates[2]*1.40, labels = 'IGD+IGM+ B cells', cex = 0.7)
    
    #------------ plot CD19+ with gate PB
    plotDens(fs.bcells[[i]],c(fluorochrome.chans["CD38"],fluorochrome.chans["CD27"]),main="B cells")  
    lines(PB_poly[[i]]@boundaries,type="l",lwd=2)
    text(median(PB_poly[[i]]@boundaries[,1]), median(PB_poly[[i]]@boundaries[,2]), labels = 'Plasmablasts', cex = 1)
    
    #----- plot CD20+ with gate CD10+CD27-,CD20+Cd27+
    plot(1, axes=F, xlab="",ylab="",col="white")
    plot(1, axes=F, xlab="",ylab="",col="white")
    
    plotDens(fs.CD20pos[[i]],c(fluorochrome.chans["CD10"],fluorochrome.chans["CD27"]),main="CD20+")
    lines(CD10neg_poly[[i]]@boundaries,type="l",lwd=2)
    lines(CD10pos_poly[[i]]@boundaries,type="l",lwd=2)
    text(median(CD10neg_poly[[i]]@boundaries[,1]), median(CD10neg_poly[[i]]@boundaries[,2]), labels = 'CD10-', cex = 1)
    text(median(CD10pos_poly[[i]]@boundaries[,1]), median(CD10pos_poly[[i]]@boundaries[,2]), labels = 'CD10+', cex = 1)
    
    # ------ plot transitional pop with trans gate
    plotDens(fs.igm[[i]], c(fluorochrome.chans["CD10"],fluorochrome.chans["CD38"]),main="IgM+IgD+CD27-")
    lines(trans_poly[[i]]@boundaries,type="l",lwd=2)
    text(median(trans_poly[[i]]@boundaries[,1]), median(trans_poly[[i]]@boundaries[,2]), labels = 'Immature transition B cells', cex = 1)
    
    # ------ plot CD10- with gate naive pop
    plot(1, axes=F, xlab="",ylab="",col="white")
    plot(1, axes=F, xlab="",ylab="",col="white")
    plot(1, axes=F, xlab="",ylab="",col="white")
    
    plotDens(fs.cd10neg[[i]],channels=c(fluorochrome.chans["CD27"],fluorochrome.chans["IgD"]),main="CD10-")
    abline(v=quad.3[[i]]@gates[1],h=quad.3[[i]]@gates[2],lwd=2)
    text(quad.3[[i]]@gates[1]*0.80, quad.3[[i]]@gates[2]*1.20, labels = 'Naive B cells', cex = 0.6)
    text(quad.3[[i]]@gates[1]*0.80, quad.3[[i]]@gates[2]*0.80, labels = 'Atypical B cells', cex = 0.6)
    text(quad.3[[i]]@gates[1]*1.20, quad.3[[i]]@gates[2]*0.80, labels = 'Unswitched memory B cells', cex = 0.6)
    text(quad.3[[i]]@gates[1]*1.20, quad.3[[i]]@gates[2]*1.20, labels = 'Switched memory B cells', cex = 0.6)
    
    
    dev.off()
    return(i)
  },mc.cores=n_cores)
  
}

###########################################################################################
#################### Functions related to Plots manual vs automated pop ###################
###########################################################################################

#--------  Function to import and set up the dataset of manual counts of B cells
pre_prop_manual_data <- function(){
  manual_data_xls <- read_excel("/home/rstudio/data/Gambia validation cohort (Aug2018)/Bcell_Gambia_confcohort/B cells Manual counts.xls")
  write.csv(manual_data_xls,"/home/rstudio/manual_data.csv")
  manual_data<- read.csv(file="/home/rstudio/manual_data.csv",header=TRUE,sep=",")
  #------- riassegnamo nomi colonne del dataset manuale
  pop_names_manual <- colnames(manual_data)
  pop_names_manual[which(pop_names_manual=="..1")]<- "Sample"
  pop_names_manual[which(pop_names_manual=="Time.exclusion")]<- "Clean_0"
  pop_names_manual[which(pop_names_manual=="All.cells")]<- "Size"
  pop_names_manual[which(pop_names_manual=="Live.cells")] <- "Live"
  pop_names_manual[which(pop_names_manual=="Blasts")] <- "Blast"
  pop_names_manual[which(pop_names_manual=="Non.Granulocytes")] <- "NonGran"
  pop_names_manual[which(pop_names_manual=="CD19..B.cells")] <- "Bcells" # CD19+ Bcells
  pop_names_manual[which(pop_names_manual=="CD19.CD20.")] <- "CD19.20pos" # CD19+CD20+
  pop_names_manual[which(pop_names_manual=="CD10..CD27.")] <- "CD10neg.27" # CD10- CD27-
  pop_names_manual[which(pop_names_manual=="CD27.CD10.")] <- "CD10pos.27" # CD27-CD10+
  pop_names_manual[which(pop_names_manual=="Plasma.cells")] <- "PC"
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells")] <- "IGquad4" # IgD+ IgM+ B cells
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells.1")] <- "IGquad2"# IgD+ IgM- B cells
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells.2")] <- "IGquad3" # IgD- IgM+ B cells
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells.3")] <- "IGquad1" # IgD- IgM- B cells
  pop_names_manual[which(pop_names_manual=="Plasmablasts")] <- "PB"
  pop_names_manual[which(pop_names_manual=="Plasma.cells")] <- "PC"
  pop_names_manual[which(pop_names_manual=="CD19.CD20..1")] <- "CD19.20neg"
  colnames(manual_data) <- pop_names_manual
  #------ adesso convertiamo le conte in percentuali delle colonne delle popolazioni
  # manual_data_perc <- manual_data[,2:length(colnames(manual_data))]/manual_data[,"Live cells"] # da sei se vogliamo includere le pop solo a partire dalle live,da 1 se tutte
  #--- creiamo dataframe che associa a ogni nome di una pop la sua parente
  list_child_vs_parent<-lapply(1:length(getNodes(gs)), function(i){
    splitted_string <-str_split(getNodes(gs)[i],"/")
    child_pop <- splitted_string[[1]][length(splitted_string[[1]])]
    parent_pop <- splitted_string[[1]][length(splitted_string[[1]])-1]
    return(list(child_pop=child_pop,parent_pop=parent_pop))
  })
  
  list_all_pop_freq <- lapply(1:length(list_child_vs_parent), function(i){
    child_pop_x <- list_child_vs_parent[[i]]$child_pop
    
    if(child_pop_x %in% colnames(manual_data)){
      parent_pop_x <- list_child_vs_parent[[i]]$parent_pop
      if(parent_pop_x %in% colnames(manual_data)){
        pop_x_freq <- manual_data[,child_pop_x]/manual_data[,parent_pop_x]
        
        return(list(pop_x_freq=pop_x_freq,name=child_pop_x))
      }
    }
  })
  indx<-which(list_all_pop_freq=="NULL")
  list_all_pop_freq<-list_all_pop_freq[-indx]
  
  names<-sapply(1:length(list_all_pop_freq),function(i){
    return(list_all_pop_freq[[i]]$name)
  })
  list_all_pop_freq_only<-lapply(1:length(list_all_pop_freq),function(i){
    return(list(pop=list_all_pop_freq[[i]]$pop_x_freq))
  })
  df <- as.data.frame(list_all_pop_freq_only)
  colnames(df)<-names
  manual_data_freq_on_parents <- df
  samples<-manual_data$Sample
  manual_data_freq_on_parents<-cbind(samples,manual_data_freq_on_parents)
  return(manual_data_freq_on_parents)
}

#------ function to make Correlation plot auto vs manual pop
make_corr_plot<-function(color_pop=NULL,remove_p=NULL){
  if(is.null(remove_p)==FALSE){
    indx<-grep(remove_p,final_df$Population)
    pops<-strsplit(remove_p,"|",fixed=TRUE)[[1]]
    if("CD56 dim CD16+ NK" %in% pops){
      ind<-which(final_df$Population=="CD56 dim CD16+ NK")
      indx<-c(indx,ind)
    }
    if("HLADR+ CD14+" %in% pops){
      ind<-which(final_df$Population=="HLADR+ CD14+")
      indx<-c(indx,ind)
    }
    final_df<-final_df[-indx,]
  }
  #------ Correlation spearman test
  results <- cor.test(final_df$`Counts/Parent_counts`,final_df$freq_manual_gating)
  r_value <- round(results$estimate,3)
  #----- grafich freq_man vs freq_auto
  # Correlation plot with all pops colored
  corr_plot <- ggplot(final_df,aes(x=final_df$`Counts/Parent_counts`,y=final_df$freq_manual_gating)) + geom_point(aes(colour=final_df$Population))  +
    annotate("text",x=0.2,y=1,label=sprintf("rs: %s",r_value),size=5) + labs(x="Automated_frequencies (% Parents)",y="Manual_frequencies (% Parents)") + theme(legend.title=element_blank())
  # Version in which I color only the population of interest
  # select events to be coloured in v
  if(is.null(color_pop)==FALSE){
    v<- final_df$Population==color_pop
    if(color_pop %in% unique(final_df[,"Population"])==FALSE){
      stop("the population selected is not present,please select a correct population name")
    } 
    corr_plot <- ggplot(final_df,aes(x=final_df$`Counts/Parent_counts`,y=final_df$freq_manual_gating)) + geom_point(aes(colour=v)) + scale_color_manual(values = c("black","red"),labels=c("Other pops",color_pop))+  
      annotate("text",x=0.2,y=1,label=sprintf("rs: %s",r_value),size=5) + labs(x="Automated_frequencies (% Parents)",y="Manual_frequencies (% Parents)") + theme(legend.title=element_blank())
  }
  return(corr_plot)
}


make_corr_plot_counts<-function(color_pop=NULL,remove_p=NULL){
  if(is.null(remove_p)==FALSE){
    indx<-grep(remove_p,final_df_counts_only$Population)
    pops<-strsplit(remove_p,"|",fixed=TRUE)[[1]]
    if("CD56 dim CD16+ NK" %in% pops){
      ind<-which(final_df_counts_only$Population=="CD56 dim CD16+ NK")
      indx<-c(indx,ind)
    }
    if("HLADR+ CD14+" %in% pops){
      ind<-which(final_df_counts_only$Population=="HLADR+ CD14+")
      indx<-c(indx,ind)
    }
    final_df_counts_only<-final_df_counts_only[-indx,]
  }
  #------ spearman
  results <- cor.test(final_df_counts_only$Count,final_df_counts_only$counts_manual_gating)
  r_value <- round(results$estimate,3)
  #----- plots count manual vs counts auto
  corr_plot <- ggplot(final_df_counts_only,aes(x=final_df_counts_only$Count,y=final_df_counts_only$counts_manual_gating)) + geom_point(aes(colour=final_df_counts_only$Population))  +
    labs(x="Automated_counts",y="Manual_counts") + theme(legend.title=element_blank())
  if(is.null(color_pop)==FALSE){
    v<- final_df_counts_only$Population==color_pop
    if(color_pop %in% unique(final_df_counts_only[,"Population"])==FALSE){
      stop("the population selected is not present,please select a correct population name")
    } 
    corr_plot <- ggplot(final_df_counts_only,aes(x=final_df_counts_only$Count,y=final_df_counts_only$counts_manual_gating)) + geom_point(aes(colour=v)) + scale_color_manual(values = c("black","red"),labels=c("Other pops",color_pop))+  
      labs(x="Automated_counts",y="Manual_counts") + theme(legend.title=element_blank())
  }
  return(corr_plot)
}

#------------ Box plot ------------------
generate_box_plots <- function(remove_p=NULL){

  #----- non coupled boxplot
  ordered_pop<-factor(final_df$Population,levels = unique(final_df$Population))
  
  boxplot_automated_population <- ggplot(final_df,aes(x=ordered_pop,y=final_df$`Counts/Parent_counts`)) + geom_boxplot(aes(colour=ordered_pop)) + 
    theme(legend.title=element_blank()) + labs(x="Populations",y="Count/parents") + 
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
  
  boxplot_manual_population <- ggplot(final_df, aes(x=ordered_pop, y = final_df$freq_manual_gating)) + geom_boxplot(aes(colour=ordered_pop)) + 
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
    theme(legend.title=element_blank()) + labs(x="Populations",y="Count/parents")
  #-------- face wrap version of the coupled box plot
  # melt reorganize the dataframe based on id.vars names
  df_melt<-melt(final_df,id.vars=c("name","Population","Parent","Count","ParentCount"))
  # change level name using revalue of plyr
  df_melt$variable <- revalue(df_melt$variable,c("freq_manual_gating"="Manual_Pop", "Counts/Parent_counts"="Automated_Pop"))
  # scale manual color regolulates the contour and it does not work with boxplot,you must use scale_fill_manual if you want to change the color of the boxplot.
  if(is.null(remove_p)==FALSE){
    indx<-grep(remove_p,df_melt$Population)
    df_melt<-df_melt[-indx,]
  }
  boxplots_coupled_plot <-ggplot(df_melt, aes(x=factor(variable),y=value,fill=factor(variable))) + geom_boxplot() +   
    facet_wrap(~Population) + labs(title="automated_vs_manual_boxplot") + scale_fill_manual(values = c("red","blue")) +
    theme(legend.title=element_blank()) + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) 
  #------ non face wrap version of the coupled boxplot
  boxplots_coupled_plot_JIJ_paper<- ggplot(data = df_melt, aes(x=Population, y=value)) + geom_boxplot(aes(fill=variable)) + 
    theme(legend.title=element_blank()) + labs(x="Population_name") + labs(y="frequencies_values")
  return(boxplots_coupled_plot_JIJ_paper)
}

# -------make dataframe combined between manual and automated freq
make_final_dataframe<-function(){
  pop_names_tab <- unique(tab_stat_counts_freq[,"Population"])
  
  list_df_pop_x <- list()
  for(name in pop_names_tab){
    if(name %in% colnames(manual_data_freq_on_parents)){
      tab_stat_counts_freq_pop_x<-tab_stat_counts_freq[tab_stat_counts_freq["Population"]==name,]
      row.names(tab_stat_counts_freq_pop_x)<-NULL
      manual_data_freq_on_parents_pop_x<-manual_data_freq_on_parents[c(name,"samples")]
      samples_names<-tab_stat_counts_freq_pop_x$name
      s<-as.character(manual_data_freq_on_parents_pop_x[["samples"]]) %in% samples_names
      
      manual_data_freq_on_parents_pop_x<-manual_data_freq_on_parents_pop_x[s,]
      colnames(manual_data_freq_on_parents_pop_x)[1]<-"freq_manual_gating"
      ind<-order(manual_data_freq_on_parents_pop_x$samples)
      # order the manual dataframe 
      v1<-tab_stat_counts_freq_pop_x$name
      v2 <- manual_data_freq_on_parents_pop_x$samples
      ind_v<-c()
      for (n in v1){
        ind<-which(v2==n)
        ind_v<-c(ind_v,ind)
      }
      manual_data_freq_on_parents_pop_x<-manual_data_freq_on_parents_pop_x[ind_v,]
      #--- fuse the two dataframes manual and automated
      df_automated_plus_freq_manual_pop_x <- cbind(tab_stat_counts_freq_pop_x,manual_data_freq_on_parents_pop_x)
      df_automated_plus_freq_manual_pop_x<-df_automated_plus_freq_manual_pop_x[,1:length(colnames(df_automated_plus_freq_manual_pop_x))-1]
      list_df_pop_x[[name]] <- df_automated_plus_freq_manual_pop_x
    }
  }
  final_df <- rbindlist(list_df_pop_x)
  final_df<- as.data.frame(final_df)
  return(final_df)
}


#---- function to generate stats table of counts with frequencies on the parents pops
pre_prop_automated_data <- function(i){
  tab_stat_counts <-getPopStats(gs)
  tab_stat_counts<-as.data.frame(tab_stat_counts)
  samples_names<-unique(tab_stat_counts$name)
  list_dataframes_sample<-list()
  for(name in samples_names){
    tab_stat_counts_sample <- tab_stat_counts[tab_stat_counts["name"]==name,]
    list_dataframes_sample[[name]] <- tab_stat_counts_sample
  }
  tab_stat_counts_freq <- as.data.frame(rbindlist(list_dataframes_sample))
  
  tab_stat_counts_freq[,"Counts/Parent_counts"] <- round(tab_stat_counts_freq["Count"]/tab_stat_counts_freq["ParentCount"],3)
  return(tab_stat_counts_freq)
}

#---- function to import only the manual counts of Bcells gambia(not the freq)

pre_prop_manual_counts_Bcells_Gambia <- function(){
  manual_data_xls <- read_excel("/home/rstudio/data/Gambia validation cohort (Aug2018)/Bcell_Gambia_confcohort/B cells Manual counts.xls")
  write.csv(manual_data_xls,"/home/rstudio/manual_data.csv")
  manual_data<- read.csv(file="/home/rstudio/manual_data.csv",header=TRUE,sep=",")
  #------- reassign column names of manual dataset
  pop_names_manual <- colnames(manual_data)
  pop_names_manual[which(pop_names_manual=="..1")]<- "Sample"
  pop_names_manual[which(pop_names_manual=="Time.exclusion")]<- "Clean_0"
  pop_names_manual[which(pop_names_manual=="All.cells")]<- "Size"
  pop_names_manual[which(pop_names_manual=="Live.cells")] <- "Live"
  pop_names_manual[which(pop_names_manual=="Blasts")] <- "Blast"
  pop_names_manual[which(pop_names_manual=="Non.Granulocytes")] <- "NonGran"
  pop_names_manual[which(pop_names_manual=="CD19..B.cells")] <- "Bcells" # CD19+ Bcells
  pop_names_manual[which(pop_names_manual=="CD19.CD20.")] <- "CD19.20pos" # CD19+CD20+
  pop_names_manual[which(pop_names_manual=="CD10..CD27.")] <- "CD10neg.27" # CD10- CD27-
  pop_names_manual[which(pop_names_manual=="CD27.CD10.")] <- "CD10pos.27" # CD27-CD10+
  pop_names_manual[which(pop_names_manual=="Plasma.cells")] <- "PC"
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells")] <- "IGquad4" # IgD+ IgM+ B cells
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells.1")] <- "IGquad2"# IgD+ IgM- B cells
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells.2")] <- "IGquad3" # IgD- IgM+ B cells
  pop_names_manual[which(pop_names_manual=="IgD..IgM..B.cells.3")] <- "IGquad1" # IgD- IgM- B cells
  pop_names_manual[which(pop_names_manual=="Plasmablasts")] <- "PB"
  pop_names_manual[which(pop_names_manual=="Plasma.cells")] <- "PC"
  pop_names_manual[which(pop_names_manual=="CD19.CD20..1")] <- "CD19.20neg"
  colnames(manual_data) <- pop_names_manual


  manual_data<-manual_data[ind_v,]
  manual_data<-manual_data[,2:length(colnames(manual_data))]
  return(manual_data)
}

pre_prop_automated_data_counts <- function(i){
  tab_stat_counts <-getPopStats(gs)
  tab_stat_counts<-as.data.frame(tab_stat_counts)

  return(tab_stat_counts)
}

make_final_dataframe_counts_only<-function(){
  pop_names_tab <- unique(tab_stat_counts[,"Population"])
  
  list_df_pop_x <- list()
  for(name in pop_names_tab){
    if(name %in% colnames(manual_data)){
      tab_stat_counts_pop_x<-tab_stat_counts[tab_stat_counts["Population"]==name,]
      row.names(tab_stat_counts_pop_x)<-NULL
      manual_data_pop_x<-manual_data[c(name,"Sample")]
      samples_names<-tab_stat_counts_pop_x$name
      s<-as.character(manual_data_pop_x[["Sample"]]) %in% samples_names
      
      manual_data_pop_x<-manual_data_pop_x[s,]
      colnames(manual_data_pop_x)[1]<-"counts_manual_gating"
      ind<-order(manual_data_pop_x$Sample)
      # order the manual dataframe 
      v1<-tab_stat_counts_pop_x$name
      v2 <- manual_data_pop_x$Sample
      ind_v<-c()
      for (n in v1){
        ind<-which(v2==n)
        ind_v<-c(ind_v,ind)
      }
      manual_data_pop_x<-manual_data_pop_x[ind_v,]
      #--- fuse the two dataframes manual and automated
      df_automated_plus_manual_pop_x <- cbind(tab_stat_counts_pop_x,manual_data_pop_x)
      df_automated_plus_manual_pop_x<-df_automated_plus_manual_pop_x[,1:length(colnames(df_automated_plus_manual_pop_x))-1]
      list_df_pop_x[[name]] <- df_automated_plus_manual_pop_x
    }
  }
  final_df_counts_only <- rbindlist(list_df_pop_x)
  final_df_counts_only<- as.data.frame(final_df_counts_only)
  return(final_df_counts_only)
}
################################################# Additional functions #################3
#------------------------ check fs counts --------------

# function to check the number of cells in a population of the flowSet for each sample
# the populations with a number of cells below the defined threshold should not be analyzed.
# Therefore this function return the indices of the sample to be removed from the gs
check_fs_counts<-function(fs,threshold=20000,counts=F){
  vec_cells_count_all_samples<-sapply(1:length(fs),function(i){
    cells_count_sample_x<-nrow(fs[[i]]@exprs)
    return(cells_count_sample_x)
  })
  if(counts==F){
    inds<-which(vec_cells_count_all_samples<threshold)
    return(inds)
  }
  if(counts==T){
    return(vec_cells_count_all_samples)
  }  
  
}

#------------------------ check final counts and frequencies ---------------

# Fnuction that scan the .csv files generated by the algorithm to generate the final unique files of the frequencies and counts

check_final_results<-function(files_path,check_Rym_data=F,path_Rym_data=NULL,control="H12",remove_samples=NULL){
  #---------- pre-processsing -----------------
  CSVfiles_path_list <- list.files(path = files_path, recursive = F, pattern = ".csv", full.names = T)
  print(CSVfiles_path_list)
  list_splitted<-strsplit(CSVfiles_path_list,"/")
  vector_names_all_csv_files<-sapply(1:length(list_splitted), function(i){
    name_csv_file<-tail(list_splitted[[i]],1)
  })
  
  # to be sure that vector_names_csv_file contains only csv files
  inds<-grep(".csv",vector_names_all_csv_files)
  vector_names_all_csv_files<-vector_names_all_csv_files[inds]
  #------------------ select names of main counts and frequecies file
  names_main_files<-vector_names_all_csv_files[grep("All|-|:|comp",vector_names_all_csv_files)]
  if(length(names_main_files)==0){ # it means there is no distiction between main file and separate files, so the pops files become the main files.
    names_main_files<-vector_names_all_csv_files[grep("pops_",vector_names_all_csv_files)]
  }
  # be sure that that this vector does not contain the name of previous final counts or freqs file
  ind_final<-grep("final",names_main_files)
  if(length(ind_final)>0){
    names_main_files<-names_main_files[-ind_final]
  }
  name_main_counts_file<-names_main_files[grep("counts",names_main_files)]
  name_main_freqs_file<-names_main_files[grep("freqs",names_main_files)]
  #--------- select file names that refer to samples analyzed separately from the main group
  string<-paste(names_main_files,collapse = "|") # in case there are more than one file
  inds<-grep(string,vector_names_all_csv_files)
  vector_names_no_main_files<-vector_names_all_csv_files[-inds]
  vector_names_no_main_files_counts<-vector_names_no_main_files[grep("counts",vector_names_no_main_files)]
  vector_names_no_main_files_freqs<-vector_names_no_main_files[grep("freqs",vector_names_no_main_files)]
  # be sure that that these vectors does not contain the name of previous final counts or freqs file
  ind_final<-grep("final",vector_names_no_main_files_counts)
  if(length(ind_final)>0){
    vector_names_no_main_files_counts<-vector_names_no_main_files_counts[-ind_final]
  }
  ind_final<-grep("final",vector_names_no_main_files_freqs)
  if(length(ind_final)>0){
    vector_names_no_main_files_freqs<-vector_names_no_main_files_freqs[-ind_final]
  }
  # be sure that that these vectors does not contain the name of previous sample control counts or freqs files
  ind_control_samples<-grep("sample_control",vector_names_no_main_files_counts)
  if(length(ind_control_samples)>0){
    vector_names_no_main_files_counts<-vector_names_no_main_files_counts[-ind_control_samples]
  }
  ind_control_samples<-grep("sample_control",vector_names_no_main_files_freqs)
  if(length(ind_control_samples)>0){
    vector_names_no_main_files_freqs<-vector_names_no_main_files_freqs[-ind_control_samples]
  }
  # paste information
  print(paste0("name_main_counts_file:",name_main_counts_file))
  print(paste0("name_main_freqs_file:",name_main_freqs_file))
  total_samples_analyzed_separately_counts<-length(vector_names_no_main_files_counts)
  total_samples_analyzed_separately_freqs<-length(vector_names_no_main_files_freqs)
  print(paste0("total_samples_analyzed_separately_counts:",total_samples_analyzed_separately_counts))
  print(paste0("total_samples_analyzed_separately_freqs:",total_samples_analyzed_separately_freqs))
  if(total_samples_analyzed_separately_counts!=total_samples_analyzed_separately_freqs){
    stop("something goes wrong: total_samples_analyzed_separately_counts!=total_samples_analyzed_separately_freqs is True")
  }

  #------------------------ counts file generation ---------------------
  # read main files counts
  string<-paste(name_main_counts_file,collapse = "|")
  inds<-grep(string,CSVfiles_path_list)
  path_main_counts_files<-CSVfiles_path_list[inds]
  list_main_counts_files<-lapply(1:length(path_main_counts_files), function(i){
    file_main_counts_file_x<-read.csv(path_main_counts_files[[i]])
    return(file_main_counts_file_x)
  })
  if(length(list_main_counts_files)>1){
    single_main_counts_file<-as.data.frame(rbindlist(list_main_counts_files))
  }else{
    single_main_counts_file<-list_main_counts_files[[1]]
  }
  unique_samples_name<-unique(single_main_counts_file$name)
  total_samples_analyzed_main_files_counts<-length(unique_samples_name)
  print(paste0("total_samples_analyzed_main_files_counts:",total_samples_analyzed_main_files_counts))
  # remove the names of the sample analyzed separately from the main counts file
  all_names_df_vector<-single_main_counts_file$name
  if(total_samples_analyzed_separately_counts>0){
    list_splitted<-strsplit(vector_names_no_main_files_counts,"_")
    vec_temp<-sapply(1:length(list_splitted),function(i){
      splitted_part<-tail(list_splitted[[i]],1)
      return(splitted_part)
    })
    list_splitted<-strsplit(vec_temp,".csv")
    vec_samples_name<-sapply(1:length(list_splitted),function(i){
      sample_name<-tail(list_splitted[[i]],1)
      sample_name<-paste0(sample_name,"_")
      return(sample_name)
    })
    string<-paste(vec_samples_name,collapse = "|")
    print(string)
    inds<-grep(string,all_names_df_vector)
    single_main_counts_file<-single_main_counts_file[-inds,]
  
  unique_samples_name<-unique(single_main_counts_file$name)
  samples_analyzed_main_files_counts<-length(unique_samples_name)
  print(paste0("samples_analyzed_main_files_counts:",samples_analyzed_main_files_counts))
  tot<-samples_analyzed_main_files_counts+total_samples_analyzed_separately_counts
  # find the samples analyzed seperately and not present in main files 
  main_files_samples_names_unique<-unique(all_names_df_vector)
  v_temp<-sapply(1:length(vec_samples_name), function(i){
    samples_name_x<-vec_samples_name[[i]]
    ind<-grep(samples_name_x,main_files_samples_names_unique)
    if(length(ind)!=0){
      return(1)
    }else{
      return(0)
    }
  })
  inds<-which(v_temp==0)
  name_samples_analyzed_separately_not_in_main_files<-vec_samples_name[inds]
  number_samples_analyzed_separately_not_in_main_files<-length(inds)
  print(paste0("number_samples_analyzed_separately_not_in_main_files:",number_samples_analyzed_separately_not_in_main_files))
  # check correctness of operations
  print(name_main_counts_file)
  if(tot!=(total_samples_analyzed_main_files_counts+number_samples_analyzed_separately_not_in_main_files)){
    stop("something goes wrong:tot!=(total_samples_analyzed_main_files_counts+number_samples_analyzed_separately_not_in_main_files) is True")
  }
  # read no main files counts
  string<-paste(vector_names_no_main_files_counts,collapse = "|")
  inds<-grep(string,CSVfiles_path_list)
  path_no_main_counts_files<-CSVfiles_path_list[inds]
  list_no_main_counts_files<-lapply(1:length(path_no_main_counts_files), function(i){
    file_no_main_counts_file_x<-read.csv(path_no_main_counts_files[[i]])
    return(file_no_main_counts_file_x)
  })
  single_no_main_counts_file<-as.data.frame(rbindlist(list_no_main_counts_files))
  n1<-nrow(single_no_main_counts_file)
  n2<-nrow(single_main_counts_file)
  n_tot<-n1+n2
  final_counts_df<-rbind(single_main_counts_file,single_no_main_counts_file)
  n_final<-nrow(final_counts_df)
  if(n_final!=n_tot){
    stop("error during final counts file generation: n_final!=n_tot is True")
  }
  }else{
    final_counts_df<-single_main_counts_file 
  }
  # order dataframe 
  final_counts_df<-final_counts_df[order(final_counts_df$name),]
  # remove X column
  v<-colnames(final_counts_df)
  ind<-grep("X",v)
  v_less_x_column<-v[-ind]
  final_counts_df<-final_counts_df[,v_less_x_column]

  # remove specific populations counts
  final_counts_df$Population<-as.character(final_counts_df$Population)
  final_counts_df$Parent<-as.character(final_counts_df$Parent)
  inds<-grep("Margin|Clean_1",final_counts_df$Population)
  final_counts_df<-final_counts_df[-inds,]
  # rename clean_0 with Time
  inds<-grep("Clean_0",final_counts_df$Population)
  final_counts_df$Population[inds]<-"Time"
  inds<-grep("Clean_0",final_counts_df$Parent)
  final_counts_df$Parent[inds]<-"Time"
  # separate positive control samples counts or not
  if(control=="None"){
    write.csv(final_counts_df,file = sprintf("%s/final_counts_file.csv",files_path),row.names = F)
  }else if(control=="compensation"){
    write.csv(final_counts_df,file = sprintf("%s/sample_control_counts.csv",files_path),row.names = F)
  }else{
    string_sample_x<-sprintf("%s$|%s_|%s\\.",control,control,control)
    ind<-grep(string_sample_x,final_counts_df$name)
    sample_control_counts<-final_counts_df[ind,]
    # export the final counts file
    write.csv(final_counts_df,file = sprintf("%s/final_counts_file.csv",files_path),row.names = F)
    write.csv(sample_control_counts,file = sprintf("%s/sample_control_counts.csv",files_path),row.names = F)
  }
  


  #------------------------ freqs file generation ---------------------
  # read main files freqs
  string<-paste(name_main_freqs_file,collapse = "|")
  inds<-grep(string,CSVfiles_path_list)
  path_main_freqs_files<-CSVfiles_path_list[inds]
  list_main_freqs_files<-lapply(1:length(path_main_freqs_files), function(i){
    file_main_freqs_file_x<-read.csv(path_main_freqs_files[[i]])
    column_names<-as.character(file_main_freqs_file_x$X)
    transposed_df<-as.data.frame(t(file_main_freqs_file_x))
    colnames(transposed_df)<-column_names
    transposed_df <- cbind(rownames(transposed_df), data.frame(transposed_df, row.names=NULL))
    colnames(transposed_df)[2:length(colnames(transposed_df))]<-column_names
    colnames(transposed_df)[1]<-"samples_name"
    return(transposed_df[-1,])
  })
  if(length(list_main_freqs_files)>1){
    single_main_freqs_file<-as.data.frame(rbindlist(list_main_freqs_files))
  }else{
    single_main_freqs_file<-list_main_freqs_files[[1]]
  }
  #print(single_main_freqs_file)
  total_samples_analyzed_main_files_freqs<-nrow(single_main_freqs_file)
  print(paste0("total_samples_analyzed_main_files_freqs:",total_samples_analyzed_main_files_freqs))
  # check correctness of operations
  if(total_samples_analyzed_main_files_freqs!=total_samples_analyzed_main_files_counts){
    stop("something goes wrong:total_samples_analyzed_main_files_freqs!=total_samples_analyzed_main_files_counts is True")
  }
  # remove the names of the sample analyzed separately from the main freqs file
  samples_name<-single_main_freqs_file$samples_name
  if(total_samples_analyzed_separately_counts>0){
  list_splitted<-strsplit(vector_names_no_main_files_freqs,"_")
  vec_temp<-sapply(1:length(list_splitted),function(i){
    splitted_part<-tail(list_splitted[[i]],1)
    return(splitted_part)
  })
  list_splitted<-strsplit(vec_temp,".csv")
  vec_samples_name<-sapply(1:length(list_splitted),function(i){
    sample_name<-tail(list_splitted[[i]],1)
    sample_name<-paste0(sample_name,"_")
    return(sample_name)
  })
  string<-paste(vec_samples_name,collapse = "|")
  inds<-grep(string,samples_name)
  single_main_freqs_file<-single_main_freqs_file[-inds,]
  samples_analyzed_main_files_freqs<-nrow(single_main_freqs_file)
  print(paste0("samples_analyzed_main_files_freqs:",samples_analyzed_main_files_freqs))
  tot<-samples_analyzed_main_files_freqs+total_samples_analyzed_separately_freqs
  # check correctness of operations
  if(tot!=(total_samples_analyzed_main_files_freqs+number_samples_analyzed_separately_not_in_main_files)){
    stop("something goes wrong:tot!=(total_samples_analyzed_main_files_freqs+number_samples_analyzed_separately_not_in_main_files) is True")
  }
  # read no main files freqs
  string<-paste(vector_names_no_main_files_freqs,collapse = "|")
  inds<-grep(string,CSVfiles_path_list)
  path_no_main_freqs_files<-CSVfiles_path_list[inds]
  list_no_main_freqs_files<-lapply(1:length(path_no_main_freqs_files), function(i){
    file_no_main_freqs_file_x<-read.csv(path_no_main_freqs_files[[i]])
    column_names<-as.character(file_no_main_freqs_file_x$X)
    transposed_df<-as.data.frame(t(file_no_main_freqs_file_x))
    colnames(transposed_df)<-column_names
    transposed_df <- cbind(rownames(transposed_df), data.frame(transposed_df, row.names=NULL))
    colnames(transposed_df)[2:length(colnames(transposed_df))]<-column_names
    colnames(transposed_df)[1]<-"samples_name"
    return(transposed_df[-1,])
  })
  single_no_main_freqs_file<-as.data.frame(rbindlist(list_no_main_freqs_files))
  n1<-nrow(single_no_main_freqs_file)
  n2<-nrow(single_main_freqs_file)
  n_tot<-n1+n2
  final_freqs_df<-rbind(single_main_freqs_file,single_no_main_freqs_file)
  n_final<-nrow(final_freqs_df)
  if(n_final!=n_tot){
    stop("error during final freqs file generation: n_final!=n_tot is True")
  }}else{
    final_freqs_df<-single_main_freqs_file
  }
  # order dataframe 
  final_freqs_df<-final_freqs_df[order(final_freqs_df$samples_name),]
  # remove specific populations frequencies
  inds<-grep("Margin|Clean_1|Beads",colnames(final_freqs_df))
  final_freqs_df<-final_freqs_df[,-inds]
  # rename clean_0 with Time
  inds<-grep("Clean_0",colnames(final_freqs_df))
  colnames(final_freqs_df)[inds]<-"Time"
  # format freqs df correclty
  samples_name<-final_freqs_df[,1]
  final_freqs_df<-final_freqs_df[,-1]
  final_freqs_df[]<-lapply(final_freqs_df,as.character)
  final_freqs_df[]<-lapply(final_freqs_df,as.numeric)
  final_freqs_df<-round(final_freqs_df,2)
  final_freqs_df<-cbind(samples_name,final_freqs_df)
  # separate positive control samples frequencies
  if(control=="None"){
    write.csv(final_freqs_df,file = sprintf("%s/final_freqs_file.csv",files_path),row.names = F)
  }else if(control=="compensation"){
    write.csv(final_freqs_df,file = sprintf("%s/sample_control_freqs.csv",files_path),row.names = F)
  }else{
    string_sample_x<-sprintf("%s$|%s_|%s\\.",control,control,control)
    ind<-grep(string_sample_x,final_freqs_df$samples_name)
    sample_control_freqs<-final_freqs_df[ind,]
    # export the final freqs file
    write.csv(final_freqs_df,file = sprintf("%s/final_freqs_file.csv",files_path),row.names = F)
    write.csv(sample_control_freqs,file = sprintf("%s/sample_control_freqs.csv",files_path),row.names = F)
  }

  # ------------------- check that both freqs and counts files are correct ----------------------
  if(check_Rym_data){
    # They shoud report the results of the same number of samples reported in Rym excel file except the files flagged by Rym at the beginning(first criteria)
    Rym_csv_file_counts<-read.csv(path_Rym_data)
    samples_name_Rym<-Rym_csv_file_counts[,1]
    samples_name_analysis<-as.character(final_freqs_df$samples_name)
    len_1<-length(samples_name_Rym)
    len_2<-length(samples_name_analysis)
    print(paste0("length_samples_name_Rym:",len_1))
    print(paste0("length_samples_name_analysis:",len_2))
    temp_split<-strsplit(samples_name_analysis,"_")
    vec_s<-sapply(1:length(temp_split), function(i){
        string_sample_x<-temp_split[[i]][3]
        string_sample_x<-sprintf("%s$|%s_|%s\\.",string_sample_x,string_sample_x,string_sample_x)
        return(string_sample_x)
      })
    string<-paste(vec_s,collapse = "|")
    samples_not_analyzed<-grep(string,samples_name_Rym,invert = T,value = T) # strings samples not analyzed
    if(length(samples_not_analyzed)==0){
      samples_not_analyzed<-"Null"
    }
    print(paste0("samples_not_analyzed:",samples_not_analyzed))
  }
  #------------------- remove samples from the counts and frequencies files -----------
  # Remove the samples flagged after evaluation of the 2nd,third and fourth Criteria
  # Removing of the control sample in case it is in the main samples plate
  if(is.null(remove_samples)==F){
   inds_c<-grep(remove_samples,final_counts_df$name)
   inds_f<-grep(remove_samples,final_freqs_df$samples_name)
   matched_samples<-grep(remove_samples,final_freqs_df$samples_name,value = T)
   if((length(inds_f)==0)||(length(inds_c)==0)){
     stop(sprintf("no matching for %s in final_counts_df or final_freqs_df",remove_samples))
   }
   print(paste0("inds_f:",inds_f))
   print(paste0("matched_samples_to_remove:",matched_samples))
   final_counts_df<-final_counts_df[-inds_c,]
   final_freqs_df<-final_freqs_df[-inds_f,]
   write.csv(final_counts_df,file = sprintf("%s/final_counts_file.csv",files_path),row.names = F)
   write.csv(final_freqs_df,file = sprintf("%s/final_freqs_file.csv",files_path),row.names = F)

  }  

}

#read.csv(file="/home/rstudio/Results_Gambia_Samples_myeloid/myeloid_results/4th_batch/plate_1_samples/Stats/final_counts_file.csv")

#------------------------ Generating final document description file ---------------------------------

# function to generate the file that describe how the analysis of the samples was performed,indicating the samples excluded,and the criteria of exclusion
# it takes as input the file of the counts or the frequencies

# As selection criteria:
# - we remove the samples that Rym has flagged because of benchmarch problem (The machine has produced evident poor results): 1st criteria
# - we remove the samples that both me and Rym have identified as poor quality:second criteria
# - we select samples above certain threshold of all cells gate: third criteria 
# - we select samples above certain threshold of live gate (forth criteria)

generate_final_description_file<-function(files_path,path_Rym_data,files_path_all_plates,samples_2nd_crit=NULL,live_threshold="counts",threshold_method="mean_sd",manual_treshold_counts_all_cells=15000,
                                          manual_treshold_freqs=0.40,manual_treshold_counts_live=5000,flag_control_samples=NULL){
  #---------------- import counts and freqs file of a specific plate of a specific batch ----------------------
  path_final_counts_file <- list.files(path = files_path, recursive = F, pattern = "final_counts_file.csv", full.names = T)
  path_final_freqs_file <- list.files(path = files_path, recursive = F, pattern = "final_freqs_file.csv", full.names = T)
  if((length(path_final_counts_file)>1) ||(length(path_final_freqs_file)>1)){
    stop("something goes wrong:(length(path_final_counts_file)>1) ||(length(path_final_freqs_file)>1) is True")
  }
  final_counts_file<-read.csv(file = path_final_counts_file)
  final_freqs_file<-read.csv(file = path_final_freqs_file)
  #-------------------- import counts and freqs of all plates of all batches ---------------------
  files_path_all_counts_files<- files_path_all_plates
  path_all_final_counts_files<- list.files(path = files_path_all_counts_files, recursive = T, pattern = "final_counts_file.csv", full.names = T)
  print(path_all_final_counts_files)
  all_df_counts_all_cells<-lapply(1:length(path_all_final_counts_files),function(i){
    final_counts_file<-read.csv(path_all_final_counts_files[i])
    df_counts_plate_x_all_cells<-final_counts_file[final_counts_file$Population=="All cells",]
    return(df_counts_plate_x_all_cells)
  })
  combined_df_counts_all_cells<-rbindlist(all_df_counts_all_cells)
  all_df_counts_Live<-lapply(1:length(path_all_final_counts_files),function(i){
    final_counts_file<-read.csv(path_all_final_counts_files[i])
    df_counts_plate_x_live<-final_counts_file[final_counts_file$Population=="Live cells",]
    return(df_counts_plate_x_live)
  })
  combined_df_counts_live<-rbindlist(all_df_counts_Live)
  path_all_final_freqs_files<- list.files(path = files_path_all_counts_files, recursive = T, pattern = "final_freqs_file.csv", full.names = T)
  all_df_freqs_live<-lapply(1:length(path_all_final_freqs_files),function(i){
    final_freqs_file<-read.csv(path_all_final_freqs_files[i])
    df_freqs_plate_x_all_cells<-final_freqs_file[,which(colnames(final_freqs_file)=="Live.cells")]
    return(df_freqs_plate_x_all_cells)
  })
  combined_df_freqs_live<-unlist(all_df_freqs_live)
  print(paste0("total_number_of_samples_analyzed_by_the_automated_pipeline:",nrow(combined_df_counts_all_cells)))
  #------------------ import Rym counts (just to get the original number of samples) ---------------
  Rym_data<-read.csv(file = path_Rym_data)
  print(paste0("original_number_of_samples_current_batch:",nrow(Rym_data)))
  #----------------- Generation of the final description file ------------------
  unique_samples_name<-as.character(unique(Rym_data[,1]))
  status_column<-rep("Normal",length(unique_samples_name))
  description_column<-rep(" ",length(unique_samples_name))
  criteria_column<-rep(" ",length(unique_samples_name))
  df<-as.data.frame(cbind(unique_samples_name,status_column,description_column,criteria_column))
  colnames(df)<- c("Samples","Status","Description","Criteria")
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  inds<-grep("control",df$Samples)
  if(length(inds)>0){
    df<-df[-inds,] 
  }

  #------------------------ first criteria flagging ----------------------
  # Samples flagged by Rym before the analysis, impossible to be analyzed by 
  # the automated pipeline because there are too low counts. Very poor quality samples.
  temp_split<-strsplit(as.character(unique(final_counts_file$name)),"_")
  vec_s<-sapply(1:length(temp_split), function(i){
    string_sample_x<-temp_split[[i]][3]
    string_sample_x<-sprintf("%s$|%s_|%s\\.",string_sample_x,string_sample_x,string_sample_x)
    return(string_sample_x)
  })
  string<-paste(vec_s,collapse = "|") # to indicate that the number of the samples should be present only one time
  inds<-grep(string,df$Samples,invert = T)
  if(length(inds)!=0){
    df$Status[inds]="Flagged"
    df$Description[inds]="Impossible to analyze,too low counts"
    df$Criteria[inds]="First Criteria"
  }
  #------------------------ second criteria flagging ----------------------
  # Samples flagged by Rym and me after the analysis, because the distribution 
  # of the cells in some main gates is weird(e.g. cells near the debris). Problably poor quality samples.
  if(is.null(samples_2nd_crit)==F){
    samples_selected<-samples_2nd_crit
    inds<-grep(samples_selected,df$Samples)
    df$Status[inds]="Flagged"
    df$Description[inds]="Poor quality"
    df$Criteria[inds]="Second Criteria"
  }
  #------------------ analysis of the counts in the all cells gates (third criteria) ------------
  # Samples flagged after the analysis because the counts of the all cells gate are inferior to the threshold.
  # Threshold=mean-3*standard deviation
  All_cells_counts<-combined_df_counts_all_cells
  # set of the threshold (sd around the mean)
  if(threshold_method=="mean_sd"){
    mean_counts<-mean(All_cells_counts$Count)
    sd_counts<-sd(All_cells_counts$Count)
    threshold<-mean_counts-3*sd_counts
    print(paste0("mean_counts_all_cells:",mean_counts))
  }else if(threshold_method=="MAD"){
    # set of the threshold (absolute deviations around the median)
    median_counts<-median(sort(All_cells_counts$Count))
    absolute_deviations<-abs(All_cells_counts$Count-median_counts)
    MAD<-median(sort(absolute_deviations))
    threshold<-median_counts-3*MAD
    print(paste0("median_counts_all_cells:",median_counts))
  }else if(threshold_method=="manual"){
    threshold<-manual_treshold_counts_all_cells
  }
  
  print(paste0("treshold_all_cells_counts:",threshold))
  # find the samples with counts inferior to the threshold
  All_cells_counts<-final_counts_file[final_counts_file$Population=="All cells",]
  inds<-which(All_cells_counts$Count<threshold)
  if(length(inds)>0){
    samples_name_flagged<-as.character(All_cells_counts$name[inds])
    temp_split<-strsplit(samples_name_flagged,"_")
    vec_s<-sapply(1:length(temp_split), function(i){
      splitted_part<-temp_split[[i]][3]
      sample_string<-paste0("_",splitted_part)
      sample_string<-sprintf("%s$|%s_|%s\\.",sample_string,sample_string,sample_string)
      return(sample_string)
    })
    string<-paste(vec_s,collapse = "|")
    inds<-grep(string,df$Samples) # indices of samples flagged for the Third Criteria
    samples_status<-df$Status
    samples_Criteria<-df$Criteria
    # the samples flagged with First Criteria cannot be associated with any another Flag,because they are not analyzed.
    # the samples flagged with the Second Criteria are not evaluated for the third and fourth criteria
    inds_1<-which(samples_status!="Flagged") # indices of samples not previously flagged
    inds_2<-intersect(inds,inds_1) # indices of the samples not flagged previously associated with third criteria
    if(length(inds_2)>0){
      df$Status[inds_2]="Flagged"
      df$Description[inds_2]="Counts lower than threshold: All cells"
      df$Criteria[inds_2]="Third Criteria"
    }
  }




  #------------------ analysis of the counts in the live gates (fourth criteria) ------------
  if(live_threshold=="counts"){
    # Samples flagged after the analysis because the counts of the live cells gate are inferior to the threshold.
    Live_cells_counts<-combined_df_counts_live
    if(threshold_method=="mean_sd"){
      # set of the threshold (sd around the mean)
      mean_counts<-mean(Live_cells_counts$Count)
      sd_counts<-sd(Live_cells_counts$Count)
      threshold<-mean_counts-3*sd_counts
      print(paste0("mean_Live_cells:",mean_counts))
    }else if(threshold_method=="MAD"){
      # set of the threshold (absolute deviations around the median)
      median_counts<-median(sort(Live_cells_counts$Count))
      absolute_deviations<-abs(Live_cells_counts$Count-median_counts)
      MAD<-median(sort(absolute_deviations))
      threshold<-median_counts-3*MAD
      print(paste0("median_counts_Live_cells:",median_counts))
    }else if(threshold_method=="manual"){
      threshold<-manual_treshold_counts_live
    }

    print(paste0("treshold_Live_cells_counts:",threshold))
    
    
    # find the samples with counts inferior to the threshold
    Live_cells_counts<-final_counts_file[final_counts_file$Population=="Live cells",]
    inds<-which(Live_cells_counts$Count<threshold)
    if(length(inds)>0){
      samples_name_flagged<-as.character(Live_cells_counts$name[inds])
      temp_split<-strsplit(samples_name_flagged,"_")
      vec_s<-sapply(1:length(temp_split), function(i){
        splitted_part<-temp_split[[i]][3]
        sample_string<-paste0("_",splitted_part)
        sample_string<-sprintf("%s$|%s_|%s\\.",sample_string,sample_string,sample_string)
        return(sample_string)
      })
      string<-paste(vec_s,collapse = "|")
      inds<-grep(string,df$Samples) # indices of samples flagged for the Third Criteria
      samples_status<-df$Status
      samples_Criteria<-df$Criteria
      # the samples flagged with First Criteria cannot be associated with any another Flag,because they are not analyzed.
      # the samples flagged with the Second Criteria are not evaluated for the third and fourth criteria
      inds_1<-which(samples_status!="Flagged") # indices of samples not previously flagged
      inds_2<-intersect(inds,inds_1) # indices of the samples not flagged previously associated with the fourth criteria
      if(length(inds_2)>0){ 
        df$Status[inds_2]=paste0("Flagged")
        df$Description[inds_2]=paste0("Counts lower than threshold: Live cells")
        df$Criteria[inds_2]=paste0("Fourth Criteria")
      }
    }
  }else if(live_threshold=="freqs"){
    #------------------ analysis of the freqs in the live gates (fourth criteria) ------------
    Live_cells_freqs<-combined_df_freqs_live
    if(threshold_method=="mean_sd"){
      # set of the threshold (sd around the mean)
      mean_freqs<-mean(Live_cells_freqs)
      sd_freqs<-sd(Live_cells_freqs)
      threshold<-mean_freqs-3*sd_freqs
      print(paste0("mean_freqs_Live_cells:",mean_freqs))
    }else if(threshold_method=="MAD"){
      # set of the threshold (absolute deviations around the median)
      median_freqs<-median(sort(Live_cells_freqs))
      absolute_deviations<-abs(Live_cells_freqs-median_freqs)
      MAD<-median(sort(absolute_deviations))
      threshold<-median_freqs-3*MAD
      print(paste0("median_freqs_Live_cells:",median_freqs))
    }else if(threshold_method=="manual"){
      threshold<-manual_treshold_freqs
    }
    print(paste0("treshold_Live_cells_freqs:",threshold))
    # find the samples with freqs inferior to the threshold
    Live_cells_freqs<-final_freqs_file[,which(colnames(final_freqs_file)=="Live.cells")]
    inds<-which(Live_cells_freqs<threshold)
    if(length(inds)>0){
      samples_name_flagged<-as.character(final_freqs_file$samples_name[inds])
      temp_split<-strsplit(samples_name_flagged,"_")
      vec_s<-sapply(1:length(temp_split), function(i){
        splitted_part<-temp_split[[i]][3]
        sample_string<-paste0("_",splitted_part)
        sample_string<-sprintf("%s$|%s_|%s\\.",sample_string,sample_string,sample_string)
        return(sample_string)
      })
      string<-paste(vec_s,collapse = "|")
      inds<-grep(string,df$Samples) # indices of samples flagged for the Third Criteria
      samples_status<-df$Status
      samples_Criteria<-df$Criteria
      # the samples flagged with First Criteria cannot be associated with any another Flag,because they are not analyzed.
      # the samples flagged with the Second Criteria are not evaluated for the third and fourth criteria
      inds_1<-which(samples_status!="Flagged") # indices of samples not previously flagged
      inds_2<-intersect(inds,inds_1) # indices of the samples not flagged previously associated with the fourth criteria
      if(length(inds_2)>0){
        df$Status[inds_2]=paste0("Flagged")
        df$Description[inds_2]=paste0("freqs lower than threshold: Live cells")
        df$Criteria[inds_2]=paste0("Fourth Criteria")
      }
    }
  }
  #--------------------------------- flag control samples --------------------
  if(is.null(flag_control_samples)==F){
    ind_comp<-grep("comp",flag_control_samples) 
    if(length(ind_comp)>0){ # it means that the control sample is not in the samples plate but in the compensation plate
      row_df<-data.frame(flag_control_samples,"Flagged","Control sample in compensation plate","None")
      names(row_df)<-colnames(df)
      df<-rbind(df,row_df)
    }else{
      inds<-grep(flag_control_samples,df$Samples)
      if(length(inds)>0){
        df$Status[inds]<-paste0("Flagged")
        df$Description[inds]<-paste0("Control sample")
        df$Criteria[inds]<-paste0("None")
        
      }else{
        warning(sprintf("%s not found",flag_control_samples))
      }
    }
  }


  

  #print(df)
  # export the final description file
  write.csv(df,file = sprintf("%s/final_description_file.csv",files_path),row.names = F)
  #--------------- evaluation of the normality of the counts -------
  # shapiro_all_cells_counts<-shapiro.test(All_cells_counts$Count)
  # shapiro_live_cells_counts<-shapiro.test(Live_cells_counts$Count)
  # shapiro_all_cells_freqs<-shapiro.test(All_cells_freqs)
  # shapiro_live_cells_freqs<-shapiro.test(Live_cells_freqs)
  # print(shapiro_all_cells_counts)
  # print(shapiro_live_cells_counts)
  # print(shapiro_all_cells_freqs)
  # print(shapiro_live_cells_freqs)
}

#--------------------------------- generating final samples control file ------------------
final_sample_control<-function(files_path,comparison_plots=F,remove_pops=NULL,legend=F,cleaning_v=T,select_batch="All_batches"){
  #--------------------- generate counts control file------------------------
  path_sample_control_counts_files <- list.files(path = files_path, recursive = T, pattern = "sample_control_counts.csv", full.names = T)
  #------------ exclude sample controls maybe present inside the plates folder
  inds<-grep("Plate_[0-9]_results|plate_[0-9]_results",path_sample_control_counts_files)
  if(length(inds)>0){
    path_sample_control_counts_files<-path_sample_control_counts_files[-inds]
  }
  #----------------- choose between cleaning or no cleaning controls 
  if(cleaning_v==F){ # I want the no-cleaning version controls
    inds<-grep("no_cleaning",path_sample_control_counts_files)
    path_sample_control_counts_files<-path_sample_control_counts_files[inds]
  }else{ # I want only the cleaning version controls
    inds<-grep("no_cleaning",path_sample_control_counts_files)
    if(length(inds)>0){
      path_sample_control_counts_files<-path_sample_control_counts_files[-inds]
    }
    inds<-grep("Plate_single",path_sample_control_counts_files)
    if(length(inds)>0){
      path_sample_control_counts_files<-path_sample_control_counts_files[-inds]
    }
  }
  #print(path_sample_control_counts_files)
  #--------- remove  FCS files samples plate 2 
  ind<-grep("FCS files samples plate 2",path_sample_control_counts_files)
  if(length(ind)>0){
    path_sample_control_counts_files <- path_sample_control_counts_files[-ind] # exclude the data that refer to the FCS samples 2 folder of plate 5 in Batch 3 (they are just a continuation of the samples in FCS samples 1)
  }
  print(path_sample_control_counts_files)
  #---------------- generate single controls file 
  all_df_control_counts<-lapply(1:length(path_sample_control_counts_files),function(i){
    control_counts_file<-read.csv(path_sample_control_counts_files[i])
    # print(path_sample_control_counts_files[i])
    path<-strsplit(path_sample_control_counts_files[i],"/")[[1]]
    ind<-grep("1st_batch|Batch_1",path)
    if(length(ind)>0){ # the 1st batch has only 1 plate
      batch<-grep("batch|Batch",path,value = T)
      plate<-"plate_single"
    }else{
      plate<-grep("plate_[0-9]|Plate_[0-9]",path,value = T)
      plate<-strsplit(plate,"results")[[1]]
      plate<-strsplit(plate,"samples")[[1]]
      batch<-grep("batch|Batch",path,value = T)
      if(length(batch)==0){
        batch<-"None"
      }
    }
    current_sample_name<-as.character(unique(control_counts_file$name))
    current_sample_name<-strsplit(current_sample_name,".fcs")[[1]][1]
    current_sample_name<-strsplit(current_sample_name,"_")[[1]][3]
    new_sample_name<-sprintf("%s_%s_%s",current_sample_name,plate,batch)
    new_column_name<-rep(new_sample_name,length(control_counts_file$name))
    control_counts_file<-cbind(control_counts_file,new_column_name)
    ind<-which(colnames(control_counts_file)=="name")
    control_counts_file<-control_counts_file[,-ind]
    n_col<-length(colnames(control_counts_file))
    control_counts_file<-control_counts_file[,c(n_col,1:(n_col-1))] # move last column to first position
    colnames(control_counts_file)[1]<-"samples_name"
    # Remove Time if necessary (remeber that Time is also called Clean_0 population)
    ind<-which(as.character(control_counts_file$Population)=="Time")
    if(length(ind)>0){
      control_counts_file<-control_counts_file[-ind,]
    }
    return(control_counts_file)
  })
  single_control_counts_file<-as.data.frame(rbindlist(all_df_control_counts))
  single_control_counts_file <- data.frame(lapply(single_control_counts_file, as.character), stringsAsFactors=FALSE)
  single_control_counts_file$Count<-as.numeric(single_control_counts_file$Count)

  #------------------------------------ generate freqs controls file ----------------------
  path_sample_control_freqs_files <- list.files(path = files_path, recursive = T, pattern = "sample_control_freqs.csv", full.names = T)
  #------------ exclude sample controls maybe present inside the plates folder
  inds<-grep("Plate_[0-9]_results|plate_[0-9]_results",path_sample_control_freqs_files)
  if(length(inds)>0){
    path_sample_control_freqs_files<-path_sample_control_freqs_files[-inds]
  }
  #----------------- choose between cleaning or no cleaning controls
  if(cleaning_v==F){ # I want the no-cleaning version controls
    inds<-grep("no_cleaning",path_sample_control_freqs_files)
    path_sample_control_freqs_files<-path_sample_control_freqs_files[inds]
  }else{ # I want only the cleaning version controls
    inds<-grep("no_cleaning",path_sample_control_freqs_files)
    if(length(inds)>0){
      path_sample_control_freqs_files<-path_sample_control_freqs_files[-inds]
    }
    inds<-grep("Plate_single",path_sample_control_freqs_files)
    if(length(inds)>0){
    path_sample_control_freqs_files<-path_sample_control_freqs_files[-inds]
    }
  }
  #print(path_sample_control_freqs_files)
  #--------- remove  FCS files samples plate 2 
  ind<-grep("FCS files samples plate 2",path_sample_control_freqs_files)
  if(length(ind)>0){
    path_sample_control_freqs_files <- path_sample_control_freqs_files[-ind] # exclude the data that refer to the FCS samples 2 folder of plate 5 in Batch 3 (they are just a continuation of the samples in FCS samples 1)
  }
  #---------------- generate single controls file
  #print(path_sample_control_freqs_files)
  all_df_control_freqs<-lapply(1:length(path_sample_control_freqs_files),function(i){
    control_freqs_file<-read.csv(path_sample_control_freqs_files[i],check.names = FALSE) # check.names to avoid that read.csv generates column names with ".." between words 
    path<-strsplit(path_sample_control_freqs_files[i],"/")[[1]]
    ind<-grep("1st_batch|Batch_1",path)
    if(length(ind)>0){ # the 1st batch has only 1 plate
      batch<-grep("batch|Batch",path,value = T)
      plate<-"plate_single"
    }else{
      plate<-grep("plate_[0-9]|Plate_[0-9]",path,value = T)
      plate<-strsplit(plate,"results")[[1]]
      plate<-strsplit(plate,"samples")[[1]]
      batch<-grep("batch|Batch",path,value = T)
      if(length(batch)==0){
        batch<-"None"
      }
    }
    current_sample_name<-as.character(unique(control_freqs_file$samples_name))
    current_sample_name<-strsplit(current_sample_name,".fcs")[[1]][1]
    current_sample_name<-strsplit(current_sample_name,"_")[[1]][3]
    new_sample_name<-sprintf("%s_%s_%s",current_sample_name,plate,batch)
    new_column_name<-rep(new_sample_name,length(control_freqs_file$samples_name))
    control_freqs_file<-cbind(control_freqs_file,new_column_name)
    ind<-which(colnames(control_freqs_file)=="samples_name")
    control_freqs_file<-control_freqs_file[,-ind]
    n_col<-length(colnames(control_freqs_file))
    control_freqs_file<-control_freqs_file[,c(n_col,1:(n_col-1))]
    colnames(control_freqs_file)[1]<-"samples_name"
    
    # remove Time pop counts if necessary
    ind<-which(colnames(control_freqs_file)=="Time")
    if(length(ind)>0){
      control_freqs_file<-control_freqs_file[,-ind]
    }
    
    return(control_freqs_file)
  })
  single_control_freqs_file<-as.data.frame(rbindlist(all_df_control_freqs))

  if(is.null(remove_pops)==F){
    inds<-grep(remove_pops,single_control_counts_file$Population)
    if(length(inds)==0){
      stop("selected pops to remove not found")
    }
    single_control_counts_file<-single_control_counts_file[-inds,]
    inds<-grep(remove_pops,colnames(single_control_freqs_file))
    single_control_freqs_file<-single_control_freqs_file[,-inds]
  }
  print(single_control_counts_file)
  # check correctness of operations
  if(nrow(single_control_freqs_file)>16){
    ind<-grep("Bcells",path_sample_control_freqs_files[1])
    if(length(ind)>0){
      if(nrow(single_control_freqs_file)>17)
        # Bcells panel :one control for each plate,so 1 ctrl batch 1, 6 ctrl batch 2, 5 ctrl batch 3, 5 ctrls batch 4 total is 17 controls
        stop("Something goes wrong: nrow(single_control_freqs_file) has more than 17 entries ")
    }else{
      # myeloid panel :one control for each plate,so 1 ctrl batch 1, 5 ctrl batch 2, 5 ctrl batch 3, 5 ctrls batch 4 total is 16 controls
      stop("Something goes wrong: nrow(single_control_freqs_file) has more than 16 entries ")
      
    }
  }
  #------------- make comparison plot control samples ----------
  if(comparison_plots==T){
    inds_2ndbatch<-grep("2nd_batch|Batch_2",single_control_counts_file$samples_name)
    # print(single_control_counts_file$samples_name)
    # print(inds_2ndbatch)
    inds_3rdbatch<-grep("3rd_batch|Batch_3",single_control_counts_file$samples_name)
    # #------------------ plots batch 2
    if(select_batch=="Batch_2"){
      single_control_counts_file_batch_2<-single_control_counts_file[inds_2ndbatch,]
      inds<-grep("Plate_1_2nd_batch|Plate_2_Batch_2",single_control_counts_file_batch_2$samples_name)
      x_coord_line_plot<-c(1:length(inds))
      x_coord_line_plot<-as.integer(rep(x_coord_line_plot,6))
      x_coord_labels<-c(1:length(inds))
      pops<-single_control_counts_file$Population[inds]
      color_pops<-colorRampPalette(c("red", "green","blue","yellow"))(length(inds))
      ind<-grep("Bcells",path_sample_control_freqs_files[1])
      if(length(ind)>0){ # code for the Bcells control comparison
        color_samples<-c("red","pink","green","blue","yellow","orange")
        single_control_counts_file_batch_2<-cbind(single_control_counts_file_batch_2,x_coord_line_plot)
        single_control_counts_file_batch_2$samples_name<-factor(single_control_counts_file_batch_2$samples_name)
        print(levels(single_control_counts_file_batch_2$samples_name))
        levels(single_control_counts_file_batch_2$samples_name)<-levels(single_control_counts_file_batch_2$samples_name)[c(2,1,3,4,5,6)]
        print(levels(single_control_counts_file_batch_2$samples_name))
        levels(single_control_counts_file_batch_2$samples_name)<-c("Plate_1","Plate_2","Plate_3","Plate_4","Plate_5","Plate_6")
        print(levels(single_control_counts_file_batch_2$samples_name))
      }else{ # code for the myeloid control comparison
        color_samples<-c("red", "green","blue","yellow","orange")
        single_control_counts_file_batch_2<-cbind(single_control_counts_file_batch_2,x_coord_line_plot)
        single_control_counts_file_batch_2$samples_name<-factor(single_control_counts_file_batch_2$samples_name)
        print(levels(single_control_counts_file_batch_2$samples_name))
        levels(single_control_counts_file_batch_2$samples_name)<-levels(single_control_counts_file_batch_2$samples_name)[c(2,1,3,4,5)]
        print(levels(single_control_counts_file_batch_2$samples_name))
        levels(single_control_counts_file_batch_2$samples_name)<-c("Plate_1","Plate_2","Plate_3","Plate_4","Plate_6")
        print(levels(single_control_counts_file_batch_2$samples_name))
      }

      #Note: group=1,because we have only each group is composed by only 1 element (a population for each group, a total of 31 groups for each plot)
      plot_batch_2_line_plot_plus_ribbon<-ggplot(single_control_counts_file_batch_2) + geom_line(aes(x=x_coord_line_plot,y=Count,color=factor(x_coord_line_plot),group=1),size=1) + scale_x_continuous(breaks=x_coord_labels,labels=pops) + geom_ribbon(aes(x=x_coord_line_plot,ymin=0,ymax=Count),alpha=0.3)+facet_grid(samples_name~.)+
        theme(legend.title=element_blank()) + labs(x="Population") + labs(y="Count") + scale_color_manual(values=color_pops,labels=pops) + theme(legend.position = "none") + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.x = element_text(angle = 90))
      plot_batch_2_line_plot_plus_ribbon_v2<-ggplot(single_control_counts_file_batch_2) + geom_line(aes(x=x_coord_line_plot,y=Count,group=1),size=0.7) + scale_x_continuous(breaks=x_coord_labels,labels=pops) + geom_ribbon(aes(x=x_coord_line_plot,ymin=0,ymax=Count,fill=samples_name),alpha=0.5)+facet_grid(samples_name~.,switch = "both")+
        theme(legend.title=element_blank()) + labs(x="Population") + labs(y="Count") + scale_fill_manual(values=color_samples,guide=FALSE) + theme(legend.position = "none") + theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(position = "right")
      # switch="bot" to put move the labels of the samples on the left.
      show(plot_batch_2_line_plot_plus_ribbon_v2)
    }

     
    # #------------------ plots batch 3
    if(select_batch=="Batch_3"){
      single_control_counts_file_batch_3<-single_control_counts_file[inds_3rdbatch,]
      inds<-grep("Plate_1_3rd_batch",single_control_counts_file_batch_3$samples_name)
      x_coord_line_plot<-c(1:length(inds))
      x_coord_line_plot<-as.integer(rep(x_coord_line_plot,5))
      x_coord_labels<-c(1:length(inds))
      pops<-single_control_counts_file$Population[inds]
      color_pops<-colorRampPalette(c("red", "green","blue","yellow"))(length(inds))
      color_samples<-c("red", "green","blue","yellow","orange")
      single_control_counts_file_batch_3<-cbind(single_control_counts_file_batch_3,x_coord_line_plot)
      single_control_counts_file_batch_3$samples_name<-factor(single_control_counts_file_batch_3$samples_name)
      levels(single_control_counts_file_batch_3$samples_name)<-c("Plate_1","Plate_2","Plate_3","Plate_4","Plate_5")
      #Note: group=1,because we have only each group is composed by only 1 element (a population for each group, a total of 31 groups for each plot)
      plot_batch_3_line_plot_plus_ribbon<-ggplot(single_control_counts_file_batch_3) + geom_line(aes(x=x_coord_line_plot,y=Count,color=factor(x_coord_line_plot),group=1),size=1) + scale_x_continuous(breaks=x_coord_labels,labels=pops) + geom_ribbon(aes(x=x_coord_line_plot,ymin=0,ymax=Count),alpha=0.5)+facet_grid(samples_name~.)+
        theme(legend.title=element_blank()) + labs(x="Population") + labs(y="Count") + scale_color_manual(values=color_pops,labels=pops) + theme(legend.position = "none") + theme(axis.text.x = element_text(size = 7)) + theme(axis.text.x = element_text(angle = 90))
      plot_batch_3_line_plot_plus_ribbon_v2<-ggplot(single_control_counts_file_batch_3) + geom_line(aes(x=x_coord_line_plot,y=Count,group=1),size=0.7) + scale_x_continuous(breaks=x_coord_labels,labels=pops) + geom_ribbon(aes(x=x_coord_line_plot,ymin=0,ymax=Count,fill=samples_name),alpha=0.5)+facet_grid(samples_name~.,switch = "both")+
        theme(legend.title=element_blank()) + labs(x="Population") + labs(y="Count") + scale_fill_manual(values=color_samples,guide=FALSE) + theme(legend.position = "none") + theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(position = "right")
     
      # show(plot_batch_3_line_plot_plus_ribbon_v2)
    }

    #------------------------- plots all batches
    if(select_batch=="All_batches"){
      single_control_counts_file$samples_name<-factor(single_control_counts_file$samples_name)
      inds<-grep("Plate_1_3rd_batch|Plate_1_Batch_3",single_control_counts_file$samples_name)
      n_plates<-length(unique(single_control_counts_file$samples_name))
      x_coord_line_plot<-c(1:length(inds))
      x_coord_line_plot<-as.integer(rep(x_coord_line_plot,n_plates))
      x_coord_labels<-c(1:length(inds))
      pops<-single_control_counts_file$Population[inds]
      #print(length(x_coord_line_plot))
      print(levels(single_control_counts_file$samples_name))
      ind<-grep("Bcells",path_sample_control_freqs_files[1])
      if(length(ind)>0){
        single_control_counts_file$samples_name<-factor(single_control_counts_file$samples_name,levels=levels(single_control_counts_file$samples_name)[c(12,11,6,13,2,8,7,5,4,1,17,16,15,14,3,10,9)])
      }else{
        single_control_counts_file$samples_name<-factor(single_control_counts_file$samples_name,levels=levels(single_control_counts_file$samples_name)[c(12,11,6,13,2,8,7,5,4,1,16,15,14,3,10,9)])
      }
      print(levels(single_control_counts_file$samples_name))
      # inds_prova<-grep("2nd_batch",single_control_counts_file$samples_name)
      # print(single_control_counts_file[inds_prova,])
      if(length(ind)>0){
        levels(single_control_counts_file$samples_name)<-c("Plate5b4","Plate4b4","Plate3b4","Plate2b4","Plate1b4","Plate5b3","Plate4b3","Plate3b3","Plate2b3","Plate1b3","Plate6b2","Plate5b2","Plate4b2","Plate3b2","Plate2b2","Plate1b2","Plate1b1")
      }else{
        levels(single_control_counts_file$samples_name)<-c("Plate5b4","Plate4b4","Plate3b4","Plate2b4","Plate1b4","Plate5b3","Plate4b3","Plate3b3","Plate2b3","Plate1b3","Plate6b2","Plate4b2","Plate3b2","Plate2b2","Plate1b2","Plate1b1")
      }
      print(levels(single_control_counts_file$samples_name))
      if(length(ind)>0){
        color_samples<-c("red", "green","blue","yellow","orange","red", "green","blue","yellow","orange","red","pink","green","blue","yellow","orange","red")
      }else{
        color_samples<-c("red", "green","blue","yellow","orange","red", "green","blue","yellow","orange","red","green","blue","yellow","orange","red")
      }
      print(colnames(single_control_counts_file))
      plot_all_batches_line_plot_plus_ribbon_v2<-ggplot(single_control_counts_file) + geom_line(aes(x=x_coord_line_plot,y=Count,group=1),size=0.7) + scale_x_continuous(breaks=x_coord_labels,labels=pops) + geom_ribbon(aes(x=x_coord_line_plot,ymin=0,ymax=Count,fill=samples_name),alpha=0.5)+facet_grid(samples_name~.,switch = "both")+
         theme(legend.title=element_blank()) + labs(x="Population") + labs(y="Count") + scale_fill_manual(values=color_samples,guide=FALSE) + theme(legend.position = "none") + theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(position = "right")
      show(plot_all_batches_line_plot_plus_ribbon_v2)
      # inds<-grep("b2",single_control_counts_file$samples_name)
      # print(single_control_counts_file[inds,])
      }
      # calculate standard deviations
      print(single_control_counts_file)
      unique_pops<-unique(single_control_counts_file$Population)
      single_control_counts_file$ParentCount<-as.numeric(single_control_counts_file$ParentCount)
      vec_results<-rep(0,length(unique_pops))
      i<-0
      for(p in unique_pops){
        i<-i+1
        print(sprintf("analysis of pop %s",p))
        inds<-which(single_control_counts_file$Population==p)
        single_control_counts_file_pop_x<-single_control_counts_file[inds,]
        # out<-kruskal.test(single_control_counts_file_pop_x$Count,single_control_counts_file_pop_x$samples_name)
        # value<-out$p.value
        value_counts<-sd(single_control_counts_file_pop_x$Count)
        value_freqs<-sd(single_control_counts_file_pop_x$Count/single_control_counts_file_pop_x$ParentCount)
        vec_results[i]<-value_freqs
      }
      print(vec_results)
      average_sd<-mean(vec_results)
      print(average_sd)
    }
  
}

#----------------------- other useful things ------------------------

interactive_boxplots_all_counts<-function(df,df_outliers){
  # #----- make boxplot ------
  # with outlier
  # myeloid panel
  # p<<-ggplot(df,aes(x=Population,y=Count)) + geom_boxplot()  + theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) +
  # ggtitle("boxplot pops myeloid panel:1379 samples analyzed,average of 49 outliers (3.5%) ")
  # b cell panel
  p<-ggplot() + theme(axis.text.x = element_text(face="bold",size=7,angle = 90)) +
    geom_boxplot(data=df,aes(x=Population,y=Count), outlier.shape = NA,na.rm = T) +
    geom_point_interactive(data=df_outliers,aes(x=n_points,y=value,tooltip=n_points,onclick=v_onclick_string),na.rm=T,size=0.5)
  x<-girafe(ggobj = p)
  print(x)
  #show(p)
}

#-------------------- function to build the df of the outliers -----------------------
build_df_outliers<-function(df){
  p_temp<-ggplot(df,aes(x=Population,y=Count)) + geom_boxplot()
  out<-ggplot_build(p_temp)$data[[1]]
  list_outliers<-as.data.frame(out)$outliers
  pops<-unique(df$Population)
  pops<-levels(pops)
  names(list_outliers)<-pops
  dfoutliers<-plyr::ldply(list_outliers, rbind)
  dfoutliers<-transpose(dfoutliers)
  colnames(dfoutliers)<-dfoutliers[1,]
  dfoutliers<-dfoutliers[-1,]
  row.names(dfoutliers)<-NULL
  dfoutliers<-reshape2::melt(dfoutliers,measure.vars=colnames(dfoutliers))
  unique_pops<-as.character(unique(dfoutliers$variable))
  #-------- we add the onclick strings on the df outliers
  list_all_df<-list()
  p<-0
  for(pop in unique_pops){
    p<-p+1
    inds<-grep(pop,dfoutliers$variable)
    dfoutliers_current_pop<-dfoutliers[inds,]
    n_points<-rep(p,nrow(dfoutliers_current_pop))
    dfoutliers_current_pop<-cbind(dfoutliers_current_pop,n_points)
    value_outliers_current_pop<-dfoutliers_current_pop$value
    v_onclick_string<-c()
    for(value in value_outliers_current_pop){
      if(is.na(value)==F){
        inds<-grep(pop,df$Population)
        df_current_pop<-df[inds,]
        regex_value<-sprintf("^%s$",value)
        ind<-grep(regex_value,as.character(df_current_pop$Count))
        onclick_string<-df_current_pop$onclick[ind]
      }else{
        onclick_string<-"None"
      }
      # print(ind)
      # print(onclick_string)
      v_onclick_string<-append(v_onclick_string,onclick_string)
    }
    
    dfoutliers_current_pop_plus_onclick<-cbind(dfoutliers_current_pop,v_onclick_string)
    list_all_df[[pop]]<-dfoutliers_current_pop_plus_onclick
  }
  dfoutliers<-as.data.frame(rbindlist(list_all_df,fill = T))
  dfoutliers$value<-as.numeric(dfoutliers$value)
  dfoutliers$v_onclick_string<-as.character(dfoutliers$v_onclick_string)
  dfoutliers$v_onclick_string<-sprintf("grid::grid.raster(readPNG(%s))",df_outliers$v_onclick_string)
  return(dfoutliers)
}

# function to remove uninterested pop and to add the onclick column
pre_process_df<-function(list_path_all_results_folder,df){
  #-------- remove pops we don't want to plot ---------
  pops_to_remove<-c("CD64-","Time","Singlets","Beads")
  string<-paste0(pops_to_remove,collapse = "|")
  inds_pops_to_remove<-grep(string,df$Population)
  df<-df[-inds_pops_to_remove,]
  #----------- set correct class for each variable -------
  df$Population<-factor(df$Population)
  df$cells_ul_blood<-as.numeric(df$cells_ul_blood)
  df$Count<-as.numeric(df$Count)
  #------------ find path plots -------------------------
  inds<-grep("Plots",list_path_all_results_folder)
  path_all_plots<<-list_path_all_results_folder[inds]
  unique_samples<-unique(df[,1])
  stringsplitted<-strsplit(unique_samples,".fcs")
  vec_unique_name<-c()
  for(i in 1:length(stringsplitted)){
    name<-stringsplitted[[i]][1]
    vec_unique_name<-append(vec_unique_name,name)
  }
  list_all_dfs<-list()
  i<-0
  for(name_sample in vec_unique_name){
    i<-i+1
    path_plot_current_sample<-grep(name_sample,path_all_plots,value = T)
    inds<-grep(name_sample,df[,1])
    df_current_sample<-df[inds,]
    if(length(path_plot_current_sample)>1){
      batch_current_sample<-unique(df_current_sample$Batch)
      plate_current_sample<-unique(df_current_sample$Plate)
      string<-paste(batch_current_sample,plate_current_sample,sep="|")
      onclick<-c()
      for(s in string){
        path_plot_current_sample_s<-grep(s,path_plot_current_sample,value = T)
        onclick_s<-rep(path_plot_current_sample_s,22)
        onclick<-append(onclick,onclick_s)
      }
      onclick<-paste0(onclick,sep = ",")
      df_current_sample<-cbind(df_current_sample,onclick)
    }else{
      onclick<-rep(path_plot_current_sample,nrow(df_current_sample))
      df_current_sample<-cbind(df_current_sample,onclick)
    }
    
    list_all_dfs[[i]]<-df_current_sample
  }
  df<-as.data.frame(rbindlist(list_all_dfs))
  df$onclick<-as.character(df$onclick)
  return(df)
}

# function to generate CLR data
generate_CLR_data<-function(fs,cellpop,name_pop,channels,filters=F,markers=NA,density=F,path.output="./"){
  name_pop_1<-strsplit(name_pop,"_")[[1]][1]
  path.output<-paste0(path.output,sprintf("/CLR-%s",name_pop_1))
  suppressWarnings ( dir.create (path.output,recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output,"/data/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output,"/plots/"),recursive = T ))
  suppressWarnings ( dir.create ( paste0(path.output,"/clrs/"),recursive = T ))
  output<-sapply(1:length(fs),function(i){
    if(class(fs)=="flowSet"){
      frame<-fs[[i]]
    }else{ # fs is a list of cell population objects
      frame<-fs[[i]]@flow.frame
    }
    ind_1<-grep(channels[1],frame@parameters@data$name)
    if(length(ind_1)==0){
      print("channels[1] not found, checking by markers")
      if(is.na(markers[1])==T){
        stop("markers[1] argument is empty")
      }
      ind_1_marker<-grep(markers[1],as.vector(frame@parameters@data$desc))
      print(ind_1_marker)
      channels[1]<-frame@parameters@data$name[ind_1_marker]
    }
    ind_2<-grep(channels[2],frame@parameters@data$name)
    if(length(ind_2)==0){
      print("channels[2] not found, checking markers")
      if(is.na(markers[2])==T){
        stop("markers[2] argument is empty")
      }
      ind_2_marker<-grep(markers[2],as.vector(frame@parameters@data$desc))
      channels[2]<-frame@parameters@data$name[ind_2_marker]
    }
    print(sprintf("channels after checking:%s",paste0(channels,collapse = ",")))
    print(sprintf("analysis populations for sample %s",identifier(frame)))
    # generation of the events assignments csv file
    list_gateAssignments<-list()
    CL<-c()
    for(p in 1:length(cellpop)){
      print(sprintf("extract gate pop %s",strsplit(name_pop,"_")[[1]][p]))
      cellpop_p<-cellpop[[p]]
      cellpop_i<-cellpop_p[[i]]
      gateAssignments <-1:nrow(frame) %in% cellpop_i@index
      list_gateAssignments[[p]]<-gateAssignments
    }

    cl_temp<-do.call("cbind",list_gateAssignments)
    CL <- cbind(CL,  cl_temp)
    pops_name<-strsplit(name_pop,"_")[[1]]
    colnames(CL) <- pops_name
    if(class(fs)!="flowSet"){
      df_data<-exprs(frame)[,channels]
      inds_na<-which(is.na(df_data)==T)
      CL<-CL[-inds_na,]
    }
    write.table(CL,file = paste0(paste0(path.output,"/clrs/"),
                                   gsub(identifier(frame),pattern = ".fcs",replacement = ".csv")),row.names=FALSE, sep=",")
    # generation of the single bivariate plots
    png(filename=paste0(paste0(path.output,"/plots/"),gsub(identifier(frame),pattern = ".fcs",replacement = ".png")),width = 800,height = 800)
    if(density==F){
      plotDens(frame, channels)
      if(is.list(filters)==F){
        for(p in 1:length(cellpop)){
          cellpop_p<-cellpop[[p]]
          cellpop_i<-cellpop_p[[i]]
          lines(cellpop_i@filter,type="l", lwd=2)
        }
      }else{
        for(p in 1:length(filters)){
          filter_p<-filters[[p]]
          filter_i<-filter_p[[i]]
          lines(filter_i@boundaries,type="l", lwd=2)
          #points(exprs(frame)[c(1,30),channels],pch=".",col="red",cex=9)
        } 
      }
    }else{
      plot.density(frame,channels)
      for(p in 1:length(cellpop)){
        cellpop_p<-cellpop[[p]]
        cellpop_i<-cellpop_p[[i]]
        lines(cellpop_i@filter,type="l", lwd=2)
      }
    }
    dev.off()
    # generation of the expression data of the bivariate plots
    df_data<-exprs(frame)[,channels]
    if(class(fs)!="flowSet"){
      inds_na<-which(is.na(df_data)==T)
      df_data<-df_data[-inds_na,]
    }
    write.table(df_data,file=paste0(paste0(path.output,"/data/"), gsub(identifier(frame),pattern = ".fcs",replacement = ".csv")), row.names=FALSE, sep=",")
      
  })
  
  
}


#------------------- fix markers names desc column slot -----------------------
# path_fcs_files_fixing indicates the path of fcs files to be used to fix the description column of current fsc files.
fix_markers_names<-function(fs_to_fix,path_fcs_files_fixing){
  # The flowFrame of Bcells samples Gambia presents the desc column with all NA. 
  # We assume that the names of the markers are equal to the PNG samples ones
  # I take only a png sample
  files_path_PNG_cells <-  paste0(path_fcs_files_fixing)
  all_files_list_paths <- list.files(files_path_PNG_cells, pattern = ".fcs",full.names = TRUE)
  f_png <- read.FCS(all_files_list_paths[1])
  fs<-fs_to_fix
  list_flowFrames_modified <- lapply(1:length(fs),function(i){
    f <- fs[[i]]
    v_all_chans <- as.vector(pData(f@parameters)$name)
    v_no_fl_chans_indx <- grep("FSC+|SSC+|Time+",v_all_chans) 
    if (anyNA(pData(f@parameters)$desc[-v_no_fl_chans_indx])){
      
      all_channel_names_PNG <- as.vector(f_png@parameters@data$name)
      v_no_fl_chans_indx_PNG <- grep("FSC+|SSC+|Time+",all_channel_names_PNG)
      channel_names_PNG <- all_channel_names_PNG[-v_no_fl_chans_indx_PNG]
      all_marker_names_PNG <- as.vector(f_png@parameters@data$desc)
      marker_names_PNG <- all_marker_names_PNG[-v_no_fl_chans_indx_PNG]
      df_channel_vs_marker <- data.frame(channel_names_PNG,marker_names_PNG,stringsAsFactors = FALSE) 
      
      all_chans_names <- as.vector(pData(f@parameters)$name)
      for(index in (1:length(df_channel_vs_marker$channel_names))){
        chan <- df_channel_vs_marker$channel_names[index] 
        marker <- df_channel_vs_marker$marker_names[index] 
        chan_indx <- grep(chan,all_chans_names) 
        pData(f@parameters)$desc[chan_indx] <- marker 
      }
    }
    # if no NA in the desc column for the fluorescent channel,no fix is needed
    else{
      return(f)
    }
    if (anyNA(pData(f@parameters)$desc[-v_no_fl_chans_indx])){
      stop("fix markers names of desc column failed,please check the flowFrame")
    }
    return(f)
  })
  
  fs <- as(list_flowFrames_modified,"flowSet")
  return(fs)
}

# make dataframe to store the thresholds calculated by flowDensity (Myeloid panel)
make_tresholds_df_myeloid<-function(gs){
  #----- make empty threholds dataframe 
  columns_name_myeloid<-c("All cells","B cells","Basophils","CD11b+CD16+ Mature Neutrophils","CD11b+CD16- Granulocytes","CD11b-CD16+ Immature Neutrophils 2",
                          "CD11-CD16- Immature Neutrophils 1","CD3+ T cells","CD3-","CD45+CD66- (Non Granulocytes)","CD56 Hi NK","CD56- dim CD16- NK","CD56- dim CD16+ NK",
                          "CD56+CD16+ NKT cells","CD56-CD16+ NK","CD56-CD16- cells","CD64+","CD64-","Classical Monocytes","Time","Granulocytes",
                          "HLADR+CD14+ Monocytes","HLADR+ CD14-","HLADR-CD14-","Live cells","Non classical Monocytes","Singlets","gd T cells","gd- T cells","mDC",
                          "pDC","root")
  m <- matrix(0, ncol = length(columns_name_myeloid) , nrow = length(gs))
  df_thresholds<-data.frame(m) # it will contain the thresholds of all channels for all samples
  colnames(df_thresholds)<-columns_name_myeloid
  row.names(df_thresholds)<-sampleNames(gs)
  return(df_thresholds)
}

# make dataframe to store the thresholds calculated by flowDensity (Bcells panel)
make_tresholds_df_bcells<-function(gs){
  #----- make empty threholds dataframe 
  columns_name_Bcells<-c("All cells","Atypical B cells","Blasts","CD10+","CD10-","CD19+ B cells",
                         "CD19+CD20+","CD19+CD20-","Time","Granulocytes","IGD-IGM- B cells","IGD+IGM- B cells","IGD-IGM+ B cells",
                         "IGD+IGM+ B cells","IGM","Immature transition B cells","Live cells","Naive B cells","Non granulocytes","Plasmablasts",
                         "Singlets","Switched memory B cells-","Unswitched memory B cells","plasma cells","root")
  m <- matrix(0, ncol = length(columns_name_Bcells) , nrow = length(gs))
  df_thresholds<-data.frame(m) # it will contain the thresholds of all channels for all samples
  colnames(df_thresholds)<-columns_name_Bcells
  row.names(df_thresholds)<-sampleNames(gs)
  return(df_thresholds)
}

