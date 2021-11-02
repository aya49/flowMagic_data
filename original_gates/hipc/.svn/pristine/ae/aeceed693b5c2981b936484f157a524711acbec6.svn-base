############# import libraries ##################
library(ggplot2)
library(ggiraph)
library(png)
library(grid)
library(data.table)
library(dplyr)
library(htmlwidgets)

############################################# define functions ##############################

import_all_csv_files_same_type<-function(list_path_all_results_folder,type,panel,ML_like="None"){
  
  inds_counts_files<-grep("final_counts",list_path_all_results_folder)
  vec_all_counts_files_results<-list_path_all_results_folder[inds_counts_files]
  inds_freqs_files<-grep("final_freqs",list_path_all_results_folder)
  vec_all_freqs_files_results<-list_path_all_results_folder[inds_freqs_files]
  if(type=="counts"){
    list_all_df<-list()
    c<-0
    for(f in vec_all_counts_files_results){
      c<-c+1
      stringsplitted<-strsplit(f,"/")[[1]]
      ind_batch<-grep("Batch|batch",stringsplitted)
      batch<-stringsplitted[ind_batch]
      inds<-grep("batches",batch) # esclude strings containing batches
      batch<-batch[-inds]
      ind_plate<-grep("Plate|plate",stringsplitted)
      if(length(ind_plate)!=0){
        plate<-stringsplitted[ind_plate]
        if(length(plate)>1){ # case of 3rd batch plate 5 (because there two subplate) 
          plate<-plate[1]
        }
      }else{
        plate<-"None"
      }
      df <- read.csv(f)
      ind_X<-grep("^X$",colnames(df)) 
      if(length(ind_X)>0){ # one column contains the row names
        df<-df[,-ind_X]
      }
      
      # add column of cell/ul blood
      inds<-which(df$Population=="Beads")
      inds_1<-which(df$Population=="Time") # first population of each sample
      if(panel=="myeloid"){
        n_populations<-31 # 32 populations in each sample included time,so +31, or 31 populations included Singlets, so +30
        if(length(inds_1)==0){
          inds_1<-which(df$Population=="Singlets") # first population of each sample is Singlets (Time removed for this plate)
          n_populations<-30
        }
      }else if(panel=="Bcells"){
        n_populations<-24 # 25 populations in each sample included time,so +24, or 24 populations included Singlets, so +23
        if(length(inds_1)==0){
          inds_1<-which(df$Population=="Singlets") # first population of each sample is Singlets (Time removed for this plate)
          n_populations<-23
        }
      }else{
        stop("error input: incorrect panel value")
      }
      final_column_cells_ulbloods<-c()
      s<-0
      for(indstime in inds_1){
        s<-s+1 # counter of the sample
        counts_all_pops_sample_x<-df$Count[indstime:(indstime+n_populations)] 
        
        ind_bead_current_sample<-inds[s]
        counts_beads_sample_x<-df$Count[ind_bead_current_sample]
        temp_counts<-counts_all_pops_sample_x/counts_beads_sample_x
        C<-480000
        D<-200
        dilution<-112.5
        temp_counts<-temp_counts*C
        temp_counts<-temp_counts/D
        final_counts<-temp_counts/dilution
        final_column_cells_ulbloods<-append(final_column_cells_ulbloods,final_counts)
      }
      final_column_cells_ulbloods<-format(final_column_cells_ulbloods,scientific = F)
      final_column_cells_ulbloods<-as.numeric(final_column_cells_ulbloods)
      final_column_cells_ulbloods<-round(final_column_cells_ulbloods,2)
      
      df<-cbind(df,final_column_cells_ulbloods)
      
      ind<-grep("final_column_cells_ulbloods",colnames(df))
      colnames(df)[ind]<-"cells_ul_blood"
      # add plate and batch
      vec_batch<-rep(batch,nrow(df))
      df<-cbind(df,vec_batch)
      
      vec_plate<-rep(plate,nrow(df))
      df<-cbind(df,vec_plate)
      inds<-grep("vec_batch|vec_plate",colnames(df))
      colnames(df)[inds]<-c("Batch","Plate")
      
      list_all_df[[c]]<-df
    }
    df<-as.data.frame(rbindlist(list_all_df))
    # we can also represent the counts in a Machine learning like structure (row=samples,cols=features).
    if(ML_like!="None"){
      unique_samples<-unique(df$name)
      list_df_reshaped<-list()
      list_df_reshaped<-mclapply(1:length(unique_samples),function(i){
        current_sample_name<-unique_samples[i]
        inds_current_sample<-grep(current_sample_name,df$name)
        df_current_sample<-df[inds_current_sample,]
        if(ML_like=="counts"){
          df_current_sample<-df_current_sample[,c("name","Population","Count")]
          df_current_sample_reshaped<-reshape(df_current_sample,idvar = "name",timevar = "Population",direction = "wide")
        }else if(ML_like=="ul_blood"){
          df_current_sample<-df_current_sample[,c("name","Population","cells_ul_blood")]
          df_current_sample_reshaped<-reshape(df_current_sample,idvar = "name",timevar = "Population",direction = "wide")
        }
        return(df_current_sample_reshaped)
      },mc.cores = 40)
      df<-as.data.frame(rbindlist(list_df_reshaped,fill = T))
    }
  }else{
    list_all_df<-list()
    c<-0
    for(f in vec_all_freqs_files_results){
      c<-c+1
      stringsplitted<-strsplit(f,"/")[[1]]
      ind_batch<-grep("Batch|batch",stringsplitted)
      batch<-stringsplitted[ind_batch]
      inds<-grep("batches",batch) # esclude strings containing batches
      batch<-batch[-inds]
      ind_plate<-grep("Plate|plate",stringsplitted)
      if(length(ind_plate)!=0){
        plate<-stringsplitted[ind_plate]
        if(length(plate)>1){ # case of 3rd batch plate 5 (because there two subplate) 
          plate<-plate[1]
        }
      }else{
        plate<-"None"
      }
      df <- read.csv(f)
      ind_X<-grep("^X$",colnames(df))
      if(length(ind_X)>0){
        df<-df[,-ind_X]
      }
      # round the frequencies
      df[,2:ncol(df)]<-round(df[,2:ncol(df)],2)
      # add plate and batch info
      vec_batch<-rep(batch,nrow(df))
      df<-cbind(df,vec_batch)
      vec_plate<-rep(plate,nrow(df))
      df<-cbind(df,vec_plate)
      inds<-grep("vec_batch|vec_plate",colnames(df))
      colnames(df)[inds]<-c("Batch","Plate")
      list_all_df[[c]]<-df
    }
    df<-as.data.frame(rbindlist(list_all_df,fill = T))
  }
  df[, ] <- lapply(df[, ], as.character)
  return(df)
}



plot_results <- function(alldata,showoutliers,outlierdata,popstoplot,samplestoplotdata, shiny){
  
  newdata <- alldata[1,]
  #unique(alldata$population)
  for(i in 1:length(popstoplot)){
    newdata <- rbind(newdata,alldata[which(alldata$population==popstoplot[i]),])
  }
  alldata <- unique(newdata)
  
  if(showoutliers==T && length(outlierdata[,1]) > 0){ #if there are outliers to plot...
    
    if(!is.null(samplestoplotdata)){  #...and user-specified samples to plot, add multiple sets of geom_point_interactive.
      r <- ggplot() +
        geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=population,
                                                  tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                                "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
        geom_point_interactive(position=position_jitter(width=0.45,height=.45),size=.7,data=outlierdata,
                               aes(x = population, y = proportion,colour=population,
                                   tooltip=paste0(round(proportion,digits=2),"%, ",file),data_id=file))+
        geom_point_interactive(position=position_jitter(width=0.1,height=.1),data=samplestoplotdata,size=0.7,
                               aes(x = population, y = proportion,color=samplestoplotdata$file,fill="black",
                                   tooltip=paste0(round(proportion,digits=2),"%, ",file),data_id=file))+
        scale_y_continuous(name="Proportion % (parent population)",limits = c(0,100),breaks=seq(0,100,by=10),expand=c(0,0))+
        xlab(label = "Cell Population")+
        theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
              axis.text.y = element_text(hjust = 1, size=10,color=1),
              legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
      if (!shiny){
        r <- girafe(ggobj=r)
        r <- girafe_options(r, opts_toolbar(saveaspng = TRUE),
                            #specifies properties of outlier dots when cursor hovers over them
                            opts_hover(css = "stroke:black;stroke-width:5pt"))
      }
    }else{  #...and no user-specified samples to plot, only add outlier set of geom_point_interactive.
      r <- ggplot() +
        geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=population,
                                                  tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                                "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
        geom_point_interactive(position=position_jitter(width=0.45,height=.45),size=.7,data=outlierdata,
                               aes(x = population, y = proportion,colour=population,
                                   tooltip=paste0(round(proportion,digits=2),"%, ",file),data_id=file))+
        scale_y_continuous(name="Proportion % (parent population)",limits = c(0,100),breaks=seq(0,100,by=10),expand=c(0,0))+
        xlab(label = "Cell Population")+
        theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
              axis.text.y = element_text(hjust = 1, size=10,color=1),
              legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
      if (!shiny){
        r <- girafe(ggobj=r)
        r <- girafe_options(r, opts_toolbar(saveaspng = TRUE),
                            #specifies properties of outlier dots when cursor hovers over them
                            opts_hover(css = "stroke:black;stroke-width:5pt;r:3pt"))
      }
    }
  }else{  #if there are no outlier points...
    if(length(samplestoplotdata[,1])>0){  #...but user-specified samples to plot
      r <- ggplot() +
        geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=population,
                                                  tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                                "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
        geom_point_interactive(position=position_jitter(width=0.45,height=.45),size=.7,data=samplestoplotdata,
                               aes(x = population, y = proportion,colour=file,
                                   tooltip=paste0(round(proportion,digits=2),"%, ",file),data_id=file))+
        scale_y_continuous(name="Proportion % (parent population)",limits = c(0,100),breaks=seq(0,100,by=10),expand=c(0,0))+
        xlab(label = "Cell Population")+
        theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
              axis.text.y = element_text(hjust = 1, size=10,color=1),
              legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
      if (!shiny){
        r <- girafe(ggobj=r)
        r <- girafe_options(r, opts_toolbar(saveaspng = TRUE),
                            #specifies properties of outlier dots when cursor hovers over them
                            opts_hover(css = "stroke:black;stroke-width:5pt;r:3pt"))
      }
    }else{  #...and no user-specified samples to plot, no need to add geom_point_interactive.
      r <- ggplot() +
        geom_boxplot_interactive(data=alldata,aes(x = population, y = proportion, colour=population,
                                                  tooltip=paste(population,", median = ",round(as.numeric(median),digits=2),
                                                                "%, mean = ",round(as.numeric(mean),digits=2),"%")),outlier.shape=NA,na.rm=T)+
        scale_y_continuous(name="Proportion % (parent population)",limits = c(0,100),breaks=seq(0,100,by=10),expand=c(0,0))+
        xlab(label = "Cell Population")+
        theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10),
              axis.text.y = element_text(hjust = 1, size=10,color=1),
              legend.position="none",plot.margin=unit(c(0.8,0,0,0.85),"cm"))
      if (!shiny)
        r <- girafe(ggobj=r)
    }
  }
  return(r)
}



#---------- function of merhousnugh to generate the interactive boxplot ----------------------
boxplot.svg <- function(gs,sdconst=3,popstoplot,rows.to.avoid=NULL,row.name.no=2,percentage=F
                        ,showoutliers=T,samplestoplot=NULL,shiny=T,tooltip.id =c("file.name","file+pop")[1])
{
  # 1) If gs is a gatingSet, then popstoplot are nodes to be plotted from the gs
  # 2) If gs is a matrix of props or data.frame of props, then rows correspond to samples and columns to populations. 
  # And popstoplot corresponds to rownames of gs
  # 3) If gs is a vector of path to multiple csv file, then they will be read and merged together, so make sure they have same number of rows, and are consistent.
  # And popstoplot corresponds to rownames of each of these csv files
  # And rows.to.avoid, excludes the rows that are not related to samples, and should not be included in the boxplot
  # And row.name.no will be passed to read.csv to grab the rownames from the csv file.
  #Rebecca (190808): added samplestoplot <- vector of samples user would like to be plotted as individual points over boxplot,
  #and showoutliers <- if T, will plot outliers,
  #and percentage <- if True, assumes cell proportions are already in units of % and will not multiply by 100%. Set this to False if proportions are in decimal.
  
  if (class(gs)=="GatingSet")
  {
    sampleNames(gs) <- paste(unlist(lapply(sampleNames(gs), function(x) unlist(strsplit(x,".fcs"))[1])),".fcs",sep="")
    stats <- gs_pop_get_count_fast(gs,format="wide",statistic="freq",path="auto")
    tmp <- t(stats)
    tmp <- as.matrix(tmp[,popstoplot])
    colnames(tmp)<-popstoplot
  }else if (class(gs)=="matrix" |class(gs)=="data.frame")
  {
    tmp<- as.matrix(gs[,popstoplot])
    colnames(tmp)<- popstoplot
  }else if (class(gs)=="character" & length(grep(gs,pattern = ".csv"))==length(gs))
  {
    if (is.null(rows.to.avoid))
      tmp <- as.matrix(do.call(rbind,lapply(gs, function(path) read.csv(path, check.names=F,row.names=row.name.no)))[,popstoplot])
    else
      tmp <- as.matrix(do.call(rbind,lapply(gs, function(path) read.csv(path, check.names=F,row.names=row.name.no)[-rows.to.avoid,]))[,popstoplot])
    colnames(tmp)<- popstoplot
  }else{
    stop("Unsupported input as gs, check the documentation of boxplot.svg")
  }
  
  propmatrix <- stack(as.data.frame(tmp))
  propmatrix$file <- rownames(tmp)
  names(propmatrix) <- c("proportion","population","file")
  if(percentage){
    propmatrix$proportion <- round(propmatrix$proportion,digits = 5)
  }else{
    propmatrix$proportion <- round(propmatrix$proportion*100,digits = 5)
  }
  if (tooltip.id=="file+pop")
    propmatrix$file <- paste (propmatrix$file, propmatrix$population,sep="_p-")
  #remove empty rows
  propmatrix <-propmatrix[order(propmatrix$population,propmatrix$file,propmatrix$proportion),]
  delete <- which(propmatrix$population == "")
  if(length(delete) >0){
    propmatrix <- propmatrix[-delete,]
  }
  populationNumber <- length(unique(propmatrix$population))
  propmatrix$median <- ""
  propmatrix$mean <- ""
  
  outlierpropmatrix <- propmatrix[1,]
  outlierpropmatrix[1,] <- NA
  idx <- 1
  for(i in unique(propmatrix$population)){
    if(i != ""){
      indices <- which(propmatrix$population==i)
      mean <- mean(x=as.numeric(propmatrix$proportion[indices]))
      propmatrix$mean[indices] <- mean
      median <- median(as.numeric(propmatrix$proportion[indices]))
      propmatrix$median[indices] <- median
      sd <- sd(as.numeric(propmatrix$proportion[indices]))
      upplim <- mean + sd*sdconst
      lowlim <- mean - sd*sdconst
      
      outliers <- which(as.numeric(propmatrix$proportion[indices]) > upplim)
      outliers <- append(outliers,values=which(as.numeric(propmatrix$proportion[indices]) < lowlim))
      if(length(outliers) >0){
        for(k in outliers){
          outlierpropmatrix[idx,] <- propmatrix[indices[k],]
          idx <- idx+1      
        }
      }
    }
  }
  
  if(outlierpropmatrix$population[1] == "" || is.na(outlierpropmatrix$population[1])){
    outlierpropmatrix <- outlierpropmatrix[-which(outlierpropmatrix$population==""),]
  }
  
  if(!is.null(samplestoplot)){
    samplestoplotpropmatrix <- propmatrix[which(unlist(lapply(propmatrix$file, function(f) 
      unlist(strsplit(f,"_p-"))[1]))==samplestoplot[1]),]
    
    for(i in 2:length(samplestoplot)){
      samplestoplotpropmatrix <- rbind(samplestoplotpropmatrix,propmatrix[which(unlist(lapply(propmatrix$file, function(f) 
        unlist(strsplit(f,"_p-"))[1]))==samplestoplot[i]),])
    }
    samplestoplotpropmatrix <- unique(samplestoplotpropmatrix)
  }else{
    samplestoplotpropmatrix <- NULL
  }
  
  library(ggplot2)
  library(ggiraph)
  boxplot <- plot_results(alldata=propmatrix,showoutliers=showoutliers,outlierdata=outlierpropmatrix,
                          popstoplot=popstoplot,samplestoplotdata=samplestoplotpropmatrix,shiny = shiny)
  boxplot
}
# # function to make the df counts with pops as columns and row.names equals to the samples names
# reshape_df<-function(df){
#   df[name]
# }

# function to pre_process the df of freqs
pre_process_df_freqs<-function(df_freqs,remove_pops="None",renamed_pops){
  inds_col_remove<-grep("samples_name|Batch|Plate|name",colnames(df_freqs))
  samples_names<-df_freqs[,1]
  df_freqs<-df_freqs[,-inds_col_remove]
  df_freqs[, ] <- lapply(df_freqs[, ], as.numeric)
  # rename the columns
  colnames(df_freqs)<-renamed_pops
  # we remove the pops we don't want to plot
  inds_pops_to_remove<-grep(remove_pops,colnames(df_freqs))
  df_freqs<-df_freqs[,-inds_pops_to_remove]
  rownames(df_freqs)<-samples_names
  return(df_freqs)
}
############################################################ execute functions #######################################
#-------- list paths Bcells and myeloid files -----------------------
list_path_all_Bcells_results_folder<-list.files("/home/rstudio/data/Renamed_Bcells_all_batches_results_S3_system",full.names = T,recursive=T)
list_path_all_myeloid_results_folder<-list.files("/home/rstudio/data/Renamed_Myeloid_all_batches_results_s3",full.names = T,recursive=T)

# # import counts as data frame
# df_Bcells<-import_all_csv_files_same_type(list_path_all_Bcells_results_folder,type="counts",panel = "Bcells",ML_like = "None")
# df_myeloid<-import_all_csv_files_same_type(list_path_all_myeloid_results_folder,type="counts",panel = "myeloid",ML_like = "None")
# df_2<-pre_process_df(list_path_all_Bcells_results_folder,df)
# df_Bcells$Population[1:25]
# 
# 
# 
# # generate interactive boxplot with my code
# inds<-grep("Batch_1",df_2$Batch)
# df_2<-df_2[inds,]
# df_outliers<-build_df_outliers(df_2)
# interactive_boxplots_all_counts(df_2,df_outliers)
# # generate interactive boxplot with merhnoush code


# imports freqs as data frames
df_freqs_Bcells<-import_all_csv_files_same_type(list_path_all_Bcells_results_folder,type="freqs",panel = "Bcells",ML_like = "None")
df_freqs_myeloid<-import_all_csv_files_same_type(list_path_all_myeloid_results_folder,type="freqs",panel = "myeloid",ML_like = "None")
nrow(df_freqs_Bcells)
nrow(df_freqs_myeloid)
df_counts_Bcells<-import_all_csv_files_same_type(list_path_all_Bcells_results_folder,type="counts",panel = "Bcells",ML_like = "counts")
df_counts_Myeloid<-import_all_csv_files_same_type(list_path_all_myeloid_results_folder,type="counts",panel = "myeloid",ML_like = "counts")

nrow(df_counts_Bcells)
nrow(df_counts_Myeloid)
columns_name_myeloid<-c("All cells","B cells","Basophils","CD11b+CD16+ Mature Neutrophils","CD11b+CD16- Granulocytes","CD11b-CD16+ Immature Neutrophils 2",
                       "CD11-CD16- Immature Neutrophils 1","CD3+ T cells","CD3-","CD45+CD66- (Non Granulocytes)","CD56 Hi NK","CD56- dim CD16- NK","CD56- dim CD16+ NK",
                       "CD56+CD16+ NKT cells","CD56-CD16+ NK","CD56-CD16- cells","CD64+","CD64-","Classical Monocytes","Time","Granulocytes",
                       "HLADR+CD14+ Monocytes","HLADR+ CD14-","HLADR-CD14-","Live cells","Non classical Monocytes","Singlets","gd T cells","gd- T cells","mDC",
                       "pDC","root")

columns_name_Bcells<-c("All cells","Atypical B cells","Blasts","CD10+","CD10-","CD19+ B cells",
                        "CD19+CD20+","CD19+CD20-","Time","Granulocytes","IGD-IGM- B cells","IGD+IGM- B cells","IGD-IGM+ B cells",
                        "IGD+IGM+ B cells","IGM","Immature transition B cells","Live cells","Naive B cells","Non granulocytes","Plasmablasts",
                        "Singlets","Switched memory B cells-","Unswitched memory B cells","plasma cells","root")


df_freqs_Bcells<-pre_process_df_freqs(df_freqs_Bcells,remove_pops="Time|root|Singlets|IGM|All cells",renamed_pops = columns_name_Bcells)
df_freqs_myeloid<-pre_process_df_freqs(df_freqs_myeloid,remove_pops="Time|root|Singlets|All cells",renamed_pops = columns_name_myeloid)

df_counts_Bcells<-pre_process_df_freqs(df_counts_Bcells,remove_pops="Time|root|Singlets|IGM|All cells",renamed_pops = columns_name_Bcells)
df_counts_myeloid<-pre_process_df_freqs(df_counts_Myeloid,remove_pops="Time|root|Singlets|All cells",renamed_pops = columns_name_myeloid)


# generate interactive boxplot with merhnoush code
bp_bcells <-boxplot.svg(df_freqs_Bcells,sdconst=3,popstoplot=colnames(df_freqs_Bcells),rows.to.avoid=NULL,row.name.no=2,percentage=F
            ,showoutliers=T,samplestoplot=NULL,shiny=F)

bp_myeloid <-boxplot.svg(df_freqs_myeloid,sdconst=3,popstoplot=colnames(df_freqs_myeloid),rows.to.avoid=NULL,row.name.no=2,percentage=F
                 ,showoutliers=T,samplestoplot=NULL,shiny=F)

saveWidget(widget = bp_bcells, file="/home/rstudio/results/interactive_boxplots/boxplot_bcell_panel.html", selfcontained = TRUE, libdir = ,
           background = "white")
saveWidget(widget = bp_myeloid, file="/home/rstudio/results/interactive_boxplots/boxplot_myeloid_panel.html", selfcontained = TRUE, libdir = ,
           background = "white")

bp_bcells_counts <-boxplot.svg(df_counts_Bcells,sdconst=3,popstoplot=colnames(df_freqs_Bcells),rows.to.avoid=NULL,row.name.no=2,percentage=T
                        ,showoutliers=T,samplestoplot=NULL,shiny=F)

bp_myeloid_counts <-boxplot.svg(df_counts_myeloid,sdconst=3,popstoplot=colnames(df_freqs_myeloid),rows.to.avoid=NULL,row.name.no=2,percentage=T
                         ,showoutliers=T,samplestoplot=NULL,shiny=F)

saveWidget(widget = bp_bcells_counts, file="/home/rstudio/results/interactive_boxplots/boxplot_bcell_panel_counts.html", selfcontained = TRUE, libdir = ,
           background = "white")
saveWidget(widget = bp_myeloid_counts, file="/home/rstudio/results/interactive_boxplots/boxplot_myeloid_panel_counts.html", selfcontained = TRUE, libdir = ,
           background = "white")


