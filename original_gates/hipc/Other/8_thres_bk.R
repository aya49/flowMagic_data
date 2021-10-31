# load data and functions from sebastiano's code
load("~/code/flowTypeFIlter/test_data/HIPC environment.RData")
install.packages("microbenchmark")
library(microbenchmark)
library("flowType", lib.loc="/usr/local/lib/R/site-library")
library(flowCore)
library(stringr)
library(data.table)
library(flowTypeFilterC)
library(parallel)

# Change fs_Bcells or fs_myeloid everytime to execute the function on a different set of fcs file (fcs file of each plate)
#fs_Bcells<-import_fs("/home/rstudio/data/Tobias Kollmann/flowType_data/Bcell_panel/Batch_1/batch_1/fcs_files")
#fs_myeloid<-import_fs("/home/rstudio/data/Tobias Kollmann/flowType_data/Myeloid_panel/Batch_1/fcs_files")

# in my docker the path are differents,change again this in your code.
fs_Bcells<-import_fs("/home/rstudio/data/flowType_data/Bcell_panel/Batch_1/batch_1/fcs_files")
fs_myeloid<-import_fs("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_1/fcs_files")

# simplified excution of flowTypeFilter
execute_flowType_filter <- function(test_filter, fs, path_df_thresholds, path_output, panel, n_cores){
  if(panel=="Bcells"){
    # import df_trehsholds current plate
    df_thresholds<-import_thresholds(path_df_thresholds=path_df_thresholds)
    
    #---------- get flowType results--------------
    list_results<-mclapply(1:length(fs),function(i){
      f<-fs[[i]]
      name_current_sample<-identifier(f)
      print(name_current_sample)
      
      #---- get info current flowFrame
      df_f<-pData(parameters(f))
      MarkerNames<-as.vector(df_f$desc)
      inds<-which(is.na(MarkerNames)==T)
      MarkerNames[inds]<-df_f$name[inds]
      
      #---- get thresholds of interest
      ind<-which(df_thresholds$samples_name==name_current_sample)
      
      thr_non_gran_CD66<-as.numeric(strsplit(as.character(df_thresholds$Non.granulocytes[ind]),",")[[1]][1])
      thr_pc_CD38<-as.numeric(strsplit(as.character(df_thresholds$plasma.cells[ind]),",")[[1]][1])
      thr_pc_CD138<-as.numeric(strsplit(as.character(df_thresholds$plasma.cells[ind]),",")[[1]][2])
      thr_CD20<-as.numeric(strsplit(as.character(df_thresholds$CD19.CD20..1[ind]),",")[[1]][2])
      thr_IgD<-as.numeric(strsplit(as.character(df_thresholds$IGD.IGM..B.cells[ind]),",")[[1]][1])
      thr_IgM<-as.numeric(strsplit(as.character(df_thresholds$IGD.IGM..B.cells[ind]),",")[[1]][2])
      thr_CD27<-as.numeric(strsplit(as.character(df_thresholds$Plasmablasts[ind]),",")[[1]][2])
      thr_CD10<-as.numeric(strsplit(as.character(df_thresholds$CD10..1[ind]),",")[[1]][1])
      thr_CD27_2<-as.numeric(strsplit(as.character(df_thresholds$Naive.B.cells[ind]),",")[[1]][1])
      thr_IgD_2<-as.numeric(strsplit(as.character(df_thresholds$IGD.IGM..B.cells[ind]),",")[[1]][2])
      
      ind_CD66<-get_indices_marker(df_f,"CD66")
      ind_CD38<-get_indices_marker(df_f,"CD38")
      ind_CD138<-get_indices_marker(df_f,"CD138")
      ind_CD20<-get_indices_marker(df_f,"CD20")
      ind_IgD<-get_indices_marker(df_f,"IgD")
      ind_IgM<-get_indices_marker(df_f,"IgM")
      ind_CD27<-get_indices_marker(df_f,"CD27")
      ind_CD10<-get_indices_marker(df_f,"CD10")
      
      #---set arguments value for flowType filter (8 thresholds)
      list_thresholds<-list(thr_non_gran_CD66, thr_pc_CD38, thr_pc_CD138, thr_CD20, c(thr_IgD,thr_IgD_2), thr_IgM, c(thr_CD27,thr_CD27_2), thr_CD10)
      PropMarkers<<-c(ind_CD66, ind_CD38, ind_CD138, ind_CD20, ind_IgD, ind_IgM, ind_CD27, ind_CD10) # only indices markers with simple thresholds
      PartitionsPerMarker_vec<-rep(2,8)
      
      # execute flowType filter or flowType
      print("execution of flowTypefilter")
      if(test_filter == TRUE){
        print("execution of flowTypefilter")
        output_flowtype_filter<-flowTypeFilterC::flowTypeFilter(Frame=f,
                                                                PropMarkers=PropMarkers, 
                                                                Methods='Thresholds', # A bug here. you kept "Filters" option. However,if you remove all filters from the code, you cannot use the Filters option.
                                                                Thresholds=list_thresholds,
                                                                MarkerNames=MarkerNames,
                                                                verbose = F,cores = 10,  # only flowTypeFilter can use muliple cores to generate the populations
                                                                MaxMarkersPerPop=8, 
                                                                PartitionsPerMarker = PartitionsPerMarker_vec) 
        
        
      }else if(test_filter == F){
        # I removed the cores option from here,that was another bug. You can use parallelizzation (mutiple cores) with flowType 
        # only with flowTypeFilter,flowType has no parallelization option
        # Note: there are two parallezzizations here,one at the level of samples (how many samples to analyze in parallel)
        # and one at the level of populations generation (how many populations to generate in parallel in one single sample)
        # flowType cannnot parallezize the populations generation for each single sample (second level of parallelizzaion).
        output_flowtype_filter<-flowType::flowType(Frame=f,
                                                                PropMarkers=PropMarkers, 
                                                                Methods='Thresholds',
                                                                Thresholds=list_thresholds,
                                                                MarkerNames=MarkerNames,
                                                                verbose = F,
                                                                MaxMarkersPerPop=8, # why you changed this? the lowest pops have 6  markers
                                                                PartitionsPerMarker = PartitionsPerMarker_vec) 
      }

      # the deepest pop has 8 markers
      
      print("decodification of the phenotypes codes")
      phenotype.names<<-unlist(lapply(output_flowtype_filter@PhenoCodes,function(x){return(flowTypeFilterC::decodePhenotype(x,c("CD66","CD38","CD138","CD20","IgD","IgM","CD27","CD10"),output_flowtype_filter@PartitionsPerMarker,c(rep(FALSE,8))))}))
      phenotype.names[1]<-"root_live"
      cell_counts<-output_flowtype_filter@CellFreqs
      Pheno_codes<-output_flowtype_filter@PhenoCodes
      return(list(phenotype_names=phenotype.names,cell_counts=cell_counts,sample_name=name_current_sample,Pheno_codes=Pheno_codes))
    },mc.cores = n_cores)
    
  }else if(panel=="myeloid"){
    # import df_trehsholds current plate
    df_thresholds<-import_thresholds(path_df_thresholds=path_df_thresholds)
    
    #---------- get flowType results--------------
    list_results<-mclapply(1:length(fs), function(i){
      f<-fs[[i]]
      name_current_sample<-identifier(f)
      print(name_current_sample)
      #---- get info current flowFrame
      df_f<-pData(parameters(f))
      MarkerNames<-as.vector(df_f$desc)
      inds<-which(is.na(MarkerNames)==T)
      MarkerNames[inds]<-df_f$name[inds]
      
      # #---- get thresholds of interest
      ind<-which(df_thresholds$samples_name==name_current_sample)
      
      thr_CD64<-as.numeric(df_thresholds$CD64.[ind])
      thr_mature_neutrophils_CD11B<-as.numeric(strsplit(as.character(df_thresholds$CD11b.CD16..Mature.Neutrophils[ind]),",")[[1]][1])
      thr_mature_neutrophils_CD16<-as.numeric(strsplit(as.character(df_thresholds$CD11b.CD16..Mature.Neutrophils[ind]),",")[[1]][2])
      thr_CD3<-as.numeric(df_thresholds$CD3.[ind])
      thr_gd<-as.numeric(strsplit(as.character(df_thresholds$gd.T.cells[ind]),",")[[1]][2])
      thr_CD16_NKT<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..NKT.cells[ind]),",")[[1]][1])
      thr_CD56NKT<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..NKT.cells[ind]),",")[[1]][2])
      thr_CD16_NK<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..cells[ind]),",")[[1]][1])
      thr_CD56_NK<-as.numeric(strsplit(as.character(df_thresholds$CD56.CD16..cells[ind]),",")[[1]][2])
      thr_CD123<-as.numeric(strsplit(as.character(df_thresholds$B.cells[ind]),",")[[1]][1])
      thr_CD11c<-as.numeric(strsplit(as.character(df_thresholds$B.cells[ind]),",")[[1]][2])
      thr_CD16_monocytes<-as.numeric(strsplit(as.character(df_thresholds$Classical.Monocytes[ind]),",")[[1]][2])
      
      
      ind_CD64<-get_indices_marker(df_f,"CD64",channel_name = "Alexa Fluor 700-A")
      ind_CD11B<-get_indices_marker(df_f,"CD11b",channel_name = "BV786-A")
      ind_CD16<-get_indices_marker(df_f,"CD16",channel_name = "FITC-A")
      if(length(ind_CD16)>1){
        ind_CD16<-ind_CD16[2]
      }
      ind_CD3<-get_indices_marker(df_f,"CD3",channel_name = "PE-CF594-A")
      ind_gd<-get_indices_marker(df_f,"gd",channel_name = "PE-A")
      if(length(ind_gd)==0){
        ind_gd<-get_indices_marker(df_f,"gd TCR")
      }
      ind_CD56<-get_indices_marker(df_f,"CD56",channel_name = "BV650-A")
      ind_CD123<-get_indices_marker(df_f,"CD123",channel_name = "PE-Cy7-A")
      ind_CD11c<-get_indices_marker(df_f,"CD11c",channel_name = "APC-A")
    
      #---set arguments value for flowType filter
      list_thresholds<-list(thr_CD64,thr_mature_neutrophils_CD11B,c(thr_mature_neutrophils_CD16,thr_CD16_NKT,thr_CD16_NK,thr_CD16_monocytes),thr_CD3,thr_gd,c(thr_CD56NKT,thr_CD56_NK),thr_CD123,thr_CD11c)
      PropMarkers<-c(ind_CD64,ind_CD11B,ind_CD16,ind_CD3,ind_gd,ind_CD56,ind_CD123,ind_CD11c) # only indices markers with simple thresholds
      PartitionsPerMarker_vec<-rep(2, 8)
      
      # execute flowType or flowTypeFilter
      # I fixed the same things here
      if (test_filter == TRUE){
        print("execution of flowTypefilter")
        output_flowtype <- flowTypeFilterC::flowTypeFilter(Frame=f,
                                                                PropMarkers=PropMarkers,
                                                                Methods="Thresholds",
                                                                Thresholds=list_thresholds,
                                                                MarkerNames=MarkerNames,
                                                                verbose = T,
                                                                cores = 10,
                                                                MaxMarkersPerPop=8,
                                                                PartitionsPerMarker = PartitionsPerMarker_vec) 
        # the deepest pop has 8 markers
      }else{
        print("execution of flowType")
        output_flowtype <- flowType::flowType(Frame=f,
                                                                PropMarkers=PropMarkers,
                                                                Methods="Thresholds",
                                                                Thresholds=list_thresholds,
                                                                MarkerNames=MarkerNames,
                                                                verbose = T,
                                                                MaxMarkersPerPop=8,
                                                                PartitionsPerMarker = PartitionsPerMarker_vec) 
        
        # the deepest pop has 8 markers
      }
      
      
      
      # the output_flowtype_filter variable should contain everything. The counts,the phenocodes,the thresholds and filters.
      print("decodification of the phenotypes codes")
      phenotype.names<-unlist(lapply(output_flowtype_filter@PhenoCodes,function(x){return(flowTypeFilterC::decodePhenotype(x,c("CD64","CD11B","CD16","CD3","gd","CD56","CD123","CD11c"),output_flowtype_filter@PartitionsPerMarker,c(rep(FALSE,8))))}))
      phenotype.names[1]<-"root_live"
      cell_counts<-output_flowtype_filter@CellFreqs
      Pheno_codes<-output_flowtype_filter@PhenoCodes
      return(list(phenotype_names=phenotype.names,cell_counts=cell_counts,sample_name=name_current_sample,Pheno_codes=Pheno_codes))
    },mc.cores = n_cores)
    
  }
  
  vec_samples_name<<-sapply(1:length(list_results),function(i){
    sample_name_i<-list_results[[i]]$sample_name
    return(sample_name_i)
  })
  names(list_results)<-vec_samples_name
  #---- get final df ---------
  df_flowType_results<-get_df_from_results(list_results)
  df_flowType_results[, ] <- lapply(df_flowType_results[, ], as.character)
  
  return(df_flowType_results)
}

# excute flowType

# change the path to excecute benchmarking on diffrent batch
# root of my file dic the root /home/rstudio/data/Tobias Kollmann/flowType_data/...

# Bcells batch 1, 35 samples

#path_df_thresholds ="/home/rstudio/data/Tobias Kollmann/flowType_data/Bcell_panel/Batch_1/batch_1/df_thresholds" 
path_df_thresholds ="/home/rstudio/data/flowType_data/Bcell_panel/Batch_1/batch_1/df_thresholds" 

# I tested also this one and it works 
mbm_b1 <- microbenchmark("flowTypeFilter" = {execute_flowType_filter(test_filter = TRUE, 
                                                                     fs = fs_Bcells,
                                                                  path_df_thresholds,
                                                                  panel = "Bcells",
                                                                  n_cores=10)}, 
                      "flowType" = {execute_flowType_filter(test_filter = FALSE,
                                                            fs = fs_Bcells,
                                                            path_df_thresholds,
                                                            panel = "Bcells",
                                                            n_cores=10)}, 
                      times = 10L)

# Sebastiano test
start<-Sys.time()
execute_flowType_filter(test_filter = T, 
                        fs = fs_Bcells,
                        path_df_thresholds,
                        panel = "Bcells",
                        n_cores=1)
end<-Sys.time()
time_elapsed<-end-start
print(time_elapsed)
start<-Sys.time()
execute_flowType_filter(test_filter = F, 
                        fs = fs_Bcells,
                        path_df_thresholds,
                        panel = "Bcells",
                        n_cores=1)
end<-Sys.time()
time_elapsed<-end-start
print(time_elapsed)