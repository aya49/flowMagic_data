execute_flowType_filter<-function(fs,path_df_thresholds,path_filter_objects,path_output,panel,n_cores=1){
  if(panel=="Bcells"){
    # import df_trehsholds current plate
    df_thresholds<-import_thresholds(path_df_thresholds=path_df_thresholds)
    # import filters objects current plate
    list_filters_bcells<-import_filters(path_filter_objects=path_filter_objects,"Bcells")
    list_filters_blast<-import_filters(path_filter_objects=path_filter_objects,"blast")
    list_filters_trans<-import_filters(path_filter_objects=path_filter_objects,"trans")
    
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
      #----- get filters df
      ind<-which(names(list_filters_bcells)==name_current_sample)
      df_filters_bcells<-list_filters_bcells[[ind]]
      df_filters_blast<-list_filters_blast[[ind]]
      df_filters_trans<-list_filters_trans[[ind]]
      # print(df_filters_bcells)
      # print(df_filters_blast)
      # print(df_filters_trans)
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
      
      # print(thr_non_gran_CD66)
      # print(thr_pc_CD138)
      # print(thr_CD20)
      # print(thr_IgD)
      # print(thr_IgM)
      # print(thr_CD27)
      # print(thr_CD10)
      # print(thr_CD27_2)
      # print(thr_IgD_2)
      
      # print(ind_CD66)
      # print(ind_CD38)
      # print(ind_CD20)
      # print(ind_CD138)
      # print(ind_IgD)
      # print(ind_IgM)
      # print(ind_CD27)
      # print(ind_CD10)
      #---set arguments value for flowType filter
      list_thresholds<-list(thr_non_gran_CD66,thr_pc_CD38,thr_pc_CD138,thr_CD20,c(thr_IgD,thr_IgD_2),thr_IgM,c(thr_CD27,thr_CD27_2),thr_CD10,df_filters_bcells,df_filters_blast,df_filters_trans)
      PropMarkers<<-c(ind_CD66,ind_CD38,ind_CD138,ind_CD20,ind_IgD,ind_IgM,ind_CD27,ind_CD10) # only indices markers with simple thresholds
      filter_markers_bcells<-get_markers_filter_df(df_f,df_filters_bcells)
      filter_markers_blast<-get_markers_filter_df(df_f,df_filters_blast)
      filter_markers_trans<-get_markers_filter_df(df_f,df_filters_trans)
      PartitionsPerMarker_vec<-rep(2,11)
      # execute flowType filter
      print("execution of flowTypefilter")
      output_flowtype_filter<-flowTypeFilterC::flowTypeFilter(Frame=f,PropMarkers=PropMarkers,Methods='Filters',Thresholds=list_thresholds,MarkerNames=MarkerNames,verbose = F,cores = 10,MaxMarkersPerPop=6,PartitionsPerMarker = PartitionsPerMarker_vec) # the deepest pop has 6 markers
      print("decodification of the phenotypes codes")
      phenotype.names<<-unlist(lapply(output_flowtype_filter@PhenoCodes,function(x){return(flowTypeFilterC::decodePhenotype(x,c("CD66","CD38","CD138","CD20","IgD","IgM","CD27","CD10",filter_markers_bcells,filter_markers_blast,filter_markers_trans),output_flowtype_filter@PartitionsPerMarker,c(rep(FALSE,8),T,T,T)))}))
      phenotype.names[1]<-"root_live"
      cell_counts<-output_flowtype_filter@CellFreqs
      Pheno_codes<-output_flowtype_filter@PhenoCodes
      return(list(phenotype_names=phenotype.names,cell_counts=cell_counts,sample_name=name_current_sample,Pheno_codes=Pheno_codes))
    },mc.cores = n_cores)
  }else if(panel=="myeloid"){
    # import df_trehsholds current plate
    df_thresholds<-import_thresholds(path_df_thresholds=path_df_thresholds)
    # import filters objects current plate
    list_filters_basophils<-import_filters(path_filter_objects=path_filter_objects,"basophils")
    list_filters_CD66negCD45pos<-import_filters(path_filter_objects=path_filter_objects,"CD66negCD45pos")
    list_filters_DRn_14n<-import_filters(path_filter_objects=path_filter_objects,"DRn_14n")
    list_filters_DRp_14n<-import_filters(path_filter_objects=path_filter_objects,"DRp_14n")
    list_filters_granulocytes<-import_filters(path_filter_objects=path_filter_objects,"granulocytes")
    list_filters_monocytes<-import_filters(path_filter_objects=path_filter_objects,"monocytes")
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
      #----- get filters df
      ind<-which(names(list_filters_basophils)==name_current_sample)
      df_filters_basophils<-list_filters_basophils[[ind]]
      df_filters_CD66negCD45pos<-list_filters_CD66negCD45pos[[ind]]
      df_filters_DRn_14n<-list_filters_DRn_14n[[ind]]
      df_filters_DRp_14n<-list_filters_DRp_14n[[ind]]
      df_filters_granulocytes<-list_filters_granulocytes[[ind]]
      df_filters_monocytes<-list_filters_monocytes[[ind]]
      
      # print(df_filters_basophils)
      # print(df_filters_CD66negCD45pos)
      # print(df_filters_DRn_14n)
      # print(df_filters_DRp_14n)
      # print(df_filters_granulocytes)
      # print(df_filters_monocytes)
      
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
      
      # print(thr_CD64)
      # print(thr_mature_neutrophils_CD11B)
      # print(thr_mature_neutrophils_CD16)
      # print(thr_CD3)
      # print(thr_gd)
      # print(thr_CD16_NKT)
      # print(thr_CD56NKT)
      # print(thr_CD16_NK)
      # print(thr_CD56_NK)
      # print(thr_CD123)
      # print(thr_CD11c)
      # print(thr_CD16_monocytes)
      
      print(ind_CD64)
      print(ind_CD11B)
      print(ind_CD16)
      print(ind_CD3)
      print(ind_gd)
      print(ind_CD56)
      print(ind_CD123)
      print(ind_CD11c)
      #---set arguments value for flowType filter
      list_thresholds<-list(thr_CD64,thr_mature_neutrophils_CD11B,c(thr_mature_neutrophils_CD16,thr_CD16_NKT,thr_CD16_NK,thr_CD16_monocytes),thr_CD3,thr_gd,c(thr_CD56NKT,thr_CD56_NK),thr_CD123,thr_CD11c,
                            df_filters_basophils,df_filters_CD66negCD45pos,df_filters_DRn_14n,df_filters_DRp_14n,df_filters_granulocytes,df_filters_monocytes)
      PropMarkers<-c(ind_CD64,ind_CD11B,ind_CD16,ind_CD3,ind_gd,ind_CD56,ind_CD123,ind_CD11c) # only indices markers with simple thresholds
      PartitionsPerMarker_vec<-rep(2,14)
      
      # execute flowType filter
      print("execution of flowTypefilter")
      output_flowtype_filter<-flowTypeFilterC::flowTypeFilter(Frame=f,PropMarkers=PropMarkers,Methods='Filters',Thresholds=list_thresholds,MarkerNames=MarkerNames,verbose = T,cores = 10,MaxMarkersPerPop=5,PartitionsPerMarker = PartitionsPerMarker_vec) # the deepest pop has 5 markers
      # the output_flowtype_filter variable should contain everything. The counts,the phenocodes,the thresholds and filters.
      print("decodification of the phenotypes codes")
      phenotype.names<-unlist(lapply(output_flowtype_filter@PhenoCodes,function(x){return(flowTypeFilterC::decodePhenotype(x,c("CD64","CD11B","CD16","CD3","gd","CD56","CD123","CD11c","basophils","NonGran","HLDR-CD14-","HLDR+CD14-","granulocytes","monocytes"),output_flowtype_filter@PartitionsPerMarker,c(rep(FALSE,8),T,T,T,T,T,T)))}))
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
  write.csv(df_flowType_results,path_output,row.names = F)
  return(df_flowType_results)
}

