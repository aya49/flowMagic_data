import_fs<-function(path_fcs_files){
  paths<-list.files(path_fcs_files,recursive = T,full.names = T,pattern = ".fcs")
  names_samples<-list.files(path_fcs_files,recursive = T,full.names = F,pattern = ".fcs")
  frames<-mclapply(paths, read.FCS,mc.cores = 8)
  fs<-as(frames, "flowSet")
  sampleNames(fs)<-names_samples
  return(fs)
}


import_thresholds<-function(path_df_thresholds){
  paths<-list.files(path_df_thresholds,recursive = T,full.names = T)
  df<-read.csv(paths,header = T)
  colnames(df)[1]<-"samples_name"
  return(df)
}

import_filters<-function(path_filter_objects,pop){
  paths<-list.files(path_filter_objects,recursive = T,full.names = T)
  ind<-grep(pop,paths)
  list_filter<-readRDS(paths[ind])
  return(list_filter)
}

get_indices_marker<-function(df_f,marker_name,channel_name="None",mode=1){
  ind<-which(df_f$desc==marker_name)
  if((length(ind)==0)||(length(ind)>1)){
    ind<-which(df_f$name==channel_name)
  }
  if(mode==2){
    ind<-which(df_f$desc==marker_name)
    if(length(ind)==0){
      ind<-which(df_f$name==marker_name)
    }
  }
  return(as.integer(ind))
}

get_markers_filter_df<-function(df_f,df_filters){
  markers<-colnames(df_filters)
  ind_marker_1<-get_indices_marker(df_f,markers[1],mode = 2)
  ind_marker_2<-get_indices_marker(df_f,markers[2],mode = 2)
  marker_1<-as.character(df_f$desc[ind_marker_1])
  marker_2<-as.character(df_f$desc[ind_marker_2])
  if(is.na(marker_1)==T){
    marker_1<-as.character(df_f$name[ind_marker_1])
  }
  if(is.na(marker_2)==T){
    marker_2<-as.character(df_f$name[ind_marker_2])
  }
  marker_1<-str_remove(marker_1,"-")
  marker_2<-str_remove(marker_2,"-")
  filter_markers<-paste0(marker_2,marker_1)
  return(filter_markers)
}

get_df_from_results<-function(list_results){
  list_df<-list()
  samples_names<-names(list_results)
  for(i in 1:length(list_results)){
    name_current_sample<-samples_names[i]
    df_i<-as.data.frame(list_results[[i]][c(1,2,4)]) # for each sample, we save counts,pheno_names and pheno_codes
    samples_name<-rep(name_current_sample,nrow(df_i))
    df_i<-cbind(samples_name,df_i)
    list_df[[i]]=df_i
  }
  df_flowType_results<-as.data.frame(rbindlist(list_df))
  return(df_flowType_results)
}



# function to combine all the flowtype results 
combine_all_flowtype_dfs<-function(paths_all_flowType_results,
                                   type_data="great_LOQ",type_visit="V2",n_cores=12,
                                   path_samples_info){
  list_all_df<-list()
  # These are the subjects IDs with the best quality regarding the antibody concentration (reported on the server)
  subject_IDs_great_LOQ<-"G014K|G016D|G018C|G026E|G028G|G033F|G053B|G063J|G067C|G068A|G073K|G078D|G089J|G093G|G094J|G102K|G113F|G118B|G124C|G134F|G144E|G150K|G152D|G154G|G159F|G174A|G175E|G181C|G183G|G200J|G208B|G216J|G233G|G241A|G260G|G264A|G290B|G306H|G308A|G311F|G312J|G346D|G366B|G388G|G404K|G411J|G425G|G433D|G453H|G458E|G468K|G483J|G492K|G499D|G501C|G502A|G511G|G521E|G530H|G532K|G533C|G537J|G542C|G556D|G557E|G567B|G572E|G575G|G576C|G582J|G587H|G616F|G617K|G625J|G634G|G650H|G655F|G659D|G667A|G685K|G688E|G691G|G700D|G703E|G704C|G713J"
  info_df<-read.csv(file = path_samples_info)
  inds<-grep("V2|Visit 2",info_df$Visit.Num)
  visitID_v2<-as.character(info_df$Visit.ID[inds])
  string_visit_v2<-paste0(visitID_v2,collapse = "|")
  list_all_df<-mclapply(1:length(paths_all_flowType_results),function(i){
    current_df<-read.csv(paths_all_flowType_results[i],colClasses = "character")
    if(type_data=="great_LOQ"){
      inds<-grep(subject_IDs_great_LOQ,current_df$samples_name)
      if(length(inds)!=0){
        current_df<-current_df[inds,]
      }else{
        return(i)
      }
    }
    if(type_visit=="V2"){
      inds<-grep(string_visit_v2,current_df$samples_name)
      if(length(inds)!=0){
        current_df<-current_df[inds,]
      }else{
        return(i)
      }
    }
    ind<-grep("Batch_2/Plate_5/",paths_all_flowType_results[i])
    if(length(ind)!=0){
      ind<-grep("FCT_G216J_4L42.fcs",current_df[,1])
      current_df<-current_df[ind,]
    }
    return(current_df)
  },mc.cores = n_cores )
  vec_classes<-sapply(list_all_df, class)
  inds<-which(vec_classes!="data.frame")
  if(length(inds)!=0){
    list_all_df<-list_all_df[-inds]
  }
  df<-as.data.frame(rbindlist(list_all_df))
  return(df)
}
