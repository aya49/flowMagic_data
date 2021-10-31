################################## define functions #####################################################
renaming_files<-function(path_files,path_names,exclude_files="None",exclude_names="None",type_of_data){
  # --- import data ----------
  list_files_to_modify<-list.files(path_files)
  new_naming_info<-read.csv(file=path_names, header = TRUE)
  #------ pre-processing--------
  if(exclude_files!="None"){
    inds<-grep(exclude_files,list_files_to_modify)
    if(length(inds!=0)){
      list_files_to_modify<-list_files_to_modify[-inds]
    }
  }
  if(exclude_names!="None"){
    old_names<-new_naming_info[,1]
    inds<-grep(exclude_names,old_names)
    if(length(inds)!=0){
      new_naming_info<-new_naming_info[-inds,]
    }
  }
  visit_id<-as.character(new_naming_info[,3])
  subject_id<-as.character(new_naming_info[,2])
  old_names<-as.character(new_naming_info[,1])
  print(paste0("Number_samples_new_names:",length(visit_id)))
  print(paste0("Number_samples_to_modify:",length(list_files_to_modify)))
  print("############## list files to modify ###################")
  print(list_files_to_modify)
  print("############## list old names reported on the map ###################")
  print(old_names)
  #--- check duplicates ------
  for(id in visit_id){
    ind<-grep(id,visit_id)
    if(length(ind)>1){
      print("there are duplications of visit id")
      print(id)
      visit_id[ind[2]]<-paste0(visit_id[ind[2]],"_dup")
    }
  }
  #-------- renaming---------
  setwd(path_files)
  if(type_of_data!="Stats"){
    v<-sapply(1:length(old_names),function(i){
      str<-strsplit(old_names[i],".fcs")[[1]][1]
      ind_1<-grep(str,list_files_to_modify)
      if(length(ind_1)>1){
        str<-paste0(str,"_")
        ind_1<-grep(str,list_files_to_modify)
      }
      # print(i)
      # print(str)
      # print(ind_1)
      # print(list_files_to_modify[ind_1])
      # print(sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i]))
      #----- check duplicates ------
      
      if(length(ind_1)!=0){
        if((type_of_data=="Cleaning")||(type_of_data=="Plot")){
          #print(length(list_files_to_modify[ind_1]))
          file.rename(sprintf("%s/%s",path_files,list_files_to_modify[ind_1]),sprintf("FCT_%s_%s.png",subject_id[i],visit_id[i]))
        }else if(type_of_data=="WSP"){
          file.rename(sprintf("%s/%s",path_files,list_files_to_modify[ind_1]),sprintf("FCT_%s_%s.wsp",subject_id[i],visit_id[i]))
        }else if(type_of_data=="raw_files"){
          file.rename(sprintf("%s/%s",path_files,list_files_to_modify[ind_1]),sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i]))
        }
      }
      if(type_of_data=="WSP"){
        str_2<-strsplit(str,"_")[[1]][3]
        str_2<-paste0(str_2,"_AutomatedWSP")
        ind_2<-grep(str_2,list_files_to_modify)
        print(length(list_files_to_modify[ind_2]))
        if(length(ind_2)!=0){
          file.rename(sprintf("%s/%s",path_files,list_files_to_modify[ind_2]),sprintf("FCT_%s_%s_AutomatedWSP.wsp",subject_id[i],visit_id[i]))
        }
      }
      
    })
    #-------- print final results---------
    list_files_modified<-list.files(path_files)
    print("############## list files modified ###################")
    print(list_files_modified)
  }else{
    # import files
    counts_csv<-read.csv(file=paste0(path_files,"/final_counts_file.csv"), header = TRUE)
    counts_csv[,1]<-as.character(counts_csv[,1])
    freqs_csv<-read.csv(file=paste0(path_files,"/final_freqs_file.csv"), header = TRUE)
    freqs_csv[,1]<-as.character(freqs_csv[,1])
    description_csv<-read.csv(file=paste0(path_files,"/final_description_file.csv"), header = TRUE)
    description_csv[,1]<-as.character(description_csv[,1])
    for (i in 1:length(old_names)){
      str<-strsplit(old_names[i],".fcs")[[1]][1]
      ind_3<-grep(str,freqs_csv[,1])
      print(str)
      if(length(ind_3)>1){
        str<-paste0(str,"_")
        ind_3<-grep(str,freqs_csv[,1])
        if(length(ind_3)!=0){
          print(subject_id[i])
          freqs_csv[ind_3,1]<-sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i])
        }
        ind_1<-grep(str,counts_csv[,1])
        if(length(ind_1)!=0){
          counts_csv[ind_1,1]<-sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i])
        }
        ind_2<-grep(str,description_csv[,1])
        if(length(ind_2)!=0){
          description_csv[ind_2,1]<-sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i])
        }
      }else{
        if(length(ind_3)!=0){
          freqs_csv[ind_3,1]<-sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i])
        }
        ind_1<-grep(str,counts_csv[,1])
        if(length(ind_1)!=0){
          counts_csv[ind_1,1]<-sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i])
        }
        ind_2<-grep(str,description_csv[,1])
        if(length(ind_2)!=0){
          description_csv[ind_2,1]<-sprintf("FCT_%s_%s.fcs",subject_id[i],visit_id[i])
        }
      }




    }
    write.csv(counts_csv,file = sprintf("%s/final_counts_file.csv",path_files),row.names = F)
    write.csv(description_csv,file = sprintf("%s/final_description_file.csv",path_files),row.names = F)
    write.csv(freqs_csv,file = sprintf("%s/final_freqs_file.csv",path_files),row.names = F)
    print("######################## final renamed files #############################")
    print(description_csv)
  }


}

################################################ Execute renaming Bcells #########################
##################################################################################################

#----------- Validation data-------------
renaming_files(path_files="/home/rstudio/data/Renamed_Validation_data_Gambia_PNG/Flow_cytometry_data/Bcells_panel/Plate_2",
               path_names="/home/rstudio/data/Renamed_Validation_data_Gambia_PNG/Plate_maps/Plate_2_map.csv",
               type_of_data = "raw_files")

#------------------- Batch 1 ----------------------

# Cleaning 
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_1/Batch_1_files/Plate_1_map.csv",
               exclude_files = "Tube",
               type_of_data = "Cleaning")

# Cleaning replaced samples
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Gambia_Bcells_final/Batch_1/Batch_1_replaced_6_samples/Plate_1_replaced_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_1/Batch_1_files/Plate_1_map.csv",
               exclude_files = "Tube",
               type_of_data = "WSP")

# WSP replaced samples
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Gambia_Bcells_final/Batch_1/Batch_1_replaced_6_samples/Plate_1_replaced_map.csv",
               type_of_data = "WSP")

# Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_1/Batch_1_files/Plate_1_map.csv",
               exclude_files = "Tube",
               type_of_data = "Plot")

# Plot replaced samples
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Gambia_Bcells_final/Batch_1/Batch_1_replaced_6_samples/Plate_1_replaced_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_1/Batch_1_files/Plate_1_map.csv",
               exclude_files = "Tube",
               type_of_data = "Stats")
# Stats replaced samples

renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_1/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Gambia_Bcells_final/Batch_1/Batch_1_replaced_6_samples/Plate_1_replaced_map.csv",
               type_of_data = "Stats")

# Raw files

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_1/Batch_1_files/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_1/Batch_1_files/Plate_1_map.csv",
               exclude_files = "Tube",
               exclude_names="A8",
               type_of_data = "raw_files")

# Raw files replaced samples
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_1/Batch_1_replaced_6_samples/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Gambia_Bcells_final/Batch_1/Batch_1_replaced_6_samples/Plate_1_replaced_map.csv",
               type_of_data = "raw_files")

#-------------------------- Batch 2 --------------------------
#---- Plate 1
# Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_1/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_1/Plate_1_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_1/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_1/Plate_1_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_1/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_1/Plate_1_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_1/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_1/Plate_1_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_1/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_1/Plate_1_map.csv",
               type_of_data = "Stats")
#--- Plate 2
# Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_2/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_2/Plate_2_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_2/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_2/Plate_2_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_2/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_2/Plate_2_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_2/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_2/Plate_2_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_2/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_2/Plate_2_map.csv",
               type_of_data = "Stats")

#---- Plate 3
# Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_3/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_3/Plate_3_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_3/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_3/Plate_3_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_3/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_3/Plate_3_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_3/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_3/Plate_3_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_3/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_3/Plate_3_map.csv",
               type_of_data = "Stats")
#---- Plate 4
# Raw files 
renaming_files(path_files="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_4/FCS_files",
               path_names="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_4/Plate_4_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/rstudio/data/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_4/Cleaning",
               path_names="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_4/Plate_4_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/rstudio/data/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_4/WSP",
               path_names="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_4/Plate_4_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/rstudio/data/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_4/Plots",
               path_names="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_4/Plate_4_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/rstudio/data/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_4/Stats",
               path_names="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_4/Plate_4_map.csv",
               type_of_data = "Stats")
#----Plate 5
# Raw files 
renaming_files(path_files="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_5/FCS_files",
               path_names="/home/rstudio/data/Renamed_Gambia_Bcells_final/Batch_2/Plate_5/Plate_5_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_5/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_5/Plate_5_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_5/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_5/Plate_5_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_5/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_5/Plate_5_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_5/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_5/Plate_5_map.csv",
               type_of_data = "Stats")
#----Plate 6
#Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_6/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_6/Plate_6_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_6/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_6/Plate_6_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_6/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_6/Plate_6_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_6/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_6/Plate_6_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_2/Plate_6/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_2/Plate_6/Plate_6_map.csv",
               type_of_data = "Stats")
#---------------------------- Batch 3 ------------------------
#---- Plate 1
# Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_1/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_1/Plate_1_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_1/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_1/Plate_1_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_1/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_1/Plate_1_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_1/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_1/Plate_1_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_1/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_1/Plate_1_map.csv",
               type_of_data = "Stats")
#------ Plate 2
# Raw files P
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_2/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_2/Plate_2_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_2/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_2/Plate_2_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_2/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_2/Plate_2_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_2/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_2/Plate_2_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_2/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_2/Plate_2_map.csv",
               type_of_data = "Stats")

#-----Plate 3
# Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_3/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_3/Plate_3_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_3/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_3/Plate_3_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_3/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_3/Plate_3_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_3/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_3/Plate_3_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_3/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_3/Plate_3_map.csv",
               type_of_data = "Stats")
#-----Plate 4
# Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_4/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_4/Plate_4_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_4/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_4/Plate_4_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_4/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_4/Plate_4_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_4/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_4/Plate_4_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_4/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_4/Plate_4_map.csv",
               type_of_data = "Stats")
#----Plate 5
#Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/FCS_files/plate_5",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_5_map.csv",
               type_of_data = "raw_files")

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/FCS_files/plate_6",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_6_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_5/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_5_map.csv",
               type_of_data = "Cleaning")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_6/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_6_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_5/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_5_map.csv",
               type_of_data = "WSP")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_6/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_6_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_5/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_5_map.csv",
               type_of_data = "Plot")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_6/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_6_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_5/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_5_map.csv",
               type_of_data = "Stats")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_3/Plate_5/plate_6/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_3/Plate_5/Plate_5_plate_6_map.csv",
               type_of_data = "Stats")


#----------------------------- Batch 4 -----------------------
#---Plate 1
#Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_1/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_1/Plate_1_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_1/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_1/Plate_1_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_1/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_1/Plate_1_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_1/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_1/Plate_1_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_1/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_1/Plate_1_map.csv",
               type_of_data = "Stats")
#----- Plate 2
#Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_2/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_2/Plate_2_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_2/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_2/Plate_2_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_2/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_2/Plate_2_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_2/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_2/Plate_2_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_2/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_2/Plate_2_map.csv",
               type_of_data = "Stats")

#----Plate 3
#  Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_3/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_3/Plate_3_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_3/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_3/Plate_3_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_3/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_3/Plate_3_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_3/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_3/Plate_3_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_3/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_3/Plate_3_map.csv",
               type_of_data = "Stats")
#-----Plate 4
# Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_4/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_4/Plate_4_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_4/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_4/Plate_4_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_4/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_4/Plate_4_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_4/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_4/Plate_4_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_4/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_4/Plate_4_map.csv",
               type_of_data = "Stats")
#---- Plate 5
#Raw files 
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_5/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_5/Plate_5_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_5/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_5/Plate_5_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_5/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_5/Plate_5_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_5/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_5/Plate_5_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Bcells_all_batches_results_S3_system/Batch_4/Plate_5/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_Bcells_final/Batch_4/Plate_5/Plate_5_map.csv",
               type_of_data = "Stats")


################################################ Execute renaming Myeloid panel #########################
########################################################################################################

#------------------- Batch 1 ----------------------
# Raw files

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_1st_Batch/Myeloid_samples_fcsfiles",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_1st_Batch/Plate_1_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/1st_batch/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_1st_Batch/Plate_1_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/1st_batch/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_1st_Batch/Plate_1_map.csv",
               type_of_data = "WSP")
#Plot
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/1st_batch/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_1st_Batch/Plate_1_map.csv",
               type_of_data = "Plot")
#Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/1st_batch/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_1st_Batch/Plate_1_map.csv",
               type_of_data = "Stats")


#------------------- Batch 2 ----------------------

# Raw files plate 1

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_1/Flowjo files/96 Well plate1/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_1/Plate_1_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_1_results/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_1/Plate_1_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_1_results/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_1/Plate_1_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_1_results/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_1/Plate_1_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_1_results/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_1/Plate_1_map.csv",
               type_of_data = "Stats")

# Raw files plate 2

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_2/Flowjo files/96 Well plate1/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_2/Plate_2_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_2_results/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_2/Plate_2_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_2_results/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_2/Plate_2_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_2_results/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_2/Plate_2_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_2_results/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_2/Plate_2_map.csv",
               type_of_data = "Stats")
# Raw files plate 3

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_3/Flowjo files/96 Well plate1/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_3/Plate_3_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_3_results/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_3/Plate_3_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_3_results/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_3/Plate_3_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_3_results/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_3/Plate_3_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_3_results/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_3/Plate_3_map.csv",
               type_of_data = "Stats")

# Raw files plate 4

renaming_files(path_files="/home/rstudio/data/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_4/Flowjo files/96 Well plate1/",
               path_names="/home/rstudio/data/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_4/Plate_4_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/rstudio/data/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_4_results/Cleaning",
               path_names="/home/rstudio/data/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_4/Plate_4_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/rstudio/data/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_4_results/WSP",
               path_names="/home/rstudio/data/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_4/Plate_4_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/rstudio/data/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_4_results/Plots",
               path_names="/home/rstudio/data/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_4/Plate_4_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/rstudio/data/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_4_results/Stats",
               path_names="/home/rstudio/data/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_4/Plate_4_map.csv",
               type_of_data = "Stats")
# Raw files plate 6

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_6/Flowjo files/96 Well plate1/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_6/Plate_6_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_6_results/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_6/Plate_6_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_6_results/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_6/Plate_6_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_6_results/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_6/Plate_6_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/2nd_batch/plate_6_results/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 2nd Batch (G022-G285)/Myeloid_Gambia_(G022-G285)/Plate_6/Plate_6_map.csv",
               type_of_data = "Stats")
#------------------- Batch 3 ----------------------
# Raw files plate 1

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_1/Flowjo files/96 Well plate1/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_1/Plate_1_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_1_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_1/Plate_1_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_1_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_1/Plate_1_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_1_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_1/Plate_1_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_1_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_1/Plate_1_map.csv",
               type_of_data = "Stats")

# Raw files plate 2

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_2/Flowjo files/FCS files_samples/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_2/Plate_2_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_2_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_2/Plate_2_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_2_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_2/Plate_2_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_2_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_2/Plate_2_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_2_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_2/Plate_2_map.csv",
               type_of_data = "Stats")
# Raw files plate 3

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_3/Flowjo files/Samples files/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_3/Plate_3_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_3_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_3/Plate_3_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_3_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_3/Plate_3_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_3_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_3/Plate_3_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_3_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_3/Plate_3_map.csv",
               type_of_data = "Stats")
# Raw files plate 4

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_4/Flowjo files/96 Well plate1/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel//Plate_4/Plate_4_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_4_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_4/Plate_4_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_4_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_4/Plate_4_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_4_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_4/Plate_4_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_4_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_4/Plate_4_map.csv",
               type_of_data = "Stats")

# Raw files plate 5

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel//Plate_5/FCS files/FCS files samples plate 1",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel//Plate_5/Plate_5_map.csv",
               type_of_data = "raw_files")
renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel//Plate_5/FCS files/FCS files Samples  plate 2",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel//Plate_5/Plate_5_map_2.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 1/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map.csv",
               type_of_data = "Cleaning")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 2/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map_2.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 1/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map.csv",
               type_of_data = "WSP")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 2/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map_2.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 1/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map.csv",
               type_of_data = "Plot")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 2/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map_2.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 1/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map.csv",
               type_of_data = "Stats")
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/3rd_batch/Plate_5_samples/FCS files samples plate 2/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia 3rd Batch (G286-528)/Myeloid Panel/Plate_5/Plate_5_map_2.csv",
               type_of_data = "Stats")

#------------------- Batch 4 ----------------------
# Raw files plate 1

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_1/FCS_files/",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_1/Plate_1_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_1_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_1/Plate_1_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_1_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_1/Plate_1_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_1_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_1/Plate_1_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_1_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_1/Plate_1_map.csv",
               type_of_data = "Stats")

# Raw files plate 2

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_2/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_2/Plate_2_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_2_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_2/Plate_2_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_2_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_2/Plate_2_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_2_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_2/Plate_2_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_2_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_2/Plate_2_map.csv",
               type_of_data = "Stats")


# Raw files plate 3

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_3/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_3/Plate_3_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_3_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_3/Plate_3_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_3_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_3/Plate_3_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_3_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_3/Plate_3_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_3_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_3/Plate_3_map.csv",
               type_of_data = "Stats")
# Raw files plate 4

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_4/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_4/Plate_4_map.csv",
               type_of_data = "raw_files")

# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_4_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_4/Plate_4_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_4_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_4/Plate_4_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_4_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_4/Plate_4_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_4_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_4/Plate_4_map.csv",
               type_of_data = "Stats")
# Raw files plate 5

renaming_files(path_files="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_5/FCS_files",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_5/Plate_5_map.csv",
               type_of_data = "raw_files")
# Cleaning
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_5_samples/Cleaning",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_5/Plate_5_map.csv",
               type_of_data = "Cleaning")
# WSP
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_5_samples/WSP",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_5/Plate_5_map.csv",
               type_of_data = "WSP")
# Plots
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_5_samples/Plots",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_5/Plate_5_map.csv",
               type_of_data = "Plot")
# Stats
renaming_files(path_files="/home/smontante/Desktop/Brinkman_lab_projects/Second_project_Gambia/Renamed_Myeloid_all_batches_results_s3/4th_batch/plate_5_samples/Stats",
               path_names="/mnt/f/FCS data/Tobias Kollmann/Renamed_Gambia_samples_myeloid_final/Gambia_Batch_4(529-711)/Plate_5/Plate_5_map.csv",
               type_of_data = "Stats")
