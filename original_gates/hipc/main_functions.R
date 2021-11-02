#--------Executing all gating functions myeloid panel -----
main_myeloid<-function(n_sample="All",path_comp_matrix,path_dir_files,path_fcs_files_fixing=NULL,
               plots=T,flowjo_wsp=F,name_output="Default",n_cores=1,
               pre_processing=T,export_gs=F,export_thresholds=T,density_axis=F,
               path_output_clr="/home/rstudio/results/CLR_result",export_fs_clean=F,
               export_fcs_live_cells=F,export_CLR=F){
  set.seed(40,kind = "L'Ecuyer-CMRG") 
  #----- preprocessing
  output_import_files <- import_files_v2(n_sample=n_sample,files_path=path_dir_files,n_cores=n_cores)
  fs <- output_import_files$fs
  
  FCSfiles_path_list <- output_import_files$FCSfiles_path_list
  path.output <- create_path_output_dir("Myeloid_Gambia")
  output_find_markers <- find_markers_indices(fs=fs)
  scat.chans <-output_find_markers$scat.chans
  fluorochrome.chans <-output_find_markers$fluorochrome.chans
  Time.channel <-output_find_markers$Time.channel
  # fix fs if desc column has only NA
  test_na<-is.na(fs[[1]]@parameters@data$desc)
  if(all(test_na)==T){
    marker_names<-names(fluorochrome.chans)
    indices_marker_names<-as.vector(fluorochrome.chans)
    list_frames<-lapply(1:length(fs),function(i){
      fs[[i]]@parameters@data$desc[indices_marker_names]<-marker_names
      return(fs[[i]])
    })
    fs<-as(list_frames,"flowSet")
  }
  # generate empty gatingSet from fs
  gs <-get_GatingSet(fs,FCSfiles_path_list)
  if(pre_processing==T){
    # --- gating dei marginal events
    output<-margin_gating(fs,gs)
    fs.marg<-output$fs_marg
    gs<-output$gs
    # ---- Compensation e trasformation
    comp<-import_comp_matrix_v2(path_comp_matrix=path_comp_matrix)
    output<-comp_and_transform(gs,comp)
    fs.marg<-output$fs.marg
    gs<-output$gs
    # -------cleaning
    output<-cleaning(fs.marg,gs,n_cores,path.output)
    fs.clean<-output$fs.clean
    gs <-output$gs
    clean.inds <- output$clean.inds
    if(export_fs_clean==T){
      write.flowSet(fs.clean,outdir = "/home/rstudio/results/shiny_app_files/Cleaned_files")
      stop()
    }
  }else{ # cleaning already executed
    fs<-getData(gs,"root")
    fs.clean<-fs
  }
  gc()
  # ---- make df tresholds
  df_thresholds<-make_tresholds_df_myeloid(gs)
  #------ Singlets
  output_2<-gating_Time_to_singlets(fs.clean,gs,
                                    pre_processing = pre_processing,n_cores,
                                    scat.chans)
  fs.sngl<- output_2$fs.sngl
  gs<-output_2$gs
  singlets <- output_2$list_singlets
  singlets_pop <- output_2$list_singlets_pop
  
  gc()
  if(export_CLR==T){
    generate_CLR_data(fs.clean,list(singlets_pop),name_pop = "Singlets",channels=c("FSC-A","FSC-H"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  #------- Size e beads
  output_3<-gating_singlets_to_Size_beads(fs.sngl,gs,n_cores,scat.chans)
  fs.size<-output_3$fs.size
  gs<-output_3$gs
  beads_poly <- output_3$beads_poly_list
  size_poly <- output_3$size_poly_list
  beads <-output_3$list_beads
  list_size <-output_3$list_size
  
  gc()
  if(export_CLR==T){
    generate_CLR_data(fs.sngl,list(list_size,beads),name_pop = "Size_beads",channels=c("FSC-A","SSC-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  #---- Live
  output_4<-gating_size_to_live(fs.size,gs,n_cores,fluorochrome.chans,scat.chans)
  fs.live<-output_4$fs.live
  gs<-output_4$gs
  live_poly <- output_4$live_poly_list
  list_live <- output_4$list_live
  
  if(export_fcs_live_cells==T){
    # export live cells root for flowTypeFilter
    write.flowSet(fs.live,outdir = "/home/rstudio/results/flowType_files/fcs_files")
  }
  gc()
  if(export_CLR==T){
    generate_CLR_data(fs.size,list(list_live),name_pop = "Live",channels=c("APC-eF780-A","SSC-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  #------ to Granulocytes e non gran
  output_5<-gating_live_to_granulocytes(fs.live,gs,n_cores,fluorochrome.chans)
  fs.gran<-output_5$fs.gran
  fs.66n45p<-output_5$fs.66n45p
  gs<-output_5$gs
  gran_poly<- output_5$lista_poly_gran
  cd66_45_poly <- output_5$lista_poly_cd45
  lista_CD66negCD45pos <-output_5$lista_CD66negCD45pos
  lista_granulocytes <-output_5$lista_granulocytes
  
  filter_CD66negCD45pos<-sapply(1:length(lista_CD66negCD45pos), function(i){
    filter_i<-lista_CD66negCD45pos[[i]]@filter
    return(filter_i)
  })
  filter_granulocytes<-sapply(1:length(lista_granulocytes), function(i){
    filter_i<-lista_granulocytes[[i]]@filter
    return(filter_i)
  })
  saveRDS(filter_CD66negCD45pos,"/home/rstudio/results/flowType_files/filter_objects/filter_CD66negCD45pos.rds")
  saveRDS(filter_granulocytes,"/home/rstudio/results/flowType_files/filter_objects/filter_granulocytes.rds")
  if(export_CLR==T){
    generate_CLR_data(fs.live,list(lista_CD66negCD45pos,lista_granulocytes),name_pop = "CD66-CD45+_CD66+CD45+",channels=c("BV711-A","V450-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  gc()
  # ------ Gating to monocytes (14-) e 14+
  output_6<-gating_ngran_to_DRpos(fs.66n45p,gs,fluorochrome.chans)
  fs.mono<-output_6$fs.mono
  fs.drn.14n<-output_6$fs.DRn.14n
  fs.drp.14n<-output_6$fs.DRp.14n
  gs<-output_6$gs
  list_monocytes<-output_6$list_monocytes
  list_DRp_14n <- output_6$list_DRp_14n
  list_DRn_14n <- output_6$list_DRn_14n
  # export filters
  filter_monocytes<-sapply(1:length(list_monocytes), function(i){
    filter_i<-list_monocytes[[i]]@filter
    return(filter_i)
  })
  filter_DRp_14n<-sapply(1:length(list_DRp_14n), function(i){
    filter_i<-list_DRp_14n[[i]]@filter
    return(filter_i)
  })
  filter_DRn_14n<-sapply(1:length(list_DRn_14n), function(i){
    filter_i<-list_DRn_14n[[i]]@filter
    return(filter_i)
  })
  saveRDS(filter_monocytes,"/home/rstudio/results/flowType_files/filter_objects/filter_monocytes.rds")
  saveRDS(filter_DRp_14n,"/home/rstudio/results/flowType_files/filter_objects/filter_DRp_14n.rds")
  saveRDS(filter_DRn_14n,"/home/rstudio/results/flowType_files/filter_objects/filter_DRn_14n.rds")
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.66n45p,list(list_monocytes,list_DRp_14n,list_DRn_14n),name_pop = "HLADR+CD14+_HLADR+CD14-_HLADR-CD14-",channels=c("BV605-A","V500-A"),filters=F,markers = c("HLA-DR|HLADR","CD14"),density=density_axis,path.output=path_output_clr)
  }
  
  #  ------ Gating to classical mono
  output_7<-gating_mono_to_classic_mono(fs.mono,gs,fluorochrome.chans)
  gs<-output_7$gs
  fs.classmono <- output_7$fs.classmono
  CD16mon_gate <- output_7$list_CD16mon_gate
  list_classicalmono<-output_7$list_classicalmono
  list_nonclassicalmono<-output_7$list_nonclassicalmono
  # export thresholds
  thresholds_classicalmono<-sapply(1:length(list_classicalmono), function(i){
    thr_i<-paste0(list_classicalmono[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="Classical Monocytes")
  df_thresholds[,inds_col]<-thresholds_classicalmono
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.mono,list(list_classicalmono,list_nonclassicalmono),name_pop = "CD16-_CD16+",channels=c("V500-A","FITC-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  
  # ------- Gating to B_pdc_mdc
  output_8<-gating_drp_14n_to_B_pdc_mdc(fs.drp.14n,gs,fluorochrome.chans,n_cores)
  gs<-output_8$gs
  mdc_poly <- output_8$lista_mdc_poly
  pdc_poly <- output_8$lista_pdc_poly
  bcell_poly <- output_8$lista_bcell_poly
  list_Bcells<-output_8$list_Bcells
  list_mDC<-output_8$list_mDC
  list_pDC<-output_8$list_pDC
  # export thresholds
  thresholds_Bcells<-sapply(1:length(list_Bcells), function(i){
    thr_i<-paste0(list_Bcells[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="B cells")
  df_thresholds[,inds_col]<-thresholds_Bcells
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.drp.14n,list(list_Bcells,list_mDC,list_pDC),name_pop = "CD123-CD11c-_CD123CD11c+_CD123+CD11c-",channels=c("PE-Cy7-A","APC-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  
  #------- Gating to t_gdt_cd3n
  output_9<-gating_drn_14n_to_T_gdt_cd3(fs.drn.14n,gs,fluorochrome.chans,scat.chans)
  fs.tcell<-output_9$fs.tcell
  fs.cd3n <- output_9$fs.cd3n
  gs<-output_9$gs
  list_CD3_gate <- output_9$list_CD3_gate
  list_gdt_poly <- output_9$list_gdt_poly
  list_gdneg_poly <- output_9$list_gdneg_poly
  list_CD3Tcells<-output_9$list_CD3Tcells
  list_CD3neg<-output_9$list_CD3neg
  # export threholds
  thresholds_CD3Tcells<-sapply(1:length(list_CD3Tcells), function(i){
    thr_i<-paste0(list_CD3Tcells[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  list_gdTcells<-output_9$list_gdTcells
  list_gdnegcells<-output_9$list_gdnegcells
  thresholds_gdTcells<-sapply(1:length(list_gdTcells), function(i){
    thr_i<-paste0(list_gdTcells[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  thresholds_cd3<-sapply(1:length(list_CD3_gate), function(i){
    thr_i<-list_CD3_gate[[i]]
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD3+ T cells")
  df_thresholds[,inds_col]<-thresholds_CD3Tcells
  inds_col<-which(colnames(df_thresholds)=="gd T cells")
  df_thresholds[,inds_col]<-thresholds_gdTcells
  inds_col<-which(colnames(df_thresholds)=="CD3-")
  df_thresholds[,inds_col]<-thresholds_cd3
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.drn.14n,list(list_CD3Tcells,list_CD3neg),name_pop = "CD3+_CD3-",channels=c("PE-CF594-A","SSC-A"),filters=F,density=density_axis,path.output=path_output_clr)
    generate_CLR_data(list_CD3Tcells,list(list_gdTcells,list_gdnegcells),name_pop = "gd+_gd-",channels=c("PE-CF594-A","PE-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  
  #------- gating to NKT e NK
  output_10<-gating_nkt_nk(fs.tcell,fs.cd3n,gs,fluorochrome.chans)
  gs<-output_10$gs
  list_all_nkt_related_pops <- output_10$list_all_nkt_related_pops
  list_CD3.16p56p <- output_10$list_CD3.16p56p
  list_CD3.16n56n <- output_10$list_CD3.16n56n
  list_CD3.16n56p <- output_10$list_CD3.16n56p
  list_CD3.16p56n <- output_10$list_CD3.16p56n
  
  fs.nk.56n16n <- output_10$fs.nk.56n16n
  list_nk_56n16n <- output_10$list_nk_56n16n
  list_NK_temp <- output_10$list_NK_temp
  list_nk_56hi <- output_10$list_nk_56hi
  list_nk.56n16p <- output_10$list_nk.56n16p
  list_nk.56dim16n <- output_10$list_nk.56dim16n
  list_nk.56p16p <- output_10$list_nk.56p16p
  
  # export thresholds
  thresholds_nkt<-sapply(1:length(list_all_nkt_related_pops), function(i){
    thr_i<-paste0(list_all_nkt_related_pops[[i]]$CD3.16p56p@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD56+CD16+ NKT cells")
  df_thresholds[,inds_col]<-thresholds_nkt
  thresholds_nk_56n16n<-sapply(1:length(list_nk_56n16n), function(i){
    thr_i<-paste0(list_nk_56n16n[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD56-CD16- cells")
  df_thresholds[,inds_col]<-thresholds_nk_56n16n
  
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.tcell,list(list_CD3.16p56p,list_CD3.16n56n,list_CD3.16n56p,list_CD3.16p56n),name_pop = "CD16+CD56+_CD16-CD56-_CD16-CD56+_CD16+CD56-",channels=c("FITC-A","BV650-A"),filters=F,density=density_axis,path.output=path_output_clr)
    generate_CLR_data(fs.cd3n,list(list_nk_56n16n,list_nk.56n16p,list_nk_56hi,list_nk.56dim16n,list_nk.56p16p),name_pop = "NKCD16-CD56-_CD16+CD56-_CD16-CD56++_CD16-CD56dim_CD56+CD16+",channels=c("FITC-A","BV650-A"),filters=F,density=density_axis,path.output=path_output_clr)
    
  }
  
  # ------ gating Baso
  output_11 <- gating_nk_56n16n_to_Baso(fs.nk.56n16n,gs,fluorochrome.chans)
  gs<-output_11$gs
  list_basophils<-output_11$list_basophils
  # export filters
  filter_basophils<-lapply(1:length(list_basophils), function(i){
    filter_i<-list_basophils[[i]]@filter
    return(filter_i)
  })
  names(filter_basophils)<-sampleNames(gs)
  saveRDS(filter_basophils,"/home/rstudio/results/flowType_files/filter_objects/filter_basophils.rds")
  baso_poly <- output_11$list_basophil_poly
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.nk.56n16n,list(list_basophils),name_pop = "CD123-HLADR-",channels=c("PE-Cy7-A","BV605-A"),filters=F,markers = c("CD123","HLA-DR|HLADR"),density=density_axis,path.output=path_output_clr)
  }
  
  #---- gating  to mature e immature granulocytes
  output_12 <- gating_gran_to_mgran_im_gran(fs.gran,gs,fluorochrome.chans)
  gs<-output_12$gs
  fs.11p16p<-output_12$fs.11p16p
  CD16granulo_gate <- output_12$list_CD16granulo_gate
  CD11bgranulo_gate <- output_12$list_CD11bgranulo_gate
  list_maturegranulo<-output_12$list_maturegranulo
  list_immaturegranulo<-output_12$list_immaturegranulo
  list_CD11bposCD16neg.granulo<-output_12$list_CD11bposCD16neg.granulo
  list_CD11bnegCD16pos.neutro<-output_12$list_CD11bnegCD16pos.neutro
  
  # fill the df thresholds
  thresholds_maturegranulo<-sapply(1:length(list_maturegranulo), function(i){
    thr_i<-paste0(list_maturegranulo[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD11b+CD16+ Mature Neutrophils")
  df_thresholds[,inds_col]<-thresholds_maturegranulo
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.gran,list(list_maturegranulo,list_immaturegranulo,list_CD11bposCD16neg.granulo,list_CD11bnegCD16pos.neutro),name_pop = "CD11b+CD16+_CD11b-CD16-_CD11b+CD16-_CD11b-CD16+",channels=c("BV786-A","FITC-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  
  #------------ Gating to CD64+
  output_prova_13 <- gating_mneutro_to_CD64pos(fs.11p16p,gs,fluorochrome.chans)
  gs<-output_prova_13$gs
  CD64_gate <- output_prova_13$list_CD64_gate
  # fill the df threholds
  thresholds_CD64_gate<-sapply(1:length(CD64_gate), function(i){
    thr_i<-CD64_gate[[i]]
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD64+")
  df_thresholds[,inds_col]<-thresholds_CD64_gate
  #---- Generate final Data
  if(flowjo_wsp==T){
    generate_final_data(path.output,name_output = name_output,gs)
  }
  #---- Generates plots of each Gating step for each sample
  if(plots==T){
    generate_plots_myeloid(gs,n_cores,fs.marg,fs.clean,singlets,fs.sngl,
                                    beads_poly,size_poly,fs.size,live_poly,
                                    fs.live,cd66_45_poly,gran_poly,fs.11p16p,CD64_gate,
                                    fs.gran,CD11bgranulo_gate,fs.66n45p,list_DRn_14n,list_DRp_14n,
                                    list_monocytes,fs.drn.14n,list_CD3_gate,fs.tcell,
                                    list_gdt_poly,list_gdneg_poly,list_all_nkt_related_pops,
                                    fs.nk.56n16n,baso_poly,fs.cd3n,list_nk_56n16n,list_NK_temp,
                                    list_nk_56hi,list_mDC,list_pDC,pdc_poly,mdc_poly,bcell_poly,
                                    fs.mono,CD16mon_gate,path.output,Time.channel,fluorochrome.chans,
                           scat.chans,clean.inds,CD16granulo_gate,fs)
  }
  if(export_gs==T){
    setwd("/home/rstudio/results/shiny_app_files/GatingSet")
    n_samples<-length(gs)
    for(i in 1:n_samples){
      gs[[i]]@data[[i]]@description$`$FIL`<<- sampleNames(gs)[i]
    }
    flowWorkspace::save_gs(gs,path = "/home/rstudio/results/shiny_app_files/GatingSet/GatingSet_files",overwrite = T)
  }
  if(export_thresholds==T){
    write.csv(df_thresholds,"/home/rstudio/results/flowType_files/df_thresholds/df_thresholds.csv")
  }
  return(gs)
}


#--------Executing all gating functions Bcells panel -----
main_bcells<-function(n_sample="All",path_comp_matrix,path_dir_files,path_fcs_files_fixing=NULL,plots=T,
               flowjo_wsp=F,n_cores=1,name_output="Default",pre_processing=T,export_gs=F,export_thresholds=T,
               density_axis=F,path_output_clr="/home/rstudio/results/CLR_result",
               export_fs_clean=F,export_fcs_live_cells=F,export_CLR=F){
  set.seed(40)
  #----- preprocessing
  output_import_files <- import_files_v2(n_sample=n_sample,files_path=path_dir_files,n_cores=n_cores)
  fs <- output_import_files$fs
  
  FCSfiles_path_list <- output_import_files$FCSfiles_path_list
  path.output <- create_path_output_dir("Bcells_Gambia")
  if(is.null(path_fcs_files_fixing)==FALSE){
    fs<-fix_markers_names(fs_to_fix=fs,path_fcs_files_fixing=path_fcs_files_fixing)
  }
  output_find_markers <- find_markers_indices_bcells_v2(fs=fs)
  scat.chans <-output_find_markers$scat.chans
  fluorochrome.chans <-output_find_markers$fluorochrome.chans
  # generate empty gatingSet from fs
  gs<-get_GatingSet(fs,FCSfiles_path_list)
  # code to skip the cleaning step
  if(pre_processing==T){
    # ---- gating dei marginal events
    output<-margin_gating(fs,gs)
    fs.marg<-output$fs.marg
    gs<-output$gs
    # ----- compensation and trasformation
    comp<-import_comp_matrix_v2(path_comp_matrix=path_comp_matrix)
    out<-comp_and_transform(gs,comp)
    gs<-out$gs
    fs.marg<-out$fs.marg
    #----- cleaning
    output <- cleaning(fs.marg,gs,n_cores,path.output)
    fs.clean <- output$fs.clean
    gs <- output$gs
    clean.inds<-output$clean.inds
    if(export_fs_clean==T){
      write.flowSet(fs.clean,outdir = "/home/rstudio/results/shiny_app_files/Cleaned_files")
      stop()
    }
  }else{
    fs<-getData(gs,"root")
    fs.clean<-fs
  }
  # ------ make df tresholds
  df_thresholds<-make_tresholds_df_bcells(gs)
  #----- gating time to singlets
  output_2 <- gating_time_to_singlets(fs.clean,gs,pre_processing,n_cores,scat.chans)
  gs <- output_2$gs
  fs.sngl <- output_2$fs.sngl
  singlets <- output_2$list_singlets
  singlets_pop <- output_2$list_singlets_pop
  
  # fill the df thresholds
  list_rot.gate <- output_2$list_rot.gate
  thresholds_singlets<-sapply(1:length(list_rot.gate), function(i){
    thr_i<-list_rot.gate[[i]]
    return(thr_i)
  })
  inds_col<-grep("Singlets",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_singlets
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.clean,list(singlets_pop),name_pop = "Singlets",channels=c("FSC-A","FSC-H"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  #------- Size e beads
  output_3 <- gating_singlets_to_size_beads(fs.sngl,gs,n_cores,scat.chans)
  fs.size <- output_3$fs.size
  gs <- output_3$gs
  size_poly <- output_3$size_poly_list
  beads_poly <- output_3$beads_poly_list
  beads <-output_3$list_beads
  # fill the df thresholds
  list_size <- output_3$list_size
  thresholds_size<-sapply(1:length(list_size), function(i){
    thr_i<-paste0(list_size[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("All cells",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_size
  if(export_CLR==T){
    generate_CLR_data(fs.sngl,list(list_size,beads),name_pop = "Size_beads",channels=c("FSC-A","SSC-A"),filters=list(size_poly,beads_poly),density=density_axis,path.output=path_output_clr)
  }

  #---- Live
  output_4 <- gating_size_to_live(fs.size,gs,n_cores,fluorochrome.chans,scat.chans)
  fs.live<-output_4$fs.live
  gs<-output_4$gs
  live_poly <- output_4$lista_poly_live
  # fill the df thresholds
  list_live <- output_4$list_live
  thresholds_live<-sapply(1:length(list_live), function(i){
    thr_i<-paste0(list_live[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("Live cells",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_live
  if(export_fcs_live_cells==T){
    # export live cells root for flowTypeFilter
    write.flowSet(fs.live,outdir = "/home/rstudio/results/flowType_files/fcs_files")
  }
  if(export_CLR==T){
    generate_CLR_data(fs.size,list(list_live),name_pop = "Live",channels=c("APC-eF780-A","SSC-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  #------ to Granulocytes e non gran
  output_5 <- gating_live_to_non_Granulocytes(fs.live,gs,n_cores,fluorochrome.chans)
  gs <- output_5$gs
  fs.non_gran <- output_5$fs.non_gran
  gran_poly <- output_5$lista_poly_gran
  non_gran_poly <- output_5$lista_poly_non_gran
  # fill the df thresholds
  list_non_gran <- output_5$list_non_gran
  thresholds_non_gran<-sapply(1:length(list_non_gran), function(i){
    thr_i<-paste0(list_non_gran[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("Non granulocytes",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_non_gran
  list_gran <- output_5$list_gran
  thresholds_gran<-sapply(1:length(list_gran), function(i){
    thr_i<-paste0(list_gran[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("Granulocytes",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_gran
  # export filter list
  filter_non_gran<-sapply(1:length(list_non_gran), function(i){
    filter_i<-list_non_gran[[i]]@filter
    return(filter_i)
  })
  names(filter_non_gran)<-sampleNames(gs)
  saveRDS(filter_non_gran,"/home/rstudio/results/flowType_files/filter_objects/filter_non_gran.rds")
  if(export_CLR==T){
    generate_CLR_data(fs.live,list(list_non_gran,list_gran),name_pop = "CD14CD66-_CD14CD66+",channels=c("BV711-A","V500-A"),filters=F,density=density_axis,path.output=path_output_clr)
  }
  #--- Bcells(CD19+)
  output_6 <-gating_nonGran_to_Bcells(fs.non_gran,gs,fluorochrome.chans,scat.chans)
  gs<-output_6$gs
  fs.bcells<-output_6$fs.bcells
  Bcells_poly<-output_6$list_poly_Bcells
  # fill the df thresholds
  list_Bcells <- output_6$list_Bcells
  list_NotBcells <- output_6$list_NotBcells
  
  thresholds_Bcells<-sapply(1:length(list_Bcells), function(i){
    thr_i<-paste0(list_Bcells[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD19+ B cells")
  df_thresholds[,inds_col]<-thresholds_Bcells
  # export filter list
  filter_Bcells<-sapply(1:length(list_Bcells), function(i){
    filter_i<-list_Bcells[[i]]@filter
    return(filter_i)
  })
  names(filter_Bcells)<-sampleNames(gs)
  saveRDS(filter_Bcells,"/home/rstudio/results/flowType_files/filter_objects/filter_Bcells.rds")

  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.non_gran,list(list_Bcells,list_NotBcells),name_pop = "CD19+_CD19-",channels=c("APC-A","SSC-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  
  # #--- bcells to PC
  output_7 <- gating_bcells_to_pc(fs.bcells,gs,fluorochrome.chans)
  gs<-output_7$gs
  fs.pc<-output_7$fs.pc
  pc_poly<-output_7$list_poly_pc
  # fill the df thresholds
  list_pc <- output_7$list_pc
  thresholds_pc<-sapply(1:length(list_pc), function(i){
    thr_i<-paste0(list_pc[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("plasma cells",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_pc
  # export filter list
  filter_pc<-sapply(1:length(list_pc), function(i){
    filter_i<-list_pc[[i]]@filter
    return(filter_i)
  })
  names(filter_pc)<-sampleNames(gs)
  saveRDS(filter_pc,"/home/rstudio/results/flowType_files/filter_objects/filter_pc.rds")

  if(export_CLR==T){
    generate_CLR_data(fs.bcells,list(list_pc),name_pop = "CD38+CD138+",channels=c("PE-Cy5-A","V450-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  
  # #--- bcells to PB
  output_8<-gating_bcells_to_PB(fs.bcells,gs,fluorochrome.chans)
  gs <- output_8$gs
  fs.PB <- output_8$fs.PB
  PB_poly<-output_8$list_poly_PB
  
  # fill the df thresholds
  list_PB<- output_8$list_PB
  thresholds_PB<-sapply(1:length(list_PB), function(i){
    thr_i<-paste0(list_PB[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("Plasmablasts",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_PB
  # export filter list
  filter_PB<-sapply(1:length(list_PB), function(i){
    filter_i<-list_PB[[i]]@filter
    return(filter_i)
  })
  names(filter_PB)<-sampleNames(gs)
  saveRDS(filter_PB,"/home/rstudio/results/flowType_files/filter_objects/filter_PB.rds")

  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs.bcells,list(list_PB),name_pop = "CD38+CD27+",channels=c("PE-Cy5-A","PE-Cy7-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  
  # --- Bcells to CD19/20
  output_9<-gating_bcells_to_CD20p(fs.bcells,gs,fluorochrome.chans)
  gs<-output_9$gs
  fs.CD20pos<-output_9$fs.CD20pos
  fs.CD20neg<-output_9$fs.CD20neg
  CD20pos_poly<-output_9$list_CD20_pos_poly
  CD20neg_poly<-output_9$list_CD20_neg_poly
  
  # fill the df thresholds
  list_cd20.neg<- output_9$list_cd20.neg
  list_cd20.pos<- output_9$list_cd20.pos
  
  thresholds_cd20.neg<-sapply(1:length(list_cd20.neg), function(i){
    thr_i<-paste0(list_cd20.neg[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD19+CD20-")
  df_thresholds[,inds_col]<-thresholds_cd20.neg
  # export filter list
  filter_cd20_neg<-sapply(1:length(list_cd20.neg), function(i){
    filter_i<-list_cd20.neg[[i]]@filter
    return(filter_i)
  })
  names(filter_cd20_neg)<-sampleNames(gs)
  saveRDS(filter_cd20_neg,"/home/rstudio/results/flowType_files/filter_objects/filter_cd20_neg.rds")
  # export CLR if needed
  if(export_CLR==T){
    generate_CLR_data(fs.bcells,list(list_cd20.neg,list_cd20.pos),name_pop = "CD19CD20-_CD19CD20+",channels=c("APC-A","BV786-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  
  # ------- Bcells to igd/igm
  output_10<-gating_bcells_to_igm_igd(fs.bcells,gs,fluorochrome.chans)
  gs<-output_10$gs
  quad1<-output_10$list_quad1_pop
  quad2<-output_10$list_quad.2
  quad3<-output_10$list_quad.3
  quad4<-output_10$list_quad.4
  
  fs.igm<-output_10$fs.igm
  fs.IGDnegIGMpos<-getData(gs,"IGD-IGM+ B cells")
  fs.IGDnegIGMneg<-getData(gs,"IGD-IGM- B cells")
  fs.IGDposIGMpos<-getData(gs,"IGD+IGM+ B cells")
  fs.IGDposIGMneg<-getData(gs,"IGD+IGM- B cells")
  
  # fill the df thresholds
  thresholds_quad1<-sapply(1:length(quad1), function(i){
    thr_i<-paste0(quad1[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="IGD-IGM- B cells")
  df_thresholds[,inds_col]<-thresholds_quad1

  
  if(export_CLR==T){
    generate_CLR_data(fs.bcells,list(quad1,quad2,quad3,quad4),name_pop = "IgD-IgM-_IgD+IgM-_IgD-IgM+_IgD+IgM+",channels=c("PE-CF594-A","FITC-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  
  #------- cd20 to cd10neg
  output_11<-gating_cd20_to_cd10neg(fs.CD20pos,gs,fluorochrome.chans)
  gs<-output_11$gs
  fs.cd10pos<-output_11$fs.cd10pos
  fs.cd10neg<-output_11$fs.cd10neg
  CD10pos_poly<-output_11$list_CD10pos_poly
  CD10neg_poly <-output_11$list_CD10neg_poly
  list_cd10.neg<-output_11$list_cd10.neg
  # fill the df thresholds
  list_cd10.neg<- output_11$list_cd10.neg
  list_cd10.pos<- output_11$list_cd10.pos
  
  thresholds_cd10.neg<-sapply(1:length(list_cd10.neg), function(i){
    thr_i<-paste0(list_cd10.neg[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-which(colnames(df_thresholds)=="CD10-")
  df_thresholds[,inds_col]<-thresholds_cd10.neg
  # export filter list
  filter_cd10_neg<<-sapply(1:length(list_cd10.neg), function(i){
    filter_i<-list_cd10.neg[[i]]@filter
    return(filter_i)
  })
  names(filter_cd10_neg)<-sampleNames(gs)
  saveRDS(filter_cd10_neg,"/home/rstudio/results/flowType_files/filter_objects/filter_cd10_neg.rds")

  # export CLR if needed
  if(export_CLR==T){
    generate_CLR_data(fs.CD20pos,list(list_cd10.neg,list_cd10.pos),name_pop = "CD27CD10-_CD27CD10+",channels=c("Alexa Fluor 700-A","PE-Cy7-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }

  #----- cd10neg to Naive
  output_12<-gating_cd10n_to_naive(fs.cd10neg,gs,fluorochrome.chans)
  gs<-output_12$gs
  fs.naive<-getData(gs,"Naive B cells") 
  naive_poly<-output_12$list_quad3_poly
  quad.3<-output_12$list_quad.3
  quad.1<-output_12$list_quad.1
  quad.2<-output_12$list_quad.2
  quad.4<-output_12$list_quad.4
  
  # fill the df thresholds
  thresholds_quad.3<-sapply(1:length(quad.3), function(i){
    thr_i<-paste0(quad.3[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("Naive B cells",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_quad.3

  
  if(export_CLR==T){
    generate_CLR_data(fs.cd10neg,list(quad.1,quad.2,quad.3,quad.4),name_pop = "CD27-IgD-_CD27+IgD-_CD27-IgD+_CD27+IgD+",channels=c("PE-Cy7-A","PE-CF594-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  
  # ---- to trans
  output_13<-gating_IGM_IGD_CD27_to_tran(fs.igm,quad.3,gs,fluorochrome.chans)
  gs<-output_13$gs
  trans_poly<-output_13$list_trans_poly
  fs.trans<-getData(gs,"Immature transition B cells")
  list_trans<-output_13$list_trans
  list_f_cd27.igd<-output_13$list_f_cd27.igd
  fs_cd27.igd<-as(list_f_cd27.igd,"flowSet")
  sampleNames(fs_cd27.igd)<-sampleNames(gs)
  # fill the df thresholds
  thresholds_list_trans<-sapply(1:length(list_trans), function(i){
    thr_i<-paste0(list_trans[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("Immature transition B cells",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_list_trans
  # export filter list
  filter_trans<-sapply(1:length(list_trans), function(i){
    filter_i<-list_trans[[i]]@filter
    return(filter_i)
  })
  names(filter_trans)<-sampleNames(gs)
  saveRDS(filter_trans,"/home/rstudio/results/flowType_files/filter_objects/filter_trans.rds")
  
  # export clr if needed
  if(export_CLR==T){
    generate_CLR_data(fs_cd27.igd,list(list_trans),name_pop = "CD10-CD38+",channels=c("Alexa Fluor 700-A","PE-Cy5-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  
  # ----- from live to blast
  output_14<-gating_live_to_blast(fs.live,gs,fluorochrome.chans,scat.chans)
  gs<-output_14$gs
  blast_poly<-output_14$list_poly_blast
  fs.blast<-getData(gs,"Blasts")
  
  # fill the df thresholds
  list_blast<-output_14$list_blast
  thresholds_list_blast<-sapply(1:length(list_blast), function(i){
    thr_i<-paste0(list_blast[[i]]@gates,collapse = ",")
    return(thr_i)
  })
  inds_col<-grep("Blasts",colnames(df_thresholds))
  df_thresholds[,inds_col]<-thresholds_list_blast
  # export filter list
  filter_blast<-sapply(1:length(list_blast), function(i){
    filter_i<-list_blast[[i]]@filter
    return(filter_i)
  })
  names(filter_blast)<-sampleNames(gs)
  saveRDS(filter_blast,"/home/rstudio/results/flowType_files/filter_objects/filter_blast.rds")
  if(export_CLR==T){
    generate_CLR_data(fs.live,list(list_blast),name_pop = "CD34+SSCA-",channels=c("PE-A","SSC-A"),filters=F,density = density_axis,path.output=path_output_clr)
  }
  #--- generate final data
  if(flowjo_wsp==T){
    print("---- generate wsp")
    generate_final_data(path.output,name_output = name_output,gs)
  }
  #----- generate plots
  if(plots==T){
    print("---- generate plots")
    generate_plots_Bcells(gs,n_cores,path.output,fs.marg,clean.inds,fs.clean,singlets,fs.sngl,beads_poly,size_poly,
                          fs.size,live_poly,fs.live,gran_poly,non_gran_poly,fs.non_gran,Bcells_poly,
                          fs.bcells,pc_poly,blast_poly,CD20pos_poly,CD20neg_poly,quad1,PB_poly,
                          fs.CD20pos,CD10neg_poly,CD10pos_poly,fs.igm,trans_poly,fs.cd10neg,quad.3,
                          fluorochrome.chans,scat.chans)
  }
  if(export_gs==T){
    print("---- export gs")
    setwd("/home/rstudio/results/shiny_app_files/GatingSet")
    outFile<-"GatingSet.wsp"
    n_samples<-length(gs)
    for(i in 1:n_samples){
      gs[[i]]@data[[i]]@description$`$FIL`<<- sampleNames(gs)[i]
    }
    #gatingset_to_flowjo(gs, outFile=outFile)
    save_gs(gs,path = "/home/rstudio/results/shiny_app_files/GatingSet/GatingSet_files",overwrite = T)
  }
  if(export_thresholds==T){
    print("---- export thresholds")
    write.csv(df_thresholds,"/home/rstudio/results/flowType_files/df_thresholds/df_thresholds.csv")
  }
  return(gs)
}


