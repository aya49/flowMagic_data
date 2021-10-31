########################## import libraries ###########################
library("flowCore")
library("flowDensity")
library("ggcyto")
library(flowCut)
library(flowWorkspace)
library('MASS')
library(sp) # to make SpatialPolygon object
library(rgeos) # to use the function gIntersect,gDifference ecc..,for polygon comparison
library('flowPeaks') # Note : you need the library libgsl-dev
library(stringr)
#library(cytoUtils)
library(CytoML) # to use GatingSetFlowjo
library(data.table) # to use rbindlist
library(reshape2) # to use facet_wrap
library(plyr) # to use revalue
library(readxl)
library(parallel) # to enable parallization
library(RColorBrewer)

#########################################################################################################
########################################### Definitions of  functions ###################################
#########################################################################################################

#------------------------ function import fcs files -----
# fcs file that contains cleaned live cells (flowType input)

import_fcs_file<-function(path_fcs_file){
  f<-read.FCS(path_fcs_file)
  return(f)
}


#-------  function to execute gating of flowtype phenotypes --------------
# f= flowFrame
gating_flowtype_pheno<-function(f,channel_1="None",channel_2="None",phenotype_pheno="None",
                                gate_channel_1="None",gate_channel_2="None",name_gated_pop="None",
                                position_ch_1="pos",position_ch_2="pos",filter_pop=F,filter_object="none",
                                min_y=1,max_x=2,min_x){
  if(filter_pop==F){
    # we take the gates from the df_thresholds file
    gate_channel_1 <- gate_channel_1
    gate_channel_2 <- gate_channel_2
    if(is.na(gate_channel_1)==T){
      position_1<-NA
    }else{
      if(position_ch_1=="pos"){
        position_1<-T
      }else{
        position_1<-F
      }
    }
    if(is.na(gate_channel_2)==T){
      position_2<-NA
    }else{
      if(position_ch_2=="pos"){
        position_2<-T
      }else{
        position_2<-F
      }
    }
    
    # we gate the population
    print(channel_1)
    print(channel_2)
    print(position_1)
    print(position_2)
    print(gate_channel_1)
    print(gate_channel_2)
    cell_pop <- flowDensity(f, channels = c(channel_1, channel_2), position = c(position_1,position_2), gates = c(gate_channel_1, gate_channel_2))
  }else if(filter_pop==T){
    cell_pop <- flowDensity(f,channels = c(channel_1, channel_2), position = c(F,F),filter=filter_object)
    
  }
  par(mar=c(5,5,5,5))
  plotDens(f, channels = c(channel_1, channel_2),  main = phenotype_pheno, cex.lab = 3, cex.axis = 1, cex.main=3,
           xlim = c(min_x,max_x),ylim=c(min_y,4))

  lines(cell_pop@filter,lwd=2)
  text(mean(cell_pop@filter[,1]), mean(cell_pop@filter[,2]), labels = name_gated_pop, cex = 2)
  f_gated<-getflowFrame(cell_pop)
  return(f_gated)
}








#########################################################################################################
########################################### execute functions ###########################################
#########################################################################################################

########### myeloid  panel################

#-------- positive group ----------------

f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_4/fcs_files/FCT_G181C_3U65.fcs")
#f@parameters@data$desc[c(15,19,12,7)]<-c("CD66","CD45","CD11b","CD11c")

# first gate
#  it is a filter:we load the filter object from the flow type folder
filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_4/filter_objects/filter_CD66negCD45pos.rds")
filter_df<-filter_list[[79]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])

f_non_gran<-gating_flowtype_pheno(f,channel_1="BV711-A",channel_2="V450-A",phenotype_pheno="All cells",
                                   gate_channel_1=channel_1_gate_value,gate_channel_2 = channel_2_gate_value,
                                   position_ch_1="neg",position_ch_2="pos",name_gated_pop="CD66-CD45+(NonGran)",
                                   filter_pop = T,filter_object = filter_df,max_x = 4,min_y = -1,min_x = -1)
# second gate
filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_4/filter_objects/filter_monocytes.rds")
filter_df<-filter_list[[79]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])

f_mono<-gating_flowtype_pheno(f_non_gran,channel_1="eF605-A",channel_2="V500-A",phenotype_pheno="CD66-CD45+(NonGran)",
                                   gate_channel_1=channel_1_gate_value,gate_channel_2 = channel_2_gate_value,
                                    position_ch_1="neg",position_ch_2="pos",
                                   name_gated_pop="HLADR-CD14+(Monocytes)",filter_pop=T,filter_object = filter_df,
                              max_x = 4,min_y = -1,min_x = -1)

# third gate

f_CD11b_pos<-gating_flowtype_pheno(f_mono,channel_1="BV786-A",channel_2="APC-A",phenotype_pheno="CD66-CD45+HLADR-CD14+",
                                   gate_channel_1=2.07,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                   name_gated_pop="CD11b+",max_x = 3.5,min_y = -1,min_x=-0.5)
#-------- negative group ------------------
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_6/fcs_files/FCT_G264A_5K58.fcs")
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_3/fcs_files/FCT_G124C_7A78.fcs")

# first gate
#  it is a filter:we load the filter object from the flow type folder
filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_6/filter_objects/filter_CD66negCD45pos.rds")
filter_df<-filter_list[[52]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])

filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_3/filter_objects/filter_CD66negCD45pos.rds")
filter_df<-filter_list[[58]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])


f_non_gran<-gating_flowtype_pheno(f,channel_1="BV711-A",channel_2="V450-A",phenotype_pheno="All cells",
                                  gate_channel_1=channel_1_gate_value,gate_channel_2 = channel_2_gate_value,
                                  position_ch_1="neg",position_ch_2="pos",name_gated_pop="CD66-CD45+(NonGran)",
                                  filter_pop = T,filter_object = filter_df,min_y = -1,max_x=4,min_x = -1)
# second gate
filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_4/filter_objects/filter_monocytes.rds")
filter_df<-filter_list[[52]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])
filter_list<-readRDS("/home/rstudio/data/flowType_data/Myeloid_panel/Batch_2/Plate_3/filter_objects/filter_monocytes.rds")
filter_df<-filter_list[[58]]
channel_1_gate_value<-max(filter_df[,1])
channel_2_gate_value<-max(filter_df[,2])

f_mono<-gating_flowtype_pheno(f_non_gran,channel_1="eF605-A",channel_2="V500-A",phenotype_pheno="CD66-CD45+(NonGran)",
                              gate_channel_1=channel_1_gate_value,gate_channel_2 = channel_2_gate_value,
                              position_ch_1="neg",position_ch_2="pos",
                              name_gated_pop="HLADR-CD14+(Monocytes)",filter_pop=T,filter_object = filter_df,
                              min_y = -1,max_x = 4,min_x = -1)

# third gate

f_CD11b_pos<-gating_flowtype_pheno(f_mono,channel_1="BV786-A",channel_2="APC-A",phenotype_pheno="CD66-CD45+HLADR-CD14+",
                                   gate_channel_1=2.39,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                   name_gated_pop="CD11b+",min_y = -1,max_x = 3.5,min_x = -0.5)



########### Bcell panel ################
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_4/fcs_files/FCT_G688E_5U39.fcs")
#f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_1/fcs_files/FCT_G530H_8K01.fcs")
#f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_4/Plate_1/fcs_files/FCT_G542C_5S57.fcs")
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_2/Plate_4/fcs_files/FCT_G152D_2S62.fcs")

#---------- positive sample ---------

# first gate

f_IgMpos<-gating_flowtype_pheno(f,channel_1="FITC-A",channel_2="PE-Cy5-A",phenotype_pheno="All cells",
                                 gate_channel_1=2.39,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="IgM+",max_x = 4,min_y = -0.2,min_x = 0)

# second gate

f_CD38pos<-gating_flowtype_pheno(f_IgMpos,channel_1="PE-Cy5-A",channel_2="FITC-A",phenotype_pheno="IgM+",
                                 gate_channel_1=1.61,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="CD38+",max_x = 4,min_y = 2,min_x = 0)

# third gate 
f_138pos<-gating_flowtype_pheno(f_CD38pos,channel_1="V450-A",channel_2="FITC-A",phenotype_pheno="IgM+CD38+",
                                 gate_channel_1=1.67,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="CD138+",max_x = 4,min_y = 2,min_x = 0)

#---------- negative sample  ----------
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_5/plate_5/fcs_files/FCT_G483J_7F23.fcs")
f<-import_fcs_file(path_fcs_file="/home/rstudio/data/flowType_data/Bcell_panel/Batch_3/Plate_4/fcs_files/FCT_G458E_2R11.fcs")

# first gate

f_IgMpos<-gating_flowtype_pheno(f,channel_1="FITC-A",channel_2="PE-Cy5-A",phenotype_pheno="All cells",
                                gate_channel_1=2.28,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                name_gated_pop="IgM+",max_x = 4,min_y = 0,min_x = 0)

# second gate

f_CD38pos<-gating_flowtype_pheno(f_IgMpos,channel_1="PE-Cy5-A",channel_2="FITC-A",phenotype_pheno="IgM+",
                                 gate_channel_1=1.59,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                 name_gated_pop="CD38+",max_x = 4,min_y = 2,min_x = 0)

# third gate 
f_138pos<-gating_flowtype_pheno(f_CD38pos,channel_1="V450-A",channel_2="FITC-A",phenotype_pheno="IgM+CD38+",
                                gate_channel_1=0.88,gate_channel_2 = NA,position_ch_1="pos",position_ch_2="none",
                                name_gated_pop="CD138+",max_x = 4,min_y = 2,min_x = -0.4)

