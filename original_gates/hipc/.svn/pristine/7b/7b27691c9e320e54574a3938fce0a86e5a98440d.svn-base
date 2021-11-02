########################## importing libraries and functions ########################
library("flowCore")
library("flowDensity")
library("ggcyto")
library("devtools") # you need devtools to use the function install_github that you need to install flowCut
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
library(MASS)
source("/home/rstudio/Code_Bcells_Myeloid_cells/main_functions.R")
source("/home/rstudio/Code_Bcells_Myeloid_cells/Utils.R")
source("/home/rstudio/Code_Bcells_Myeloid_cells/Gating_functions_myeloid.R")
source("/home/rstudio/Code_Bcells_Myeloid_cells/flowPrep.R")

####################################### Execute functions ##################################
# Execute gating myeloid panel
main_myeloid("All",path_comp_matrix ="/home/rstudio/data/Renamed_Validation_data_Gambia_PNG/Flow_cytometry_data/Myeloid_panel/Compensation_Matrix/Comp_matrix.csv",
     path_dir_files = "/home/rstudio/data/Renamed_Validation_data_Gambia_PNG/Flow_cytometry_data/Myeloid_panel/Plate_2",
     flowjo_wsp=T,
     name_output="All",
     n_cores=1,
     pre_processing = T,
     plots = T,
     export_gs = F, export_fs_clean = F,export_fcs_live_cells = T,export_thresholds = T,export_CLR = F,
     path_output_clr = "/home/rstudio/results/CLR_result")



main_myeloid(87,path_comp_matrix ="/home/rstudio/data/Gambia_samples_myeloid_final/Gambia_1st_Batch/Myeloid_comp.csv",
     path_dir_files = "/home/rstudio/results/shiny_app_files/Cleaned_files",
     flowjo_wsp=F,
     name_output="All",
     n_cores=1,
     pre_processing = F,
     plots = F,
     export_gs = F)


list.files("/home/rstudio/data/Validation_data_Gambia_PNG/Flow_cytometry_data/Myeloid_panel/Plate_2",recursive = T,full.names = T)

# checking results gating 

check_final_results(files_path ="/home/rstudio/results/Results_Renamed_validation_samples_Myeloid/Plate_2/Stats", check_Rym_data = F,
                    control="None")

# Generating final document description file
generate_final_description_file(files_path ="/home/rstudio/results/Results_Renamed_validation_samples_Myeloid/Plate_1/Stats",
                                path_Rym_data = "/home/rstudio/data/Validation_data_Gambia_PNG/Plate_maps/Plate_1_map.csv",
                                files_path_all_plates = "/home/rstudio/results/Results_Renamed_validation_samples_Myeloid",
                                flag_control_samples = "E10")

# note: the reference tresholds of 3rd and 4th criteria 
# of validation cohort are the same of main cohort
#--- main cohort
# we esclude the controls from the following calculations:
# total number of samples Batch 1: 41 no ctrl (41 + 1 no_info D7)
# total number of samples Batch 2: 48,96,95,95,95 ---> 429 no ctrl
# total number of samples Batch 3: 96,96,96,96,96,2 --> 482 no ctrl
# total number of samples Batch 4: 96,95,96,80(79 + 1 no_info G9),94 ---> 460 no ctrl
# 429+482+41+460 = 1412
# D7 and G9 don't exist

#Generating plot comparison control samples
final_sample_control(files_path = "/home/rstudio/results/Results_validation_samples_myeloid/",comparison_plots = T,legend = F,cleaning_v = T,remove_pops = "Singlets|All cells|Beads")



