
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
source("/home/rstudio/Code_Bcells_Myeloid_cells/Gating_functions_bcells.R")
source("/home/rstudio/Code_Bcells_Myeloid_cells/flowPrep.R")

####################################### Execute functions ##################################

# Executing gating bcells panel
main_bcells("All",path_comp_matrix ="/home/rstudio/data/Renamed_Validation_data_Gambia_PNG/Flow_cytometry_data/Bcells_panel/Compensation_matrix/Comp_matrix.csv",
     path_dir_files = "/home/rstudio/data/Renamed_Validation_data_Gambia_PNG/Flow_cytometry_data/Bcells_panel/Plate_2",
     path_fcs_files_fixing="/home/rstudio/data/PNG B cell panel (Extended)/PNG B cells expended panel 01032018/",
     flowjo_wsp=T,
     name_output="All",
     n_cores=1,
     pre_processing = T,
     plots = T,
     export_gs = F, export_fs_clean = F,export_fcs_live_cells = T,export_thresholds = T,export_CLR = F,
     path_output_clr = "/home/rstudio/results/CLR_result")

# D3 to be excluded
# "/home/rstudio/data/PNG B cell panel (Extended)/PNG B cells expended panel 01032018/"
# only on cleaned files
main_bcells("All",path_comp_matrix ="/home/rstudio/data/Gambia_Bcells_final/Batch_1/Batch_1_replaced_6_samples/Compensation Matrix.csv",
     path_dir_files = "/home/rstudio/results/shiny_app_files/Cleaned_files",
     path_fcs_files_fixing="/home/rstudio/data/PNG B cell panel (Extended)/PNG B cells expended panel 01032018/",
     flowjo_wsp=F,
     name_output="All",
     n_cores=1,
     pre_processing = F,
     plots = F,
     export_gs = T)

path<-list.files("/home/rstudio/data/Validation_data_Gambia_PNG/Flow_cytometry_data/Bcells_panel/Plate_1",recursive = T,full.names = T)
f_test<-read.FCS(path[1])
#No compensation file in plate 3 !!!!!! use compensation of plate 2
#No compensation file in plate 2!!!!!! use compensation of plate 1
# Plate 2 no compensation file, I use compensation matrix of plate 1


# checking results gating 

check_final_results(files_path ="/home/rstudio/results/Results_Renamed_validation_samples_Bcells/Plate_2/Stats", 
                    check_Rym_data = F,control = "none")


# Generating final document description file 

generate_final_description_file(files_path ="/home/rstudio/results/Results_validation_samples_Bcells/Plate_1/Stats",
                                path_Rym_data = "/home/rstudio/data/Validation_data_Gambia_PNG/Plate_maps/Plate_1_map.csv",
                                files_path_all_plates = "/home/rstudio/results/Results_validation_samples_Bcells",
                                flag_control_samples = "H12")
# we esclude the controls from the following calculations:
# total number of samples Batch 1: 42 (36+5 replaced +1 rep,so it is actualy 41) no ctrl 
# total number of samples Batch 2: 48,96,95,95,95(94 rep + 1),95 no ctrl ---> 524 (or 429)
# total number of samples Batch 3: 96,96,96,96,96,2 no ctrl ----> 482
# total number of samples Batch 4: 96,95,96,79,94 no ctrl ---> 460
# 42+524+482+460 = 1508
# 41+429+482+460 = 1412

# Generating plot comparison control samples 
final_sample_control(files_path = "/home/rstudio/Results_Gambia_Samples_Bcells",comparison_plots = T,
                     legend = F,cleaning_v = T,remove_pops = "Singlets|All cells|Beads|Live cells|Non granulocytes|Granulocytes",select_batch = "All_batches")



