

### python_script_example.py

####### ----------------------------------------------------------

### Let's import some libraries
### ***
### ***
### ***

### ***

import numpy as np
import pandas as pd

import matplotlib.path as mplPath


from sklearn import metrics

import scipy.sparse

import time
import random
import math
import sys
import os

### For matching file names
import fnmatch

### For copying the files
import shutil


### To read MMOS outputs
import json

### F1 score
from sklearn.metrics import f1_score

import argparse





####### ----------------------------------------------------------

print("Code started ... ")

####### ----------------------------------------------------------

########################### Functions

####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------



### ***



def func_sim_matrix_pixel(events_label_df, user_id, consider_label_zero=True):
    '''
    This function finds similarity matrix (co-occurance matrix) for each user clustering
    '''
    df_l_t = events_label_df
    usr_id = user_id

    sim_m = np.zeros((len(df_l_t["pixel_label"]), len(df_l_t["pixel_label"]) ), dtype=np.uint8)

    l_label =  df_l_t["pixel_label"].unique()

    for l_l in l_label:
        print(l_l, end="\r")
        if l_l ==0:
            if  consider_label_zero==True:
                for i in df_l_t[df_l_t["pixel_label"]==l_l].index:

                    for j in df_l_t[df_l_t["pixel_label"]==l_l].index:
                        if j >= i:
                            sim_m[i, j]= int(1)
                            sim_m[j, i]= int(1)
        else:
            for i in df_l_t[df_l_t["pixel_label"]==l_l].index:


                for j in df_l_t[df_l_t["pixel_label"]==l_l].index:
                    if j >= i:
                        sim_m[i, j]= int(1)
                        sim_m[j, i]= int(1)



### Uncomment below if you want to save output
#    sparse_matrix = scipy.sparse.csc_matrix(sim_m)
#    scipy.sparse.save_npz('similarity_matrix_pixel_usr'+str(usr_id)+'.npz', sparse_matrix)

    return sim_m



####### ----------------------------------------------------------


def func_pixel_coor(res_x_FCS, res_y_FCS, x_events, y_events):
    '''
    This function recieves the resolution of pixels (number of pixels in each direction), and events x and y,
    and calculates the centre of each pixels, and form a meshgrid based on them.
    
    It returns the info below:
    
    pixel centre: x_pix_c ,    y_pix_c   -  their sizes are len(res_x_FCS) and len(res_y_FCS)
    pixel edges:  x_pix_edges, y_pix_edges   -  their sizes are len(res_x_FCS)+1 and len(res_y_FCS)+1

    pixel meshgrid: x_pixel_coor, y_pixel_coor, Grid_pixel_coor    - their size are len(res_x_FCS)*len(res_y_FCS)
    
    '''
        
    res_x_FCS = res_x_FCS
    res_y_FCS = res_y_FCS
    x_events = x_events
    y_events = y_events
    
    # make the heatmap of events to find edges

    heatmap_events, xedges, yedges = np.histogram2d(x_events, y_events, bins=(res_x_FCS, res_y_FCS))

    
    
    # Here, we define the centre of each bin/pixel using its edges information from heatmap_events.
    # let's find x and y of the centre of each pixel by averaging the edges:

    x_pix_c =[]
    x_pix_edges =[]

    for i in np.arange(len(xedges)):
        if i>0:
            x_a = (xedges[i]+xedges[i-1])/2
            x_pix_c.append(x_a)

            x_pix_edges.append([xedges[i-1],xedges[i]])


    y_pix_c =[]
    y_pix_edges =[]

    for j in np.arange(len(yedges)):
        if j>0:
            y_a = (yedges[j]+yedges[j-1])/2
            y_pix_c.append(y_a)

            y_pix_edges.append([yedges[j-1],yedges[j]])
            
            
            
    # Let's form a grid coordinates using the centre of bins/pixels from heatmap_events:


    x_pixel_coor, y_pixel_coor = np.meshgrid(x_pix_c, y_pix_c) # make a canvas with coordinates
    x_pixel_coor, y_pixel_coor = x_pixel_coor.flatten(), y_pixel_coor.flatten()

    Grid_pixel_coor = np.vstack((x_pixel_coor,y_pixel_coor)).T
    
    
    return x_pix_c, y_pix_c, x_pix_edges, y_pix_edges, x_pixel_coor, y_pixel_coor, Grid_pixel_coor




####### ----------------------------------------------------------



def func_df_pixels_events_info(res_x_FCS, res_y_FCS, x_events, y_events):
    '''
    This function returns a dataframe containing information about pixels and index of events inside each pixel etc.
    
    It calls the "func_pixel_coor()" to calculate pixels.
     
    '''
    
    res_x_FCS = res_x_FCS
    res_y_FCS = res_y_FCS
    x_events = x_events
    y_events = y_events
    
    
    ### Call "func_pixel_coor()" function to find pixel info. Otherwise one need to calculate this beforehand
    ### and give the variables as input
    x_pix_c, y_pix_c, x_pix_edges, y_pix_edges, x_pixel_coor, y_pixel_coor, Grid_pixel_coor = \
    func_pixel_coor(res_x_FCS, res_y_FCS, x_events, y_events)
    
    
    heatmap_events, xedges, yedges = np.histogram2d(x_events, y_events, bins=(res_x_FCS, res_y_FCS))

    
    ### Define an initial dataframe for pixel info:
    df_pixel_info = pd.DataFrame({"x_pixel_coor":x_pixel_coor, "y_pixel_coor":y_pixel_coor, "events_density":heatmap_events.T.flatten() })

    
    ### Define dataframe for pixel information:
    df_pixel_edges_info = pd.DataFrame({"x_pixel_coor":x_pix_c, "y_pixel_coor":y_pix_c, "x_pix_edges":x_pix_edges, "y_pix_edges":y_pix_edges })



    
    # To be able to follow pixels, we form a meshgrid with "index" values of pixels:
    ### Make meshgrid with "index" of x and y for pixels

    x_p_index, y_p_index =  np.meshgrid(np.arange(len(x_pix_c) ), np.arange(len(y_pix_c) ) )
    x_p_index, y_p_index = x_p_index.flatten(), y_p_index.flatten()

    Grid_p_index = np.vstack((x_p_index, y_p_index)).T


    x_edges_pixel_list = []
    y_edges_pixel_list = []

    events_index_pixel = []

    events_points = np.vstack((x_events,y_events)).T

    df_e = pd.DataFrame({"x":events_points[:,0], "y":events_points[:,1]})

    events_index = list(np.arange(len(events_points)))

    for i_G in np.arange(len(Grid_p_index)):
    #    print(i_G, end="\r")

        [x_p_edge_1, x_p_edge_2] = df_pixel_edges_info["x_pix_edges"][Grid_p_index[i_G][0]]
        [y_p_edge_1, y_p_edge_2] = df_pixel_edges_info["y_pix_edges"][Grid_p_index[i_G][1]]


        x_p_e = [x_p_edge_1, x_p_edge_1, x_p_edge_2, x_p_edge_2, x_p_edge_1]
        y_p_e = [y_p_edge_1, y_p_edge_2, y_p_edge_1, y_p_edge_2, y_p_edge_1]

        pixel_vertices = np.vstack((x_p_e,y_p_e)).T



        e_i_p_temp = []

        min_x_edge = min(x_p_edge_1, x_p_edge_2)
        max_x_edge = max(x_p_edge_1, x_p_edge_2)

        min_y_edge = min(y_p_edge_1, y_p_edge_2)
        max_y_edge = max(y_p_edge_1, y_p_edge_2)


        if df_pixel_info["events_density"][i_G]>0:

            events_check_p = (df_e["x"]>= min_x_edge) & (df_e["x"]<= max_x_edge) & (df_e["y"]>= min_y_edge) & (df_e["y"]<= max_y_edge)

            e_i_p_temp = np.array(x_events.index[events_check_p] )




        events_index_pixel.append(e_i_p_temp)
        x_edges_pixel_list.append([x_p_edge_1, x_p_edge_2])
        y_edges_pixel_list.append([y_p_edge_1, y_p_edge_2])

        
        
    # Let's undate the dataframe for pixel info:
    
    df_pixel_info["events_index_pixel"] = events_index_pixel
    df_pixel_info["x_pix_edges"] = x_edges_pixel_list
    df_pixel_info["y_pix_edges"] = y_edges_pixel_list

    
    return df_pixel_info



####### ----------------------------------------------------------

### ***

def extract_gates(user_all_gates, u_id):
    '''
    This function recieves vertices from all gates of a user (with user id of u_id),
    and extarct each gates x and y coordinates.
    
    user_all_gates is a dataframe that contains "x", "y", and "gate" (e.x. G1, G2, ...) info from a user.
    
    output is a dataframe thatn contains "gate_vertices" and "gate_name" (i.e. user_id + gate info: e.x. 1_G1, 1_G2 ...)
    
    '''
    
    user_temp = user_all_gates
    gate_names_unique = user_temp["gate"].unique()
    gate_list = []
    gate_list_name = []
    for i in np.arange(len(gate_names_unique)):
        gate_t = user_temp[user_temp["gate"]==gate_names_unique[i]]

# ### do not repeat the last point
#        gate_list.append(gate_t[['x','y']][:-1])

        gate_list.append(gate_t[['x','y']])
        gate_list_name.append(str(u_id)+"_"+str(gate_names_unique[i]))
        
    output = pd.DataFrame({"gate_vertices":gate_list, "gate_name": gate_list_name})
    
    return output
  

####### ----------------------------------------------------------

### ***

def func_label_events(user_gates_info, x_events, y_events):
    
    '''
    This function labels events for a user considering gates.
    
    Info pixels information is given rather than  x_events, y_events, it will label pixels (i.e. pixels are events then).
    
    '''
    user_gates_temp = user_gates_info
    
    x_events = x_events
    y_events = y_events

    events_l_temp = np.zeros(len(x_events))

    l_n = 0

    for i_gn in np.arange(len(user_gates_temp["gate_name"])):

        l_n = l_n+1

        # Extract the point values that define the perimeter of the polygon
        x_pol = user_gates_temp["gate_vertices"][i_gn]["x"]
        y_pol = user_gates_temp["gate_vertices"][i_gn]["y"]

        poly_vertices = np.vstack((x_pol,y_pol)).T

        events_points = np.vstack((x_events,y_events)).T

        # Check to see which points are inside the polygon
        bbPath = mplPath.Path(poly_vertices)
        events_check  = bbPath.contains_points(events_points)

        events_l_temp = events_l_temp + np.array(events_check).astype(int) * l_n
        
        
    return events_l_temp




####### ----------------------------------------------------------

### This is a function to get gates for each user from MMOS output, and scale them back.
### One can turn off scaling by setting "do_scaling==False"

def func_get_gates_user(df_results_MMOS, f_p_minmax, do_scaling=True):


    df_results_MMOS = df_results_MMOS

    ### scaling back info f_p_minmax should be [minX, minY, maxX, maxY]
    f_p_minmax = f_p_minmax
    
    if do_scaling == True:

        min_x_s = f_p_minmax[0]
        min_y_s = f_p_minmax[1]

        max_x_s = f_p_minmax[2]
        max_y_s = f_p_minmax[3]

    if do_scaling == False:

        min_x_s = 0
        min_y_s = 0

        max_x_s = 1
        max_y_s = 1

        
        
    n_pol_user_f_p = df_results_MMOS.shape[0]

    x_u_f_p = []
    y_u_f_p = []

    g_u_f_p = []

    for i_pol in np.arange(n_pol_user_f_p):


        df_results_MMOS_v = pd.DataFrame(df_results_MMOS["polygons"][i_pol])

        g_name_temp = "G"+str(int(i_pol+1) )

        for i in np.arange(df_results_MMOS_v["vertices"].shape[0]):
            
            x_scale_back_temp = df_results_MMOS_v["vertices"][i][0]*(max_x_s-min_x_s) + min_x_s
            y_scale_back_temp = df_results_MMOS_v["vertices"][i][1]*(max_y_s-min_y_s) + min_y_s
            
            ### Keep track of first point to close polygon
            first_x_u_f_p = df_results_MMOS_v["vertices"][0][0]*(max_x_s-min_x_s) + min_x_s
            first_y_u_f_p = df_results_MMOS_v["vertices"][0][1]*(max_y_s-min_y_s) + min_y_s
            first_g_u_f_p = g_name_temp

            x_u_f_p.append(x_scale_back_temp)
            y_u_f_p.append(y_scale_back_temp)
            g_u_f_p.append(g_name_temp)

        x_u_f_p.append(first_x_u_f_p)
        y_u_f_p.append(first_y_u_f_p)
        g_u_f_p.append(first_g_u_f_p)



    df_as_user_temp = pd.DataFrame({"x":x_u_f_p, "y":y_u_f_p, "gate":g_u_f_p   })
    
    return df_as_user_temp



####### ----------------------------------------------------------



####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------






####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### MAIN CODE   ------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------

### Read FCM Export Files

### Define the name of folders

folder_of_FCM_Export_Files = "/data/user_data_by_ccp"

folder_of_gold_standards = "/data/golden_samples"

folder_of_gold_standards_events = "/data/data"

folder_of_gold_standards_list = "/code/support_files"
folder_for_outputs = "/data/users_performance_outputs"

### This defines the type of average for F1-measure
average_type='micro'

print(sys.argv[1])
file_to_read = int(sys.argv[1])

####### MAKE SURE THIS IS ADJUSTED - See the main loop for the place we should change
###gold_standard_labels_file_name = "labels_"+gold_standard_events_file_name


####### ----------------------------------------------------------

### Here we make a list of FCM Export Files

list_of_FCM_Export_Files = []
full_list = []
### Find the current directory where the script is running.
### Other folders containing users and events info should be under this directory.

for file in os.listdir(folder_of_FCM_Export_Files):
    full_list.append(file)


for file in os.listdir(folder_of_FCM_Export_Files):
    if fnmatch.fnmatch(file, full_list[file_to_read]):
         list_of_FCM_Export_Files.append(file)

print(list_of_FCM_Export_Files[0])
####### ----------------------------------------------------------



### Let's read the list of gold standards files

### "golds.csv" is the actual list. "gold_standards_list.csv" is for testing
#gs_list_file_name = "golds.csv"
gs_list_file_name = "goldstandard_CCP.csv"


### the path includes file name too
path_gold_standards_list = folder_of_gold_standards_list+"/"+gs_list_file_name

df_gold_standards_list = pd.read_csv(path_gold_standards_list, engine='python')
print(df_gold_standards_list)
### If there is no header read below
# df_gold_standards_list = pd.read_csv(path_gold_standards_list, header=None)

column_name_df_gold_standards_list = df_gold_standards_list.columns[0]




####### ----------------------------------------------------------



### Here, we use the list of FCM Export Files, go to each one,
### find the user for gold reference and do calculation
### and save an output for each batch of FCM Export Files


for fcm_file_temp in list_of_FCM_Export_Files:
    FCM_Export_Files_name_temp = fcm_file_temp

##############################---------------------------------------------------
    
    ### Let's read the FCM Export File

    ### the path includes file name too
    path_FCM_Export_Files = folder_of_FCM_Export_Files+"/"+FCM_Export_Files_name_temp

    FCM_export_temp = pd.read_csv(path_FCM_Export_Files, names=["filepath","score","result","minmax", "player", "created_at"], header=0, sep="','", dtype={"score" :np.float64})
    FCM_export_temp["filepath"] = FCM_export_temp["filepath"].str[1:]
    ### Now this should be applied to "player" columns
    # FCM_export_temp["minmax"] = FCM_export_temp["minmax"].str[:-1]
    FCM_export_temp["result"] = FCM_export_temp["result"].apply(json.loads)
    FCM_export_temp["minmax"] = FCM_export_temp["minmax"].apply(json.loads)
    # ### NEW for player - UserID
    # FCM_export_temp["player"] = FCM_export_temp["player"].str[:-1]
    FCM_export_temp["player"] = FCM_export_temp["player"].astype(int)
    ### NEW 2 - After adding time : "created_at"
    FCM_export_temp["created_at"] = FCM_export_temp["created_at"].str[:-1]
    ### Convert string to timestamp for "created_at"
    # FCM_export_temp["created_at"] =pd.to_datetime( FCM_export_temp["created_at"])

    
    
##############################---------------------------------------------------


    ### Below the loop for each batch of FCM Export Files starts ....
    

    df_users_performance = pd.DataFrame(columns= ["player", "filepath", "our_score_F1", "our_score_ARI", "MMOS_score", "created_at" ])

    ### Let's check if a "filepath" corresponds to a gold standards

    for i_f_p in np.arange(FCM_export_temp.shape[0]):

        if (i_f_p%10000) == 0:
            print("i_f_p: ", i_f_p)

        ### Here we get the plot name to be analyzed
        plot_name_to_be_analyzed = FCM_export_temp["filepath"][i_f_p]

        gold_standard_events_file_name = plot_name_to_be_analyzed

    ##############################---------------------------------------------------
        ### FCM_export_temp_one_user can come here
        # Let's get the FCM export data for one row (i.e. one user)

        FCM_export_temp_one_user = pd.DataFrame()
        FCM_export_temp_one_user = FCM_export_temp_one_user.append(FCM_export_temp.iloc[i_f_p])
        FCM_export_temp_one_user.reset_index(drop=True, inplace=True)


    ##############################---------------------------------------------------

        ### Here is some relevant info to be saved
        player_to_save = FCM_export_temp_one_user["player"][0].astype(int)
        filepath_to_save = FCM_export_temp_one_user["filepath"][0]
        ### our_score_to_save = TO BE DEFINED AT THE END
        MMOS_score_to_save = FCM_export_temp_one_user["score"][0]
        created_at_to_save = FCM_export_temp_one_user["created_at"][0]


    ##############################---------------------------------------------------
        ### If this is true, the code will go ahead to calculate "our_score_to_save"
        if plot_name_to_be_analyzed in list(df_gold_standards_list[column_name_df_gold_standards_list]):

    ##############################---------------------------------------------------

            ### Let's get the information of gates etc. for this user

            ### DataFrame to keep gates
            gate_info_all = pd.DataFrame( columns=["user_id", "gates"])

            ### define an array of user ids
            ### We only alayze one user at a time here
            only_one_user = 1
            user_id = np.arange(1, only_one_user+1)

            for u_i in user_id:

                ### results for a particular user
                results_FCM_export_temp = FCM_export_temp_one_user["result"][0]
                df_results_FCM_export_temp = pd.DataFrame(results_FCM_export_temp)
                ### min max for scaling vertices
                #f_p_minmax = FCM_export_temp["minmax"][0]
                #BELOW IS THE CORRECT min max
                f_p_minmax = FCM_export_temp_one_user["minmax"][0]

                gates_u_i = func_get_gates_user(df_results_FCM_export_temp, f_p_minmax, do_scaling= True)
                gate_info_all = gate_info_all.append({"user_id":u_i, "gates": gates_u_i}, ignore_index=True)

    ##############################---------------------------------------------------

            ### Let's read the corresponding events for this "filepath"
            ### ***
            events_csv = pd.read_csv(folder_of_gold_standards_events+'/'+gold_standard_events_file_name)

            ### ***

            x_events = events_csv.iloc[:,0]
            y_events = events_csv.iloc[:,1]


    ##############################---------------------------------------------------

            ### Now, let's get the events label for this user
            ### ***

            df_user_gates_label_events = pd.DataFrame(columns = ["user_id", "events_label"])

            for u_i in user_id:

                index_u_i = gate_info_all[gate_info_all["user_id"]==u_i].index[0]
                user_temp = gate_info_all["gates"][index_u_i]
                user_gates_temp = extract_gates(user_temp, u_i)

                label_temp = func_label_events(user_gates_temp, x_events, y_events)

                df_user_gates_label_events = df_user_gates_label_events.append({"user_id":u_i, "events_label": label_temp}, ignore_index=True)


    ##############################---------------------------------------------------

            ### Let's read the gold standard

            #gold_standard_events_file_name = plot_name_to_be_analyzed
            #gold_standard_labels_file_name = '265-00 CNTRP_Phase2 ONE_01 YVR_RG 2016-01-15 101.csv'

            gold_standard_labels_file_name = gold_standard_events_file_name
            print(gold_standard_labels_file_name)

            df_gold_standard_labels = pd.read_csv(folder_of_gold_standards+'/'+gold_standard_labels_file_name)


            ### Let's get the columns names and replace them with numbers

            list_of_gsl_columns = list(df_gold_standard_labels.columns.values)

            columns_name_numbers = np.arange(1, len(list_of_gsl_columns)+1 )

            for i in np.arange(len(list_of_gsl_columns)):
                df_gold_standard_labels.rename(columns= {list_of_gsl_columns[i]:columns_name_numbers[i]}, inplace= True)


            ### Let's convert True False to numbers of new columns names
            list_of_gsl_columns_new = list(df_gold_standard_labels.columns.values)

            for i in list_of_gsl_columns_new:
                df_gold_standard_labels[i] = df_gold_standard_labels[i].astype(int)*i


            # Let's make a new data frame with only one column "labels" that contains the corresponding label numbers

            df_gold_standard_labels_numbers = pd.DataFrame(columns= ["labels"])
            df_gold_standard_labels_numbers["labels"] = df_gold_standard_labels.sum(axis=1)


    ##############################---------------------------------------------------

            ### Let's calculate the F1 Score

            true_gold_standard = list(df_gold_standard_labels_numbers["labels"])
            pred_user = list(df_user_gates_label_events["events_label"][0])


            ### Define what type of averaging for F-measure:

            ### average_type='micro'
            F1_measure_gold = f1_score(true_gold_standard, pred_user, average=average_type)

            
            
    ##############################---------------------------------------------------

    
            ### Let's also calculated adjusted Rand Index
            ARI_gold = metrics.adjusted_rand_score(true_gold_standard, pred_user)

    
    ##############################---------------------------------------------------

    
    
            
    ##############################---------------------------------------------------


            ### This is for the data frame to be saved.
            our_score_to_save_1 = F1_measure_gold
            our_score_to_save_2 = ARI_gold

            ### Let's add the new information here.
            df_users_performance = df_users_performance.append({"player":player_to_save , "filepath":filepath_to_save , "our_score_F1":our_score_to_save_1 , "our_score_ARI":our_score_to_save_2 , "MMOS_score":MMOS_score_to_save, "created_at":created_at_to_save }, ignore_index=True)


    ##############################---------------------------------------------------

    ### Let's save the data for this batch of FCM Export Files
    df_users_performance.to_csv(folder_for_outputs+"/"+FCM_Export_Files_name_temp+'_users_performace.csv', index=False)

    ##############################---------------------------------------------------


    ##############################---------------------------------------------------


    ##############################---------------------------------------------------


    ##############################---------------------------------------------------


    ##############################---------------------------------------------------


    ##############################---------------------------------------------------
    ##############################---------------------------------------------------
    ##############################---------------------------------------------------










####### ----------------------------------------------------------

print("Code ended ...")


####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------













