

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
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering

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

# For background label:

# from statistics import mode

### To read MMOS outputs
import json

### F1 score
from sklearn.metrics import f1_score


import matplotlib.pyplot as plt
import argparse

### ADDED 11 AUG 2020 - to use statistics.mode
import statistics

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

#### NEW FUNCTIONS FOR BACKGROUND PROCESS

### Functions below try to find the label that corresponds to background label.



def func_background_label_1(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c):

# Method 1: Check pixels at 4 corners, and choose the label which t least 3 of those pixels belong to.

    col_name_to_analyze = col_name_to_analyze
    df_pixel_info_example = df_pixel_info_example
    df_pixel_final_labels_hybrid_all_n_c = df_pixel_final_labels_hybrid_all_n_c
    

    max_x_pixel_coor = max(df_pixel_info_example["x_pixel_coor"])
    min_x_pixel_coor = min(df_pixel_info_example["x_pixel_coor"])

    max_y_pixel_coor = max(df_pixel_info_example["y_pixel_coor"])
    min_y_pixel_coor = min(df_pixel_info_example["y_pixel_coor"])



    index_bottom_left_corner = df_pixel_info_example.index[(df_pixel_info_example["x_pixel_coor"]== min_x_pixel_coor) & (df_pixel_info_example["y_pixel_coor"]== min_y_pixel_coor) ].tolist()
    index_bottom_right_corner = df_pixel_info_example.index[(df_pixel_info_example["x_pixel_coor"]== min_x_pixel_coor) & (df_pixel_info_example["y_pixel_coor"]== max_y_pixel_coor) ].tolist()

    index_top_left_corner = df_pixel_info_example.index[(df_pixel_info_example["x_pixel_coor"]== max_x_pixel_coor) & (df_pixel_info_example["y_pixel_coor"]== min_y_pixel_coor) ].tolist()
    index_top_right_corner = df_pixel_info_example.index[(df_pixel_info_example["x_pixel_coor"]== max_x_pixel_coor) & (df_pixel_info_example["y_pixel_coor"]== max_y_pixel_coor) ].tolist()


    label_b_l_c = df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze][index_bottom_left_corner]
    label_b_r_c = df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze][index_bottom_right_corner]

    label_t_l_c = df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze][index_top_left_corner]
    label_t_r_c = df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze][index_top_right_corner]

    list_of_corner_labels = list(label_b_l_c)+list(label_b_r_c)+list(label_t_l_c)+list(label_t_r_c)


    df_background_check_m1 = pd.DataFrame(columns= ["label", "label_repeat"])

    list_label_temp = []
    list_label_repeat_temp = []

    for l_v in set(df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze]):

        n_c_r = 0

        for j in list_of_corner_labels:
            if j == l_v:
                n_c_r = n_c_r +1

        list_label_temp.append(l_v)
        list_label_repeat_temp.append(n_c_r)


    df_background_check_m1["label"] = list_label_temp

    df_background_check_m1["label_repeat"] = list_label_repeat_temp

    df_background_check_m1.sort_values("label_repeat", ascending=False, inplace=True)
    df_background_check_m1 = df_background_check_m1.reset_index(drop=True)

    ### ADDED 17 JULY 2020
    print("df_background_check_m1: ", df_background_check_m1)
    

    if df_background_check_m1["label_repeat"][0]==4:

        background_label_m1 = df_background_check_m1["label"][0]
        m1_has_output = True
        # background check warning
        b_ch_m1_warining  = False

    if df_background_check_m1["label_repeat"][0]==3:
        background_label_m1 = df_background_check_m1["label"][0]
        m1_has_output = True
        # background check warning
        b_ch_m1_warining  = True

    if df_background_check_m1["label_repeat"][0]<3:
        background_label_m1 = 999
        m1_has_output = False
        # background check warning
        b_ch_m1_warining  = True

    ### ADDED 28 July 2020
    df_background_check_m1_only2corners = 999
    if df_background_check_m1["label_repeat"][0]==2:
        background_label_m1 = 999
        m1_has_output = True
        # background check warning
        b_ch_m1_warining  = True
        df_background_check_m1_only2corners = df_background_check_m1["label"][0]


    ### ADDED 28 July 2020  
    return background_label_m1, m1_has_output , df_background_check_m1_only2corners

    #return background_label_m1, m1_has_output




####### ----------------------------------------------------------



def func_background_label_2(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c):

    # Method 2 : Find the difference between pixel_coor for each label, and choose the label with maximum diff (most separated)
    ### Let's form a new dataframe with pixels and it's labels

    col_name_to_analyze = col_name_to_analyze
    df_pixel_info_example = df_pixel_info_example
    df_pixel_final_labels_hybrid_all_n_c = df_pixel_final_labels_hybrid_all_n_c

    
    df_background_m2 = pd.concat([df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c], axis=1)

    df_background_check_m2 = pd.DataFrame(columns= ["label", "x_diff", "y_diff"])

    #col_name_to_analyze = "pfl_n_c_3"

    x_diff_list = []
    y_diff_list = []

    ### ADDED 28 July 2020
    x_diff_times_y_diff_list = []
    ###

    list_label_temp = []

    for l_v in set(df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze]):

        list_label_temp.append(l_v)

        df_check_m2_temp = df_background_m2[["x_pixel_coor", "y_pixel_coor", col_name_to_analyze]][df_background_m2[col_name_to_analyze]==l_v]

        x_pix_max_temp = max(df_check_m2_temp["x_pixel_coor"])
        x_pix_min_temp = min(df_check_m2_temp["x_pixel_coor"])

        x_diff_list.append(x_pix_max_temp - x_pix_min_temp)

        y_pix_max_temp = max(df_check_m2_temp["y_pixel_coor"])
        y_pix_min_temp = min(df_check_m2_temp["y_pixel_coor"])

        y_diff_list.append(y_pix_max_temp - y_pix_min_temp)

        ### ADDED 28 July 2020
        x_diff_times_y_diff_temp = (x_pix_max_temp - x_pix_min_temp)*(y_pix_max_temp - y_pix_min_temp)
        x_diff_times_y_diff_list.append(x_diff_times_y_diff_temp)
        ####
        


    df_background_check_m2["label"] = list_label_temp

    df_background_check_m2["x_diff"] = x_diff_list
    df_background_check_m2["y_diff"] = y_diff_list

    ### ADDED 28 July 2020
    df_background_check_m2["x_times_y"] = x_diff_times_y_diff_list
    ###


    df_background_check_m2.sort_values("x_diff", ascending=False, inplace=True)
    df_background_check_m2 = df_background_check_m2.reset_index(drop=True)
    background_label_m2_x = df_background_check_m2["label"][0]

    df_background_check_m2.sort_values("y_diff", ascending=False, inplace=True)
    df_background_check_m2 = df_background_check_m2.reset_index(drop=True)
    background_label_m2_y = df_background_check_m2["label"][0]

    ### ADDED 28 July 2020
    df_background_check_m2.sort_values("x_times_y", ascending=False, inplace=True)
    df_background_check_m2 = df_background_check_m2.reset_index(drop=True)
    background_label_m2_x_times_y = df_background_check_m2["label"][0]
    ###


    ### ADDED 17 July 2020
    print("df_background_check_m2: ", df_background_check_m2)

    if background_label_m2_x == background_label_m2_y:

        background_label_m2 = background_label_m2_x
        m2_has_output = True
    ### ADDED 28 July 2020
    elif background_label_m2_x_times_y == background_label_m2_x:
        background_label_m2 = background_label_m2_x
        m2_has_output = True
    elif background_label_m2_x_times_y == background_label_m2_y:
        background_label_m2 = background_label_m2_y
        m2_has_output = True   
    ###
    else:
        background_label_m2 = 999
        m2_has_output = False

    return background_label_m2, m2_has_output




####### ----------------------------------------------------------


def func_background_label_3(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c):

    # Method 3 : Count the number of events and choose the label with minimum number of events

    ### Let's form a new dataframe with pixels and it's labels
    
    col_name_to_analyze = col_name_to_analyze
    df_pixel_info_example = df_pixel_info_example
    df_pixel_final_labels_hybrid_all_n_c = df_pixel_final_labels_hybrid_all_n_c
    

    df_background_m3 = pd.concat([df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c], axis=1)

    df_background_check_m3 = pd.DataFrame(columns= ["label", "n_events"])

    list_label_temp = []
    list_label_n_events_temp = []

    for l_v in set(df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze]):

        list_label_temp.append(l_v)

        df_check_m3_temp = df_background_m3[["x_pixel_coor", "y_pixel_coor", "events_index_pixel" , col_name_to_analyze]][df_background_m3[col_name_to_analyze]==l_v]

        n_l_e_c = 0
        index_list_temp = df_check_m3_temp.index.values.tolist()
        for i_x in index_list_temp:
            n_l_e_c = n_l_e_c + len(df_check_m3_temp["events_index_pixel"][i_x])

        list_label_n_events_temp.append(n_l_e_c)


    df_background_check_m3["label"] = list_label_temp
    df_background_check_m3["n_events"] = list_label_n_events_temp


    df_background_check_m3.sort_values("n_events", ascending=True, inplace=True)
    df_background_check_m3 = df_background_check_m3.reset_index(drop=True)

    background_label_m3 = df_background_check_m3["label"][0]
    m3_has_output = True
    
    
    return background_label_m3, m3_has_output
    



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

####### ----------------------------------------------------------

euro_unicode = "\u20A0"
#user_filename = user_filename.replace("/", euro_unicode)


### This is the plot name to be analyzed - name as given back to us by CCP.
### For some files this can be different than the name we sent them.
events_csv = pd.read_csv('/code/support_files/files_to_analyze.csv',engine='python')
print(sys.argv)
file_index = sys.argv[1]
file_index = int(file_index)
plot_name_to_be_analyzed = events_csv.loc[file_index-1]['x']

#plot_name_to_be_analyzed = "FCT_G596J_9F34_CD19_CD20-.csv"


### This folder contains separated users
folder_of_users_to_be_analyzed = "/data/users"

### This folder is for users already analyzed
folder_of_users_already_analyzed = "/mount/users_analyzed"


### This folder contains files of events for plots
folder_of_events_to_be_analyzed = "/data/data"

### This folder if for results
folder_for_results = "/data/results"


### Below is to set whether we want to consider subfolders
### and whether we want to copy analyzed users, and remove them.

subfolders_active = False
copy_users_already_analyzed = False
remove_users_already_analyzed = False





### Let's read the filematching file.
### This is given by Jerome group because they changed the names.
'''
filename_for_csv_of_filename_matching = "filename_matching.csv"
folder_of_filename_matching = "extra_files"
df_filename_matching = pd.read_csv(folder_of_filename_matching+'/'+filename_for_csv_of_filename_matching, header = None)
#df_filename_matching = pd.read_csv(folder_of_events_to_be_analyzed+'/'+filename_for_csv_of_filename_matching, header = None)

df_filename_matching.rename(columns= {0:"updated_name", 1:"old_name"}, inplace = True)
'''




#plot_name_to_be_analyzed = "COP8/COP-FB-025_Tcell_Tcell_COP-11-09-12_D07_CCRCD45RA.csv"

### Let's define the name to be used to read users
### Adjustment for "/" in names to euro_unicode
if "/" in plot_name_to_be_analyzed:
    plot_name_for_users_folder = plot_name_to_be_analyzed.replace("/", euro_unicode)
    
else:
    plot_name_for_users_folder = plot_name_to_be_analyzed

plot_name_for_users = plot_name_for_users_folder
print(plot_name_for_users)

### Let's define the name to be used for reading events
plot_name_for_events = plot_name_to_be_analyzed

print(plot_name_for_events)

'''
### Check whether the plot name is in the filename matching csv
if sum(df_filename_matching["updated_name"]==plot_name_for_events)>0:
    print("adjusting the name of files ...")
    ### Check for uniqueness of updated names
    if len(df_filename_matching[df_filename_matching["updated_name"]==plot_name_for_events].index)>1:
        print("More than one file name was found in the updated file names ...")
        sys.exit()
    
    index_of_plot_name = df_filename_matching[df_filename_matching["updated_name"]==plot_name_for_events].index[0]
    
    ### Introducing these two names just in case to be used later ...
    plot_name_for_events_updated = df_filename_matching["updated_name"][index_of_plot_name]
    plot_name_for_events_old = df_filename_matching["old_name"][index_of_plot_name]
    ### Here, we change the name to old one if there was a match
    plot_name_for_events = plot_name_for_events_old

    
print(plot_name_for_events)
''' 
 

####### ----------------------------------------------------------

### If we want to consider subfolders with plot name (without .csv)

'''
if subfolders_active == True:
    
    ### Adjust the path for reading users
    sub_folder_of_users_to_be_analyzed = plot_name_for_users[:-4]
    folder_of_users_to_be_analyzed = folder_of_users_to_be_analyzed+"/"+ sub_folder_of_users_to_be_analyzed
    
    ### Adjust the path for saving results
    sub_folder_for_results = plot_name_for_users[:-4]
    folder_for_results = folder_for_results+"/"+sub_folder_for_results
    
    # Create target directory & all intermediate directories if don't exists
    dirName= folder_for_results

    if not os.path.exists(dirName):
        os.makedirs(dirName)
        print("Directory " , dirName ,  " Created ")
    else:
        print("Directory " , dirName ,  " already exists")
    
    
    ### Adjust the path for moving/copying the users already analyzed
    sub_folder_of_users_already_analyzed = plot_name_for_users[:-4]
    folder_of_users_already_analyzed = folder_of_users_already_analyzed+"/"+sub_folder_of_users_already_analyzed

    # Create target directory & all intermediate directories if don't exists
    dirName2= folder_of_users_already_analyzed

    if not os.path.exists(dirName2):
        os.makedirs(dirName2)
        print("Directory " , dirName2 ,  " Created ")
    else:
        print("Directory " , dirName2 ,  " already exists")
'''

    


####### ----------------------------------------------------------

### Here we make a list of users filenames

list_of_users_files = []

### Find the current directory where the script is running.
### Other folders containing users and events info should be under this directory.
for file in os.listdir(folder_of_users_to_be_analyzed+"/"+plot_name_for_users_folder):
    list_of_users_files.append(file)

#print(list_of_users_files)
####### ----------------------------------------------------------



### Let's set a limit for number of users before analyzing
# limit_n_users = 200
# if (len(list_of_users_files)+1) > limit_n_users:
#     sys.exit()


### DataFrame to keep gates
gate_info_all = pd.DataFrame( columns=["user_id", "gates"])

### define an array of user ids
user_id = np.arange(1, len(list_of_users_files)+1)

### MAXIMUM NUMBER OF USERS TO READ
### Put a maximum for number of users that we read
### set the check to False if do not need it (i.e. all available users will be read)
check_max_n_u_read = True
max_number_of_users_to_read = 1000

if check_max_n_u_read == True:
    if len(user_id) > max_number_of_users_to_read:
        user_id = np.arange(1, max_number_of_users_to_read+1)  
    
### BELOW IS NUMBER OF USERS TO BE INCLUDED FOR FINAL CLUSTERING
### set the check to False if do not need it (i.e. all available users will be included)
check_max_n_u_include = True
max_number_of_users_to_include = 500


### THIS IS NEW - FOR WEIGHTED ANALYSIS
### WE WANT TO FOLLOW USERS - user id vs. "player"
df_user_player_info = pd.DataFrame( columns=["user_id", "player"] )


for u_i in user_id:
    plot_name_for_users_temp = plot_name_for_users_folder+"_u"+str(u_i)+".csv"
    print(plot_name_for_users_temp)
    
    path = folder_of_users_to_be_analyzed+"/"+plot_name_for_users_folder+"/"+plot_name_for_users_temp
    
    FCM_export_temp = pd.read_csv(path, names=["filepath","score","result","minmax", "player", "created_at"], header=0, sep="','", dtype={"score" :np.float64})
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
    
    
    
    
    ### results for a particular user
    results_FCM_export_temp = FCM_export_temp["result"][0]
    df_results_FCM_export_temp = pd.DataFrame(results_FCM_export_temp)
    ### min max for scaling vertices
    f_p_minmax = FCM_export_temp["minmax"][0]

    gates_u_i = func_get_gates_user(df_results_FCM_export_temp, f_p_minmax, do_scaling= True)
    gate_info_all = gate_info_all.append({"user_id":u_i, "gates": gates_u_i}, ignore_index=True)
    
    
    
    ### THIS IS NEW - FOR WEIGHTED ANALYSIS
    ### WE WANT TO FOLLOW USERS - user id vs. "player"
    player_number_of_this_user = FCM_export_temp["player"][0]
    df_user_player_info = df_user_player_info.append({"user_id":u_i, "player": player_number_of_this_user}, ignore_index=True)





####### ----------------------------------------------------------


### Let's copy users that are already analyzed to "users_analyzed"

'''
if copy_users_already_analyzed == True:

    for u_i in user_id:
        plot_name_for_users_temp = plot_name_for_users_folder+"_u"+str(u_i)+".csv"
        print(plot_name_for_users_temp)

        ### path of sourceto be copied - includes file name
        path_src_with_filename = folder_of_users_to_be_analyzed+"/"+plot_name_for_users_temp
        ### path of destination
        path_dst = folder_of_users_already_analyzed

        ### Let's copy
        shutil.copy2(path_src_with_filename, path_dst)
 


####### ----------------------------------------------------------


### Let's remove the users already analyzed
if remove_users_already_analyzed == True:
    for u_i in user_id:
        plot_name_for_users_temp = plot_name_for_users[:-4]+"_u"+str(u_i)+".csv"
        print(plot_name_for_users_temp)

        ### path of sourceto be copied - includes file name
        path_src_with_filename = folder_of_users_to_be_analyzed+"/"+plot_name_for_users_temp

        os.remove(path_src_with_filename)
    
    ### if there is a subfolder, this will remove it
    if subfolders_active == True:
        os.rmdir(folder_of_users_to_be_analyzed)
'''     
    



####### ----------------------------------------------------------

####### ----------------------------------------------------------
####### EVENTS - DOCKER - PIXEL RESOLUTION - SPECTRAL ANALYSIS OPTION   ------------------------------------------
####### ----------------------------------------------------------



fcs_csv_file_name = plot_name_for_events
fcs_csv_file_folder = folder_of_events_to_be_analyzed
events_csv = pd.read_csv(fcs_csv_file_folder+'/'+fcs_csv_file_name+'.csv')


######### ---------------------------------------- ###
######### --- DEALING WITH NA VALUES IN EVENTS --- ###
######### ------------ Added: 15 July 2020 --------###
######### ---------------------------------------- ###

### column name for events x and y
x_col_name = events_csv.columns.values[0]
y_col_name = events_csv.columns.values[1]

### index of the all events
all_events_csv_index = list(events_csv.index.values)

### Let's find the index of NA for each x and y
x_NA_index = list(events_csv[events_csv[x_col_name].isnull()].index)
y_NA_index = list(events_csv[events_csv[y_col_name].isnull()].index)

### Check to see if we have NA values, and if so, replace them ..
total_number_of_rows_with_NA_events = len(x_NA_index) + len(y_NA_index)

if total_number_of_rows_with_NA_events > 0:
    
    print("NA values found for events. They will be replaced ...")

    ### Let's find index with no NA
    non_NA_index =  list(set(all_events_csv_index) - set(x_NA_index) - set(y_NA_index))
    non_NA_index.sort()

    ### Below we find index if EITHER of x or y has NA values.
    index_to_change= list(set(y_NA_index + x_NA_index))
    index_to_change.sort()


    if len(non_NA_index) == 0:
        print("All events had NA values ...")
        sys.exit()

    ### We want to find the first events with no NA values.
    ### We will use this event to replace all NAs.
    ### So, we choose the first event in non NA list to repeat for NA ones
    index_event_to_repeat =  non_NA_index[0]   
    x_event_to_repeat = events_csv[x_col_name][index_event_to_repeat]
    y_event_to_repeat = events_csv[y_col_name][index_event_to_repeat]



    ### Now, Let's replace the NA values
    events_csv[x_col_name][index_to_change] = x_event_to_repeat
    events_csv[y_col_name][index_to_change] = y_event_to_repeat


   
    
    





### ***
### Here we define x and y of events

x_events = events_csv.iloc[:,0]
y_events = events_csv.iloc[:,1]


####### ----------------------------------------------------------

docker_path = '/mnt/app_docker_Hybrid'




### Here we define resolution for binning the FCS data.

res_x_FCS_example = 50
res_y_FCS_example = 50


### If this is set to True, the code will do ARI calculation to assign a weight to each user and do calculations
### TO BE ADDED LATER
weighted_analysis = False



### We should decide about assign_labels for Spectral Clustering to be "discretize" or "kmeans" - Default is kmeans

### THIS HAS BEEN CHANGED 15 July 2020
#choose_assign_labels="kmeans"
choose_assign_labels="discretize"
print("Spectral Clustering parameter: ", choose_assign_labels)

### If we want to save output separately for each n_clusters set below as True
# save_n_c_output_separately = False



####### ----------------------------------------------------------

### Here we define resolution for binning the FCS data.

# res_x_FCS_example = 100
# res_y_FCS_example = 100

### Heatmap for events and edges of bins

heatmap_events_example, xedges_example, yedges_example = np.histogram2d(x_events, y_events, bins=(res_x_FCS_example, res_y_FCS_example))


####### ----------------------------------------------------------

### Calculate some information about pixels

x_pix_c_example, y_pix_c_example, x_pix_edges_example, y_pix_edges_example, x_pixel_coor_example, y_pixel_coor_example, Grid_pixel_coor_example = \
func_pixel_coor(res_x_FCS_example, res_y_FCS_example, x_events, y_events)


####### ----------------------------------------------------------

### Form a dataframe for pixels information

df_pixel_info_example = func_df_pixels_events_info(res_x_FCS_example, res_y_FCS_example, x_events, y_events)


####### ----------------------------------------------------------


### Form a dataframe having information about labels of pixels for each user based on their gates

df_user_gates_label_pixel = pd.DataFrame(columns = ["user_id", "pixel_label"])

for u_i in user_id:

    index_u_i = gate_info_all[gate_info_all["user_id"]==u_i].index[0]
    user_temp = gate_info_all["gates"][index_u_i]
    user_gates_temp = extract_gates(user_temp, u_i)
    
    label_temp = func_label_events(user_gates_temp, x_pixel_coor_example, y_pixel_coor_example)
    
    df_user_gates_label_pixel = df_user_gates_label_pixel.append({"user_id":u_i, "pixel_label": label_temp}, ignore_index=True)


####### ----------------------------------------------------------


### Let's have a look at the distribution of number of gates from users:

n_gate_users_list = []

for i in np.arange(len(df_user_gates_label_pixel["user_id"])):
    n_c_t = len(np.unique(df_user_gates_label_pixel["pixel_label"][i])) - 1
    n_gate_users_list.append(n_c_t)
    
    
### Average number of users gate - Which one to choose?
### We choose median because in general it is more robust to outliers.

#n_gates_users_avg = np.floor(np.mean(n_gate_users_list) + 0.5)
#n_gates_users_avg = np.floor(np.median(n_gate_users_list) + 0.5)

### ADDED 11 AUG 2020
### Below, we "push" the number of gates to be higher by using 0.7 rather than 0.5
### We find mean, median and mode and choose the maximum between them as the number of gates from users
n_g_users_1 = np.floor(np.mean(n_gate_users_list) + 0.7)
n_g_users_2 = np.floor(np.median(n_gate_users_list) + 0.7)
n_g_users_3 = np.floor(statistics.mode(n_gate_users_list))
n_gates_users_avg = max([n_g_users_1, n_g_users_2, n_g_users_3])




### Let's make dataframe from user_is and n_gates, so we can use later to filter/include/exclude users

df_user_id_n_gates = pd.DataFrame({"user_id":df_user_gates_label_pixel["user_id"], "n_gates": n_gate_users_list})
 


####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### -----------  SANITY CHECK 1 --------------------------
####### ----------------------------------------------------------


### Let's read the clusterinfo.csv to do a sanity check
### we want to see whether avg number of clusters for users match fully gated examples etc.

folder_of_clusterinfo_sanity = "/data/pre_clusters"
filename_of_clusterinfo_sanity = plot_name_for_users_folder
clusterinfo_sanity = pd.read_csv(folder_of_clusterinfo_sanity+"/"+filename_of_clusterinfo_sanity+".csv")



####### ----------------------------------------------------------


analyze_this_file = False

### f_p_name is the same as "filepath" from MMOS results.
### It should be the same as
''''
f_p_name = plot_name_to_be_analyzed

if f_p_name in clusterinfo_sanity["filepath"]:
    i_f_p_name = clusterinfo_sanity[clusterinfo_sanity["filepath"]==f_p_name].index[0]
    if clusterinfo_sanity["flowMeans"][i_f_p_name] == n_gates_users_avg:
        analyze_this_file=True
    if clusterinfo_sanity["flowPeaks"][i_f_p_name] == n_gates_users_avg:
        analyze_this_file=True
    if clusterinfo_sanity["truth"][i_f_p_name] == n_gates_users_avg:
        analyze_this_file=True    

print(clusterinfo_sanity["truth"][0])
print(n_gates_users_avg)
if clusterinfo_sanity["truth"][0] == n_gates_users_avg:
        analyze_this_file=True
'''
### If analyze_this_file is True, the code will continue
### Otherwise, it will  exit
### Make below as comment to analyze it anyway

# if analyze_this_file == False:
#     sys.exit()


####### ----------------------------------------------------------




####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### --------  Find Number of Gates to Check Against ----------
####### ----------------------------------------------------------

### Let's read the clusterinfo.csv to do a sanity check
### we want to see whether avg number of clusters for users match fully gated examples etc.
#it's the same as sanity check 1
'''
folder_of_clusterinfo_sanity = "extra_files"
filename_of_clusterinfo_sanity ='clustersinfo.csv'
clusterinfo_sanity = pd.read_csv(folder_of_clusterinfo_sanity+"/"+filename_of_clusterinfo_sanity)
'''

#clusterinfo_sanity = pd.read_csv('clustersinfo.csv')


####### ----------------------------------------------------------




### Let's fins what is the true number of gates for this file based on other flowWHATEVER package we have

f_p_name = plot_name_to_be_analyzed
'''
if f_p_name in clusterinfo_sanity["filepath"]:
    i_f_p_name = clusterinfo_sanity[clusterinfo_sanity["filepath"]==f_p_name].index[0]
    ### CHOOSE WHICH ONE TO USE
    true_n_gates_to_ckeck_against = clusterinfo_sanity["clusters"][i_f_p_name]
#    true_n_gates_to_ckeck_against = clusterinfo_sanity["clusters_excluded"][i_f_p_name]
'''
if "truth" in list(clusterinfo_sanity.columns.values):
    true_n_gates_to_ckeck_against = clusterinfo_sanity["truth"][0]
else:
    if "flowPeaks" in list(clusterinfo_sanity.columns.values):
        true_n_gates_to_ckeck_against = clusterinfo_sanity["flowPeaks"][0]
    else:
        true_n_gates_to_ckeck_against = clusterinfo_sanity["clusters"][0]

'''
if not (f_p_name in clusterinfo_sanity["filepath"]):
    ### THIS IS JUST FOR NOW.
    ### WE SHOULD NOT HAVE THIS HERE.
    ### Delete = 2 and use = 12345
    true_n_gates_to_ckeck_against_temp = 2
#    true_n_gates_to_ckeck_against_temp = 12345
    true_n_gates_to_ckeck_against = true_n_gates_to_ckeck_against_temp
'''
 
 ### If we can not find the true number of gates, we exit
# if true_n_gates_to_ckeck_against == 12345:
#     sys.exit()
    


####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### --------  MAIN LOOP --------------------------------------
####### ----------------------------------------------------------


######################################################
############ COPY IMPORTANT ORIGINAL DATAFRAMES ######
######################################################

### Before going forward, let's make sure we get a copy of the original important dataframes etc.
### This is because for cases that we "excluded" users, we delete some of the users.
### If we want to do a calculation using all users again, we can user these "..._original" cases.

gate_info_all_original = gate_info_all.copy()

df_user_gates_label_pixel_original = df_user_gates_label_pixel.copy()

user_id_original = user_id

####### ----------------------------------------------------------


####### Now, we define a nested loop. One loop is to consider type of run (i.e. "all", "excluded", ...)
####### and another loop inside above loop to run for cases with various weights.



### Here, we define cases with and without excluding users.
### "all" considers all users, "excluded" excludes users with unmatched number of gates
### If a new type of run is added, an extension for column naming of output should be defined like "ex", "L10", etc.
### See    ### EXTENTION FOR COLUMN NAME below
### We also need to define which user to include or exclude
### See INCLUDE/EXCLUDE USERS below

list_type_of_run = ["all", "excluded", "limit10", "limit20", "limit30", "limit50", "limit100", "limit111"]

special_type_of_run_for_testing = "limit111"


### What weights to be considered
### NOTE ###
### PLEASE NOTE THAT THESE NAMES SHOULD BE THE SAME AS WE HAVE IN THE "average_player_performace_results.csv" file
### EITHER ALL OF THEM OR SOME OF THEM CAN BE CHOSEN
### NOTE ###
list_of_various_weights = [ 'No_weight', 'F1_weight', 'ARI_weight', 'MMOS_weight' ]
#list_of_various_weights = [ 'No_weight', 'F1_weight' ]

### This flag becomes True if number of users after exclusion is less than n_flag_for_excluded
n_flag_for_excluded = 30
flag_for_excluded= False

### We set this below in the loop to make sure we reset it everytime for various limit cases
#flag_for_limit = False


### If we want to do weighted analysis and we could not find a performance history for a player
### we use below value to assign a weight
score_for_player_with_no_performance_history = 0.5


### These are dataframes that contain everything at the end to be saved
df_all_to_be_saved_pixels = pd.DataFrame()
df_all_to_be_saved_events = pd.DataFrame()

###
### MAIN LOOP
###
for type_of_run in list_type_of_run:
    
    ### type of run for current loop
    type_of_run_temp = type_of_run
    print("\n", "TYPE OF RUN IS: ", type_of_run)
    #print("TYPE OF RUN IS: ", type_of_run)
    
    ### Let's set some important dataframes to their "original" as a copy was made above
    ### So, everytime, we start with all users, and then limit/exclude if wanted
    gate_info_all = gate_info_all_original.copy()
    df_user_gates_label_pixel = df_user_gates_label_pixel_original.copy()
    user_id = user_id_original

    
    
    ### EXTENTION FOR COLUMN NAME
    ### Here, for each type_of_run, we define an extention for column name to be considered in out final output
    if type_of_run_temp == "all":
        type_of_run_extention_name = "_all"
    if type_of_run_temp == "excluded":
        type_of_run_extention_name = "_ex"
        
    if "limit" in type_of_run_temp:
        n_limit_temp = int(type_of_run_temp[5:])
        type_of_run_extention_name = "_L"+str(n_limit_temp)
    ### FOR EXAMPLE:
    # if type_of_run_temp == "limit10":
    #     type_of_run_extention_name = "_L10"
    
    
    ######---------------------------------
    
    ######## HERE WE MAKE CHANGES RELATED TO WHICH USERS TO INCLUDE
    ######## For a new type of run, we need to define how to include/exclude users
    ########  --------     INCLUDE/EXCLUDE USERS ----------
    #################-----------------------------------------

    
    ### For "all" we do not need to do anything. We should use the original dataframes
    if type_of_run_temp == "all":
        ### To make sure we are working with original dataframes
        gate_info_all = gate_info_all_original.copy()
        df_user_gates_label_pixel = df_user_gates_label_pixel_original.copy()
        user_id = user_id_original

    ######---------------------------------
    
    ### This flag becomes True if number of users are not enough for "limit**"
    ### For example, for "limit50" in "list_type_of_run" above, if we do not have 50 users
    ### then "flag_for_limit" becomes True
    flag_for_limit = False
    
    ### This is for cases that we put a limit on the number of users
    if "limit" in type_of_run_temp:
        

        ### Find what is the limit number; i.e. 10, 20, 30 ...
        n_limit_temp = int(type_of_run_temp[5:])

        ### if n_limit_temp is larger than number of users we have, we include everything we have.
        ### so, basically it will be similar to the case of "all" in type_of_run
        if n_limit_temp >  len( user_id_original):
            ### Choose the first NN users for "limitNN"
            gate_info_all = gate_info_all_original.copy()
            df_user_gates_label_pixel = df_user_gates_label_pixel_original.copy()
            user_id = user_id_original

            flag_for_limit = True

        ### If we have more users than n_limit_temp
        ### Choose the first n_limit_temp users for "limitNN"
        if n_limit_temp <=  len( user_id_original):
            gate_info_all = gate_info_all_original.head(n_limit_temp).copy()
            df_user_gates_label_pixel = df_user_gates_label_pixel_original.head(n_limit_temp).copy()
            user_id = user_id_original[:n_limit_temp]

        
    ######---------------------------------
        
        
    
    ### If we need to exclude users, we do it here
    ### If we want to do any other type of exclusion (like consider only 10, 20, or 40 users etc.)
    ### we can ad "excluded_2" in a similar manner
    if type_of_run_temp == "excluded":

        ######################---------------------------------------------------

        ### Let's find the list of users to include

        df_users_to_include = df_user_id_n_gates[df_user_id_n_gates["n_gates"] == true_n_gates_to_ckeck_against]
        list_of_users_to_include = list(df_users_to_include["user_id"])
        number_of_users_to_include = len(list_of_users_to_include)
        
        ### "n_flag_for_excluded" is set above at the beginning of the loop
        if number_of_users_to_include < n_flag_for_excluded:
            flag_for_excluded= True
            

        ### Let's find the list to exclude
        list_of_users_to_exclude =  list(set( list(df_user_id_n_gates["user_id"]) ) - set(list_of_users_to_include) )

        ######################---------------------------------------------------

        ### Let's update the"df_user_gates_label_pixel" and "gate_info_all"
        ### and exclude users with unmatched number of gates

        for u_i_ex in list_of_users_to_exclude:

            index_gate_info_all_to_drop = gate_info_all[gate_info_all["user_id"]== u_i_ex].index[0]
            gate_info_all.drop(index_gate_info_all_to_drop, axis= 0, inplace= True)


            index_df_user_gates_label_pixel_to_drop = df_user_gates_label_pixel[df_user_gates_label_pixel["user_id"]== u_i_ex].index[0]
            df_user_gates_label_pixel.drop(index_df_user_gates_label_pixel_to_drop, axis= 0, inplace= True)

        ######################---------------------------------------------------

        ### THIS IS VERY VERY IMPORTANT
        ### THIS IS VERY VERY IMPORTANT
        ### THIS IS VERY VERY IMPORTANT
        ### If any user is excluded, the user_id must get updated

        user_id = list(gate_info_all["user_id"] )
#
    
    print("number of users: ", len(user_id) )


    ######## -------------------------------- #####
    ######## ------ empty user_id error ----- #####
    ######## ------ Added: 15 July 2020 ----- #####
    ######## -------------------------------- #####

    #print("check for number of users error - user_id: ", len(user_id) )

    ### We put a check on number of users. If no user was found (for example due to excluding some users etc.)
    ### we will skip that case. 

    if len(user_id) == 0:
        print("There is no user. We need to skip this case ...")
        continue

    print("There are at least one user. We will consider this case ...")
    





    ########################################
    ### MAXIMUM NUMBER OF USERS TO INCLUDE
    ########################################
    ### Plese note that this is different that maximum number of user to read (which is defined at the beginning of this script when defining user_id)
    ### SEE THE BEGINNING OF THE CODE TO DEFINE IT - SAME PLACE AS user_id
    #max_number_of_users_to_include = 10
    

    if check_max_n_u_include == True:

        ### Adjust the Dataframes
        gate_info_all_max_temp = gate_info_all.head(max_number_of_users_to_include).copy()
        df_user_gates_label_pixel_max_temp =     df_user_gates_label_pixel.head(max_number_of_users_to_include).copy()
        user_id_max_temp = user_id[:max_number_of_users_to_include]
   
        ### Let's copy back the dataframes
        gate_info_all = gate_info_all_max_temp.copy()
        df_user_gates_label_pixel = df_user_gates_label_pixel_max_temp.copy()
        user_id = user_id_max_temp
    
    
    
    
    print("number of users after maximum number allowed:", len(user_id) )








    ######################---------------------------------------------------
    
    ### Let's find the superimpose of all gates based on pixel repetition:
    #print(df_user_gates_label_pixel)
    print("we are here ... ")
    

    #pixel_repeats_gates_all = np.zeros(df_user_gates_label_pixel["pixel_label"][0].shape)
    pixel_repeats_gates_all = np.zeros(int(res_x_FCS_example*res_y_FCS_example) )

    for u_id in list(df_user_gates_label_pixel["user_id"]):
        u_id_index = df_user_gates_label_pixel[df_user_gates_label_pixel["user_id"]== u_id].index[0]

        temp_rep = np.array(df_user_gates_label_pixel["pixel_label"][u_id_index]>0).astype(int)

        pixel_repeats_gates_all = pixel_repeats_gates_all + temp_rep


    ######################---------------------------------------------------

    # To be able to follow pixels, we form a meshgrid with "index" values of pixels:
    ### Make meshgrid with "index" of x and y for pixels

    x_p_index, y_p_index =  np.meshgrid(np.arange(len(x_pix_c_example) ), np.arange(len(y_pix_c_example) ) )
    x_p_index, y_p_index = x_p_index.flatten(), y_p_index.flatten()

    Grid_p_index = np.vstack((x_p_index, y_p_index)).T


    ######################---------------------------------------------------

    ### Let's update pixel info with pixel index and pixel repetitions as well:

    df_pixel_info_example["x_pix_index"] = x_p_index
    df_pixel_info_example["y_pix_index"] = y_p_index

    df_pixel_info_example["pix_repetition"] = pixel_repeats_gates_all

    ######################---------------------------------------------------

    


    ######################---------------------------------------------------

    ### WEIGHTING RELATED CALCULATIONS

    ### Let's do the same calculations with weighting for users

    ### First, let's read the player weights

    folder_of_player_weights = "/code/support_files"

    filename_of_player_weights = "average_player_performace_results.csv"

    ######################---------------------------------------------------

    ### The path also includes the file name
    path_file_player_weights_temp = folder_of_player_weights+"/"+filename_of_player_weights

    df_average_player_performace = pd.read_csv(path_file_player_weights_temp)

        
  
    ### Let's see what weight we have in the df_average_player_performace
    ### We should make sure there a unique part of name in "df_average_player_performace" column names
    ### The original column names were ['player', 'avg_score_F1', 'avg_score_ARI', 'avg_MMOS_score']
    ### For new weights, we need to add another if statement
    ### We also need to change the code below defining "NEWNAME_weight"
    
    list_of_available_weight_name = list(df_average_player_performace.columns.values)
     
    


    ######################---------------------------------------------------

    ### Let's assign each u_id the corresponding weight

    #df_weights_for_users = pd.DataFrame(columns= ["user_id", "player", "F1_weight", "ARI_weight","MMOS_weight"])
    df_weights_for_users = pd.DataFrame(columns= list_of_available_weight_name)
    
    
    ### This is moved to above - to the beginning of the loop
    #score_for_player_with_no_performance_history = 0.5

    for u_i in user_id:
        
        ### Here we find the correponding "player" for u_i
        index_u_i_in_df_user_player_info = df_user_player_info[df_user_player_info["user_id"]== u_i].index[0]
        player_number_temp = df_user_player_info["player"][index_u_i_in_df_user_player_info]
        
        
        ### Let's make a temporary dictionary to be added to df_weights_for_users at the end of each loop
        ### we add the correponding weights to this dictionary for each user and add it to df_weights_for_users at the end of each loop
        dic_temp_df_weights_for_users = {"user_id":u_i , "player":player_number_temp}


        ### Let's see if we have any score for player_number_temp in df_average_player_performace
        ### Otherwise, we assign a score of "score_for_player_with_no_performance_history" - see above
        
        if player_number_temp in list(df_average_player_performace["player"]):
            index_p_n_t = df_average_player_performace[df_average_player_performace["player"] == player_number_temp].index[0]
            
            for weight_name_temp in list_of_available_weight_name:
                weight_value_temp = df_average_player_performace[weight_name_temp][index_p_n_t]
                dic_temp_df_weights_for_users[weight_name_temp] = weight_value_temp
            
                 
        else:
            ### see above for "score_for_player_with_no_performance_history"
            for weight_name_temp in list_of_available_weight_name:
                weight_value_temp = score_for_player_with_no_performance_history
                dic_temp_df_weights_for_users[weight_name_temp] = weight_value_temp

                
                
        ### This dataframe contains the information to be used for weights in hybrid matrix
        #df_weights_for_users = df_weights_for_users.append({"user_id":u_i , "player":player_number_temp , "F1_weight":player_our_weight_temp, "ARI_weight":player_our_weight_temp_2, "MMOS_weight":player_MMOS_weight_temp }, ignore_index=True)
        #print(dic_temp_df_weights_for_users)
        df_weights_for_users = df_weights_for_users.append(dic_temp_df_weights_for_users, ignore_index=True)

    df_weights_for_users["user_id"] = df_weights_for_users["user_id"].astype(int)
    df_weights_for_users["player"] = df_weights_for_users["player"].astype(int)






    ######################---------------------------------------------------

    ### Make a dataframe containing gates name and user_id -
    ### A new gate by LABEL 0 is added for each user for pixels not included in any of the gates


    df_gates_hypergraph = pd.DataFrame(columns = ["user_id", "gate_name"])

    u_i_repeat_list = []
    gate_name_list = []


    for u_i in user_id:

        index_u_i = gate_info_all[gate_info_all["user_id"]==u_i].index[0]
        user_temp = pd.DataFrame()
        user_temp = gate_info_all["gates"][index_u_i].copy()
        #user_gates_temp = extract_gates(user_temp, u_i)
        user_gates_temp = pd.DataFrame()
        user_gates_temp = extract_gates(user_temp, u_i).copy()



        label_temp = func_label_events(user_gates_temp, x_pixel_coor_example, y_pixel_coor_example)

        unique_label_temp = np.unique(label_temp)


        ### THIS WAS NOT WORKING FOR ALL CASES BECAUSE OF SOME REASON!!!
        #user_id_list_rep = np.repeat(u_i, len(user_gates_temp["gate_name"]) +1 ) # +1 is for label 0
        user_id_list_rep = np.repeat(u_i, len(unique_label_temp) )  
        u_i_repeat_list = u_i_repeat_list + (list(user_id_list_rep))

        #print("next user ... ")
        #print("")




        # gate 0 for pixels/events not included
        gate_0_temp = str(user_gates_temp["gate_name"][0])
        gate_0_temp = gate_0_temp[:-1] + "0"

        g_n_l_temp = []
        ### If there is a background pixel
        if 0 in unique_label_temp:
            for i_g_track in np.arange(int(len(unique_label_temp))):
                g_n_i = str(u_i)+"_G"+str(i_g_track)
                g_n_l_temp.append(g_n_i) 
            gate_name_list = gate_name_list + g_n_l_temp

        ### If there is NOT a background pixel - SHOULD BE VERY VERY RARE
        if not (0 in unique_label_temp):
            print("WARNING!!! NO BACKGROUND PIXEL ... !!! ALL PIXELS GATED!!!")
            print("u_i: ", u_i)
            for i_g_track in np.arange(int(len(unique_label_temp))):
                g_n_i = str(u_i)+"_G"+str(i_g_track)
                g_n_l_temp.append(g_n_i) 
            #gate_name_list = gate_name_list + [gate_0_temp] + g_n_l_temp
            gate_name_list = gate_name_list + g_n_l_temp


        # THIS WAS NOT WORKING FOR SOME CASES!!!   
        #gate_name_list = gate_name_list + [gate_0_temp] + (list(user_gates_temp["gate_name"]) )

        #print("gate_name_list: ", gate_name_list)
        #print("next user ... ")
        #print("")



    df_gates_hypergraph["user_id"] = u_i_repeat_list
    df_gates_hypergraph["gate_name"] = gate_name_list

    #print("df_gates_hypergraph: ", df_gates_hypergraph)
    #print("df_gates_hypergraph.shape: ", df_gates_hypergraph.shape)






    ######################---------------------------------------------------
    ######################---------------------------------------------------
    ############### MAIN CALCULATION LOOP STARTS HERE---------------------
    ######################---------------------------------------------------
    ######################---------------------------------------------------


    ### Let's define the weights we want to consider. It can be one or all options
    ### See above to define it.
    #list_of_various_weights = [ 'No_weight' ]
    #list_of_various_weights = [ 'No_weight', 'F1_weight', 'ARI_weight', 'MMOS_weight' ]


    for which_weight in list_of_various_weights:

        which_weight_temp = which_weight
        print(which_weight_temp)
        
        ### Define a col name to be used to save final outputs columns
        ### "type_of_run_extention_name" is defined above
        col_name_df_all_to_be_saved_temp = which_weight_temp+type_of_run_extention_name
         
            
        print(col_name_df_all_to_be_saved_temp)



        #################################--------------------------------------------


        ### Now, let's define our matrix for Hybrid Bipartite Graph


        size_hybrid_matrix = len(df_pixel_info_example["x_pixel_coor"]) + len(df_gates_hypergraph["gate_name"])
        print("len(df_pixel_info_example[x_pixel_coor]): ", len(df_pixel_info_example["x_pixel_coor"]))
        print("len(df_gates_hypergraph[gate_name]): ", len(df_gates_hypergraph["gate_name"]) )

        ### We had to change the dtype from usual case because now we may have float values it is now default - `numpy.float64`.
        hybrid_matrix = np.zeros((size_hybrid_matrix, size_hybrid_matrix ))

        i_hyb_m =  len(df_pixel_info_example["x_pixel_coor"]) - 1

        for u_i in user_id:

        #    print(u_i, end= "\r")
            index_u_i = gate_info_all[gate_info_all["user_id"]==u_i].index[0]
            user_temp = gate_info_all["gates"][index_u_i]
            user_gates_temp = extract_gates(user_temp, u_i)


            label_temp = func_label_events(user_gates_temp, x_pixel_coor_example, y_pixel_coor_example)

            unique_label_temp = np.unique(label_temp)

            n_check = 0

            ### THIS IS FOR WEIGHTING - Some weights needs extra adjustments (like ARI_weight or No_weight)
            ### NEW WEIGHT - If a new weight is defined WE NEED TO MODIFY BELOW ACCORDINGLY - see below
            index_weight_for_u_i = df_weights_for_users[df_weights_for_users["user_id"]== u_i].index[0]
            if which_weight_temp == "F1_weight":
                weight_for_u_i = df_weights_for_users[which_weight_temp][index_weight_for_u_i]

            ### ARI ranges between -1 to 1. If ARI is used, we add +1 to weight to make the weight is always positive.
            if which_weight_temp == "ARI_weight":
                weight_for_u_i = (df_weights_for_users[which_weight_temp][index_weight_for_u_i] +1)/2 # +1 is to make ARI positive. and /2 to make it between 0-1

            if which_weight_temp == "MMOS_weight":
                weight_for_u_i = df_weights_for_users[which_weight_temp][index_weight_for_u_i]

            if which_weight_temp == "No_weight":
                weight_for_u_i = 1
            
            ### NEW WEIGHT - If a new weight is added (this is general, and it only assigns the weight as it is.)
            if not (which_weight_temp in ["No_weight", "F1_weight", "ARI_weight", "MMOS_weight"] ):
                print("this is a new weight ...")
                weight_for_u_i = df_weights_for_users[which_weight_temp][index_weight_for_u_i]


            ### Now, let's check if the weight is "EXACTLY" zero.
	    ### If weight is zero, we need to change it to a very small value rather than zero,
	    ### because Spectral Clustering does not work with not fully connected graphs.
            if weight_for_u_i == 0:
                print("Weight was zero. It will be changed to 0.000001 for Spectral Clustering.")		
                weight_for_u_i = 0.000001


            #print("u_i: ", u_i)	
            #print("label_temp: ", label_temp)		
            #print("unique_label_temp: ", unique_label_temp)		


            for u_l_t in unique_label_temp:

                i_hyb_m = i_hyb_m + 1 # we go to next hyper node (next cluster/gate)

                index_u_l_t = np.where(label_temp == u_l_t)[0]

                #if type_of_run_temp == "excluded":

                #    print("gotch yeah ...")
                #    print("u_l_t:  ", u_l_t  ,"-   len(index_u_l_t): ", len(index_u_l_t) )
                #    print("u_i: ", u_i)	
                #    print("label_temp: ", label_temp)		
                #    print("unique_label_temp: ", unique_label_temp)
                #    print("weight_for_u_i: ", weight_for_u_i)		
                
                ### if any particular function of weight_for_u_i is required, it can be defined here
		### let's define a function
                a_function_of_weight = weight_for_u_i
                ### CHANGED 16 JULY 2020
                #a_function_of_weight = (10*weight_for_u_i)**2

                #print(a_function_of_weight)
                #if a_function_of_weight < 10**(-2):
                #    print(a_function_of_weight)
                #    a_function_of_weight = 10**(-2)


                hybrid_matrix[i_hyb_m, index_u_l_t] = 1.0 * a_function_of_weight # (1+weight_for_u_i)
                hybrid_matrix[index_u_l_t, i_hyb_m] = 1.0 * a_function_of_weight # (1+weight_for_u_i)

                ### n_check after this loop should be the same as total number of pixels
                n_check = n_check + len(index_u_l_t)



        #################################--------------------------------------------


        ### We should decide about assign_labels to be "discretize" or "kmeans" - Default is kmeans


        # choose_assign_labels="kmeans"
        # #choose_assign_labels="discretize"

        ### Now we run the clustering algorithm for various n_clusters
        ### NOTE that +1 is needed for LABEL 0 (i.e. pixels that are not gates)
        ### So if average number of gates from user is 3 (n_gates_users_avg = 3), n_clusters should be 4

        ### VERY VERY VERY IMPORTANT
        # n_c_avg_users = n_gates_users_avg +1
        n_c_avg_users = true_n_gates_to_ckeck_against +1

        ### THIS IS FOR SPECIAL CASE - Here we consider the average of users rather than our input
        if type_of_run_temp ==  special_type_of_run_for_testing:
            print("Dealing with specail type of run: ", type_of_run_temp)
            n_c_avg_users = n_gates_users_avg +1

        n_cluster_list = [int(n_c_avg_users)]
        #n_cluster_list = [int(n_c_avg_users -1), int(n_c_avg_users), int(n_c_avg_users +1)]

        df_pixel_final_labels_hybrid_all_n_c = pd.DataFrame()
        final_e_df_example_all_n_c = pd.DataFrame()



        for n_clusters in n_cluster_list:
            
            print("number of clusters: n_gate +1 : ", n_clusters)

            clustering_SC = SpectralClustering(n_clusters=n_clusters,  affinity= "precomputed", assign_labels=choose_assign_labels, random_state=0).fit(hybrid_matrix)



            pixel_final_labels_hybrid = clustering_SC.labels_[:len(df_pixel_info_example["x_pixel_coor"]) ]

            df_pixel_final_labels_hybrid = pd.DataFrame({"pixel_final_label":pixel_final_labels_hybrid})




            ### After clustering the pixels, we can go back to our events and label the events inside each pixel.
            ### First, we update our dataframe that contains pixel information with the labels we found for pixels.

            ### Update the dataframe with pixel labels

            df_pixel_info_example["pixel_final_label"] = pixel_final_labels_hybrid


            events_final_index_example = []
            events_final_label_example = []

            for i in np.arange(df_pixel_info_example.shape[0]):

                events_final_index_example = events_final_index_example + (list(df_pixel_info_example["events_index_pixel"][i]) )
                events_final_label_example = events_final_label_example + list(np.zeros(len(df_pixel_info_example["events_index_pixel"][i])) + df_pixel_info_example["pixel_final_label"][i])


            final_e_df_example = pd.DataFrame({"events_index":events_final_index_example, "events_label":events_final_label_example})

            final_e_df_example = final_e_df_example.sort_values("events_index").reset_index().drop("index", axis= 1).drop("events_index", axis= 1)




            ### Below dataframes will keep the results of all n_clusters. One can also save these files
            df_pixel_final_labels_hybrid_all_n_c["pfl_n_c_"+str(n_clusters)] = pixel_final_labels_hybrid
            final_e_df_example_all_n_c["efl_n_c_"+str(n_clusters)] = final_e_df_example["events_label"].astype(np.int64)


        # #### Save output for all_n_c
        # #df_pixel_final_labels_hybrid_all_n_c.to_csv(str(docker_path)+'/'+ fcs_csv_file_name+'_pixel_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)
        # #final_e_df_example_all_n_c.to_csv(str(docker_path)+'/'+ fcs_csv_file_name+'_events_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)

        # df_pixel_final_labels_hybrid_all_n_c.to_csv(fcs_csv_file_name+'_pixel_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)
        # final_e_df_example_all_n_c.to_csv(fcs_csv_file_name+'_events_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)


        #################################--------------------------------------------

        ### ADDED 28 July 2020
        ### This is to skip the loop if finding background label was incolclusive.
        ### When it is inconclusive, it's value will change to True. See background analysis part of the code
        check_background_inconclusive = False
        ###
        

        ### HERE WE DEAL WITH BACKGROUND.
        ### We want to set the label 0 to be background
        for i in np.arange(len(df_pixel_final_labels_hybrid_all_n_c.columns)):
            col_name_p = df_pixel_final_labels_hybrid_all_n_c.columns[i]

            col_name_e = final_e_df_example_all_n_c.columns[i]


        ### Let's find the background label

            col_name_to_analyze = col_name_p
            col_name_to_analyze_events = col_name_e

            background_label_m1, m1_has_output, background_label_m1_only2corners = func_background_label_1(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c)

            background_label_m2, m2_has_output = func_background_label_2(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c)

            background_label_m3, m3_has_output = func_background_label_3(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c)

            
            ### ADDED: 17 July 2020
            print("background_label_m1: ", background_label_m1)
            print("background_label_m1_only2corners: ", background_label_m1_only2corners)
            print("background_label_m2: ", background_label_m2)
            print("background_label_m3: ", background_label_m3)




            # Check if at least two of three methods have outputs:
            if sum(np.array([(m1_has_output == True) , (m2_has_output == True) , (m3_has_output == True)]).astype(int))<=1:
                print("Bakcground label could not be found.")
                sys.exit()

            if (background_label_m1 == background_label_m2) & (background_label_m1 == background_label_m3):
                background_label_decided = background_label_m1
                print("C1 case for background")

            elif (background_label_m1 == background_label_m2) & (background_label_m1 != background_label_m3):
                background_label_decided = background_label_m1
                print("C2 case for background")

            elif (background_label_m1 != background_label_m2) & (background_label_m1 == background_label_m3):
                background_label_decided = background_label_m1
                print("C3 case for background")

            elif (background_label_m1 != background_label_m2) & (background_label_m2 == background_label_m3):
                background_label_decided = background_label_m2
                print("C4 case for background")

            ### ADDED 28 July 2020
            elif (background_label_m1_only2corners == background_label_m2):
                background_label_decided = background_label_m1_only2corners
                print("Warning. Background decided based on: (background_label_m1_only2corners == background_label_m2)")
                print("C5 case for background")
            ### Below elif checks "..._only2corners" against background_label_m3
            ### because m3 is based on the lowest number of events, it can be a bit risky that can result in the wrong label.
            ### So, we can comment it out if necessary ...
            elif (background_label_m1_only2corners == background_label_m3):
                background_label_decided = background_label_m1_only2corners
                print("Warning. Background decided based on: (background_label_m1_only2corners == background_label_m3)")
                print("C6 case for background")

            else:
                print("WARNING! Could not match methods for Background label ...")
                print("We need to skip this case ...")
                print("")
                check_background_inconclusive = True
                #sys.exit()
                #continue

            ### ADDED 28 July 2020
            ### This is to exit the loop if background label is inconclusive ...
            if check_background_inconclusive == True:
                print("inconclusive background ... skipping this case ...")
                print("")
                continue



            print("background label is: ", background_label_decided)


#             run_time_after_Background = time.time() - start_t #time.process_time() - start_t
#             print("run_time_after_Background : ", run_time_after_Background)

        ### Let's adjust the labels - We always keep the label 0 for background
        ### So, if background label is not 0, we need to switch the labels
            if background_label_decided != 0:
                df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze][df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze]==0]=1234
                df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze][df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze]==background_label_decided]=0
                df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze][df_pixel_final_labels_hybrid_all_n_c[col_name_to_analyze]==1234]= background_label_decided

                ### Check whether there is any events in background pixels, and if yes, change the label acoordingly:
                if len(final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==background_label_decided])>0:
                    print("number of background events:", len(final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==background_label_decided]) )
                    final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==0]=1234
                    final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==background_label_decided]=0
                    final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==1234]= background_label_decided

                if len(final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==background_label_decided])==0:
                    print("number of background events:", len(final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==background_label_decided]) )
                    final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==0]=1234
                    final_e_df_example_all_n_c[col_name_to_analyze_events][final_e_df_example_all_n_c[col_name_to_analyze_events]==1234]= background_label_decided

            print(" ")
         


        #################################--------------------------------------------

        ###
        ### Let's update the columns names
        for c_n in list(df_pixel_final_labels_hybrid_all_n_c.columns.values):
            
            if (type_of_run_temp == "excluded") & (flag_for_excluded==True):
                df_pixel_final_labels_hybrid_all_n_c.rename(columns= {c_n: c_n+"_"+col_name_df_all_to_be_saved_temp+"_flag"}, inplace = True)

            elif ("limit" in type_of_run_temp) & (flag_for_limit==True):
                df_pixel_final_labels_hybrid_all_n_c.rename(columns= {c_n: c_n+"_"+col_name_df_all_to_be_saved_temp+"_flag"}, inplace = True)
                
            else:
                df_pixel_final_labels_hybrid_all_n_c.rename(columns= {c_n: c_n+"_"+col_name_df_all_to_be_saved_temp}, inplace = True)
            
                                                                               

        ### Let's update the columns names
        for c_n in list(final_e_df_example_all_n_c.columns.values):
            if (type_of_run_temp == "excluded") & (flag_for_excluded==True):
                final_e_df_example_all_n_c.rename(columns= {c_n: c_n+"_"+col_name_df_all_to_be_saved_temp+"_flag"}, inplace = True)

            elif ("limit" in type_of_run_temp) & (flag_for_limit==True):
                final_e_df_example_all_n_c.rename(columns= {c_n: c_n+"_"+col_name_df_all_to_be_saved_temp+"_flag"}, inplace = True)
                
            else:
                final_e_df_example_all_n_c.rename(columns= {c_n: c_n+"_"+col_name_df_all_to_be_saved_temp}, inplace = True)

        #################################--------------------------------------------

        ### ADDED 28 July 2020
        ### This is to exit the loop if background label is inconclusive ...
        if check_background_inconclusive == True:
            print("... exiting the loop due to inconclusive background ... ")
            #print("second")
            print("")
            continue
        
        

        ### Let's add the new analysis to our final dataframe to be saved:

        df_all_to_be_saved_pixels = pd.concat([df_all_to_be_saved_pixels, df_pixel_final_labels_hybrid_all_n_c], axis = 1)
        df_all_to_be_saved_events = pd.concat([df_all_to_be_saved_events, final_e_df_example_all_n_c] , axis = 1)




        #################################--------------------------------------------

print("########## here #############")
print("########## here #############")
    

## saving the outputs
#### Let's save the output
df_all_to_be_saved_pixels.to_csv(folder_for_results+'/'+plot_name_for_users_folder+'_results_pixel'+'.csv', index=False)
df_all_to_be_saved_events.to_csv(folder_for_results+'/'+plot_name_for_users_folder+'_results'+'.csv', index=False)




######################---------------------------------------------------






####### ----------------------------------------------------------

print("Code ended ...")


####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------













