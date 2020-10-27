

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


### To read MMOS outputs
import json

### F1 score
from sklearn.metrics import f1_score


import matplotlib.pyplot as plt

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

    return background_label_m1, m1_has_output




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


    df_background_check_m2["label"] = list_label_temp

    df_background_check_m2["x_diff"] = x_diff_list
    df_background_check_m2["y_diff"] = y_diff_list


    df_background_check_m2.sort_values("x_diff", ascending=False, inplace=True)
    df_background_check_m2 = df_background_check_m2.reset_index(drop=True)
    background_label_m2_x = df_background_check_m2["label"][0]

    df_background_check_m2.sort_values("y_diff", ascending=False, inplace=True)
    df_background_check_m2 = df_background_check_m2.reset_index(drop=True)
    background_label_m2_y = df_background_check_m2["label"][0]

    if background_label_m2_x == background_label_m2_y:

        background_label_m2 = background_label_m2_x
        m2_has_output = True
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

### This is the plot name to be analyzed - name as given back to us by CCP.
### For some files this can be different than the name we sent them.
events_csv = pd.read_csv('/code/support_files/files_to_analyze.csv',engine='python')
print(sys.argv)
file_index = sys.argv[1]
file_index = int(file_index)
plot_name_to_be_analyzed = events_csv.loc[file_index-1]['x']

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
"""
filename_for_csv_of_filename_matching = "filename_matching.csv"
df_filename_matching = pd.read_csv(folder_of_events_to_be_analyzed+'/'+filename_for_csv_of_filename_matching, header = None)

df_filename_matching.rename(columns= {0:"updated_name", 1:"old_name"}, inplace = True)
"""



#plot_name_to_be_analyzed = "COP8/COP-FB-025_Tcell_Tcell_COP-11-09-12_D07_CCRCD45RA.csv"

### Let's define the name to be used to read users
### Adjustment for "/" in names to euro_unicode
if "/" in plot_name_to_be_analyzed:
    plot_name_for_users_folder = plot_name_to_be_analyzed.replace("/", euro_unicode)
    
else:
    plot_name_for_users_folder = plot_name_to_be_analyzed

plot_name_for_users = plot_name_to_be_analyzed
print(plot_name_for_users)
### Let's define the name to be used for reading events
plot_name_for_events = plot_name_to_be_analyzed

#print(plot_name_for_events)

"""
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
"""    
 
####### ----------------------------------------------------------

"""
### If we want to consider subfolders with plot name (without .csv)


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

"""
    


####### ----------------------------------------------------------

### Here we make a list of users filenames

list_of_users_files = []

### Find the current directory where the script is running.
### Other folders containing users and events info should be under this directory.
project_folder = os.getcwd()


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

for u_i in user_id:
    plot_name_for_users_temp = plot_name_for_users_folder+"_u"+str(u_i)+".csv"
    path = folder_of_users_to_be_analyzed+"/"+plot_name_for_users_folder+"/"+plot_name_for_users_temp
    FCM_export_temp = pd.read_csv(path, names=["filepath","score","result","minmax", "player", "created_at"], header=0, sep="','", dtype={"score" :np.float64},engine='python')
    FCM_export_temp["filepath"] = FCM_export_temp["filepath"].str[1:]
    ### Now this should be applied to "player” columns
    # FCM_export_temp["minmax”] = FCM_export_temp["minmax”].str[:-1]
    FCM_export_temp["result"] = FCM_export_temp["result"].apply(json.loads)
    FCM_export_temp["minmax"] = FCM_export_temp["minmax"].apply(json.loads)
    # ### NEW for player - UserID
    # FCM_export_temp["player”] = FCM_export_temp["player”].str[:-1]
    FCM_export_temp["player"] = FCM_export_temp["player"].astype(int)
    ### NEW 2 - After adding time : "created_at”
    FCM_export_temp["created_at"] = FCM_export_temp["created_at"].str[:-1]
    ### Convert string to timestamp for "created_at”
    # FCM_export_temp["created_at”] =pd.to_datetime( FCM_export_temp["created_at”])
    ### results for a particular user
    results_FCM_export_temp = FCM_export_temp["result"][0]
    df_results_FCM_export_temp = pd.DataFrame(results_FCM_export_temp)
    ### min max for scaling vertices
    f_p_minmax = FCM_export_temp["minmax"][0]
    gates_u_i = func_get_gates_user(df_results_FCM_export_temp, f_p_minmax, do_scaling= True)
    gate_info_all = gate_info_all.append({"user_id":u_i, "gates": gates_u_i}, ignore_index=True)



####### ----------------------------------------------------------

### Let's copy users that are already analyzed to "users_analyzed"


if copy_users_already_analyzed == True:

    for u_i in user_id:
        plot_name_for_users_temp = list_of_users_files[u_i-1]
        print(plot_name_for_users_temp)

        ### path of sourceto be copied - includes file name
        path_src_with_filename = folder_of_users_to_be_analyzed+"/"+plot_name_for_users_temp
        print(path_src_with_filename)        
	### path of destination
        path_dst = folder_of_users_already_analyzed+"/"+plot_name_for_users_temp
        print(path_dst)
        ### Let's move
        shutil.move(path_src_with_filename, path_dst)
 

####### ----------------------------------------------------------

####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### START MAIN ANALYSIS   ------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
fcs_csv_file_name = plot_name_for_events

fcs_csv_file_folder = folder_of_events_to_be_analyzed
events_csv = pd.read_csv(fcs_csv_file_folder+'/'+fcs_csv_file_name+'.csv')

docker_path = '/mount/temp_result'


### Here we define resolution for binning the FCS data.

res_x_FCS_example = 100
res_y_FCS_example = 100


### If this is set to True, the code will do ARI calculation to assign a weight to each user and do calculations
### TO BE ADDED LATER
weighted_analysis = False



### We should decide about assign_labels for Spectral Clustering to be "discretize" or "kmeans" - Default is kmeans

choose_assign_labels="kmeans"
#choose_assign_labels="discretize"

### If we want to save output separately for each n_clusters set below as True
save_n_c_output_separately = False



####### ----------------------------------------------------------

### ***

#events_csv = pd.read_csv(str(fcs_csv_file_folder) + '/' + str(fcs_csv_file_name) + '.csv')

### ***

x_events = events_csv.iloc[:,0]
y_events = events_csv.iloc[:,1]



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
n_gates_users_avg = np.floor(np.median(n_gate_users_list) + 0.5)
    


####### ----------------------------------------------------------

### Let's fins the superimpose of all gates based on pixel repetition:


pixel_repeats_gates_all = np.zeros(df_user_gates_label_pixel["pixel_label"][0].shape)

for i in np.arange(len(df_user_gates_label_pixel["pixel_label"])):
    
    temp_rep = np.array(df_user_gates_label_pixel["pixel_label"][i]>0).astype(int)
    
    pixel_repeats_gates_all = pixel_repeats_gates_all + temp_rep

####### ----------------------------------------------------------

# To be able to follow pixels, we form a meshgrid with "index" values of pixels:
### Make meshgrid with "index" of x and y for pixels

x_p_index, y_p_index =  np.meshgrid(np.arange(len(x_pix_c_example) ), np.arange(len(y_pix_c_example) ) )
x_p_index, y_p_index = x_p_index.flatten(), y_p_index.flatten()

Grid_p_index = np.vstack((x_p_index, y_p_index)).T



####### ----------------------------------------------------------

### Let's update pixel info with pixel index and pixel repetitions as well:

df_pixel_info_example["x_pix_index"] = x_p_index
df_pixel_info_example["y_pix_index"] = y_p_index

df_pixel_info_example["pix_repetition"] = pixel_repeats_gates_all


####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### Define Hybrid Matrix
####### ----------------------------------------------------------


### Make a dataframe containing gates name and user_id -
### A new gate by LABEL 0 is added for each user for pixels not included in any of the gates
df_gates_hypergraph = pd.DataFrame(columns = ["user_id", "gate_name"])

u_i_repeat_list = []
gate_name_list = []

for u_i in user_id:
    
    index_u_i = gate_info_all[gate_info_all["user_id"]==u_i].index[0]
    user_temp = gate_info_all["gates"][index_u_i]
    user_gates_temp = extract_gates(user_temp, u_i)
    
    
    user_id_list_rep = np.repeat(u_i, len(user_gates_temp["gate_name"]) +1 ) # +1 is for label 0
    u_i_repeat_list = u_i_repeat_list + (list(user_id_list_rep))
    
    # gate 0 for pixels/events not included
    gate_0_temp = str(user_gates_temp["gate_name"][0])
    gate_0_temp = gate_0_temp[:-1] + "0"
    
    gate_name_list = gate_name_list + [gate_0_temp] + (list(user_gates_temp["gate_name"]) )
    
    
df_gates_hypergraph["user_id"] = u_i_repeat_list
df_gates_hypergraph["gate_name"] = gate_name_list


####### ----------------------------------------------------------



### Now, let's define our matrix for Hybrid Bipartite Graph


size_hybrid_matrix = len(df_pixel_info_example["x_pixel_coor"]) + len(df_gates_hypergraph["gate_name"])

hybrid_matrix = np.zeros((size_hybrid_matrix, size_hybrid_matrix ), dtype=np.uint8)

i_hyb_m =  len(df_pixel_info_example["x_pixel_coor"]) - 1

for u_i in user_id:
    
#    print(u_i, end= "\r")
    index_u_i = gate_info_all[gate_info_all["user_id"]==u_i].index[0]
    user_temp = gate_info_all["gates"][index_u_i]
    user_gates_temp = extract_gates(user_temp, u_i)
   

    label_temp = func_label_events(user_gates_temp, x_pixel_coor_example, y_pixel_coor_example)

    unique_label_temp = np.unique(label_temp)
    
    n_check = 0
    for u_l_t in unique_label_temp:
        
        i_hyb_m = i_hyb_m + 1 # we go to next hyper node (next clauter/gate)
        
        index_u_l_t = np.where(label_temp == u_l_t)[0]
        hybrid_matrix[i_hyb_m, index_u_l_t] = 1
        hybrid_matrix[index_u_l_t, i_hyb_m] = 1
        
        n_check = n_check + len(index_u_l_t)
        



####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### Apply Spectral Clustering
####### ----------------------------------------------------------

### We should decide about assign_labels to be "discretize" or "kmeans" - Default is kmeans


# choose_assign_labels="kmeans"
# #choose_assign_labels="discretize"

### Now we run the clustering algorithm for various n_clusters
### NOTE that +1 is needed for LABEL 0 (i.e. pixels that are not gates)
### So if average number of gates from user is 3 (n_gates_users_avg = 3), n_clusters should be 4

n_c_avg_users = n_gates_users_avg +1


n_cluster_list = [int(n_c_avg_users)]
#n_cluster_list = [int(n_c_avg_users -1), int(n_c_avg_users), int(n_c_avg_users +1)]

df_pixel_final_labels_hybrid_all_n_c = pd.DataFrame()
final_e_df_example_all_n_c = pd.DataFrame()



for n_clusters in n_cluster_list:
    
    print(n_clusters)

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


    ### SAVE THE OUTPUTS
    
    if save_n_c_output_separately == True:
    
        ### Save output of pixel final labels
    #    df_pixel_final_labels_hybrid.to_csv(str(docker_path)+'/'+ fcs_csv_file_name+'_final_pixel_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_n_c_'+str(n_clusters)+'.csv', index=False)
        df_pixel_final_labels_hybrid.to_csv(fcs_csv_file_name+'_pixel_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_n_c_'+str(n_clusters)+'.csv', index=False)

        ### Save output of events final labels
    #    final_e_df_example.to_csv(str(docker_path)+'/'+ fcs_csv_file_name+'_events_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_n_c_'+str(n_clusters)+'.csv', index=False)
        final_e_df_example.to_csv(fcs_csv_file_name+'_events_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_n_c_'+str(n_clusters)+'.csv', index=False)

    
    ### Below dataframes will keep the results of all n_clusters. One can also save these files
    df_pixel_final_labels_hybrid_all_n_c["pfl_n_c_"+str(n_clusters)] = pixel_final_labels_hybrid
    final_e_df_example_all_n_c["efl_n_c_"+str(n_clusters)] = final_e_df_example["events_label"].astype(np.int64)


# #### Save output for all_n_c
# #df_pixel_final_labels_hybrid_all_n_c.to_csv(str(docker_path)+'/'+ fcs_csv_file_name+'_pixel_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)
# #final_e_df_example_all_n_c.to_csv(str(docker_path)+'/'+ fcs_csv_file_name+'_events_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)

# df_pixel_final_labels_hybrid_all_n_c.to_csv(fcs_csv_file_name+'_pixel_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)
# final_e_df_example_all_n_c.to_csv(fcs_csv_file_name+'_events_final_label_'+str(choose_assign_labels)+'_'+str(res_x_FCS_example)+'_'+str(res_x_FCS_example)+'_all_n_c.csv', index=False)

    
    



####### ----------------------------------------------------------

for i in np.arange(len(df_pixel_final_labels_hybrid_all_n_c.columns)):
    
    col_name_p = df_pixel_final_labels_hybrid_all_n_c.columns[i]
    
    col_name_e = final_e_df_example_all_n_c.columns[i]
    
    
### Let's find

    col_name_to_analyze = col_name_p
    col_name_to_analyze_events = col_name_e

    background_label_m1, m1_has_output = func_background_label_1(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c)

    background_label_m2, m2_has_output = func_background_label_2(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c)

    background_label_m3, m3_has_output = func_background_label_3(col_name_to_analyze, df_pixel_info_example, df_pixel_final_labels_hybrid_all_n_c)

    # Check if at least two of three methods have outputs:
    if sum(np.array([(m1_has_output == True) , (m2_has_output == True) , (m3_has_output == True)]).astype(int))<=1:
        print("Bakcground label could not be found.")
        sys.exit()

    if (background_label_m1 == background_label_m2) & (background_label_m1 == background_label_m3):
        background_label_decided = background_label_m1

    elif (background_label_m1 == background_label_m2) & (background_label_m1 != background_label_m3):
        background_label_decided = background_label_m1

    elif (background_label_m1 != background_label_m2) & (background_label_m1 == background_label_m3):
        background_label_decided = background_label_m1

    elif (background_label_m1 != background_label_m2) & (background_label_m2 == background_label_m3):
        background_label_decided = background_label_m2
       
    print("background label is: ", background_label_decided)
    
    

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
        
####### ----------------------------------------------------------

### saving the outputs

#df_pixel_final_labels_hybrid_all_n_c.to_csv(folder_for_results+'/'+plot_name_for_users_folder+'_results'+'.csv', index=False)
final_e_df_example_all_n_c.to_csv(folder_for_results+'/'+plot_name_for_users_folder+'_results'+'.csv', index=False)


####### ----------------------------------------------------------


####### ----------------------------------------------------------
print("Code ended ...")


####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------













