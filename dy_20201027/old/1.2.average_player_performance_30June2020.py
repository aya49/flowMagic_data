

### python_script_example.py

####### ----------------------------------------------------------

### Let's import some libraries
### ***
### ***
### ***

### ***

import numpy as np
import pandas as pd

import sys
import os


### For matching file names
import fnmatch




####### ----------------------------------------------------------

print("Code started ... ")

####### ----------------------------------------------------------



####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------
####### MAIN CODE   ------------------------------------------
####### ----------------------------------------------------------
####### ----------------------------------------------------------

####### ----------------------------------------------------------

### Folder to look for users performances from weighting users code
folder_of_users_performance_outputs = "/data/users_performance_outputs"



####### ----------------------------------------------------------

### Here we make a list of FCM Export Files

list_of_users_performance_outputs = []

### Find the current directory where the script is running.
### Other folders containing users and events info should be under this directory.
current_folder = os.getcwd()


for file in os.listdir(folder_of_users_performance_outputs):
    if fnmatch.fnmatch(file, "*"+"_users_performace.csv"):
         list_of_users_performance_outputs.append(file)

####### ----------------------------------------------------------

# Let's read all the files of users performances and combine them to one

df_users_performance_all = pd.DataFrame()

for u_p_filename_temp in list_of_users_performance_outputs:
    ### This includes the file name too
    path_u_p_filename_temp = folder_of_users_performance_outputs+"/"+u_p_filename_temp
    
    df_users_performance_temp = pd.read_csv(path_u_p_filename_temp)
    
    df_users_performance_all = df_users_performance_all.append(df_users_performance_temp)
    
    
    
df_users_performance_all.reset_index(drop=True, inplace= True)


####### ----------------------------------------------------------

### Let's loop over players and find average

df_avg_player_score = pd.DataFrame(columns=["player", "avg_score_F1", "avg_score_ARI", "avg_MMOS_score"])


for id_player in list(df_users_performance_all["player"]):
    
    df_for_player = df_users_performance_all[df_users_performance_all["player"]==id_player]
    
    avg_player_our_score_F1 = df_for_player["our_score_F1"].mean()
    avg_player_our_score_ARI = df_for_player["our_score_ARI"].mean()
    avg_player_MMOS_score = df_for_player["MMOS_score"].mean()
    
    
    df_avg_player_score = df_avg_player_score.append({"player":int(id_player) , "avg_score_F1":avg_player_our_score_F1, "avg_score_ARI":avg_player_our_score_ARI, "avg_MMOS_score":avg_player_MMOS_score }, ignore_index=True)


### Make sure "player" values are integer and not float etc.
df_avg_player_score["player"] = df_avg_player_score["player"].astype(int)

### Let's sort the data frame based on "player" befor saving
df_avg_player_score.sort_values("player", inplace= True)

### Save the output
df_avg_player_score.to_csv('/code/support_files/average_player_performace_results.csv', index=False)
 


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













