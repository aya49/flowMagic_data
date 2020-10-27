from os import path, listdir, makedirs
from csv import reader, writer, QUOTE_NONNUMERIC
from json import load, dump
import shutil
import pandas as pd

# input folder name, change this to whatever the input data folder is called
INPUT_FOLDER_NAME = "/mnt/f/Brinkman group/COVID/data/structure_test/user_data_by_ccp"
# the name of the folder the output will be in
OUTPUT_FOLDER_NAME = "/mnt/f/Brinkman group/COVID/data/structure_test/users"
# json file with how many repetitions of each filename, and any that got to 100
SAVE_FILE_NAME = "filepaths_count.json"
# Storage destination for analyzed files:
STORAGE_FOLDER_NAME = "/mnt/FCS_local2/COVID/user_data_by_ccp/"
# the unicode symbol for the euro will be used to replace any forward slash
# characters in a file path. That way, way, a csv file can be saved with the
# name as the file path (a forward slash can't be in the file name). To get the
# file path from the file name, replace any euro unicode symbols with a
# forward slash
euro_unicode = "\u20A0"

events_csv = pd.read_csv('/code/files_to_analyze.csv')
print(events_csv.loc[0]['x'])
