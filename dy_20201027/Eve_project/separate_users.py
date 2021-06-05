from os import path, listdir, makedirs
from csv import reader, writer, QUOTE_NONNUMERIC
from json import load, dump
import shutil

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

# get the folder this file is in, and the data input folder
project_folder = path.dirname(__file__)
data_folder = path.join(project_folder, INPUT_FOLDER_NAME)
output_folder = path.join(project_folder, OUTPUT_FOLDER_NAME)
storage_folder = STORAGE_FOLDER_NAME
# create a folder for the output, if it doesn't exist yet
try:
    makedirs(output_folder)
except FileExistsError:
    # folder already exists
    pass

# open the json file containing how many users have done each file
print(f"Looking for a file with saved filepath data called: {SAVE_FILE_NAME}")
try:
    analyzed_filepaths = load(open(SAVE_FILE_NAME, "r"))
    print("Found save data for analyzed file data.")
# create the json file, as it has not been created yet
except (OSError, IOError) as e:
    # no save file has been created yet
    # create a new empty dictionary to store save data in
    # "amount" is a dictionary of how many times a filepath has been analyzed
    # the key is the filepath, and the value is the number of analyzations
    # "in progess" is all filepaths that have appeared less than 100 times
    # "exact" is all filepaths that have appeared exactly 100 times
    # "over" is  all filepaths that have appeared over 100 times
    # "csv files" keeps track of which csv files have already been passed in
    # this ensures a file is not passed in more than once
    analyzed_filepaths = {"over": [],
                          "exact": [],
                          "in progress": [],
                          "amounts": {},
                          "csv files": []
                          }
    print("No save data for analyzed file data was found.")

# loop through all files in input folder
input_files = listdir(data_folder)
print(f"\nGoing through {len(input_files)} files in {data_folder}")
print("This may take some time...\n")
for file in input_files:
    # print out how many files have been tested, and how many are left
    print(f"{input_files.index(file) + 1}/{len(input_files)}", end=" ")

    # ignore files in the input folder that don't have the correct extension
    if file.endswith(".csv.part_00000"):

        # make sure this file has not been passed in before
        if file not in analyzed_filepaths['csv files']:
            print(f"Getting data from file: {file}")

            # add the csv file to the save data, to make sure the same file is
            # not run again
            analyzed_filepaths['csv files'].append(file)
            print(data_folder)
            print(file)
            print(path.join(data_folder, file))
            # open the csv file and go through the data
            with open(path.join(data_folder, file), "r") as f:
                # read the csv file
                # specify the separator and surrounding character of values
                data = reader(f, delimiter=",", quotechar="'")

                # this list will keep track of the first row, which is the
                # fieldnames (columns) of the csv file
                csv_fieldnames = []
                # loop through each row in the csv file
                for i, row in enumerate(data):
                    # don't do the first row, as it is the column names
                    if i != 0:
                        # keep track of how many people analyzed a file path
                        # the file path has already been added, so add 1 to it
                        if row[0] in analyzed_filepaths['amounts']:
                            # found the filepath again, add 1 to the counter
                            analyzed_filepaths['amounts'][row[0]] += 1

                            # if it has reached 100, add it to the exact list
                            if analyzed_filepaths['amounts'][row[0]] == 100:
                                analyzed_filepaths['in progress'].remove(
                                    row[0])
                                analyzed_filepaths['exact'].append(row[0])

                            # if it is 101, remove it from the exact list,
                            # and add it to the over list (as it is 100+)
                            if analyzed_filepaths['amounts'][row[0]] == 101:
                                analyzed_filepaths['exact'].remove(row[0])
                                analyzed_filepaths['over'].append(row[0])
                        else:
                            # a new file path is found
                            analyzed_filepaths['amounts'][row[0]] = 1
                            analyzed_filepaths['in progress'].append(row[0])

                        # insert the user number to the end of the csv file
                        # path starting at _u1, it will go up by one each time
                        # another file path is found
                        user_filename = list(row[0])
                        user_filepath = "".join(user_filename)
                        user_filename.insert((-4), ("_u" + str(
                            analyzed_filepaths['amounts'][row[0]])))
                        user_filename = "".join(user_filename)

                        # replace forward slash characters in the file path
                        # that way a csv file can be saved, where the
                        # filename of the csv file is the file path a user
                        # completed. Later, to get the file path from the
                        # saved csv's filename, replace all euro unicode
                        # symbols with a forward slash
                        user_filename = user_filename.replace("/",
                                                              euro_unicode)
								# create a folder with for the user, if it does not exist
                        user_filepath = user_filepath.replace("/",euro_unicode)
                        user_filepath = user_filepath.replace(".csv","")
                        user_output_folder = path.join(output_folder, user_filepath)
                        try:
                        	makedirs(user_output_folder)
                        except FileExistsError:
                        	# folder already exists
                        	pass
                        # save the analyzed filename as a csv file
                        user_filepath = path.join(user_output_folder, user_filename)
                        with open(user_filepath, "w") as csv_file:
                            # write the csv file for the user
                            csv_writer = writer(csv_file,
                                                quoting=QUOTE_NONNUMERIC,
                                                delimiter=",", quotechar="'")
                            # the fieldnames
                            csv_writer.writerow(csv_fieldnames)
                            # the data for this user
                            csv_writer.writerow(row)

                    # the first row is the columns
                    else:
                        csv_fieldnames = row
				# move file to storage folder...
            old_filename = path.join(data_folder, file)
            new_filename = path.join(storage_folder, file)
            shutil.move(old_filename, new_filename)  
        else:
            # this csv file has already been passed in
            print("Ignoring data from file already passed in during previous "
                  f"run: {file}")
    else:
        # file does not end with ".csv.part_00000", so it will not be checked
        print("Ignoring file which does not end in \".csv.part_00000\": "
              f"{file}")

# save the number of people who completed challenged into a json file
dump(analyzed_filepaths, open(SAVE_FILE_NAME, "w"))
print(f"\nSaved analyzed file data to the json file: {SAVE_FILE_NAME}")
