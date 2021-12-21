#!/usr/bin/python3.5

###########################################################################################
# Script name: get_mgrast_data.py
# Developed by: Haris Zafeiropoulos
# Framework: PREGO - WP2
###########################################################################################
# GOAL:
# This script aims to get all the info needed from the MG-RAST database 
# with respect to environmental metadata, project description, species and KO terms found.
###########################################################################################
# Notes:
# In this script we are making use of the multiprocessing Python library that allows us to ask
# for multiplous queries in parallel
###########################################################################################
#
# Usage: ./get_mgrast_data.py
#
###########################################################################################

# Dependencies
from  functions_mgrast import obj_from_url
import time, os, random, json, sys, datetime, re
import multiprocessing
from functions_mgrast import *

# To know how long it takes totally
start = datetime.datetime.now()
date =  datetime.date.today()
old_stdout = sys.stdout
log_file_path = "/data/databases/scripts/gathering_data/logfiles/" 
log_file_name = log_file_path + "get_mg_rast_data_"+ str(date) + ".log"
log_file = open(log_file_name,"w")
sys.stdout = log_file
print("script execution stared at:", start)

# URL to get the minimum information regarding all the MG-RAST projects that are publically available 
all_projects_url = "http://api.mg-rast.org/project?&order=name&verbosity=verbose&limit=max" 

# The 'obj_from_url' function of the mglib library, allows to get the URL needed
read_total_projects = obj_from_url(all_projects_url, auth=None, data=None, debug=False, method=None)
time.sleep(3)


# We need to be sure that we do not keep the same info multiple times due to the frequent runs of our robots
# Thus, for MG-RAST we keep the previous version as archive (REMEMBER TO ADD THIS IN THE ROBOT SCRIPT) and
# we get the whole info each time from scratch

amplicon_file = "/data/databases/mg_rast/amplicon/amplicon_metadata_retrieved.tsv"
wgs_file = "/data/databases/mg_rast/wgs/wgs_metadata_retrieved.tsv"
metabarcode_file = "/data/databases/mg_rast/metabarcode/metabarcode_metadata_retrieved.tsv"
mt_file = "/data/databases/mg_rast/mt/mt_metadata_retrieved.tsv"
unknown_file = "/data/databases/mg_rast/unknown/unknown_metadata_retrieved.tsv"
text_mining_file_name = "/data/databases/mg_rast/mgrast_metadata_for_text_mining.tsv"
projects_included_on_mgrast = "/data/databases/scripts/gathering_data/mg_rast/projects_included_in_mgrast.tsv"
retrieved_projects_file = "/data/databases/scripts/gathering_data/mg_rast/projects_retrieved_from_mgrast.tsv"

files_to_check = [amplicon_file, wgs_file, metabarcode_file, mt_file, unknown_file, text_mining_file_name, projects_included_on_mgrast, retrieved_projects_file]

for file in files_to_check:
    with open(file, "w+") as tsv:
        print("The" + file + " file has been created.")
        tsv.close()


# We make a list with the items that our function will run for
project_ids = []
for project in read_total_projects['data']:
    project_ids.append(project)
    with open(projects_included_on_mgrast, "a") as temp:
        project_id = project['id']
        temp.write(project_id + "\n")
        temp.close()

# Main code 
if __name__ == '__main__':
    
    # Set the maximum number of processes that our script is going to use
    number_processes = 3
    
    # Start a pool
    pool = multiprocessing.Pool(number_processes) 

    # Highly important line:
    # the command pool.map starts running the multiprocess task
    # with a function that will run and a variable that represents the iteration set
    missing_projects = pool.map(handle_a_project, project_ids)
    pool.close()
    pool.join()
    
    counter = 1 
    print("\n>>This is the first iteration of the script.")
    print(missing_projects)
    
    while len(missing_projects) > 0:
        counter += 1        
        if counter < 10:

            new_pool = multiprocessing.Pool(number_processes) 
            missing_projects = set(pool.map(handle_a_project, missing_projects))
            
            new_pool.close()
            new_pool.join()
            
            print("\n>>This is the " + str(counter) + "th iteration of the script. And here are the still missing projects")
            print(missing_projects)
            
        else:
            print("\n\n*** Some of the projects included in MG-RAST failed to retrieved. These are the followings. ***")
            print(missing_project)
        
    print("\n\n***All the projects have been retrieved. ***")

# write and close the log file     
end = datetime.datetime.now()
print("Script execution ended at:", end)
total_time = end - start
print("Script totally ran for :", total_time)
sys.stdout = old_stdout
log_file.close()
