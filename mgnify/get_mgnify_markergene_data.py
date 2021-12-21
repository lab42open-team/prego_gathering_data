#!/usr/bin/python3.5

########################################################################################
# script name: get_mgnify_markergene_data.py
# path on oxygen: /data/databases/scripts/gathering_data/mgnify
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL
# this script aims to get all the files needed from the MGnify database 
# with respect to amplicon studies, specifically to 16S rRNA marker gene studies.
########################################################################################
## usage: ./get_mgnify_markergene_data.py
## potential usage in the PREGO platform: as part of the mgnify_marker_gene_data.sh robot (currently located at /home/haris/prego/robots)
########################################################################################

# import libraries needed
import requests, time, datetime
import json, re, os, sys, traceback, asyncio
import concurrent.futures
from requests.exceptions import Timeout
from requests.adapters import HTTPAdapter
from requests.exceptions import ConnectionError
from mgnify_functions import get_page, load_json_file, retry, download_pages, find_missing_urls, download_url

start = datetime.datetime.now()
date =  datetime.date.today()
old_stdout = sys.stdout
log_file_path = "/data/databases/scripts/gathering_data/logfiles/" 
log_file_name = log_file_path + "mgnify_markergene_data_gathering_"+ str(date) + ".log"
log_file = open(log_file_name,"w")
sys.stdout = log_file

# keep in variables all the paths needed for this script
mgnify_wd = "/data/databases/mgnify/"
mgnify_markergene_directory = mgnify_wd + 'marker_gene_data/'
samples_directory = mgnify_markergene_directory + 'samples/'
runs_directory = mgnify_markergene_directory + 'runs/'

##########################################################################################################################
## STEP 1
## download all the "sample" pages of the "amplicon" label of the "experiment type" and keep the samples IDs 
##########################################################################################################################

# this is the main URL for all the samples from amplicon experiments.
url="https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/amplicon/samples?page=1"
# get the first page in order to get the total number of pages - do not download, just get it
data = ""
while data == "":
    samples = requests.get(url = url)
    try:
        # read the page
        data = samples.json()
    except:
        print("The first page of samples was not retrieved properly. This task is about to run again! \n")
        pass

# get the total number of pages 
number_of_pages = data['meta']['pagination']['pages']
print("Number of marker gene data pages to be retrieved: " + str(number_of_pages) + "\n")

# get all the sample urls in a list
counter = 0
sample_urls = []
while counter < number_of_pages :
    counter += 1
    url= "https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/amplicon/samples?page=" + str(counter)
    sample_urls.append(url)

# Use an executor, built thanks to the concurrent.futures package, using the number of threads you want to, in order to perform your download!
print("The samples executor is about to start. \n" + str(datetime.datetime.now()))
with concurrent.futures.ThreadPoolExecutor(max_workers = 5) as executor:
    print("The executor for downloading samples pages just started! \n")
    # Start the load operations and mark each future with its URL. HERE is where the concurrent.futures part is taking place!
    future_to_url = {executor.submit(download_pages, url, samples_directory, "_sample.json", 15): url for url in sample_urls}

# The urls that failed, they have been kept in a temp file. Read that file and keep each url as a list element.
temp_file = "/data/databases/mgnify/marker_gene_data/temp_file_with_failed_urls.tsv"
for i in range(1, number_of_pages + 1):
    item = str(i) + "_sample.json"
    file = "/data/databases/mgnify/marker_gene_data/samples/" + item
    try:
        f = open(file)
    except IOError:
        print("File " + file + " was not accessible in the first attempt! \n")
        with open(temp_file, "a") as temp:
            missing_url= "https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/amplicon/samples?page=" + str(i) + "\n"
            temp.write(missing_url)
            temp.close()
    finally:
        f.close()

# make a list with all the missing urls
failed_sample_urls = [line[:-1] for line in open(temp_file, "r").readlines()]
print("The number of the missing urls is equal to " + str(len(failed_sample_urls)) + "\n")

# for as long as there are still missing urls..
while len(failed_sample_urls) > 0 :
    print("We still have " + str(len(failed_sample_urls)) + " urls to get! \n")
    # .... perform the "retry" function. 
    retry(failed_sample_urls, samples_directory, "_sample.json")

# erase everything from the temp_file without removing it - in order to use it in the next run of the robot.
#open(temp_file, 'w').close()
# HZ: After discussed the robots architecture, it is better to remove this temp file after getting all the samples needed. 
os.remove(temp_file)

#########################################################################################################
## STEP 1.b
## count what you got!          
#########################################################################################################

# from the files returned in the previous block of code, keep all MGnify study IDs (i.e MGYS89342)
mgnify_ids_from_samples = []
# read one by one all the files that were downloaded from the previous task
counter_for_sample_files = 1
for filename in os.listdir(samples_directory):
    # get the path for this specific file
    file = samples_directory + filename
    counter_for_sample_files += 1
    # open it as a json file
    with open(file) as f:
        data = json.load(f)
    # for each entry in the 'data' part
    for datum in data['data']:
        # for each entry in this certain study
        for study in datum['relationships']['studies']['data']:
            # get the MGnify ID of the study
            mgnify_id = study['id']
            # if this ID appears for the first time, then add it on the list of studies reached
            if mgnify_id not in mgnify_ids_from_samples:
                mgnify_ids_from_samples.append(mgnify_id)

print("The number of the MGYS ids found is: " + str(len(mgnify_ids_from_samples)) + "\n")

timepoint_1 = datetime.datetime.now()
first_step_durance = timepoint_1 - start
print("The first step took: ", first_step_durance)

#########################################################################################################
## STEP 2
## download all the amplicon RUNs from MGnify          
#########################################################################################################

# this is the main URL for all the runs from amplicon experiments.
url_run="https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/amplicon/runs?page=1"

# likewise, I need to keep the number of the pages with runs
get_page_for_runs = get_page(url_run, 30)

# check if the page was open 
while True:
    # if yes, keep the number of pages with runs
    if get_page_for_runs.status_code == 200:
        page_for_runs = load_json_file(get_page_for_runs)
        number_of_run_pages = page_for_runs['meta']['pagination']['pages']
        print("The number of pages with RUNs is:" + str(number_of_run_pages) + "\n")
        break
    else:
        print("The first page of samples was not retrieved properly. This task is about to run again! \n")
        get_page_for_runs = get_page(url_run)

# make a list with all the urls with runs that we will need to download
urls_for_runs = []
counter = 0 
while counter < number_of_run_pages:
    counter += 1
    url= "https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/amplicon/runs?page=" + str(counter)
    urls_for_runs.append(url)
print("These were the urls we are about to download! \n")

# start getting the RUNs pages 
while len(urls_for_runs) > 0:
    
    url = urls_for_runs[-1]
    path_to_save = "/data/databases/mgnify/marker_gene_data/runs/"
    prefix = "_run.json" ; sleeptime = 0.5 ; sec = 30
    check = download_pages(url, path_to_save, prefix, sec)
    
    if check == "ok":
        del urls_for_runs[-1]
    else:
        urls_for_runs = urls_for_runs[-1:] + urls_for_runs[:-1]

print("Step 1 has been completed successfully \n")

#########################################################################################################
## STEP 2.b
## count what you got!          
#########################################################################################################


# from the files returned in this step (2), keep again all MGnify study IDs (i.e MGYS89342)
print("Now the MGnify study IDs are going to be kept in order to find the corresponding OTU tables later.")

mgnify_ids_from_runs = []
# for each of the pages corresponding to runs from amplicons
for filename in os.listdir(runs_directory):
    file = runs_directory + filename
    with open(file) as f:
        data = json.load(f)
    # for each entry on this specific page    
    for entry in data['data']:
        # get the MGnify ID
        mgnify_id = entry['relationships']['study']['data']['id']
        # if that id has not appeared before, then add it to the runs list
        if mgnify_id not in mgnify_ids_from_runs:
            mgnify_ids_from_runs.append(mgnify_id)

print("Number of MG IDs from samples:  " + str(len(mgnify_ids_from_samples)) + "\n")
print("Number of MG IDs from runs:  " + str(len(mgnify_ids_from_runs)) + "\n")

# Here I check whether the MGnify IDs are the same from both sources
if set(mgnify_ids_from_runs) == set(mgnify_ids_from_samples):
    mgnify_unique_ids = mgnify_ids_from_samples
    print("The MGnify IDs are identical from the Sample- and the Runs- downloads\n")
else:
    mgnify_ids = mgnify_ids_from_runs + mgnify_ids_from_samples
    mgnify_unique_ids = []
    for item in mgnify_ids:
        if item not in mgnify_unique_ids:
            mgnify_unique_ids.append(item)
    print("The MGnify IDs returned from the sample and the run query were not the same.\
          You may check the files in those two folders for differences\n")


print("Step 2 has been completed successfully \n")
timepoint_2 = datetime.datetime.now()
second_step_durance = timepoint_2 - timepoint_1
print("The first step took: ", second_step_durance)

########################################################################################################
## STEP 3
## get all the taxonomy tables of the studies you got from the previous steps     #######
########################################################################################################

# again, check for the "abundances" directory this time. if not there, then make it.
abundances_directory = mgnify_markergene_directory + 'abundances/'
processes_directory = mgnify_markergene_directory + 'processes/'
check_for_directory = os.path.isdir(processes_directory)
x = False; y = (bool(x))
if check_for_directory == y:
    os.mkdir(processes_directory)

# for each study MGnify ID of those recorded 
for study in mgnify_unique_ids:

    # read the "downloads" page of this specific study
    url="https://www.ebi.ac.uk/metagenomics/api/v1/studies/" + study + "/downloads"
    downloads = requests.get(url = url, timeout = 30) 
    files = downloads.json()
    
    # for each entry in the page
    for entry in files['data']:
        
        # check if it corresponds to a taxonomic assignments entry
        if entry['attributes']['description']['label'] == "Taxonomic assignments" :
            # if it does, download the OTU table
            otu_url = entry['links']
            otu_table = requests.get(otu_url['self'], allow_redirects = True)
            filename = abundances_directory + study + "_OTU_table.tsv"
            open(filename, 'wb').write(otu_table.content)
            
        # likewise, check whether it corresponds to a Gene Ontology annotation    
        elif entry['attributes']['description']['label'] == "Complete GO annotation" :
            processes_url = entry['links']
            processes = requests.get(processes_url['self'], allow_redirects = True)
            filename = processes_directory + study + "GO_annotations.tsv"
            open(filename, 'wb').write(processes.content)

print("Step 3 has been completed successfully \n")
timepoint_3 = datetime.datetime.now()
third_step_durance = timepoint_3 - timepoint2
print("The third step took: ", third_step_durance)

########################################################################################################
## Downloads have been completed. End of script by keeping the time it took to run & exiting the log file
########################################################################################################

print("Script execution stared at:", start)
end = datetime.datetime.now()
print("Script execution ended at:", end)
total_time = end - start
print("Script totally ran for :", total_time)

sys.stdout = old_stdout
log_file.close()
sys.exit(0)



## Comment! We now have everything needed to run Lars' script.
## The parsing for that is done in his script - we just need to fill the "runs", "samples" and "abundances" folders.


#---------------------------------------------------------------------------------------------------------------------------------------------

### Attention!
### Another way for getting the RUNs via an executor (step 2)

# # Likewise, a with statement to ensure threads are cleaned up promptly
# with concurrent.futures.ThreadPoolExecutor(max_workers = 1) as executor:
#     print("the executor for downloading run pages just started! \n")
#     # HERE is where the "concurrent.futures" part is taking place!
#     future_to_url = {executor.submit(download_pages, url, runs_directory, "_run.json", 40): url for url in urls_for_runs}
# 
# print("the executor for downloading run pages is done! \n")
# 
# # get the urls that were not retrieved in the first place
# suffix_for_runs = "_run.json"
# url_for_runs = "https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/amplicon/runs?page="
# find_missing_urls(number_of_run_pages, suffix_for_runs, runs_directory, url_for_runs)

# # read the temp file with the missing RUNs urls, make a list with those and get them all.
# temp_file = "/data/databases/mgnify/marker_gene_data/temp_file_with_failed_urls.tsv"
# with open(temp_file, "r") as temp:
#     run_urls_still_need_to_get = temp.readlines()
#     
#     while len(run_urls_still_need_to_get) > 0:
#         url = run_urls_still_need_to_get[-1]
#         url = url[:-1]
#         path_to_save = "/data/databases/mgnify/marker_gene_data/runs/"
#         prefix = "_run.json" ; sleeptime = 0.5 ; sec = 30
#         download_url(url, path_to_save, prefix, sleeptime, sec)
#         del run_urls_still_need_to_get[-1]
# print("all runs have been downloaded! \n")