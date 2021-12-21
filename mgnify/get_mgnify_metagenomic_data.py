#!/usr/bin/python3.5

########################################################################################
# script name: get_mgnify_metagenomic_data.py
# path on oxygen: /data/databases/scripts/gathering_data/mgnify
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL
# This script aims to get all the files needed from the MGnify database with respect to metagenomic studies.
# MGnify seperates "metagenomic" from "assemblies". In the latter, there are Whole Genome Shotgun sequencing reads that have been assembled to contigs.
# For the PREGO  project we decided to get just the "metagenomic" in order to have the same bias to all our data.
# Each study in the "assemblies" case, brings its own biases as it has been analyzed through a different pipeline.
########################################################################################

# import libraries needed
import requests, time, datetime, traceback
import json, re, os, sys 
import asyncio
import concurrent.futures
from mgnify_functions import get_the_job_done, load_url

# create a log file 
start = datetime.datetime.now()
print("script execution stared at:", start)
date =  datetime.date.today()
old_stdout = sys.stdout
log_file_path = "/data/databases/scripts/gathering_data/logfiles/" 
log_file_name = log_file_path + "mgnify_metagenomic_data_gathering"+ str(date) + ".log"
log_file = open(log_file_name,"w")
sys.stdout = log_file

#set the sleep time for the "get" requests
time_for_sleep = 1

# in the first place, we need to get the samples that can be found in MGnify and come from metagenomic experiments
mgnify_workdir = "/data/databases/mgnify/"
mgnify_metagenomic_workdir = mgnify_workdir + "metagenomic_data/"
mgnify_metagenomic_samples_workdir = mgnify_metagenomic_workdir + "samples/"
mgnify_metagenomic_annotations_workdir = mgnify_metagenomic_workdir + "annotations/"

# we need to get the number of pages that are related with metagenomic samples
# this is the main URL for all the samples from metagenomic experiments.
url = "https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/metagenomic/samples"

# get the first page in order to get the total number of pages - do not download, just get it
check_variable = False
while check_variable == False:
	get_page = requests.get(url, timeout = 30)
	if get_page.ok == True:
		data = get_page.json()
		check_variable = True
		print("i have the samples page!")
	else:
		print("i did not get the page asked. i am about to retry!")

# get the number of pages
number_of_pages = data['meta']['pagination']['pages']
print("the number of pages is: " + str(number_of_pages))

########################################################################################
# here is the trick!
# we need to get annotations sample - by - sample. that is because in a certain study,
# we may have samples from different biomes
# to do that, we first get all samples that are metagenomic related
# we keep the SRS id along with the MGnify STUDY id
# then we do a second query to get the MGnify ANALYSIS id
########################################################################################

# create a text file keeping sample ids, study ids, biome, feature and material
mgnify_metagenomic_ids_and_metadata_file_path = mgnify_metagenomic_workdir + "mgnify_metagenomic_data_retrieval_" + str(date) + ".txt"

# # PART A: we need to get the samples that come from the experiment type of metagenomic

# make a list with all the SRS ids that we will get their taxonomic and functional annotation
list_of_mgya_ids = []
# start reading the pages
counter = 0
annotations_list = ["Complete GO annotation", "OTUs and taxonomic assignments for LSU rRNA", "OTUs and taxonomic assignments for LSU rRNA","OTUs and taxonomic assignments for SSU rRNA", "OTUs and taxonomic assignments for SSU rRNA", "OTUs and taxonomic assignments"]

# # get all the sample pages of MGnify that come from metagenomic experiments!
# # keep those into a directory: "/data/databases/mgnify/metagenomic_data/samples"
urls_samples = []
while counter < number_of_pages:
	counter += 1
	url= "https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/metagenomic/samples?page=" + str(counter)
	urls_samples.append(url)
print("i kept all the urls for the sample pages in a list")

# # We can use a with statement to ensure threads are cleaned up promptly
# with concurrent.futures.ThreadPoolExecutor(max_workers = 4) as executor:
# 	print("i am about to start downloading the sample pages via an executor!")
# 	# Start the load operations and mark each future with its URL. HERE is where the concurrent.futures part is taking place!
# 	future_to_url = {executor.submit(download_url, url, mgnify_metagenomic_samples_workdir, "_sample.json", 0.5, 15): url for url in urls_samples}

########################################################################################
# #  PART B: in the following part, our intention is to use the sample pages we got in  
# #  the previous section in order to get the taxonomic and functional annotation files
########################################################################################

# create a file with ids and metadata (biome, feautre and material)
file_with_ids_and_metadata = open(mgnify_metagenomic_ids_and_metadata_file_path, "w")
# file_with_ids_and_metadata.write("srs\tmgys\turl\n")
file_with_ids_and_metadata.close()

# get each sample file from those we got
samples = []
for sample in os.listdir(mgnify_metagenomic_samples_workdir):
	sample_page = mgnify_metagenomic_samples_workdir + sample
	samples.append(sample_page)

# We can use a with statement to ensure threads are cleaned up promptly
with concurrent.futures.ThreadPoolExecutor(max_workers = 2) as executor:	
	# Start the load operations and mark each future with its URL. HERE is where the concurrent.futures part is taking place!
	future_to_url = {executor.submit(get_the_job_done, sample_page_under_study, list_of_mgya_ids, mgnify_metagenomic_ids_and_metadata_file_path, mgnify_metagenomic_annotations_workdir, annotations_list): sample_page_under_study for sample_page_under_study in samples}

###############################
###   INSTEAD of the second executor for the get_the_job_done function, we can also use the traditional way of a for loop.
# for sample_page_under_study in samples:
	# get_the_job_done(sample_page_under_study, list_of_mgya_ids, mgnify_metagenomic_ids_and_metadata_file_path, mgnify_metagenomic_annotations_workdir, annotations_list)
###############################

# Finally, we need to remove the relationships that are repeated.
# This is why, it is possible for a sample to be in more than one analyses but have the same metadata.
# That would lead to a bias we do not want.






# close the log file
end = datetime.datetime.now()
print("Script execution ended at:", end)
total_time = end - start
print("Script totally ran for :", total_time)
sys.stdout = old_stdout
log_file.close()



