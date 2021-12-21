#!/usr/bin/python3.5

########################################################################################
# script name: get_mgnify_metatranscriptome_data.py
# developed by: Haris Zafeiropoulos and the PREGO team
# framework: PREGO - WP2, http://prego.hcmr.gr
########################################################################################
# GOAL
# This script aims to get all the samples from metatranscriptome studies that are available on MGnify. 
# And then get all the annotation files regarding those, that can be usefull for extracting PREGO associations.
########################################################################################

########################################################################################
# latest updates
# HZ - 2020.06.23 - .....
# EP - 2020.06.23 - all output file records contain 7 columns, empty for the columns with no value
########################################################################################


import requests, traceback
import json, re, os, sys
import time, datetime, random
from mgnify_functions import first, download_samples, metatranscriptome_samples, get_download_files, download_downloads, load_url

# create a logfile for this script
start = datetime.datetime.now()
date =  datetime.date.today()
old_stdout = sys.stdout
log_file_name = "logfile_from_mgnify_metatranscriptome_downloads"+ str(date) + ".log"
log_file = open(log_file_name,"w")
sys.stdout = log_file

# present working directory
mgnify_wd = "/data/databases/mgnify/"
mgnify_metatranscriptomic_directory = mgnify_wd + "metatranscriptomic_data/"

# keep the path for the metatranscriptome samples
mgnify_metatranscriptomic_samples_directory = mgnify_metatranscriptomic_directory + "samples/"

# and the path of the dowloands that will be retrieved
mgnify_metatranscriptomic_downloads_directory = mgnify_metatranscriptomic_directory + "downloads/"
	
# here is where you can set the "sleep" time between two GET requests
time_for_sleep = random.uniform(0.5,1.5)


#####################################################################################################################
#### Part A
#### Getting the samples that come from the metatranscriptome experiment-type
#####################################################################################################################

# GOAL
# aim of the first part of this script is to find the studies that their experiment type according to MGnify is "metatranscriptome"
# and at the same time, to keep the biome of each of these studies.
# its output is a tab separated file where each line has the MGnify ID number of the study and its biome.
# e.g
# MGYS00001932	root:Environmental:Aquatic:Marine
				
# Now get all the metatranscriptome samples
metatranscriptome_samples()

# get the MGnify IDs
studies_biome = []
for filename in os.listdir(mgnify_metatranscriptomic_samples_directory):
	file = mgnify_metatranscriptomic_samples_directory + filename
	with open(file) as f:
		data = json.load(f)
	# keep in the first column the MGYS (MGnify id) number and in the second the biome
	for entry in data['data']:
		pair = entry['relationships']['studies']['data'][0]['id'], entry['relationships']['biome']['data']['id']
		if pair not in studies_biome:
			studies_biome.append(pair)

mgys_biome_path = mgnify_metatranscriptomic_directory + "mgnify_metatranscriptomic_data_retrieval.txt"
mgys_biome_create = open(mgys_biome_path, "w")
mgys_biome_create.close()

#####################################################################################################################
#### Part B
#### Get the annotation files that come from each analysis of those samples
#####################################################################################################################

# GOAL
# Aim of this part of the script is to download the metadata files for each of the studies retrieved in PART A.
# For each of those studies a directory is going to be built, on the "downloads" directory,  and all the related files will be saved in the corresponding study.
# Finally, for each study we will have something like this:
# haris@oxygen: ls metatranscriptomics/downloads/MGYS00003727
# ERP104335_GO_abundances_v4.1.tsv		ERP104335_taxonomy_abundances_SSU_v4.1.tsv	ERP104335_IPR_abundances_v4.1.tsv	MGYS00003727_downloads.json	ERP104335_taxonomy_abundances_LSU_v4.1.tsv	merged_taxa.tsv
#
# In all cases that there are both SSU and LSU analyses, a new file called "merged_taxa.tsv" is built and contains all the recorded taxa from both those analyses.

# Keep in a list all the studies we need to get their downloads 
studies = []
with open(mgys_biome_path) as f:
	for line in f:
		study = line.split("\t")[0]
		studies.append(study)

annotations_list = ["Complete GO annotation", "OTUs and taxonomic assignments for LSU rRNA", "OTUs and taxonomic assignments for LSU rRNA","OTUs and taxonomic assignments for SSU rRNA", "OTUs and taxonomic assignments for SSU rRNA", "OTUs and taxonomic assignments"]

for filename in os.listdir(mgnify_metatranscriptomic_samples_directory):
	filename = mgnify_metatranscriptomic_samples_directory + filename
	with open(filename, "r") as sample_page:
		sample_page = json.load(sample_page)
		for sample in sample_page['data']:
			sample_id = sample['id']
			print(sample_id)
			study_id = sample['relationships']['studies']['data'][0]['id']
			
			# metadata - standar
			biome = sample['attributes']['environment-biome']
			feature = sample['attributes']['environment-feature']
			material = sample['attributes']['environment-material']
			description = sample['attributes']['sample-desc']

			biome_2 = None ; feature_2 = None ; material_2 = None
			# metadata - sample metadata 
			for entry in sample['attributes']['sample-metadata']:
				if entry['key'] == "environment (biome)":
					biome_2 = entry['value']
				elif entry['key'] == "environment (feature)":
					feature_2 = entry['value']
				elif entry['key'] == "environment (material)":
					material_2 = entry['value']
				
			if biome is None and biome_2 is None:
				print("there is no biome for this sample. i have to move to the next! \n")
				continue
			# check the biome, feauture and material variables and finalize what we keep
			else:
				if biome == biome_2:
					final_biome = biome
				elif biome is not None:
					final_biome = biome
				elif biome_2 is not None:
					final_biome = biome_2
				
				if feature == feature_2:
					final_feature = feature
				elif feature is not None:
					final_feature = feature
				elif feature_2 is not None:
					final_feature = feature_2
				
				if material == material_2:
					final_material = material
				elif material is not None:
					final_material = material
				elif material_2 is not None:
					final_material = material_2
				
				# get the url for the runs 
				runs_url = sample['relationships']['runs']['links']['related']
				runs_page = load_url(runs_url, 30)
				if runs_page.status_code == 200:
					runs_page = runs_page.json()
					
					# get the analysis page for this run
					for datum in runs_page['data']:
						analysis_url = datum['relationships']['analyses']['links']['related']
						analysis_page = load_url(analysis_url, 30)
						if analysis_page.status_code == 200:
							analysis_page = analysis_page.json()
							
							# get the analysis id
							for analysis in analysis_page['data']:
								mgya_id = analysis['id']
								
								# here, we need to write the metadata found in the .txt file
								if description == None:
									new_entry_in_the_metadata_txt_file = mgya_id + "\t" + sample_id + "\t" + study_id + "\t" + final_biome + "\t" + final_feature + "\t" + final_material + "\t" + ""
								else:
									new_entry_in_the_metadata_txt_file = mgya_id + "\t" + sample_id + "\t" + study_id + "\t" + final_biome + "\t" + final_feature + "\t" + final_material + "\t" + description
								
								metadata_txt_file = open(mgys_biome_path, "a")
								metadata_txt_file.write(new_entry_in_the_metadata_txt_file)
								metadata_txt_file.write("\n")
								metadata_txt_file.close()
								
								# make a directory for this mgya id
								directory_to_save = mgnify_metatranscriptomic_downloads_directory + mgya_id
								check_for_directory = os.path.isdir(directory_to_save)
								x = False; y = (bool(x))
								if check_for_directory == y:
									os.mkdir(directory_to_save)
								
								downloads_url = analysis['relationships']['downloads']['links']['related']
								downloads_url = downloads_url + "?page_size=100"
								
								downloads_page = load_url(downloads_url, 30)
								downloads_page = downloads_page.json()
								downloads_page_number = downloads_page['meta']['pagination']['pages']
								
								if downloads_page_number == 1:
									for attribute in downloads_page['data']:
										if attribute['attributes']['description']['description'] in annotations_list:
											link = attribute['links']['self']
											get_annotation_link = load_url(link, 30)
											sleepy = random.uniform(0.3,0.6)
											time.sleep(sleepy)
											print(get_annotation_link)
											if get_annotation_link == None:
												print("we could not download this annotation file: " + link)
												continue
											else:
												file_to_save = link.split("/")[-1]
												filename = directory_to_save + "/" + file_to_save
												open(filename, 'wb').write(get_annotation_link.content)
								else:
									# counter = 1
									for page in range(1, downloads_page_number + 1):
										downloads_url = downloads_url + "?page=" + str(page) + "&page_size=100"
										downloads_page = load_url(downloads_url, 30)
										downloads_page = downloads_page.json()
						
										for attribute in downloads_page['data']:
											if attribute['attributes']['description']['description'] in annotations_list:
												link = attribute['links']['self']
												get_annotation_link = load_url(link, 30)
												print(get_annotation_link)
												sleepy = random.uniform(0.3,0.6)
												time.sleep(sleepy)
												if get_annotation_link == True:
													continue
												else:
													file_to_save = link.split("/")[-1]
													filetype = file_to_save.split(".")[-1]
													print(filetype)
													excluded = ["hdf5", "HDF5"]
													for pattern in excluded:
														if re.search(pattern, excluded):
															continue													
														else:
															filename = directory_to_save + "/" + file_to_save
															open(filename, 'wb').write(get_annotation_link.content)

print("script execution stared at:", start)
end = datetime.datetime.now()
print("Script execution ended at:", end)
total_time = end - start
print("Script totally ran for :", total_time)

sys.stdout = old_stdout
log_file.close()
