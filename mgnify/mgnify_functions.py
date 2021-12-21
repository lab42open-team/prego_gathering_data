#!/usr/bin/python3.5

########################################################################################
# script name: mgnify_functions.py
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL
# On this script you can find all the functions that are used in all the scripts that have been implemented to get data from the MGnify database.
# This script has 3 main parts; one for each of the scripts made for MGnify.
# Some functions are used in more than one scripts.
########################################################################################
# usage: functions called on the get_mgnify* scripts, by including them on the scripts
# example: from mgnify_functions import get_page, load_json_file
########################################################################################


import time, datetime, random
import json, re, os, sys, traceback
import asyncio, concurrent.futures
import requests
from requests.exceptions import Timeout
from requests.adapters import HTTPAdapter
from requests.exceptions import ConnectionError

#####################################################################################################################
#### Part A
#### Functions for getting the marker gene data
#####################################################################################################################

# just get the content of a page
def get_page(url_to_get, sec):
    page = requests.get(url = url_to_get, timeout = sec)
    return page

# read a json file
def load_json_file(file):
    return file.json()

# get urls from failed_urls list
def retry(failed_urls, directory_to_save, suffix):
    print("i am inside the retry function!")
    sec = random.uniform(0.5,1.5)
    if len(failed_urls) == 0:
        print("All urls are downloaded properly!")
        pass
    else:
        url = failed_urls[0]
        print("url under process: " + url)
        page = requests.get(url, allow_redirects = True, timeout = 130)
        time.sleep(sec)
        try:
            try_page = page.json()
        except:
            print("i had an exception...")
            pass
        if str(page.status_code) == '200':
            print("i have a new file downloaded!")
            number_of_page = url.split("=")[1]
            filename = directory_to_save + str(number_of_page) + suffix
            open(filename, 'wb').write(page.content)
            failed_urls.remove(url)

def download_pages(url, path_to_save, prefix, timeout_time):
    print("A new url is about to start downloading." + str(datetime.datetime.now()))
    # try to get a response for your url
    try:
        samples_page = requests.get(url, allow_redirects = True, timeout = timeout_time)
        sleepy = random.uniform(1.5, 3.5)
        time.sleep(sleepy)
        print(url + "\t" + "status: " + str(samples_page.status_code))
        # if everything is fine, then save the page.
        if str(samples_page.status_code) == '200':
            try_page = samples_page.json()
            suffix = url.split("=")[1]
            filename = path_to_save + suffix + prefix
            open(filename, 'wb').write(samples_page.content)
            print("the corresponding page was saved! \n")
            message = "ok"
            return message
        # if an error was returned then keep the url for the retry step.
        else:
            print("this url: " + url + "returned an error:" + str(samples_page.status_code) + " and has been kept in a temp file.")
    # if timeout is over, then also keep the url for the retry step.
    except Timeout:
        print(url + " got a TIMEOUT! \n")
    # in case that a request exceeds the configured number of maximum redirections
    except request.exceptions.TooManyRedirects:
        print("url: " + url + " ..this is a bad url! no clue how this may happened.. that is a crucial error for this script. ")
    # check for a HTTP client disconnection
    except http.client.HTTPException:
        print("a Connection Error occured while trying to get the " + url + " it was the http.client excpetion that was hit!")
    # if a connection error occurs
    except requests.exceptions.ConnectionError:
        print("a Connection Error occured while trying to get the " + url)
    # in case that any other exception takes place
    except requests.exceptions.RequestException:
        print("a kind of error that i do not understand took place.. please try again later with this url: " + url)

def find_missing_urls(number_of_expected_pages, suffix, path, url_base):
    temp_file = "/data/databases/mgnify/marker_gene_data/temp_file_with_failed_urls.tsv"
    for i in range(1, number_of_expected_pages + 1):
        item = str(i) + suffix
        file = path + item
        try:
            f = open(file)
        except IOError:
            print("File " + file + " was not accessible in the first attempt!")
            with open(temp_file, "a") as temp:
                missing_url= url_base + str(i) + "\n"
                temp.write(missing_url)
                temp.close()
        finally:
            f.close()

#####################################################################################################################
#### Part B
#### Functions for getting the metagenomic data
#####################################################################################################################

# Retrieve a single page and save its contents
def download_url(url, path_to_save, prefix, sleeptime, sec):
    samples_page = requests.get(url, allow_redirects = True, timeout = sec)
    time.sleep(sleeptime)
    if samples_page.status_code == 200:
        suffix = url.split("=")[1]
        filename = path_to_save + suffix + prefix
        open(filename, 'wb').write(samples_page.content)
        print(url + "\t" + str(samples_page.status_code))
    else:
        print("did not make it this time.. " + url + "\t" + str(samples_page.status_code))
        dowload_url(url, path_to_save, prefix, sleeptime, sec)


# Load a url - just read it
def load_url(url, sec, counter = 0):
    new_counter = counter + 1
    response = requests.get(url, timeout = sec)
    sleeping_time = random.uniform(0.5,1.5)
    time.sleep(sleeping_time)
    if response.status_code == 404:
        print("i got a 404 error. we go on with the next sample")
    elif new_counter > 50:
        print("this annotation file probably is not there.. something is going wrong with the MGnify server.")
        error = None
        return error
    elif response.status_code == 200 :
        return response
    else:
        print("i did not manage to load the url you asked for: " + url)
        print("i will try again!")
        print("the counter is: " + str(new_counter))
        load_url(url,sec, new_counter)

# for the "get_mgnify_metagenomic_data.py" script. the major function.                        
# in this function, we use as input a sample page in order to get all the annotation files that correspond to all the samples that this sample page contains
def get_the_job_done(sample_page, mylist, txt_output_file, path_to_save_annotation_files, annotations_list):    
    print(sample_page)
    # open the sample page file
    with open(sample_page) as f:
        file = json.load(f)
    # for each sample in this sample page, let the story begins....
    for entry in file['data']:
        print(" \n we handle a NEW SAMPLE!")
        # get the SRS id of the sample - a sample can be found in more than one analyses
        srs_id = entry['id']
        print(srs_id)
        # we also get the metadata needed - biome, feauture and material, but also the DESCRIPTION when available
        biome = entry['attributes']['environment-biome']
        feature = entry['attributes']['environment-feature']
        material = entry['attributes']['environment-material']
        mgnify_study_id = entry['relationships']['studies']['data'][0]['id']
        sample_description = entry['attributes']['sample-desc']
        
        biome_2 = None ; feature_2 = None ; material_2 = None
        # check the biome in the sample-metadata part, in case the biome was empty
        for record in entry['attributes']['sample-metadata']:
            if record['key'] == "environment (biome)":
                biome_2 = record['value']
            elif record['key'] == "environment (feature)":
                feature_2 = record['value']
            elif record['key'] == "environment (material)":
                material_2 = record['value']
        
        # the biome is something vital for our study; thus if there is no biome in the metadata of this sample, we need to leave it out of our study         
        if biome is None and biome_2 is None:
            print("there is no biome for this sample. i have to move to the next! \n")
            continue
        # as we might have biome, feature or/and material more than once, we finalize their value
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

            # get the RUN page's URL
            runs_of_the_sample_url = entry['relationships']['runs']['links']['related']
            
            # get the RUNs page
            get_runs_of_the_sample = load_url(runs_of_the_sample_url, 20)
            if get_runs_of_the_sample.ok == False:
                print("ERROR! I did not make it to open the run page of this sample. I need to go on with the next sample.")
                continue
            else:
                # open the RUN page
                runs_of_the_sample = get_runs_of_the_sample.json()
                # get the SRR id
                srr_id = runs_of_the_sample['data'][0]['id']
                # go to the study page - might be more than one
                study_url = "https://www.ebi.ac.uk/metagenomics/api/v1/studies/" + mgnify_study_id + "/analyses"
                get_study_of_the_sample = load_url(study_url, 30)
                if get_study_of_the_sample.ok == False:
                    print("ERROR! I did not make it to open the analyses page for this MGnify study id! I need to go on with the next sample.")
                    continue
                else:
                    study_page = get_study_of_the_sample.json()
                    number_of_pages_for_the_study = study_page['meta']['pagination']['pages']
                    
                    # parse all the pages of the study until you get the srs_id under study
                    page_of_analyses_counter = 1
                    while page_of_analyses_counter <= number_of_pages_for_the_study:
                        # really important step! 
                        page_url = "https://www.ebi.ac.uk/metagenomics/api/v1/studies/" + mgnify_study_id + "/analyses?page=" + str(page_of_analyses_counter)
                        page_of_analyses_counter += 1
                        
                        # get the specific page
                        get_page = load_url(page_url, 30)
                        if get_page.ok == False:
                            print("ERROR! Could not get the page for this analysis. I need to pass this sample.")
                            continue
                        else:
                            current_study_page = get_page.json()
                            # find the sample under study 
                            for datum in current_study_page['data']:
                                # there is a possibility that there is not the study id for a certain analysis.. Ignore it and move on.  
                                try:
                                    study_id = datum['relationships']['sample']['data']['id']
                                except KeyError:
                                    print("This entry has no study id. I cannot check whether it is the one I am looking for. I need to move on to the next entry.")
                                    continue
                                
                                # match the srs_id with the same in this page
                                if study_id == srs_id:
                                    mgnify_analysis_id = datum['attributes']['accession']
                                    mylist.append(mgnify_analysis_id)
                                    
                                    # add the new metadata in our tab seperated text file
                                    file_to_write_to = open(txt_output_file, "a")
                                    # make the final new_entry
                                    if final_feature is not None and final_material is not None:
                                        new_entry = mgnify_analysis_id + "\t" + mgnify_study_id + "\t" + study_id  + "\t" + final_biome + "\t" + final_feature + "\t" + final_material
                                        print("we have both feature and material!")
                                    elif final_feature is not None and final_material is None:
                                        new_entry = mgnify_analysis_id + "\t" + mgnify_study_id + "\t" + study_id  + "\t" + final_biome + "\t" + final_feature + "\t" + "null"
                                        print("we have only feature.. ")
                                    elif final_feature is None and final_material is not None:
                                        new_entry = mgnify_analysis_id + "\t" + mgnify_study_id + "\t" + study_id  + "\t" + final_biome + "\t" + "null" + "\t" + final_material
                                        print("we have only material!")
                                    elif final_feature is None and final_material is None:
                                        new_entry = mgnify_analysis_id + "\t" + mgnify_study_id + "\t" + study_id  + "\t" + final_biome + "\t" + "null" + "\t" + "null"
                                        print("i have only biome dude! :( ")
                                    
                                    # in case there is a description of the sample available on metadata
                                    if sample_description is not None:
                                        new_entry  = new_entry + "\t" + sample_description
                                    
                                    # WRITE your findings to the .txt file that is already opened
                                    file_to_write_to.write(new_entry)
                                    file_to_write_to.write("\n")
                                    file_to_write_to.close()
                
                                    # get the "downloads" page     
                                    downloads_of_the_analysis_url = "https://www.ebi.ac.uk/metagenomics/api/v1/analyses/" + mgnify_analysis_id + "/downloads?page_size=100"
                                    get_downloads_of_the_analysis = load_url(downloads_of_the_analysis_url, 30)
                                    # for each entry in the download page of this analysis, find the entries that have to do with taxonomic and functional annotations
                                    if get_downloads_of_the_analysis.status_code == 200:
                                        downloads_of_the_analysis_page = load_json_file(get_downloads_of_the_analysis)
                                        for download in downloads_of_the_analysis_page['data']:
                                            description = download['attributes']['description']['description']
                                            
                                            if description in annotations_list:
                                                annotation_url = download['links']['self']
                                                annotation_url_page = requests.get(annotation_url, allow_redirects = True)
                                                # create a directory for each MGYA id, where all the corresponding annotation files will be saved to
                                                directory_to_save = path_to_save_annotation_files + mgnify_analysis_id
                                                check_for_directory = os.path.isdir(directory_to_save)
                                                x = False; y = (bool(x))
                                                if check_for_directory == y:
                                                    os.mkdir(directory_to_save)
                                                # get the suffix of the file you are about to save
                                                suffix_1 = annotation_url.split(".")[-1]
                                                suffix = suffix_1.split("?")[0]
                                                filename = directory_to_save + "/" + description + "."  + suffix
                                                open(filename, 'wb').write(annotation_url_page.content)
                                                sleepy = random.uniform(0.5,1.5)
                                                time.sleep(sleepy)
                                        # when you have found the sample under study in the "study_analyses" pages, break the loop and move on to the next sample!        
                                        if study_id == srs_id:
                                            print("i found the sample in the study and i downloaded all its relative annotations!")    
                                            break
                                        
#####################################################################################################################
#### Part C
#### Functions for getting the meta-transcriptomic data
#####################################################################################################################

## Get the number of pages with samples on MGnify that contain metatranscriptome analyses

def first():
    url="https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/metatranscriptomic/samples?page=1&page_size=100"
    samples = requests.get(url = url, allow_redirects = True)
    time_for_sleep = random.uniform(0.5,1.5)
    time.sleep(time_for_sleep)

    if samples.status_code == 200:
            data = samples.json()
            number_of_pages = int(data['meta']['pagination']['pages'])
            print("I got the number of sample pages. It is equal to: " + str(number_of_pages))
            return(number_of_pages)
    else:
            print("I am in the twilight zone!")
            first()
            
## Here is how to download a specific page
def download_samples(counter):
    url= "https://www.ebi.ac.uk/metagenomics/api/v1/experiment-types/metatranscriptomic/samples?page=" + str(counter) + "&page_size=100"
    page = requests.get(url, allow_redirects = True)
    time_for_sleep = random.uniform(0.5,1.5)
    time.sleep(time_for_sleep)
    pwd = os.getcwd()
    directory_to_save = "/data/databases/mgnify/metatranscriptomic_data/samples/"

    if page.status_code == 200:
        filename = directory_to_save + str(counter) + "_sample.json"
        filepath = os.path.join(pwd, directory_to_save, filename)
        file_to_save = open(filepath, 'wb')
        file_to_save.write(page.content)
        file_to_save.close()
    else:
            download_samples(counter)

# here is the function that runs the previous for all we need
def metatranscriptome_samples():
    number_of_pages = first()
    for page in range(1, number_of_pages+1):
        print("number of page for samples is: " + str(page))
        page_number = str(page)
        download_samples(page_number)
        
        
# a function to download files and keep taxonomies for the case of metatranscriptome data
def get_download_files(download_page):
    print("i am in the get_downlaod_files function! ")
    attributes = ["Complete GO annotations",
                  "Taxonomic assignments (TSV)",
                  "Complete GO annotation",
                  "InterPro matches (TSV)",
                  "Taxonomic assignments SSU (TSV)",
                  "Taxonomic assignments LSU (TSV)"]
        
    with open(download_page) as f:
        data = json.load(f)
        
    lsu_ssu = []
    # Explore each of the elements that the download page contains
    for entry in data['data']:
        potential_download_file = entry['attributes']['description']['description']
            
        if potential_download_file in attributes:                        
            link = entry['links']['self']
            download_file = requests.get(link, allow_redirects = True)
            # Check if downloading was fine and then save all files needed        
            if download_file.status_code == 200:
                directory_to_save = download_page.split('/')[:-1]
                directory_to_save = '/'.join(directory_to_save) 

                file = link.split('/')[-1]
                filename = directory_to_save + "/" + file
                pwd = os.getcwd()
                filepath = os.path.join(pwd, directory_to_save, filename)
                
                file_to_save = open(filepath, 'wb')
                file_to_save.write(download_file.content)
                file_to_save.close()
                
                if "_SSU_" in filename:
                    lsu_ssu.append(filename)
                    
                if "_LSU_" in filename:
                    lsu_ssu.append(filename)
    # Check if we have two taxonomies for this study and if yes then create a new file with both the taxonomies 
    if len(lsu_ssu) == 2:        
        new_filename = directory_to_save + "/" + "merged_taxa.tsv"
        filepath2 = os.path.join(pwd, directory_to_save, new_filename)
                
        with open(filepath2, 'w') as outfile:
            for taxonomy in lsu_ssu:
                with open(taxonomy) as infile:
                    outfile.write(infile.read())
    # now we can go to the next study
    print("now it's time for a new study! \n")

# in order to download the files wanted, we will build a function to do that
def download_downloads(study_id, url_for_this_study):
    time_for_sleep = random.uniform(0.5,1.5)
    page = requests.get(url_for_this_study, allow_redirects = True)
    time.sleep(time_for_sleep)
    pwd = os.getcwd()
    directory_to_save = "/data/databases/mgnify/metatranscriptomic_data/downloads/" + study_id + "/"

    if page.status_code == 200:
        
        filename = directory_to_save + str(study_id) + "_downloads.json"
        filepath = os.path.join(pwd, directory_to_save, filename)
        
        file_to_save = open(filepath, 'wb')
        file_to_save.write(page.content)
        file_to_save.close()

    else:
            download_samples(study_id, url_for_this_study)
