#!/usr/bin/env python3.5

########################################################################################
# script name: functions_jgi.py
# path on oxygen: /data/databases/scripts/gathering_data/jgi
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL
# all the functions needed to get the JGI/IMG related data and metadata and extract
# all the PREGO associations types needed
########################################################################################

# to get the metadata
import time, datetime, sys
import unittest, re
from pyvirtualdisplay import Display
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import TimeoutException
from unidecode import unidecode
import mmap

# to extraxt the associations
import tagger 

# --------------------------------------------------------------------------------------------------

# Overall variables 

# Path to JGI retrieved metadata
jgi_path = "/data/databases/jgi/"

# path to look for the dictionary
dictionary_prefix = '/data/dictionary/tagger'

# import a class called "Tagger()" from the "tagger.tagger" library
prego_tagger = tagger.Tagger()

##########################################################################################################################
# Get the number of pages we need to parse
##########################################################################################################################

# Get the first page of the results in order to keep the total number of the pages
# that is because we cannot use the "all" option of the jgi/img site
def the_jgi_genome_id_function(base_url, jgi_genome_ids_file):

    # Open a display browser that.. never actually opens!
    display = Display(visible=0, size=(800, 600))
    display.start()

    # The chosen broser is the Mozilla Firefox
    driver = webdriver.Firefox()
    page = str("0")
    first_url = base_url + page + "%26sort%3DDomain%26dir%3Dasc%26c%3DDomain%26f%3D%26t%3Dtext"

    driver.get(first_url)
    try:
        check = WebDriverWait(driver, 60).until(
            EC.presence_of_element_located((By.ID, "yui-dt0-th-IMGGenomeID-liner"))
        )
        time.sleep(5)
        # With this id item i get the line that also says the total number of pages
        getting_number_of_pages = driver.find_element_by_id("yui-pg0-0-page-report").text
        print(getting_number_of_pages)
        # I keep just the number
        number_of_pages = getting_number_of_pages.split(" ")[-1]
        number_of_pages = int(number_of_pages)                    
    except NoSuchElementException as e:
        print(e)
    finally:
        driver.quit()
    
    # I make a while loop to get all the pages. But first we open the file where the IDs will be saved
    with_open_jgi_genome_ids_file = open(jgi_genome_ids_file, "w")
    
    # As i have to increase the page number but also to create a new url in each loop, i have 2 counters; one for the page number (page_counter) and one for the url (counter)
    counter = 0
    page_counter = 0
    
    while page_counter < number_of_pages:
    
        driver = webdriver.Firefox()
        # the page variable needs to be a string!
        page = str(counter)
        new_url = base_url + page + "%26sort%3DDomain%26dir%3Dasc%26c%3DDomain%26f%3D%26t%3Dtext"
        print(page, new_url)
    
        # I get the new url each time
        driver.get(new_url)
    
        try:
            print("i m in the try step")
            check = WebDriverWait(driver, 60).until(
                EC.presence_of_element_located((By.ID, "yui-dt0-th-IMGGenomeID-liner"))
            )
            print("just done check.")
            time.sleep(5)
    
            td_elements = driver.find_elements_by_tag_name("td")
            liner_elements = driver.find_elements_by_class_name("yui-dt-liner")
            
            
            ####  make an if statement for the case of SAGs!
            jgi_genome_ids_file_suffix = jgi_genome_ids_file.split("/")[-1]
            if jgi_genome_ids_file_suffix == "jgi_genome_ids_bacteria_sags.txt" or jgi_genome_ids_file_suffix == "jgi_genome_ids_archaea_sags.txt":
                
                for element in range(6,len(liner_elements),9):
                    if element != 6:
                        element += 1
                        
                    check_if_screened = liner_elements[element-2].text
                    print(check_if_screened)
                    if re.search('contamination screened', check_if_screened):
                        print("Found a match!")
                        print(liner_elements[element].text)
                        with_open_jgi_genome_ids_file.write(liner_elements[element].text)
                        with_open_jgi_genome_ids_file.write("\n")
            
            # here is the same approach for Isolates and MAGs             
            else:
                try:
                    print("length is:" + str(len(liner_elements)))
        
                    for element in range(6,len(liner_elements),9):
                        if element != 6:
                            element += 1
        
                        print(liner_elements[element].text)
                        with_open_jgi_genome_ids_file.write(liner_elements[element].text)
                        with_open_jgi_genome_ids_file.write("\n")
    
                except NoSuchElementException as e:
                    print(e)
        finally:
            driver.quit()
    
        # Increase the page_counter by 1 and the url counte by 100
        page_counter += 1
        counter += 100
        page = ''
    
    display.stop()
    with_open_jgi_genome_ids_file.close()
    
    # Now i need to remove all the first lines of the tables that were parsed.  
    with open(jgi_genome_ids_file, "r") as g:
        lines = g.readlines()
    with open(jgi_genome_ids_file, "w") as f:
        for line in lines:
            if line.strip("\n") != "IMG Genome ID":
                f.write(line)
                
    # finally, i close the .txt file with the IDs
    f.close()

##########################################################################################################################
# The following 3 functions are used to download the metadata related         
##########################################################################################################################

## Part A: With respect to the project metadata; build a function to get the metadata for a certain genome id 
def get_metadata_for_genome_id(this_jgi_genome_id, check_for_timeout):
    
    # Open a display browser that.. never actually opens!
    display = Display(visible=0, size=(800, 600))
    display.start()
    
    base_url_for_metadata = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid="
    url = base_url_for_metadata + str(this_jgi_genome_id)
    print(url)
    
    # The chosen browser is the Mozilla Firefox.
    # If I do not set a timeout valued when setting the driver, then 30seconds is the by default corresponding values;
    # that is why I got a Timeout exception when I had 50 seconds to wait in the "until" step
    driver = webdriver.Firefox()
    driver.get(url)
        
    try:
        check = WebDriverWait(driver, 25).until(
            # EC.presence_of_element_located((By.ID, "content_other"))
            EC.presence_of_element_located((By.XPATH, '//*[@id="content_other"]/form/table'))
        )
        
        if check_for_timeout == False :
            time_for_sleep = 1
        else:
            time_for_sleep = 10

        time.sleep(time_for_sleep)
    
        get_categories = driver.find_elements_by_tag_name("th")
        categories = [ get_categories[entry].text for entry in range(len(get_categories)) ]

        get_values = driver.find_elements_by_tag_name("td")
        values = [ get_values[entry].text for entry in range(len(get_values)) ]
        
        categories_to_keep = ['NCBI Taxon ID', 'Ecosystem', 'Ecosystem Category', 'Ecosystem Subtype', 'Ecosystem Type', \
            'Geographic Location', 'Habitat', 'Latitude', 'Longitude', 'Specific Ecosystem', \
            'GOLD Ecosystem	Host-associated', 'GOLD Ecosystem Category', 'GOLD Ecosystem Subtype', 'GOLD Ecosystem Type', 'GOLD Specific Ecosystem', 'Host Name']
        
        # metadata_counter = 0
        metadata_found = ''
        this_ncbi_id = ''
        
        # for category in categories:
        for item in categories_to_keep:
            
            check = False
            
            # if category in categories_to_keep:
            for category in categories:
                
                if category == item :
                
                    category_index = categories.index(category)
                    # metadata_counter += 1
                    
                    if category != "NCBI Taxon ID":                
                        value = values[category_index-2]
                        metadata_found += value + "\t"
                        check = True
                        break
                    
                    else:
                        value = values[category_index]
                        metadata_found += value + "\t"
                        this_ncbi_id = value
                        check = True
                        break
                
            if check == False:
                metadata_found = metadata_found + "\t"
                 
        print("We got the metadata of the genome: " + str(this_jgi_genome_id))
        first_function_output = [metadata_found, this_ncbi_id]
        return first_function_output
    
    except NoSuchElementException as e:
        print(e)
        
    except TimeoutException:
        print("Timeout excpetion occured!")
        check_for_timeout = True
        first_function_output = get_metadata_for_genome_id(this_jgi_genome_id, check_for_timeout)
        return first_function_output
    
    finally:
        driver.quit()
        display.stop()
    

## Part B: With respect to the KO terms
def get_KO_terms(that_jgi_genome_id, that_ncbi_id, repeat = 1):
    
    up_limit = 10
    
    # Open a display browser that.. never actually opens!
    display = Display(visible=0, size=(800, 600))
    display.start()
    
    prefix_url_for_KO_terms = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=ko&taxon_oid="
    suffix_url_for_KO_terms = "#ko=results%startIndex%3D0%26sort%3DKOID%26dir%3Dasc%26c%3DKOID%26f%3D%26t%3Dtext"
    url = prefix_url_for_KO_terms + str(that_jgi_genome_id) + suffix_url_for_KO_terms

    # The chosen broser is the Mozilla Firefox 
    driver = webdriver.Firefox()
    driver.get(url)
    try:
        check = WebDriverWait(driver, 25).until(
            EC.presence_of_element_located((By.ID, "content_other"))
            # EC.presence_of_element_located((By.XPATH, '//*[@id="content_other"]/form/table'))
        )
        
        time.sleep(1)
        
        kegg_entries = driver.find_elements_by_tag_name("td")
        
        ko_terms = []
        for entry in range(2,len(kegg_entries),5):
            ko = kegg_entries[entry].text[3:]
            ko_terms.append(ko)
        
        all_ko_terms_in_a_string = ''
        for ko in ko_terms:
            all_ko_terms_in_a_string += ko + ";"
                
        return all_ko_terms_in_a_string
                
    except NoSuchElementException as e:
        print(e)
        no_ko_terms_retrieved = 'Null'
        return no_ko_terms_retrieved
    except TimeoutException:
        print("I got a Time Excpetion erron on getting the KO terms. I am now trying again to get them.")
        repeat += 1
        if repeat < up_limit:
            all_ko_terms_in_a_string = get_KO_terms(that_jgi_genome_id, that_ncbi_id, repeat = repeat)
            return all_ko_terms_in_a_string
        else:
            all_ko_terms_in_a_string = ""
            return all_ko_terms_in_a_string
    finally:
        driver.quit()
        display.stop()

## Part C: Perform what you built in the previous 2 steps; build a function to run them to all the ids found for all JGI categories
def get_metadata_and_KO_for_each_file(file_to_read, file_to_write):
    open_file = open(file_to_read, "r")
    open_file_to_write = open(file_to_write,"a+")
    ko_terms_per_jgi_genome_tidy_file = jgi_path + "ko_terms_per_jgi_genome_tidy_format.tsv"
    ko_terms_per_jgi_genome_tidy_file = open(ko_terms_per_jgi_genome_tidy_file, "a")

    for jgi_id in open_file:
        jgi_id = str(jgi_id.rstrip("\n"))
        print("The jgi_id under study is: " + jgi_id)
        
        ####   check if jgi_id on metadata retireved and if true break
        with open(file_to_write, 'rb', 0) as file, mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
            if s.find(jgi_id.encode('ascii')) != -1:
                print("The JGI Id " + jgi_id + " has already been retrieved from previous runs of this script. \n We move to the next entry. \n")
                continue
        s.close()
        
        # --->>> run the first function to get the metadata and the ncbi id regarding this specific JGI genome Id
        check_for_time = False
        metadata_function = get_metadata_for_genome_id(jgi_id, check_for_time)
        metadata = metadata_function[0]
        ncbi_id = metadata_function[1]
        print("First function completed!")
        
        # --->>> run the function to get the KO terms
        ko_terms = get_KO_terms(jgi_id, ncbi_id)
        print("Second function has been completed!")
        
        if ko_terms != 'Null':
            
            # split the one-line KO terms in a list of those    
            list_with_KO_terms = ko_terms.split(';')
            
            # read the file with all the KO terms regarding metabolism and keep those on a list 
            metabolism_KO_terms_file = open("/data/databases/kegg_orthology/only_metabolism_KO_terms.txt", "r")
            metabolism_KO_terms = []
            for term in metabolism_KO_terms_file:
                term = term.rstrip("\n")
                metabolism_KO_terms.append(term)
            
            metabolism_KO_terms_on_genome = ''
            # check which of the KO terms of the genome are related to metabolism and keep only those
            for term in list_with_KO_terms:
                if term in metabolism_KO_terms:
                    metabolism_KO_terms_on_genome += term + ";"
                    ko_terms_per_jgi_genome_tidy_file.write(ncbi_id + "\t" + term + "\n")
            
            
            metadata_and_KO_terms_in_a_line = metadata + str(metabolism_KO_terms_on_genome[:-1])
            metadata_and_KO_terms_in_a_line = str(jgi_id) + "\t" + metadata_and_KO_terms_in_a_line
            
            # In terms of always being writable
            metadata_and_KO_terms_in_a_utf_line = metadata_and_KO_terms_in_a_line.encode('utf-8')
            metadata_and_KO_terms_in_a_utf_str_line = str(metadata_and_KO_terms_in_a_utf_line, 'utf-8') 
            metadata_and_KO_terms_in_a_utf_str_line = unidecode(metadata_and_KO_terms_in_a_utf_str_line) 

            # write in the corresponding metadata file all this info returned
            open_file_to_write.write(metadata_and_KO_terms_in_a_utf_str_line + "\n")
            
        print("\n")
        
    open_file.close()
    open_file_to_write.close()
    ko_terms_per_jgi_genome_tidy_file.close()
  

##########################################################################################################################
# Extract-the-associations step
##########################################################################################################################

# The isolates case

def single_and_multiple(file):
    
    ids_found = []
    single_entries = []
    multiple_entries = {}

    with open(file, "r") as tmp:
        
        for line in tmp:
            
            elements = line.split("\t")
            
            ncbi_id = elements[1]
            
            if ncbi_id not in ids_found:
                ids_found.append(ncbi_id)
            else:
                multiple_entries[ncbi_id] = []

    with open(file, "r") as tmp:
        for line in tmp:
    
            elements = line.split("\t")           
            ncbi_id = elements[1]
            singularity_check = True
            
            for key_id in multiple_entries.keys():
                
                if key_id == ncbi_id:
                    
                    multiple_entries[ncbi_id] += [[line]]
                    singularity_check = False
                    break
                
            if singularity_check:
                
                single_entries.append([line])
    tmp.close()
    
    return single_entries, multiple_entries
            

def tag_the_singles(list_of_singles, multi = False):
    
    # import a specific function called "load_names()" from the "Tagger()" class
    prego_tagger.load_names('%s_entities.tsv' % dictionary_prefix, '%s_names.tsv' % dictionary_prefix)
    
    # search for matches
    
    matches_overall = set()

    for genome in list_of_singles:

        envo_matches = set()
        gold_envo_matches = set()
        host_matches = set()  
        elements = genome[0].split("\t")

        for index, element in zip(range(len(elements)), elements):
            
            if index == 0:
                jgi_id = element
            
            if index == 1:                
                ncbi_id = element
                
            if index >=2 and index <= 10:
                
                envo_match = prego_tagger.get_matches(element, 0, [-27], auto_detect=False, protect_tags=False, max_tokens=10) 

                for match in envo_match:
                    for entity in match[2]:
                        envo_matches.add(entity)
                        
            if index >=11 and index <= 15:

                gold_envo_match = prego_tagger.get_matches(element, 0, [-27], auto_detect=False, protect_tags=False, max_tokens=10)
                
                for match in gold_envo_match:
                    for entity in match[2]:
                        gold_envo_matches.add(entity)
            
            if index == 16:
                org_match = prego_tagger.get_matches(element, 0, [-2], auto_detect=False, protect_tags=False, max_tokens=10)

                for match in org_match:
                    for entity in match[2]:
                        host_matches.add(entity)
        
        if multi == False: 

            new_entry = (ncbi_id, jgi_id, tuple(envo_matches), tuple(gold_envo_matches), tuple(host_matches))            
            matches_overall.add(new_entry)
              
        else:
            new_entry = (tuple(envo_matches), jgi_id, tuple(gold_envo_matches), tuple(host_matches))
            matches_overall.add(new_entry)
    
    return matches_overall    


def tag_the_multiples(dictionary_of_multiples):
    
    dictionary_for_tags = {}
    
    for org, entries in dictionary_of_multiples.items():
        
        single_matches = tag_the_singles(entries, multi = True)
        
        if org not in dictionary_for_tags.keys():
            dictionary_for_tags[org] = (single_matches)
        else:
            dictionary_for_tags[org].add(single_matches)

    return dictionary_for_tags
    

def get_kos_of_each_id(file):
    
    dictionary_with_kos = {}
    
    with open(file, "r") as tmp:
        
        for line in tmp:
            
            elements = line.split("\t")
            jgi_id = elements[0]
            kos = elements[-1][:-1]
            dictionary_with_kos[jgi_id] = kos
            
    return dictionary_with_kos
            

def get_core_ko_terms(dictionary_with_kos, threshold):
    

    # Build a dictionary with sets of KOs that the elements of each is synonym one another
    synonyms_file = '/data/databases/kegg_orthology/ko_synonym_terms.tsv'
    synonyms = []
    synonyms_sets = {}

    with open(synonyms_file, "r") as synonyms_f:
        counter = 0
        for line in synonyms_f:
            kos = line.split("\t")[0]
            ko_terms = kos.split(";")
            
            for ko in ko_terms:
                synonyms.append(ko)

                if counter not in synonyms_sets: 
                    synonyms_sets[counter] = [ko]
                else:
                    synonyms_sets[counter] += [ko]
                    
            counter += 1

    synonyms_f.close() 
    
    # Count the hits of each KO; when a KO that has synonyms occurs, count it for all     
    number_of_taxa = len(dictionary_with_kos)
    cutoff = int(threshold * number_of_taxa)
    count_ko_hits = {}
    
    for jgi_id, kos in dictionary_with_kos.items():
        
        if len(kos) != 0 :
            
            kos_list = kos.split(";")
            
            for ko in kos_list:
                
                if ko in synonyms:
                    
                    for set_kos in synonyms_sets.values():
                        
                        if ko in set_kos:
                            
                            for term in set_kos:
                                
                                if term not in count_ko_hits.keys():
                                    
                                    count_ko_hits[term] = 1
                                
                                else:
                                    count_ko_hits[term] += 1
                                
                            break
                        
                        else:
                            continue
                    
                else:
                    if ko not in count_ko_hits.keys():
                        
                        count_ko_hits[ko] = 1
                    
                    else:
                        count_ko_hits[ko] += 1
    
    core_kos = []
    non_core_kos = []
    
    for ko, count in count_ko_hits.items():
        
        if count > cutoff:
            core_kos.append(ko)
        else:
            non_core_kos.append(ko)

    return core_kos, non_core_kos
        
    
def get_the_associations(input_file, analysis_type, process_overall_threshold, output_file_knowledge):
   
    # Urls
    base_url_for_metadata = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid="
    prefix_url_for_KO_terms = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=ko&taxon_oid="
    suffix_url_for_KO_terms = "#ko=results%startIndex%3D0%26sort%3DKOID%26dir%3Dasc%26c%3DKOID%26f%3D%26t%3Dtext"
   
   
    # Keep separately the organisms with a single from those with multiple entries
    singles, multiple_ones = single_and_multiple(input_file)
    print("NCBI Ids were seperated to single and multiple entries.")
    print("The number of singles equals to: " + str(len(singles)))
    print("The number of multiples equals to: " + str(len(multiple_ones)))

    # Use the tagger to annotate envo terms and host organisms
    singles_matches = tag_the_singles(singles)
    print("The single match cases were tagged.")

    """
    FUTURE WORK: 
    The following function merges the cases where for the same NCBI ID
    we have multiple JGI occurences. 
    However, as implementd for the time being, it requires high RAM leading 
    to issues. An alternative is required. When available, an extra for loop
    for these multiple_matches, as the one that follows for single_mathces, 
    will be needed to. 
    """
    # multiples_matches = tag_the_multiples(multiple_ones)
    # print("The multiple matches cases were tagged.")

    # Get the KOs for each JGI genome Id in a dictionary
    kos_per_jgi_id = get_kos_of_each_id(input_file)
    print("KO terms per JGI Id have been retrieved.")
   
    # Find the KOs that are repeated in all taxa found
    core_kos, non_core_kos = get_core_ko_terms(kos_per_jgi_id, threshold = process_overall_threshold)
    print("The core KOs have been extracted!")
    
    # Write the associations from the singles in a .tsv file
    for entity in singles_matches:

        ncbi_id = entity[0]
        jgi_id = entity[1]

        # Keep the urls for this genome
        url = base_url_for_metadata + str(jgi_id)
        ko_url = prefix_url_for_KO_terms + str(jgi_id) + suffix_url_for_KO_terms   
      
        # In entity[2] we have the non-GOLD envo terms 
        if len(entity[2]) != 0:
      
            for term in entity[2]:
      
                with open(output_file_knowledge, "a") as temp:
                    temp.write('-2\t' + ncbi_id  + '\t-27\t' + term[1] + '\t' + 'JGI IMG' + '\t' + analysis_type + '\t' + str(3) + '\tTRUE\t' + url + '\n')
                    temp.write('-27\t' + term[1] + '\t-2\t' + ncbi_id  +'\t' + 'JGI IMG' + '\t' + analysis_type + '\t' + str(3) + '\tTRUE\t' + url + '\n')

        # In entity[3] we have the non-GOLD envo terms 
        if len(entity[3]) != 0:
         
            for term in entity[3]:
      
                with open(output_file_knowledge, "a") as temp:
                   temp.write('-2\t' + ncbi_id  + '\t-27\t' + term[1] + '\t' + 'JGI IMG' + '\t' + analysis_type + ' GOLD' + '\t' +  str(4) + '\tTRUE\t' + url + '\n')
                   temp.write('-27\t' + term[1] + '\t-2\t' + ncbi_id  + '\t' + 'JGI IMG' + '\t' + analysis_type + ' GOLD' + '\t' +  str(4) + '\tTRUE\t' + url + '\n')
           
        # Now we need to get the KOs of this genome; if there are any
        kos = kos_per_jgi_id[jgi_id]
      
        if len(kos) != 0 :
         
            kos_list = kos.split(";")

            with open(output_file_knowledge, "a") as temp:
   
                for ko in kos_list:
                    """
                    If unmute the following if, you can get only the org-proc that are not
                    in the core set of processes.
                    """
                    #if ko not in core_kos:
                    temp.write('-2\t' + ncbi_id + '\t-90\t' + 'KO:' + ko + '\t' + 'JGI IMG' + '\t' + analysis_type + '\t' + str(5) + '\tTRUE\t' + ko_url + '\n') 
                    temp.write('-90\t' + 'KO:' + ko + '\t-2\t' + ncbi_id + '\t' + 'JGI IMG' + '\t' + analysis_type + '\t' + str(5) + '\tTRUE\t' + ko_url + '\n') 
        
    return core_kos, len(singles)

