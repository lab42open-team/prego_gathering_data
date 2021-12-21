#!/usr/bin/python3.5

########################################################################################
# script name: functions_mgrast.py
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL:
# develop all the functions related to the MG-RAST resource
########################################################################################
#
# usage: add the following line in the scripts you want to use these functions
# from functions_mgrast import *
#
########################################################################################

from __future__ import print_function
import os, sys, time, json
import string, math, random
import subprocess
import requests
from urllib.request import urlopen, Request
from urllib.error import HTTPError
import tagger

###########                 PART A: Get (meta)data

## 1. Perform the queries

# Return response body from MG-RAST or Shock API
def body_from_url(url, accept, auth=None, data=None, debug=False, method=None):
    header = {'Accept': accept}
    scriptname = os.path.basename(sys.argv[0])
    header['User-Agent'] = 'mglib:' + scriptname
    if auth:
        header['Authorization'] = 'mgrast '+auth
    if data or method:
        header['Content-Type'] = 'application/json'
    if debug:
        if data:
            print("data:\t"+repr(data))
        print("header:\t"+json.dumps(header))
        print("url:\t"+url)
    try:
        print("Making request "+url, file=sys.stderr)
        req = Request(url, data, headers=header)
        if method:
            req.get_method = lambda: method
        res = urlopen(req)
    except HTTPError as error:
        if debug:
            sys.stderr.write("URL: %s\n" %url)
        try:
            eobj = json.loads(error.read().decode("utf8"))
            if 'ERROR' in eobj:
                sys.stderr.write("ERROR (%s): %s\n" %(error.code, eobj['ERROR']))
            elif 'error' in eobj:
                sys.stderr.write("ERROR (%s): %s\n" %(error.code, eobj['error'][0]))
        except:
            sys.stderr.write("ERROR (%s): %s\n" %(error.code, error.read().decode("utf8")))
        finally:
            # sys.exit(1)
            res = ''
            return res
    if not res:
        sys.stderr.write("ERROR: no results returned\n")
        # sys.exit(1)
        res = ''
        return res
    return res

# return python struct from JSON output of MG-RAST or Shock API
def obj_from_url(url, auth=None, data=None, debug=False, method=None):
    try:
        result = body_from_url(url, 'application/json', auth=auth, data=data, debug=debug, method=method)
        if type(result) != str:
            read = result.read()
    except:  # try one more time;  ConnectionResetError 
        #result = body_from_url(url, 'application/json', auth=auth, data=data, debug=debug, method=method)
        #if type(result) != str:
         #   read = result.read()
        obj = " a URL error occured! move on with the next query."
        return obj
    if type(result) == str:
        obj = ''
        print("i got nothing ! ")
        return obj
    else: 
        if result.headers["content-type"] == "application/x-download" or result.headers["content-type"] == "application/octet-stream":
            return(read)   # Watch out!
        if result.headers["content-type"][0:9] == "text/html":  # json decoder won't work
            return(read)   # Watch out!
        if result.headers["content-type"] == "application/json":  # If header is set, this should work 
            data = read.decode("utf8")
            obj = json.loads(data)
        else:
            data = read.decode("utf8")
            obj = json.loads(data)
        if obj is None:
            sys.stderr.write("ERROR: return structure not valid json format\n" + repr(data))
            obj = ''
            return obj
            #sys.exit(1)
        if len(list(obj.keys())) == 0:
            if debug:
                sys.stderr.write("URL: %s\n" %url)
            sys.stderr.write("ERROR: no data available\n")
            obj = ''
            return obj
            #sys.exit(1)
        if 'ERROR' in obj:
            sys.stderr.write("ERROR: %s\n" %obj['ERROR'])
            obj = ''
            return obj
            #sys.exit(1)
        if ('error' in obj) and obj['error']:
            if isinstance(obj['error'], basestring):
                sys.stderr.write("ERROR:\n%s\n" %obj['error'])
            else:
                sys.stderr.write("ERROR: %s\n" %obj['error'][0])
            obj = ''
            return obj
            #sys.exit(1)
        return obj

## 2. Build the pool

# A function to handle the asynchronous requests for the case of taxonomies
def asynchronous(url, clock):

    time.sleep(random.uniform(10,20))
    get_new_page = obj_from_url(url, auth = None, data = None, debug = False, method = None)
    time_in_the_loop = time.time() - clock
    
    if time_in_the_loop > 600 or type(get_new_page) == str:
        
        if type(get_new_page) == str:
            get_new_page = ''
            return get_new_page
        
        elif get_new_page['status'] == 'done':
            return get_new_page
        
        else:
            get_new_page = ''
            return get_new_page
    
    if get_new_page['status'] == 'done':
        final_page = get_new_page
        return final_page
    else:
        asynchronous(url, clock)
        return get_new_page
# returning:
# get_new_pagete

# Function to get the metadata for each sample
def get_metadata(sample_url, tik_tok):
    times = time.time() - tik_tok
    metadata = obj_from_url(sample_url, auth = None, data = None, debug = False, method = None)
    time.sleep(random.uniform(10,30))
    if metadata == '' and times <= 180.0:
        get_metadata(sample_url, tik_tok)
        return metadata
    if times > 120.0:
        metadata = ''
        return metadata
    else:
        return metadata
# returning:
# metadata

# This is our main function. It handles each project along with all its samples of the MG-RAST db
def handle_a_project(project):

    # Output directory path
    main_output_directory_path = "/data/databases/mg_rast/"

    # The URL prefixes and suffices for getting all other needed info
    url_base_for_per_sample_metadata = "http://api.mg-rast.org/metadata/export/"
    url_base_for_KO_terms_pes_sample_prefix = "http://api.mg-rast.org/profile/"
    url_base_for_KO_terms_pes_sample_suffix = "?source=KO&format=mgrast"
    url_base_for_taxonomy_per_sample_prefix = "http://api.mg-rast.org/profile/"
    url_base_for_taxonomy_per_sample_suffix = "?source=SSU&format=lca"
    text_mining_file_name = "/data/databases/mg_rast/mgrast_metadata_for_text_mining.tsv"

    # Make a list with all the metabolism - related KO terms
    file_with_metabolism_KO_terms = "/data/databases/kegg_orthology/only_metabolism_KO_terms.txt"
    list_with_metabolism_KO_terms = []
    with open(file_with_metabolism_KO_terms, 'r') as ko_file:
        for term in ko_file:
            term = term.replace('\n', '')
            list_with_metabolism_KO_terms.append(term)

    time.sleep(random.uniform(0.5,3))

    # Keep the project id  and its desciption; the latter will be saved in a special file to contribute to the text-mining PREGO module
    project_id = project['id']
    project_description = project['description']
    project_description = project_description.replace('\n',' ')

    # Build a list with all the sample ids that this project includes
    list_with_sample_ids_of_this_project = []
    for sample_entry in project['metagenomes']:
        sample_id = sample_entry['metagenome_id']
        list_with_sample_ids_of_this_project.append(sample_id)
    
    # Build a single string with all the samples of the project separated by "|" pipe        
    samples_string = ''
    for sample in list_with_sample_ids_of_this_project:
        if sample == list_with_sample_ids_of_this_project[-1]:
            samples_string += sample
        else:
            samples_string += sample + "|"

    # Write the projects information (samples, description) in a file 
    with open (text_mining_file_name, "a") as temp:
        temp.write(project_id + "\t" + samples_string + "\t" + project_description + "\n")
        temp.close()

    # Now we iterate through all samples of the project under study    
    for sample in list_with_sample_ids_of_this_project:
        
        sample_initial_time = start_1 = time.time()

        # For each sample I need to keep all the KO terms and the species related in a list
        per_sample_KO_terms = ''
        per_sample_species_taxonomies = ''
    
        # Set a random time for sleep; here is a trick! we do the "sleep" step BEFORE the query. This way, we do not risk for our first queries to start at the exact same time and get a ban!                 
        time.sleep(random.uniform(1.0,7.5))
        
        # Get the metadata for this perticular metagenome (sample)        
        sample_metadata_url = url_base_for_per_sample_metadata + sample
        time_point = time.time()
        sample_metadata = get_metadata(sample_metadata_url, time_point)
        roundtrip_1 = time.time() - start_1
        time.sleep(random.uniform(5,10))

        # In case something goes really bad in the project step, we keep the related project info to re-run the whole script just for these projects. 
        if sample_metadata == '':
            print("Sample " + str(sample_id) + " have no metadata. We move on the next sample")

        # Keep the actual metadata; we know that mixs will be there as MG-RAST does not allow data to go public without those included
        else:
            
            checking_variable = 0

            # Also check whether the sample is from a host
            sample_host = ''
            if 'env_package' in sample_metadata:
                if 'data' in sample_metadata['env_package']:
                    if 'host_common_name' in sample_metadata['env_package']['data']:
                        sample_host = sample_metadata['env_package']['data']['host_common_name']['value']

            # Now keep the envo metadata; 'mixs' MUST be present, however, we cannot be sure in such approaches, so we.. if..
            if 'mixs' in sample_metadata:
                if 'sequence_type' in sample_metadata['mixs']:
                    sample_sequence_type = sample_metadata['mixs']['sequence_type']

                    # In case this sequence type occurs for the first time, create a directory and the corresponding output files
                    sample_sequence_type = sample_sequence_type.replace(' ', '_')
                    sample_sequence_type = sample_sequence_type.lower()
                    sequence_type_path = main_output_directory_path + sample_sequence_type + '/'
                    
                    # Here is where actually I keep the file where the sample's metadata will be kept; as variable
                    sequence_type_file_name = sequence_type_path + sample_sequence_type  + "_metadata_retrieved_TEST.tsv"    
                    
                    checking_variable = 1
                    
                    if os.path.isdir(sequence_type_path) == False:
                        os.mkdir(sequence_type_path)
                    if not os.path.exists(sequence_type_file_name):
                            os.mknod(sequence_type_file_name)
                    
                if 'biome' in sample_metadata['mixs']:
                    sample_biome = sample_metadata['mixs']['biome']
                if 'feature' in sample_metadata['mixs']:
                    sample_feature = sample_metadata['mixs']['feature']
                if 'material' in sample_metadata['mixs']:    
                    sample_material = sample_metadata['mixs']['material']

            # Get the KO terms regarding this metagenome.
            # ATTENTION! In amplicon studies, we do not wait for KO to have been annotated; however, we do not need an IF statement
            # That is simply because the 'sample_KO_page' variable will be empty.
            # We choose this way as we do not know the exact number and types of the sequence types MG-RAST includes.
            if checking_variable == 1:
                
                sample_KO_url = url_base_for_KO_terms_pes_sample_prefix + sample + url_base_for_KO_terms_pes_sample_suffix
        
                start_2 = time.time()
                sample_KO_page = obj_from_url(sample_KO_url, auth=None, data=None, debug=False, method=None)
                roundtrip_2 = time.time() - start_2
                time.sleep(random.uniform(5,10))
        
                # Create a set with the uniq KO terms found    
                if type(sample_KO_page) == str:
                    per_sample_KO_terms = ''

                else:
                    
                    if sample_KO_page['status'] == 'submitted':
                        new_url = sample_KO_page['url']
                        time_at_this_point = time.time()
                        sample_KO_page = asynchronous(new_url, time_at_this_point)

                    # Now that we have the KOs terms we have to parse this data file to get only those we need to keep
                    if type(sample_KO_page) != str:
                        
                        # This is a unique dictionary for each sample
                        dictionary_for_kegg = {}
                        
                        if sample_KO_page['status'] == 'done':
        
                            # metabolic_KO_terms_counter = 0
                            if sample_KO_page != '':
                                if 'data' in sample_KO_page:
                                    if 'data' in sample_KO_page['data']:
                
                                        for entry in sample_KO_page['data']['data']:
                                            kegg_term = entry[-1][0]
                    
                                            # Check if the KO term is on the list with the metabolism - related KO terms                        
                                            if kegg_term in list_with_metabolism_KO_terms:
        
                                                # Keep the abundance of the term as well
                                                kegg_term_abundance = entry[1]
                                                
                                                ##### REPLACED 2020.09.10 with isinstance check followin
                                                # if kegg_term_abundance != 'null':
                                                if isinstance(kegg_term_abundance, int) == True or isinstance(kegg_term_abundance, float) == True:
                                                
                                                    if kegg_term in dictionary_for_kegg.keys():
                                                        
                                                        try:
                                                            dictionary_for_kegg[kegg_term] += kegg_term_abundance 
                                                        except:
                                                            print("The KO abundance value under study was not an actual number and could not be added in \
                                                                  KO's abundance. Here is the value retrieved for the KO under study: ")
                                                            print(kegg_term_abundance)
                                                        
                                                    
                                                    else:
                                                        dictionary_for_kegg[kegg_term] = kegg_term_abundance

                                                
                                                else:
                                                    continue
        
                                    else:
                                        per_sample_KO_terms = ''
                                        
                                else:
                                    per_sample_KO_terms = ''
        
                        # Check if the dictionary supposed to contain the KOs actually has some
                        if bool(dictionary_for_kegg) == True:
                            
                            length_of_dic = len(dictionary_for_kegg)
                            counter_for_dic = 1
                            for term, value in dictionary_for_kegg.items():
                                kegg_pair = term + ":" + str(value)
                                if counter_for_dic == 1:
                                    per_sample_KO_terms = kegg_pair + "|"
                                elif counter_for_dic == length_of_dic:
                                    per_sample_KO_terms = per_sample_KO_terms + kegg_pair
                                else:
                                    per_sample_KO_terms = per_sample_KO_terms + kegg_pair + "|"
                                counter_for_dic += 1
                
                # Get the taxonomies regarding this metagenome
                
                # At first, build the corresponding url
                sample_taxonomy_url = url_base_for_taxonomy_per_sample_prefix + sample + url_base_for_taxonomy_per_sample_suffix
        
                # Then get it! 
                start_3 = time.time()
                sample_taxonomies = obj_from_url(sample_taxonomy_url, auth=None, data=None, debug=False, method=None)
                roundtrip_3 = time.time() - start_3
                time.sleep(random.uniform(5,10))
        
                # Check whether the variable with the response of the taxonomies' request is empty or not
                if sample_taxonomies != '':
                    
                    # This is a unique dictionary for each sample
                    dictionary_for_taxa = {}
        
                    # Often, what we get as a response, is another url where the actual response will be delivered
                    # If this is the case, we use the "asynchronous" function to open this latter url when it is ready
                    if 'status' in sample_taxonomies:
                        if sample_taxonomies['status'] == 'submitted':
                            new_url = sample_taxonomies['url']
                            time_at_this_point = time.time()
                            sample_taxonomies = asynchronous(new_url, time_at_this_point)
        
                    # Now that we have the info, we need to get what we need.
                    if 'status' in sample_taxonomies:
                        if 'data' in sample_taxonomies:
                            if sample_taxonomies['status'] == 'done' and 'data' in sample_taxonomies['data']: 
                                for entry in sample_taxonomies['data']['data']:
                                    
                                    # Get the name of the species
                                    species = entry[0]                          
                
                                    # Get rid of the unclassified taxonomies                
                                    bad_news = "unclassified"
                                    if bad_news not in species:

                                        # Due to the format of the taxonomy, taxa found not to species level, have '-;" spaces.
                                        # We do not want to keep taxonomies that are not to the species level                                        
                                        species = species.split(";")[-1]
                                        
                                        if len(species) > 3:
                                            
                                            try:
                                                species_relative_abundance = int(entry[1])
    
                                                if species in dictionary_for_taxa.keys():
                                                    dictionary_for_taxa[species] += species_relative_abundance
    
                                                else:
                                                    dictionary_for_taxa[species] = species_relative_abundance
                                                    
                                            except Exception:
                                                print("The relative abundance of this taxon entry is not a string, a bytes-like object or a number.\
                                                      Thus, we move to the next entry.")

                            else:
                                per_sample_species_taxonomies = ''
                
                        else:
                            per_sample_species_taxonomies = ''
        
                    # Make a string with all the species, one-next-to-the-other
                    if bool(dictionary_for_taxa) == True:
                        
                        length_of_dic_tax = len(dictionary_for_taxa)
                        counter_for_dic_tax = 1
                        for species, value in dictionary_for_taxa.items():
                            
                            species_relative_abundance_pair = species + ":" + str(value)
                            
                            if counter_for_dic_tax == 1:
                                per_sample_species_taxonomies = species_relative_abundance_pair + "|"
                            elif counter_for_dic_tax == length_of_dic_tax:
                                per_sample_species_taxonomies = per_sample_species_taxonomies + species_relative_abundance_pair
                            else:
                                per_sample_species_taxonomies = per_sample_species_taxonomies + species_relative_abundance_pair + "|"
                            
                            counter_for_dic_tax += 1    
                            
                            
                # Make a list with all the elements needed for the output file
                output_list = [sample, project_id, sample_biome, sample_feature, sample_material, sample_host, per_sample_species_taxonomies, per_sample_KO_terms]
        
                # Write them in the metadata output file
                with open(sequence_type_file_name, "a") as temp:
                    for element in output_list:
                        temp.write('%s\t' % element)
                    temp.write("\n")
                    temp.close()
        
                # Clear up the KO terms and species string variables
                per_sample_KO_terms = ''
                per_sample_species_taxonomies = ''
                # sequence_type_file_name = [sequence_type_file_name]
                
                sample_roundtrip = time.time() - sample_initial_time
                # print("Total time for the sample: " + str(sample_roundtrip) + " seconds")                                                                            # HERE ENDS THE FOR LOOP FOR THE SAMPLES
            
            else:
                return project

    # Finally, keep track of the projects retrieved 
    with open("/data/databases/scripts/gathering_data/mg_rast/projects_retrieved_from_mgrast.tsv", "a") as temp:
        temp.write(project_id + "\n")
        temp.close()
# returning:
# sample_id


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

###########                 PART B: Get the associations

###                         CASE 1: MARKER-GENE DATA

# Import a class called "Tagger()" from the "tagger.tagger" library.
dictionary_prefix = '/data/dictionary/tagger'
prego_tagger = tagger.Tagger()

# Import a specific function called "load_names()" from the "Tagger()" class; to load the names from your dictionary.
prego_tagger.load_names('%s_entities.tsv' % dictionary_prefix, '%s_names.tsv' % dictionary_prefix)

# A function to check whether an environmental metadatum (e.g biome) corresponds to an ENVO term 
def find_envo_term(lexicon_term, envo_term, dictionary_of_envo_terms):
    if lexicon_term in dictionary_of_envo_terms.keys():
        envo_term = dictionary_of_envo_terms[lexicon_term]
    else:
        envo_term = ''
    return envo_term

### Part B.1: Parse the data

# A function to get the species present in each sample along with their abundances and
# the corresponding ontology terms for the environmnetal and taxonomical metadata
def match_taxa_and_envs_with_ontology_terms(input_file):

    ## Build dictionaries
        
    # Dictionary for :  { lexicon taxon name : NCBI taxonomy id}
    dictionary_of_ncbi_taxonomy_ids = {} 
    
    # Dictionary for :  { lexicon envo name : ENVO term }
    dictionary_of_envo_terms = {}
    
    # Dictionary for :  { NCBI IDs : ABUNDANCES }
    dictionary_of_taxonomy_abundances_overall = {}

    # Keep species with no corresponding NCBI taxonomy id in a set
    species_with_no_corresponding_tax_id = set()

    with open(input_file, 'r') as file_with_metadata_retrieved:
        
        # iterate through the samples of amplicon studies gathered; each 'line' is a sample
        for line in file_with_metadata_retrieved:
            
            attributes = line.split("\t")
            sample_id = attributes[0]
            taxa_found = attributes[6].split("|")#[1:]
            
            dictionary_of_taxonomy_abundances_overall[sample_id] = {}
    
            ##  STEP_1:  with respect to TAXA --> NCBI ID conversion                    
            for taxon in taxa_found:
                if taxon != '':
                    
                    taxon = taxon.split(":")
                    taxon_name = taxon[0]
                    taxon_abundance = taxon[-1]                        # if the name of the taxon the character ":" might appear..... we are doomed.. :P
    
                    # TAGGER STEP :
                    # Run the tagger to get the NCBI ids 
                    matches = prego_tagger.get_matches(taxon_name, 0, [-2], auto_detect = False, protect_tags = False, max_tokens = 10)
                    
                    # If I have just one match, then that is the corresponding NCBI id I am looking for
                    if len(matches) == 1:                                                                            
                        ncbi_taxonomy_id = matches[0][2][0][1]
                         
                    # If there is none there I add the name of the taxon in a set.
                    elif len(matches) == 0:                                                                        
                        species_with_no_corresponding_tax_id.add(taxon_name)
    
                    # If there are more than one matches, we remove any stop-words e.g unclassified and we keep the rest
                    else:
                        # Build a list of the unique ids retrieved
                        list_of_ids = []
                        # Examine ncbi ids retrieved one by one
                        for match in matches:
                            current_ncbi_taxonomy_id = match[2][0][1]
                            # Check if "unclussified"
                            if current_ncbi_taxonomy_id == '12908':   
                                continue
                            # If not, add the ncbi id in the set of ids
                            else:  
                                if current_ncbi_taxonomy_id not in list_of_ids:
                                    list_of_ids.append(current_ncbi_taxonomy_id)
    
                        # Keep the first match
                        ncbi_taxonomy_id = list_of_ids[0]
    
                    # Add ```taxon : ncbi id``` pair in the corresponding dictionary if not already there
                    if taxon_name not in dictionary_of_ncbi_taxonomy_ids.keys(): 
                        dictionary_of_ncbi_taxonomy_ids[taxon_name] = ncbi_taxonomy_id
    
                    # Overall taxonomy dictionary
                    # Fill in the overall taxonomy abundance dictionary 
                    if ncbi_taxonomy_id in dictionary_of_taxonomy_abundances_overall[sample_id]: 
                        dictionary_of_taxonomy_abundances_overall[sample_id][ncbi_taxonomy_id] += taxon_abundance
                    else:
                        dictionary_of_taxonomy_abundances_overall[sample_id][ncbi_taxonomy_id] = taxon_abundance
    
                    # Clear the variable
                    ncbi_taxonomy_id = '' 
                            
                ##  STEP_2: with respect to environmental metadata --> ENVO terms
                environmental_metadata = attributes[2:5]
                for term in environmental_metadata:
                    envo_matches = prego_tagger.get_matches(term, 0, [-27], auto_detect = False, protect_tags = False, max_tokens = 10)
                    if len(envo_matches) > 0 :
                        envo_term = envo_matches[0][2][0][1]
                        if term not in dictionary_of_envo_terms.keys():
                            dictionary_of_envo_terms[term] = envo_term
        
        file_with_metadata_retrieved.close()
        
    return dictionary_of_ncbi_taxonomy_ids, dictionary_of_envo_terms
   # returning:
   # dictionary_of_ncbi_taxonomy_ids, dictionary_of_envo_terms

# Build a dictionary for the Taxa found in each environment
def taxa_envo_in_samples(input_file, dictionary_of_ncbi_taxonomy_ids, dictionary_of_envo_terms):
   
   # Likewise, read the metadata file again
    with open(input_file, 'r') as file_with_metadata_retrieved: 
       
        # Build a dictionary for ENVO terms : { NCBI_IDs : ABUNDANCES }
        envo_taxa = {}
       
        # Again, each line is a sample
        for line in file_with_metadata_retrieved:
   
            # We split and keep the same info as before
            attributes = line.split("\t")
           
            sample_id               = attributes[0]
            sample_project          = attributes[1]
            lexicon_sample_biome    = attributes[2]
            lexicon_sample_feature  = attributes[3]
            lexicon_sample_material = attributes[4]
           
            sample_biome            = find_envo_term(lexicon_sample_biome, 'sample_biome', dictionary_of_envo_terms)
            sample_feature          = find_envo_term(lexicon_sample_feature, 'sample_feature', dictionary_of_envo_terms)
            sample_material         = find_envo_term(lexicon_sample_material, 'sample_material', dictionary_of_envo_terms)
           
           
            # But this time, the taxa found in this sample are kept in a list
            sample_taxa = attributes[6].split("|")[1:]                                                        # ATTENTION!!! this will have to change once we replace ";" with '|"
   
            # In this case, we need to map the species name retrieved from MG-RAST to a NCBI taxonomy id
            # We build a set for the unique IDs found
            sample_ncbi_id_taxa = {}
           
            # And for each name we have
            for entry in sample_taxa:
                
                try:
                    entry = entry.split(":")
                    taxon_name = entry[0]
                    taxon_abundance = float(entry[-1])
                   
                    # We are looking for the corresponding ID in the dictionary we built in the first place
                    for key, value in dictionary_of_ncbi_taxonomy_ids.items(): 
                        if taxon_name == key:
                          
                            if value in sample_ncbi_id_taxa.keys():
                                sample_ncbi_id_taxa[value] += taxon_abundance
                            else:
                                sample_ncbi_id_taxa[value] = taxon_abundance            
                except:
                    print("Could not get a value for this taxon abundance. Move on to the next one.")
   
            # We need to be sure that we have entries that have only NCBI taxonomy IDs and not words
            # Now we can make a dictionary with the sample's info and the NCBI ids found                    
            envo_taxa[sample_id] = {}
            envo_taxa[sample_id]['project'] = sample_project
            envo_taxa[sample_id]['biome'] = sample_biome
            envo_taxa[sample_id]['feature'] = sample_feature
            envo_taxa[sample_id]['material'] = sample_material

           # For now PREGO doesn't need the concanated ENVO ids biome_feature 
           # and feature_material
           
           # envo_taxa[sample_id]['biome_feature'] = sample_biome_feature
           # envo_taxa[sample_id]['feature_material'] = sample_feature_material
            envo_taxa[sample_id]['taxa'] = sample_ncbi_id_taxa
   
            """
            # Example from what we have built up to now:
            This is a key - value pair from the envo_taxa dictionary. 
            The key is a MG-RAST sample and its value is a dictionary 
            including the ENVO related metadata for this sample, 
            along with the taxa found on it with their corresponding abundances
            
            mgm4507497.3 {'feature_material': 'ENVO:00000115_ENVO:00001998', 
                'taxa': {'986': 2.0, '106592': 1.0, '389082': 1.0, '391596': 52.0, \
            '378806': 1.0, '245': 1.0, '991': 1.0, '266940': 1.0, '395019': 3.0, '996': 226.0, '119050': 14.0, '520487': 6.0, '39946': 68.0, \
            '667019': 2.0, '387094': 1.0, '189425': 1.0, 997': 1.0, '189426': 1.0, '37752': 15.0}, 
                'feature': 'ENVO:00000115', \
                'biome_feature': 'ENVO:00000446_ENVO:00000115', 
                'project': 'mgp2603', 
                'biome': 'ENVO:00000446', 
                'material': 'ENVO:00001998'}
            """
        
        file_with_metadata_retrieved.close()
    
    return envo_taxa
    # returning:
    # envo_taxa (dictionary)

### Part B.2: Functions to count co-occurences and backgroudns. Get ENVO : ORGANISMS associations

def get_envo_background(envo_taxa):
    
    # Likewise, we calculate the corresponding BACKGROUND
    # --> we are doing this step again, however, we know that the two background dictionaries should be exactly the same

    # Build a dictionary to keep those counts
    background_envo = {}
    
    for sample in envo_taxa:
        
        # For each of the levels of environmental metadata 
        envo_types = ["biome", "feature", "material"]#"biome_feature", "feature_material"] 
        for envo in envo_types:
            
            # Count how many times it appears on our data
            entity = envo_taxa[sample][envo]
            
            if entity != '':
                if entity in background_envo:
                    background_envo[entity] += 1
                else:
                    background_envo[entity] = 1
    
    return background_envo
    # returning:
    # background_envo

def get_env_org_coocurrences(envo_taxa):
    
    env_taxon_coocurrences = {}
    
    # Parse again this time through the "envo_taxa" dictionary
    for sample in envo_taxa:
        envo_types = ["biome", "feature", "material"] #"biome_feature", "feature_material"
        # The envo_types supported are controlled in function taxa_envo_in_samples.
        # Iterate through the environmental types recorded again
        for envo_type in envo_types: 
            envo = envo_taxa[sample][envo_type]
            
            # Iterate through the taxa of the sample
            for entity in envo_taxa[sample]['taxa']:
                if envo in env_taxon_coocurrences.keys():           
                    if entity in env_taxon_coocurrences[envo]: 
                        env_taxon_coocurrences[envo][entity] += 1 
                    else: 
                        env_taxon_coocurrences[envo][entity] = 1
                else: 
                    env_taxon_coocurrences[envo] = {}
                    env_taxon_coocurrences[envo][entity] = 1
    return env_taxon_coocurrences
    # returning:
    # env_taxon_coocurrences
    
def get_taxon_background(envo_taxa):
    
    taxon_background = {}
    
    for sample in envo_taxa:
        
        try:
            
            for ncbi_id in envo_taxa[sample]['taxa']:
                
                if ncbi_id not in taxon_background.keys():
                    
                    taxon_background[ncbi_id] = 1
                
                else:
                    taxon_background[ncbi_id] += 1

        except Exception:
            print("This is not a good sample: " + str(sample))
    
    return taxon_background
    # returning:
    # taxon_background


# ///////\\\\\\\\\\\\\\\//////////\\\\\\\\\\\\\\\////////////////////\\\\\\\


###                         CASE 2: WGT DATA

def tag_envs_and_filter(attribute, tokens, dict):

    matches = prego_tagger.get_matches(attribute, 0, [-27], auto_detect = False, protect_tags = False, max_tokens = 10)
    if len(matches) != 0:
        envo_term = matches[0][2][0][1]
        if envo_term not in dict.keys():
            dict[attribute] = envo_term
        return envo_term, dict
    else:
        return False

def match_wgt_metadata_with_ontology_terms(input_file):

    with open(input_file, 'r') as file_with_metadata_retrieved:
        
        # Dictionary for :  { lexicon taxon name : NCBI taxonomy id}
        dictionary_of_ncbi_taxonomy_ids = {}        
        
        dictionary_of_taxonomy_abundances_overall = {}
        
        # Keep species with no corresponding NCBI taxonomy id in a set
        species_with_no_corresponding_tax_id = set()
        
        # Create dictionaries for all types of associations we are interested in 
        # envo_kegg_terms = {}
        envo_ncbi_ids = {}
        gomf_terms_ncbi_ids = {}
        envo_gomf_terms_ncbi_ids_per_sample = {}
        dictionary_of_envo_terms = {}
        
        # Keep in mind that each line of the file is a sample; iterate for each sample
        metagenome_counter = 0
        for line in file_with_metadata_retrieved:
    
            # As we have a tab seperated file, we keep the attributes of the file in a list
            attributes = line.split("\t") 
    
            # Keep the sample and the project ids
            sample_id = attributes[0] 
            sample_project = attributes[1]
            dictionary_of_taxonomy_abundances_overall[sample_id] = {}    
            
            # Keep the environmental metadata for the sample         
            lexicon_sample_biome = attributes[2]
            lexicon_sample_feature = attributes[3]
            lexicon_sample_material = attributes[4]

            # Make a list with the ENVO terms of the environmental metadata found in the sample
            try:
                sample_biome, dictionary_of_envo_terms    = tag_envs_and_filter(lexicon_sample_biome, 'sample_biome', dictionary_of_envo_terms)
                # print(lexicon_sample_biome, sample_biome)
            except:
                sample_biome = ''
                print("No sample for lexicon attribute: " + str(lexicon_sample_biome))
            try:
                sample_feature, dictionary_of_envo_terms  = tag_envs_and_filter(lexicon_sample_feature, 'sample_feature', dictionary_of_envo_terms)
                # print(lexicon_sample_feature, sample_feature)
            except:
                sample_feature =''
                print("No feature for lexicon attribute: " + str(lexicon_sample_feature))
            try:
                sample_material, dictionary_of_envo_terms = tag_envs_and_filter(lexicon_sample_material, 'sample_material', dictionary_of_envo_terms)    
                # print(lexicon_sample_material, sample_material)
            except:
                sample_material = ''
                print("No material for lexicon attribute: " + str(lexicon_sample_material))        
            
            # Now the gomfs retrieved are kept in a list; the first element is empty, i.e. ''   
            # Keep the taxa and the gomfs found as lists
            taxa_found = attributes[6].split("|")
            sample_gomf_terms = attributes[7].split("|")  #[1:]
                       
            # This is a dictionary for the gomf terms found in the sample under study 
            sample_gomf = {}
            for gomf in sample_gomf_terms:
                if gomf != '':
                    gomf = gomf.split(":")
                    try:
                        gomf_term = gomf[0] + ":" + gomf[1]
                        gomf_abundance = float(gomf[-1])
                        if gomf_term in sample_gomf.keys():
                            sample_gomf[gomf_term] += gomf_abundance
                        else:
                            sample_gomf[gomf_term] = gomf_abundance
                    except:
                        print("\n\n\n >>>>>> gomf_term: " + gomf_term + "\n\n") 

            # Get the NCBI ID for each taxon name retrieved 
            sample_ncbi_id_taxa = {}            
            for entry in taxa_found:
                
                try:
                    taxon = entry.split(":")
                    taxon_name = taxon[0]
                    taxon_abundance = float(taxon[-1])
                except:
                    print(entry)
                    continue

                try:
                    # TAGGER STEP :
                    # Run the tagger to get the NCBI ids 
                    matches = prego_tagger.get_matches(taxon_name, 0, [-2], auto_detect = False, protect_tags = False, max_tokens = 10)

                except:
                    print(entry)
                    print("Could not get a value for this taxon abundance. Move on to the next one.")
                    continue

                # If I have just one match, then that is the corresponding NCBI id I am looking for
                if len(matches) == 1:
                    current_ncbi_taxonomy_id = matches[0][2][0][1]
                    # Check if "unclussified"
                    if current_ncbi_taxonomy_id == '12908':   
                        continue
                    else:
                        ncbi_taxonomy_id = matches[0][2][0][1]
                else:
                    for match in matches:
                        current_ncbi_taxonomy_id = match[2][0][1]
                        if current_ncbi_taxonomy_id == '12908':   
                            continue
                        else:
                            ncbi_taxonomy_id = current_ncbi_taxonomy_id
                            break


                # Add ```taxon : ncbi id``` pair in the corresponding dictionary if not already there
                if taxon_name not in dictionary_of_ncbi_taxonomy_ids.keys(): 
                    dictionary_of_ncbi_taxonomy_ids[taxon_name] = ncbi_taxonomy_id

                # Overall taxonomy dictionary
                # Fill in the overall taxonomy abundance dictionary 
                if ncbi_taxonomy_id in dictionary_of_taxonomy_abundances_overall[sample_id]: 
                    dictionary_of_taxonomy_abundances_overall[sample_id][ncbi_taxonomy_id] += taxon_abundance
                else:
                    dictionary_of_taxonomy_abundances_overall[sample_id][ncbi_taxonomy_id] = taxon_abundance

                # Clear the variable
                ncbi_taxonomy_id = ''

                # We are looking for the corresponding ID in the dictionary we built in the first place
                for key, value in dictionary_of_ncbi_taxonomy_ids.items(): 
                    if taxon_name == key:
                        
                        if value in sample_ncbi_id_taxa.keys():
                            sample_ncbi_id_taxa[value] += taxon_abundance
                        else:
                            sample_ncbi_id_taxa[value] = taxon_abundance


            # Build up the dictionary entry for this sample 
            envo_gomf_terms_ncbi_ids_per_sample[sample_id] = {}
            envo_gomf_terms_ncbi_ids_per_sample[sample_id]['project'] = sample_project
            envo_gomf_terms_ncbi_ids_per_sample[sample_id]['biome'] = sample_biome
            envo_gomf_terms_ncbi_ids_per_sample[sample_id]['feature'] = sample_feature
            envo_gomf_terms_ncbi_ids_per_sample[sample_id]['material'] = sample_material            
            envo_gomf_terms_ncbi_ids_per_sample[sample_id]['gomfs'] = sample_gomf
            envo_gomf_terms_ncbi_ids_per_sample[sample_id]['taxa'] = sample_ncbi_id_taxa
            
            # print(envo_gomf_terms_ncbi_ids_per_sample[sample_id])
            metagenome_counter += 1
            
    file_with_metadata_retrieved.close()
    
    return envo_gomf_terms_ncbi_ids_per_sample
    # returning:
    # envo_gomf_terms_ncbi_ids_per_sample

def get_backgrounds_wgs_case(envo_gomfs_terms_ncbi_ids_per_sample):
    
    # Build three dictionaries for ENVO, NCBI_IDs and gomfs correspondingly
    
    envo_backgrounds = {}
    ncbi_ids_backgrounds = {}
    gomfs_terms_backgrounds = {}
    
    # Parse the dictionary returned from the match_wgt_metadata_with_ontology_terms() function
    for sample, metadata in envo_gomfs_terms_ncbi_ids_per_sample.items():
        
        inner_env_check_1 = ''
        inner_env_check_2 = ''

        # Parse the sample's dictionary
        counter = 1
        for entry, value in metadata.items():
            # Get the environments 
            if entry == 'biome' or entry == 'feature' or entry == 'material':
                if counter == 1:
                    if value not in envo_backgrounds.keys():
                        envo_backgrounds[value] = 1
                    else:
                        envo_backgrounds[value] += 1
                    inner_env_check_1 = value
                    counter += 1 

                if counter == 2:
                    if value == inner_env_check_1:
                        print("found same envo terms!")
                        print(inner_env_check_1, value)
                        counter += 1
                        continue
                    else:
                        if value not in envo_backgrounds.keys():
                            envo_backgrounds[value] = 1
                        else:
                            envo_backgrounds[value] += 1
                        inner_env_check_2 = value
                        counter += 1

                if counter ==3 :
                    if value == inner_env_check_1 or value == inner_env_check_2:
                        continue
                    else:
                        if value not in envo_backgrounds.keys():
                            envo_backgrounds[value] = 1
                        else:
                            envo_backgrounds[value] += 1

            # Now the organisms   
            elif entry == 'taxa': 
                for taxon, abundance in metadata['taxa'].items(): 
                    if abundance > 10: 
                        if taxon not in ncbi_ids_backgrounds.keys(): 
                            ncbi_ids_backgrounds[taxon] = 1 
                        else: 
                            ncbi_ids_backgrounds[taxon] += 1 
            
            # And the processes    
            elif entry == 'gomfs': 
                for gomf, abundance in metadata['gomfs'].items(): 
                    if abundance > 20: 
                        if gomf not in gomfs_terms_backgrounds.keys(): 
                            gomfs_terms_backgrounds[gomf] = 1 
                        else: 
                            gomfs_terms_backgrounds[gomf] += 1 

    return envo_backgrounds, ncbi_ids_backgrounds, gomfs_terms_backgrounds
    # returning:
    # envo_backgrounds, ncbi_ids_backgrounds, gomfs_terms_backgrounds

def count_cooccurences(envo_gomf_terms_ncbi_ids_per_sample):
    
    envo_ncbi_ids = {} ; envo_gomf = {} ; gomf_ncbi_ids = {}
    
    for sample, metadata in envo_gomf_terms_ncbi_ids_per_sample.items():
        
        biome = metadata['biome'] ; feature = metadata['feature'] ; material = metadata['material']        
        envo_terms = set([biome, feature, material])
        for term in envo_terms:
            if term not in envo_ncbi_ids.keys():
                envo_ncbi_ids[term] = {}            
            if term not in envo_gomf.keys():
                envo_gomf[term] = {}
            
        # Parse the taxa found in the sample        
        for taxon, abundance in metadata['taxa'].items():
            
            if abundance > 10:
                
                for envo in envo_terms:
                    
                    if taxon not in envo_ncbi_ids[envo]:
                        envo_ncbi_ids[envo][taxon] = 1
                    else:
                        envo_ncbi_ids[envo][taxon] += 1
                    
        # Parse the gomfs found in the sample
        for gomf, abundance in metadata['gomfs'].items(): 
            
            if abundance > 30:
                
                for envo in envo_terms:
                    
                    if gomf not in envo_gomf[envo].keys():
                        envo_gomf[envo][gomf] = 1
                    else:
                        envo_gomf[envo][gomf] += 1

        # Parse both taxa and gomfs of the sample to get their co-occurence
        for gomf, gomf_abundance in metadata['gomfs'].items():
            
            if gomf_abundance > 30:
                if gomf not in gomf_ncbi_ids.keys():
                    gomf_ncbi_ids[gomf] = {}
                    
            else:
                continue
            
            for taxon, taxon_abundance in metadata['taxa'].items():
                if taxon_abundance > 10:
                    if taxon not in gomf_ncbi_ids[gomf].keys():
                        gomf_ncbi_ids[gomf][taxon] = 1
                    else:
                        gomf_ncbi_ids[gomf][taxon] += 1
            
    return envo_ncbi_ids, envo_gomf, gomf_ncbi_ids
    # returning:
    # envo_ncbi_ids, envo_gomf, gomf_ncbi_ids -- the latter is like:
    # {ko_1: {ncbi_id_1: 34, ncbi_id_2: 2}, ko_2: {} ....} 

