#!/usr/bin/python3.5

########################################################################################
# script name: extract_mgrast_wgs_associations.py
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL:
# this script aims to parse the (meta)data retrieved to mine associations
########################################################################################
#
# USAGE: 
# ./extract_mgrast_wgs_associations.py  --input /data/databases/mg_rast/wgs/wgs_gomfs_test.tsv\  
#           --output test_with_thresholds.tsv
#
#   OR ( in case that any of the first 3 steps have been completed )
#
# ./extract_mgrast_wgs_associations.py --eb envo_backgrounds.json --ob ncbi_ids_backgrounds.json \
#           --kb gomfs_terms_backgrounds.json --oec org_env_co_occurences.json \
#           --pec proc_env_co_occurences.json --poc proc_org_co_occurences.json
########################################################################################

import json, os, sys, getopt
import tagger
import math
from functions_mgrast import *

# Main path - working directory
scripts_path = '/data/databases/scripts/gathering_data/mg_rast'
experiments_path = '/data/experiments'

# Path to look for the dictionary
dictionary_prefix = '/data/dictionary/prego'

# Path to look for the data files we intend to use
database_prefix = '/data/databases/mg_rast'
database_prefix_wgs = database_prefix + '/wgs'

# This is a test file for the develpment of the code  
# input_file = database_prefix_wgs + '/wgs_metadata_retrieved_TEST.tsv'     

# output_file = experiments_path + "/mgrast_wgs_associations.tsv"
# output_file = scripts_path + "/output.dev"



# Build a dictionary with all info per sample. This will look like this:
# {mgm4563763.3: {'biome': ENVO_XXXXX}, 'feature : ENVO_YYYY', 'material' : ENVO_ZZZZ,  'taxa': {'5484' : 23, ....}, 'gomfss' : {} }
if len(sys.argv) == 5: 

    opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["input=", "output="])

    for i in opts:
        if i[0] == "--input" or i[0] == "i":
            if "/" in i[1][:-1]: 
                input_file = i[1]
            else:
                input_file = experiments_path + i[1]

        if i[0] == "--output" or i[0] == "o": 
            if "/" in i[1]: 
                output_file = i[1] 
            else:
                output_file = experiments_path + "/" + i[1]

    print("All steps are about to perform.\n")

    envo_taxa_gomfs_per_sample = match_wgt_metadata_with_ontology_terms(input_file)
    print("Step 1/4 : function match_wgt_metadata_with_ontology_terms() has been completed.")
    with open('envo_taxa_gomfs_per_sample.json', 'w') as fp:
        json.dump(envo_taxa_gomfs_per_sample, fp)


    # Get the background for each ENVO term
    envo_backgrounds, ncbi_ids_backgrounds, gomfs_terms_backgrounds = get_backgrounds_wgs_case(envo_taxa_gomfs_per_sample)
    print("Step 2/4 : function get_backgrounds_wgs_case() has been completed.")
    with open('envo_backgrounds.json', 'w') as fp:
        json.dump(envo_backgrounds, fp)
    with open('ncbi_ids_backgrounds.json', 'w') as fp:
        json.dump(ncbi_ids_backgrounds, fp)
    with open('gomfs_terms_backgrounds.json', 'w') as fp:
        json.dump(gomfs_terms_backgrounds, fp)


    # Count co-occurences
    org_env_co_occurences, proc_env_co_occurences, proc_org_co_occurences = count_cooccurences(envo_taxa_gomfs_per_sample)
    print("Step 3/4 : function count_cooccurences() has been completed.")
    with open('org_env_co_occurences.json', 'w') as fp:
        json.dump(org_env_co_occurences, fp)
    with open('proc_env_co_occurences.json', 'w') as fp:
        json.dump(proc_env_co_occurences, fp)
    with open('proc_org_co_occurences.json', 'w') as fp:
        json.dump(proc_org_co_occurences, fp)


else: 

    opts, args = getopt.getopt(sys.argv[1:], "e:o:k:f:g:h:", ["eb=", "ob=", "kb=", "oec=", "pec=", "poc=", "output="])

    for i in opts:

        if i[0]   == "--eb":
            print("Loading envo backgrounds...")
            with open(i[1]) as json_file:
                envo_backgrounds       = json.load(json_file)
            print("Envo backgrounds have been loaded.")
        elif i[0] == "--ob":
            print("Loading organisms backround...")
            with open(i[1]) as json_file:
                ncbi_ids_backgrounds   = json.load(json_file)
            print("Organisms background have been loaded.")
        elif i[0] == "--kb":
            print("Loading processes background...")
            with open(i[1]) as json_file:
                gomfs_terms_backgrounds = json.load(json_file)
            print("Processes background have been loaded.")
        elif i[0] == "--oec":
            print("Loading organisms-environments co-occurences...")
            with open(i[1]) as json_file:
                org_env_co_occurences  = json.load(json_file)
            print("Organisms-environments co-occurences have been loaded.")
        elif i[0] == "--pec": 
            print("Loading processes-environments co-occurrences...")
            with open(i[1]) as json_file:
                proc_env_co_occurences = json.load(json_file)
            print("Processes-environments co-occurrences have been loaded.")
        elif i[0] == "--poc":
            print("Loading processes-organisms co-occurrences...")
            with open(i[1]) as json_file:
                proc_org_co_occurences = json.load(json_file)
            print("Organisms-processes co-occurrences have been loaded.")

        elif i[0] == "--output": 
            if "/" in i[1]: 
                output_file = i[1] 
            else:
                output_file = experiments_path + "/" + i[1]


# Extract ORG - ENV associations
# As "case_1" we consider the case where an organism is our entry point
# Likewise, "case_2" is when an environment is our entry point
print("Step 4/4 : initiation of output file creation")
with open(output_file, "w+") as temp:
    
    for environment, organisms in org_env_co_occurences.items():

        background_envo = envo_backgrounds[environment]
        
        for organism, count in organisms.items():
            
            background_org = ncbi_ids_backgrounds[organism]
            
            evidence_1 = '%d of %d samples' % (count, background_org) 
            evidence_2 = '%d of %d samples' % (count, background_envo)

            if count > background_envo or count > background_org:
                print("Something is going wrong in the countings on ORG - ENV associations")
                print(organism, environment, count, background_org, background_envo)
                continue

            score_1 = 2.0*math.sqrt(float(count)/background_envo**0.2)
            score_2 = 2.0*math.sqrt(float(count)/background_org**0.2)

            # Make sure both terms are not empty 
            if organism != '' and environment != '':

                # And print the final output! again as a tab seperated file. 
                # The format requires 9 fields:
                # type 1 , id 1, type 2 , id 2 , source, evidence, score, boolean, url.
                temp.write('-2\t' + organism  + '\t-27\t' + environment + '\tMG-RAST metagenome study\t' + evidence_1 + '\t' + str(score_1) + '\tTRUE' + '\t' + '' + '\n') 
                temp.write('-27\t' + environment  + '\t-2\t' + organism + '\tMG-RAST metagenome study\t' + evidence_2 + '\t' + str(score_2) + '\tTRUE' + '\t' + '' + '\n') 
temp.close()

# Extract PROC - ENV associations
# As "case_1" we refer to process - environments associations, meaning a certain process is our entry point 
# Likewise, as "case_2" we refer to environment - processes associations, with a certan evironment as entry point
with open(output_file, "a") as temp:
    
    for environment, processes in proc_env_co_occurences.items():
           
        background_envo = envo_backgrounds[environment]
        
        for process, count in processes.items():
            
            background_proc = gomfs_terms_backgrounds[process]
            
            evidence_1 = '%d of %d samples' % (count, background_proc)
            evidence_2 = '%d of %d samples' % (count, background_envo)
            
            if count > background_envo or count > background_proc:
                print("Something is going wrong in the countings on PROC - ENV associations")
                print(process, environment)
                continue

            score_1 = 2.0 * math.sqrt(float(count) / background_envo**0.2)
            score_2 = 2.0 * math.sqrt(float(count) / background_proc**0.2)


            if process != '' and environment != '':
                # and we print the final output! again as a tab seperated file
                temp.write('-23\t' +  process + '\t-27\t' + environment + '\tMG-RAST metagenome study\t' + evidence_1 + '\t' + str(score_1) + '\tTRUE' + '\t' + '' + '\n')
                temp.write('-27\t' +  environment + '\t-23\t' + process + '\tMG-RAST metagenome study\t' + evidence_2 + '\t' + str(score_2) + '\tTRUE' + '\t' + '' + '\n')            
temp.close()


# Extract ORG - PROC associations
# Likewise, "case_1" : organism - processes, with certain taxon as entry point
# And "case_2": process - organisms, with a process as entry point
with open(output_file, "a") as temp:
    
    for process, organisms in proc_org_co_occurences.items():
        
        background_proc = gomfs_terms_backgrounds[process]
        
        for organism, count in organisms.items():
            
            org_background = ncbi_ids_backgrounds[organism]
            
            evidence_1 = '%d of %d samples' % (count, org_background) 
            evidence_2 = '%d of %d samples' % (count, background_proc)
            
            if count > background_proc or count > org_background:
                print("Something is going wrong in the countings on ORG - PROC associations")
                print(organism, process)
                continue
            
            score_1 = 2.0 * math.sqrt(float(count) / background_proc**0.2)
            score_2 = 2.0 * math.sqrt(float(count) / org_background**0.2)
            
            if process != '' and organism != '':
                # and we print the final output! again as a tab seperated file. The format requires 9 fields.
                # type 1 , id 1, type 2 , id 2 , source, evidence, score, boolean, url.
                temp.write('-2\t' +  organism + '\t-23\t' + process + '\tMG-RAST metagenome study\t' + evidence_1 + '\t' + str(score_1) + '\tTRUE' + '\t' + '' + '\n')

                temp.write('-23\t' + process + '\t-2\t' + organism + '\tMG-RAST metagenome study\t' + evidence_2 + '\t' + str(score_2) + '\tTRUE' + '\t' + '' + '\n')
temp.close()            

print("Step 4/4 : completion of output file creation")
