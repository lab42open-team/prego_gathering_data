#!/usr/bin/python3.5

########################################################################################
# script name: extract_mgrast_associations.py
# path on oxygen: /data/databases/scripts/gathering_data/mgrast
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL:
# this script aims to extract ORGANISMS -ENVIRONMENTS associations from the amplicon data
# retrieved from the MG-RAST resource.
########################################################################################
#
# usage: ./extract_mgrast_markergene_associations.py
#
########################################################################################

import json, os, sys
import tagger
import math
from functions_mgrast import *

# Paths
scripts_path = '/data/databases/scripts/gathering_data/mg_rast'
experiments_path = '/data/experiments'

database_prefix = '/data/databases/mg_rast'
database_prefix_amplicon = database_prefix + '/amplicon'

metadata_input_file = database_prefix_amplicon + '/amplicon_metadata_retrieved_TEST.tsv'   
output_file = experiments_path + "/mgrast_markergene_associations.tsv"


# Match NCBI Ids and ENVO terms to the metadata retrieved 
ncbi_taxonomy_ids_dic, envo_terms_dic = match_taxa_and_envs_with_ontology_terms(metadata_input_file)
print("Step 1/6 : function match_taxa_and_envs_with_ontology_terms() has been completed.")

# Bring together ENVO terms and NCBI Ids on each sample
envo_taxa_dic = taxa_envo_in_samples(metadata_input_file, ncbi_taxonomy_ids_dic, envo_terms_dic) 
print("Step 2/6 : function taxa_envo_in_samples() has been completed")

# Get organisms - environments cooccurences
envo_taxa_coocurrences = get_env_org_coocurrences(envo_taxa_dic) 
print("Step 3/6 : function get_env_org_coocurrences() has been completed")

# Get the background for the ENVO terms found in the samples retrieved, i.e in how many samples each ENVO terms occurs
envo_background = get_envo_background(envo_taxa_dic)
print("Step 4/6 : function get_envo_background() has been completed")

# Get the background for the NCBI Ids found; i.e in how many samples each Id occurs 
taxon_background = get_taxon_background(envo_taxa_dic)
print("Step 5/6 : function get_taxon_background() has been completed")

# Get the total number of samples retrieved
total_number_of_samples = len(envo_taxa_dic)
print(total_number_of_samples)

print("Step 6/6 : initiation of output file creation")
# Now print ENVO - ORGANISMS associations in a file
with open(output_file, "w+") as temp:
    for envo in envo_taxa_coocurrences.keys():
       
        if envo != '':
       
            for taxon in envo_taxa_coocurrences[envo]:
               
                if len(envo) > 4:
               
                    count = envo_taxa_coocurrences[envo][taxon]
                    background_envo = envo_background[envo]
                    evidence = '%d of %d samples' % (count, background_envo)
                    score = 2.0*math.sqrt(float(count)/background_envo**0.2)
                   
                    # and we print the final output! again as a tab seperated file
                    temp.write('-27\t' + envo + '\t-2\t' +  taxon  + '\tMG-RAST amplicon study\t'  + evidence + '\t' + str(score) +  '\tTRUE' + '\n')
                    temp.write('-2\t' + taxon + '\t-27\t' + envo +'\tMG-RAST amplicon study\t'  + evidence + '\t' + str(score) +  '\tTRUE' + '\n')
temp.close() 
    
print("Step 6/6 : completion of output file creation in " + output_file)
