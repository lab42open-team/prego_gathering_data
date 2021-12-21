#!/usr/bin/python3.5

########################################################################################
# script name: jgi_get_metadata_for_each_id.py
# path on oxygen: /data/databases/scripts/gathering_data/jgi
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL
# After getting the JGI genome IDs, their metadata can be retrieved
# This script aims to get the KEGG Orthology (KO) terms as well as metadata regarding
# the environment of the sample
########################################################################################

import time, datetime, sys
import unittest, re
from functions_jgi import *

# Regarding the .log file 
start = datetime.datetime.now()
date =  datetime.date.today()
old_stdout = sys.stdout
log_file_path = "/data/databases/scripts/gathering_data/logfiles/" 
log_file_name = log_file_path + "jgi_get_metadata_and_KO_terms"+ str(date) + ".log"
log_file = open(log_file_name,"w")
sys.stdout = log_file

# Open the .txt file that contains the JGI Genome IDs
jgi_path = "/data/databases/jgi/"
jgi_genome_ids_file_bacteria_isolates = jgi_path + "jgi_genome_ids_bacteria_isolates.tsv"
jgi_genome_ids_file_bacteria_sags = jgi_path + "jgi_genome_ids_bacteria_sags.tsv"
jgi_genome_ids_file_bacteria_mags = jgi_path + "jgi_genome_ids_bacteria_mags.tsv"
jgi_genome_ids_file_archaea_isolates = jgi_path + "jgi_genome_ids_archaea_isolates.tsv"
jgi_genome_ids_file_archaea_sags = jgi_path + "jgi_genome_ids_archaea_sags.tsv"
jgi_genome_ids_file_archaea_mags = jgi_path + "jgi_genome_ids_archaea_mags.tsv"

jgi_genome_ids_file_bacteria_isolates_metadata = jgi_path + "jgi_genome_ids_bacteria_isolates_metadata.tsv"
jgi_genome_ids_file_bacteria_sags_metadata = jgi_path + "jgi_genome_ids_bacteria_sags_metadata.tsv"
jgi_genome_ids_file_bacteria_mags_metadata = jgi_path + "jgi_genome_ids_bacteria_mags_metadata.tsv"
jgi_genome_ids_file_archaea_isolates_metadata = jgi_path + "jgi_genome_ids_archaea_isolates_metadata.tsv"
jgi_genome_ids_file_archaea_sags_metadata = jgi_path + "jgi_genome_ids_archaea_sags_metadata.tsv"
jgi_genome_ids_file_archaea_mags_metadata = jgi_path + "jgi_genome_ids_archaea_mags_metadata.tsv"

pairs_of_files = (
    (jgi_genome_ids_file_bacteria_isolates, jgi_genome_ids_file_bacteria_isolates_metadata), 
    (jgi_genome_ids_file_archaea_isolates, jgi_genome_ids_file_archaea_isolates_metadata), 
    (jgi_genome_ids_file_archaea_mags, jgi_genome_ids_file_archaea_mags_metadata), 
    (jgi_genome_ids_file_archaea_sags, jgi_genome_ids_file_archaea_sags_metadata), 
    (jgi_genome_ids_file_bacteria_sags, jgi_genome_ids_file_bacteria_sags_metadata),
    (jgi_genome_ids_file_bacteria_mags, jgi_genome_ids_file_bacteria_mags_metadata)
   )

# Now perform what is needed to get the metadata 
for pair in pairs_of_files:

    get_metadata_and_KO_for_each_file(pair[0], pair[1])
    
    metadata_file =  open(pair[1], 'r')
    ko_terms_per_environment_tidy_file = jgi_path + "ko_terms_per_jgi_environment_tidy_format.txt"
    ko_terms_per_environment_tidy_file = open(ko_terms_per_environment_tidy_file, "a")
    
    for line in metadata_file:
        elements = line.split('\t')
        environments = elements[2:6]
        ko_terms = elements[-1].split(';')
        
        # as we add a ';' after each KO term, we will have an empty element on the ko_term list which we need to remove
        ko_terms = ko_terms[:-1]
        for ko_term in ko_terms:
            for environment in environments:
                ko_terms_per_environment_tidy_file.write(ko_term + '\t' + environment + "\n")
    ko_terms_per_environment_tidy_file.close()            


# End with the stats and .log file
print("script execution stared at:", start)
end = datetime.datetime.now()
print("Script execution ended at:", end)
total_time = end - start
print("Script totally ran for :", total_time)

sys.stdout = old_stdout
log_file.close()

