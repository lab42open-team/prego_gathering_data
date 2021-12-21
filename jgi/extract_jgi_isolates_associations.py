#!/usr/bin/python3

########################################################################################
# script name: extract_jgi_isolates_associations.py
# path on oxygen: /data/databases/scripts/gathering_data/jgi
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL
# extract all the PREGO associations for the JGI resource
########################################################################################

from functions_jgi import *
import sys

# General variables
knowledge_path   = "/data/knowledge"
jgi_path         = "/data/databases/scripts/gathering_data/jgi"

# Output file
output_file_knowledge        = knowledge_path + "/jgi_associations.tsv"
output_file_for_core_process = jgi_path + "/metabolic_cores.tsv"

# Input files
archaea_isolates = "/data/databases/jgi/jgi_genome_ids_archaea_isolates_metadata.tsv"
bacteria_isolates = "/data/databases/jgi/jgi_genome_ids_bacteria_isolates_metadata.tsv"

archaea_sags = "/data/databases/jgi/jgi_genome_ids_archaea_sags_metadata.tsv"
bacteria_sags = "/data/databases/jgi/jgi_genome_ids_bacteria_sags_metadata.tsv"

archaea_mags = "/data/databases/jgi/jgi_genome_ids_archaea_mags_metadata.tsv"
bacteria_mags = "/data/databases/jgi/jgi_genome_ids_bacteria_mags_metadata.tsv"

list_with_input_files = [
    (bacteria_isolates, 'Isolates', 0.9), \
    (bacteria_sags, 'Single Amplified Genome', 0.9), (bacteria_mags, 'Metagenome-Assembled Genome', 0.9),\
    (archaea_isolates, 'Isolates', 0.9), (archaea_sags, 'Single Amplified Genome', 0.9), \
    (archaea_mags, 'Metagenome-Assembled Genome', 0.9)
    ]
#----------------------------------------------------------------------------------

counter = 0
for entry in list_with_input_files:

    input_metadata_file = entry[0]
    data_type = entry[1]
    threshold = entry[2]
    output1 = output_file_knowledge

    core, number_of_singles = get_the_associations(input_metadata_file, data_type, threshold, output1)

    print('Associations for the case of ' + data_type + ' have been extracted successfully!')

    file_with_metadata = input_metadata_file.split("/")[-1]
    separator = ';'
    core_kos = separator.join(core)
    resource = '-'.join(file_with_metadata.split("_")[3:5])

    if counter == 0:

        with open(output_file_for_core_process, "w+") as core_output_file:
           core_output_file.write(
               resource + '\t' + "Number of unique NCBI IDs retrieved: " +\
               str(number_of_singles) + "\t" + "Covery percentage: " + str(threshold) +\
               "\t" + core_kos + '\n')
        core_output_file.close()

    else:
        with open(output_file_for_core_process, "a") as core_output_file:
           core_output_file.write(
               resource + '\t' + "Number of unique NCBI IDs retrieved: " +\
               str(number_of_singles) + "\t" + "Covery percentage: " + str(threshold) +\
               "\t" + core_kos + '\n')
        core_output_file.close()

    counter += 1



