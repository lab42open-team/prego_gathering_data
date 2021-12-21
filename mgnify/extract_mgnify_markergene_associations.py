#!/usr/bin/python3.5

########################################################################################
# script name: extract_mgnify_markergne_associations.py
# developed by: Haris Zafeiropoulos*
# framework: PREGO - WP2
########################################################################################
# GOAL:
# Aim of this script is to parse the information retrieved (by the get_mgnify_markergene_data.py
# script) and mine associations from these data
########################################################################################
# *this script is based on the mgnify.py script that is stored at /data/experiments and 
# implemented by Lars & Haris on 1st PREGO hackathon
############################################################################################
#
# usage: ./extract_mgnify_markergene_associations.py > mgnify_markergene_associations.tsv 2> associations_log_file.txt
#
############################################################################################

# dependencies
import json
import math
import os, sys
import tagger.tagger
import matplotlib.pyplot as plt

# to call the Rscript for estimating the ML parameters and probabilities
import subprocess

###########################    here is where the coding starts!    ###########################

# here is a list with the attributes we need to get from the metadata files
# [next]: replace with a dictionary mapping the fields to the entity types the tagger should mine
attributes_to_mine = ['sample-desc', 'environment-biome', 'environment-feature', 'environment-material']

# path to look for the dictionary
dictionary_prefix = '/data/dictionary/tagger'

# import a class called "Tagger()" from the "tagger.tagger" library
prego_tagger = tagger.tagger.Tagger()

# import a specific function called "load_names()" from the "Tagger()" class
prego_tagger.load_names('%s_entities.tsv' % dictionary_prefix, '%s_names.tsv' % dictionary_prefix)
# sys.stderr.write("Tagger dictionaries loaded. \n")

# path to look for the data files we intend to use
database_prefix = '/data/databases/mgnify/marker_gene_data'

# files for statistics
parent_directory = '/data/databases/scripts/gathering_data/mgnify'
sample_entity_sources_stats_file = parent_directory + "/" +'sample_entity_sources.tsv'
taxon_sample_abundance_stats_file = parent_directory + "/" + 'taxon_sample_abundance.tsv'
distribution_for_ncbi_ids = parent_directory + "/" + 'ncbi_ids_distribution.tsv'
distribution_for_metadata_terms = parent_directory + "/" + 'metadata_terms_distribution.tsv'
background_file_1 = parent_directory + "/" + "ncbi_background.tsv"
background_file_2 = parent_directory + "/" + "metadata_background.tsv"
prob_1_file = parent_directory + "/" + "probabilities_ncbi.tsv"
prob_2_file = parent_directory + "/" + "probabilities_metadata.tsv"

#################################################################################################################################
# PART 1:    make a dictionary with SAMPLES - ATTRIBUTES - HOST (if available) (i.e {SAMPLE_ID : {(-2, '2345'): {'host-tax-id'}}})
#################################################################################################################################

# use the directory where the MGnify "run" files are deposited, as for a "for" loop
# and make a dictionary with the "run (key) : sample (value)" pairs
run_sample = {}
for filename in os.listdir(database_prefix+'/runs'):
    if filename.endswith('.json'):
        data = json.load(open(database_prefix+'/runs/'+filename, 'r'))
        for entry in  data['data']:
            run = entry['id']
            sample = entry['relationships']['sample']['data']['id']
            run_sample[run] = sample

# use the directory where the MGnify "sample" files are located, as for a "for" loop
# and make a dictionary with the "sample (key) : sample_entity_sources (value)" pairs. The "sample_entity_sources" is another dictionary.
sample_entity_sources = {}
for filename in os.listdir(database_prefix+'/samples'):
    # sys.stderr.write("*** filename under study: " + filename + "*** \n")

    if filename.endswith('.json'):
        data = json.load(open(database_prefix+'/samples/'+filename, 'r'))
        for entry in data['data']:
            sample = entry['id']
            sample_entity_sources[sample] = {}

            # remember that "entry" is the for-element that takes into account all the entries of each .json file 
            for attribute in entry['attributes']:

                # here, we set the variable values for each element of the "attributes" we examine
                value = entry['attributes'][attribute]

                # if the value of the attribute is NoneType then go to the next attribute
                if value is None:
                    continue
                else:

                    # first, we check whether the attribute we examine, is from a host species
                    # the "value" variable is the value of each attribute of attributes; changes in every in-loop 
                    if attribute == 'host-tax-id':
    
                        # here is the first TAGGER-specific structure: with "-2" the NCBI taxonomy ID is specified and we also keep the ID
                        entity = (-2, str(value))
    
                        # then we need to add the elements of the "entity" variable in our dictionary. Hence, in case that this entity has not appeared before, then it is added on the main dictionary and a structure like this is built (SAMPLE_ID: {(-2 NCBI_ID : set())})
                        if entity not in sample_entity_sources[sample]:
                            sample_entity_sources[sample][entity] = set()
    
                        # if this entity has already appeared, then the attribute under study is added and thus, we have a structure like this: {SAMPLE_ID : {(-2, '2345'): {'host-tax-id'}}}
                        sample_entity_sources[sample][entity].add(attribute)
    
                    # if attribute is not 'host-tax-id' then we need to check if it is one of the elements of our initial attributes list (line 14)
                    elif attribute in attributes_to_mine:
    
                        # if yes, then we use the TAGGER to make a list of matches of entities
                        matches = prego_tagger.get_matches(value, 0, [-21, -25,-27], auto_detect=False, protect_tags=False, max_tokens=10)
    
                        # and for each of those we do as for the 'host-tax-id' case
                        for match in matches:
                            for entity in match[2]:
                                if entity not in sample_entity_sources[sample]:
                                    sample_entity_sources[sample][entity] = set()
                                sample_entity_sources[sample][entity].add(attribute)


# Statistics time! Here we build a .tsv file with all tagger's findings for each sample. E.g
# ERS2393275    -27    ENVO:00002034	
# ERS2393275    -27    ENVO:00000137	
with open (sample_entity_sources_stats_file, 'w') as stat1:
    for sample in sample_entity_sources:
        for entity in sample_entity_sources[sample]:
            stat1.write(str(sample) + "\t" + str(entity[0]) + "\t" + str(entity[1]) + "\n")
stat1.close()

#################################################################################################################################
# PART 2:    make a dictionary with TAXA - SAMPLES - ABUNDANCES
#################################################################################################################################

# let us now move forward on the abundances data (OTU tables)
# we create a dictionary which finally will be like this: NCBI_ID : SAMPLE_ID : ABUNDANCE_ON_THE_OTU_TABLE
taxon_sample_abundance = {}
for filename in os.listdir(database_prefix+'/abundances'):
    
    if filename.endswith('.tsv'):
        with open(database_prefix+'/abundances/'+filename, 'r') as f:

            # we make a list and we read the first line of the file. Then we keep in that list only those samples that we also found on the RUNs data (present on the run_samples dictionary)
            col_sample = []
            for col in f.readline().strip().split('\t'):        # the readline() python function, reads only the first line of a file
                
                if col in run_sample:
                    col_sample.append(run_sample[col])
                else:
                    col_sample.append(None)

            # in the next lines we get the taxonomy of each OTU, which is in the first column
            # we make a dictionary called "taxlevel_name" where we put the taxonomy level and its corresponding value for each OTU
            for line in f:   
                fields = line.strip().split('\t')
                taxlevel_name = {}

                # we keep "fields[0]" as in it there is the taxonomy and we split it with the ";" as separator
                # than we make a list with all the taxa level, e.g subsubfields = [Archaea, Euryarchaeota, Candidatus_Altiarchaeales]
                # here is an example of how the taxonomy of those tables look like: "sk__Archaea;k__;p__Euryarchaeota;c__;o__Candidatus_Altiarchaeales"
                for subfield in fields[0].split(';'):
                    subsubfields = subfield.split('__')

                    # we add an entry in the dictionary called "taxlevel_name" where we have the Domain as key (e.g Bacteria) and the rest of the taxonomy as value
                    if len(subsubfields) > 1 and subsubfields[1] != '':
                        taxlevel_name[subsubfields[0]] = subsubfields[1].replace('_', ' ');

                # if we have a species-level annotation then .. 
                if 'g' in taxlevel_name and 's' in taxlevel_name:
                    taxlevel_name['s'] = taxlevel_name['g']+' '+taxlevel_name['s']

                # if not, then  we keep the lowest taxon-level, by making use of the TAGGER on the rest of the taxonomy levels
                taxa = []
                for taxlevel in ['s', 'g', 'f', 'o', 'c', 'p', 'k']:
                    if taxlevel in taxlevel_name:
                        matches = prego_tagger.get_matches(taxlevel_name[taxlevel], 0, [-2], auto_detect=False, protect_tags=False, max_tokens=5)
                        if len(matches) > 0:
                            for entity in matches[0][2]:
                                taxa.append(int(entity[1]))
                            break

                # we keep the number of samples 
                n = len(fields)-1

                # and then for each sample
                for i in range(1, n):

                    # we keep the abundance as a float for the OTU under study
                    abundance = float(fields[i])

                    # likewise, we keep the sample where the taxon and its relative abundance were found 
                    sample = col_sample[i]

                    # and if its abundance is non-zero
                    # IMPORTANT ADD: and if the corresponding sample of the run is included on MGnify 
                    if abundance > 0 and sample is not None:
                        
                        # we check through the list of taxa we made, when we were getting the taxonomy, and
                        # for each taxon level, we check if there is already in our taxon_sample_abundance dictionary, if not we add it
                        for taxon in taxa:
                            if taxon not in taxon_sample_abundance:
                                taxon_sample_abundance[taxon] = {}

                            # we do the same for the sample name 
                            if sample not in taxon_sample_abundance[taxon]:
                                taxon_sample_abundance[taxon][sample] = 0.0

                            # if we already have this combination we just sum it with the previous
                            taxon_sample_abundance[taxon][sample] += abundance


# Statistics time! Here we build a second .tsv file with all the samples ids where an NCBI taxonomy Id was found. E.g
# 85108    SRS775761
# 85108    ERS756676
with open (taxon_sample_abundance_stats_file, 'w') as stat2:
    for taxon in taxon_sample_abundance:
        for sample in taxon_sample_abundance[taxon]:
            stat2.write(str(taxon) + "\t" + str(sample) + "\n")
stat2.close()

#######################################################################################################
#    PART 3:        find the associations
#######################################################################################################

# calculate the background for better statistics
entity_background = {}
for sample in sample_entity_sources:
    for entity in sample_entity_sources[sample]:
        if entity in entity_background:
            entity_background[entity] += 1
        else:
            entity_background[entity] = 1

total_number_of_samples = len(sample_entity_sources)
unique_host_terms = [] ; unique_process_terms = [] ; unique_tissues_terms = [] ; unique_envo_terms = []
count_of_samples_for_each_ncbi = {}

for sample in sample_entity_sources:
    for entry in sample_entity_sources[sample]:
        
        entry_0 = str(entry[0])
        entry_1 = str(entry[1])
        
        if entry_0 == '-2':
            if entry_1 not in unique_host_terms:
                unique_host_terms.append(entry[1])
        elif entry_0 == '-21':
            if entry_1 not in unique_process_terms:
                unique_process_terms.append(entry[1])
        elif entry_0 == '-25':
            if entry_1 not in unique_tissues_terms:
                unique_tissues_terms.append(entry[1])
        elif entry_0 == '-27':
            if entry_1 not in unique_envo_terms:
                unique_envo_terms.append(entry[1])

total_number_of_host_terms = len(unique_host_terms)
total_number_of_tissue_terms = len(unique_tissues_terms)
total_number_of_process_terms = len(unique_process_terms)
total_number_of_envo_terms = len(unique_envo_terms)


#### the Distributions 

## The NCBI IDs -- these lead to a log normal distribution
count_samples_ncbi_ids = {}
list_for_background = []
list_with_taxa = []
for taxon in taxon_sample_abundance:
    list_with_taxa.append(taxon)
    count_of_samples_where_taxon_found = len(taxon_sample_abundance[taxon])
    list_for_background.append(count_of_samples_where_taxon_found)
    if count_of_samples_where_taxon_found in count_samples_ncbi_ids:
        count_samples_ncbi_ids[count_of_samples_where_taxon_found] += 1
    else:
        count_samples_ncbi_ids[count_of_samples_where_taxon_found] = 1

# make a file with the NCBI Id backgrounds
ncbi_id_and_its_background = {}
with open(background_file_1, "w") as background_1:
    for i in range(len(list_for_background)):
        ncbi_id_and_its_background[list_with_taxa[i]] = list_for_background[i]
        background_1.write(str(list_with_taxa[i]) + "\t" + str(list_for_background[i]) + "\n")
background_1.close()

# here is the list with the elements of the x axis of the distribution (independent variable)
number_of_ncbi_ids = []
# and here the list with the elements of the y axis
number_of_samples_where_an_id_was_found = []
# fill in the two lists described above
for key,value in count_samples_ncbi_ids.items():
    number_of_samples_where_an_id_was_found.append(key)
    number_of_ncbi_ids.append(value)    

with open (distribution_for_ncbi_ids, 'w') as distribution_file_1:
    for i in range(len(number_of_ncbi_ids)):
        distribution_file_1.write(str(number_of_ncbi_ids[i]) + "\t" + str(number_of_samples_where_an_id_was_found[i]) + "\n")
distribution_file_1.close()


## The metadata -- these lead to a 1/x distribution
metadata_term_and_its_background = {}
# for host in unique_host_terms:       ## later, check if I can have the same per data type
for sample in sample_entity_sources:
    for entry in sample_entity_sources[sample]:
        for item in entry:
            item = str(item)
            if item[0] == '-':
                continue
            else:
                term = item            
                # sys.stderr.write("term is: " + str(term) + "\n")
                if term in metadata_term_and_its_background:
                    metadata_term_and_its_background[term] += 1
                else:
                    metadata_term_and_its_background[term] = 1

# make a file with the metadata backgrounds
with open(background_file_2, "w") as background_2:
    for key,value in metadata_term_and_its_background.items():
        background_2.write(str(key) + "\t" + str(value) + "\n")
background_2.close()

# count the terms that where found in the same number of samples
counts_of_metadata_terms_found_at_number_of_samples = {}
for term, count in metadata_term_and_its_background.items():
    if count in counts_of_metadata_terms_found_at_number_of_samples:
        counts_of_metadata_terms_found_at_number_of_samples[count] += 1
    else:
        counts_of_metadata_terms_found_at_number_of_samples[count] = 1
 
# here is the list with the elements of the x axis of the distribution (independent variable)
number_of_metadata_terms = []
# and here the list with the elements of the y axis
number_of_samples_where_a_metadatum_term_was_found = []

for key,value in counts_of_metadata_terms_found_at_number_of_samples.items():
    number_of_samples_where_a_metadatum_term_was_found.append(key)
    number_of_metadata_terms.append(value)

with open(distribution_for_metadata_terms, 'w') as distribution_file_2:
    for i in range(len(number_of_metadata_terms)):
        distribution_file_2.write(str(number_of_metadata_terms[i]) + "\t" + str(number_of_samples_where_a_metadatum_term_was_found[i]) + "\n")
distribution_file_2.close()


## Maximum - likelihood estimation
# We run a R script for calculating the probabilities giving as input the distribution files we built. 
subprocess.call("/data/databases/scripts/gathering_data/mgnify/calculate_probabilities.r") 

# We now read the files that the R script ended up with and keep the NCBI ID / metadata term and their corresponding probabilities to a dictionary
ncbi_id_prob = {}
metadata_term_prob = {}

with open (prob_1_file, 'r') as read_tsv_prob_ncbi:
    for row in read_tsv_prob_ncbi:
        row_elements = row.split("\t")
        ncbi_id = row_elements[0]
        prob = float(row_elements[2][:-1])
        ncbi_id_prob[ncbi_id] = prob
read_tsv_prob_ncbi.close()
    
with open (prob_2_file, 'r') as read_tsv_prob_metadata:
    for row in read_tsv_prob_metadata:
        row_elements = row.split("\t")
        metadata_term = row_elements[0]
        prob = float(row_elements[2][:-1])
        metadata_term_prob[metadata_term] = prob
read_tsv_prob_metadata.close()


# read the probabilities dictionaries built 
for key,value in ncbi_id_prob.items():
    sys.stderr.write(key + "\t" + str(type(key)) + "\n")
    sys.stderr.write(str(value) + "\t" + str(type(value)) + "\n")
    
sys.stderr.write("\n\n\n THE SECOND DICT \n\n\n")

for key,value in metadata_term_prob.items():
    sys.stderr.write(key + "\t" + str(type(key)) + "\n")
    sys.stderr.write(str(value) + "\t" + str(type(value)) + "\n")


# get the associations 
url_prefix = "https://www.ebi.ac.uk/ebisearch/search.ebi?db=allebi&query="
for taxon in taxon_sample_abundance:
    
    # we build an entity_count dictionary where we will keep the counts that this taxon was found with each unique attribute value
    # for each taxon present in the taxon_sample_abundance dictionary we built
    # after that, in the entity_count dictionary we have something like this:
    # (-27	ENVO:00002005 ) : 51    // this number is from running the script on 2020.07.12
    entity_count = {}
    entity_samples = {}
    
    # to do this, we parse the entities of each corresponding entry in the taxon_sample_abundance dictionary (where all the samples that this taxon was found have been kept)
    for sample in taxon_sample_abundance[taxon]:

        # and for all the samples that have also been recorded in the the sample_entity_sources dictionary (where the metadata of each sample have been kept)
        if sample in sample_entity_sources:

            # we track the attributes present on the sample
            for entity in sample_entity_sources[sample]:
                
                if entity in entity_count:
                    entity_count[entity] += 1
                else:
                    entity_count[entity] = 1
                    
                if entity in entity_samples:
                    entity_samples[entity].append(sample)
                else:
                    entity_samples[entity] = [sample]

    # then, we need to get how many times the same asso
    for entity in entity_count:
        
        # count is the number of samples where the co-occurrence happens
        count = entity_count[entity]
        samples = entity_samples[entity]
        
        # build the url with the samples that lead to the association under study
        if len(samples) == 1:
            url = url_prefix + "id:(" + samples[0] + ")"
        elif len(samples) <= 10:
            ids = ''
            for i in range(len(samples) - 1):
                ids = ids + "id:(" + samples[i] + ")%20OR%20"
            ids = ids + "id:(" + samples[-1] + ")"
            url = url_prefix + ids
        else:
            ids = ''
            for i in range(9):
                ids = ids + "id:(" + samples[i] + ")%20OR%20"
            ids = ids + "id:(" + samples[10] + ")"
            url = url_prefix + ids 
        
        # metadata_background is the number of samples where each matadata term of those related with the NCBI ID under study, appears in general
        metadata_background = entity_background[entity]        

        ## we do not need the following lines when we use MLE !!!
        # # taxon_background is the number of samples where the species of the association under study, has been found
        # taxon_background = len(taxon_sample_abundance[taxon])
        # calculating the probabilities - this will need to change to add the distribution oriented information for each case.        
        # p_x = float(taxon_background/total_number_of_samples)
        # p_y = float(metadata_background/total_number_of_samples)
        
        # conversions to reach the probabilities dictionaries information 
        str_taxon = str(taxon)
        metadata_term = entity[1]
        
        p_x = ncbi_id_prob[str_taxon]
        p_y = metadata_term_prob[metadata_term]
        p_joint = float(count/total_number_of_samples)
        
        # calculate the evidence and the score
        evidence = '%d of %d samples' % (count, metadata_background)
        lars_score = 2.0*math.sqrt(float(count)/metadata_background**0.2)
        mutual_information_score = float(p_joint*math.log2(p_joint/(p_x*p_y)))
        
        # print('-2\t%d\t%d\t%s\tMGnify\t%s\t%f\tTRUE\t%s' % (taxon, entity[0], entity[1], evidence, score, url))
        print('-2\t%d\t%d\t%s\tMGnify\t%s\t%f\tTRUE\t%s' % (taxon, entity[0], entity[1],evidence, lars_score, url))
