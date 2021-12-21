#!/usr/bin/python3.5

'''
ABOUT:
Script name: get_core_met_processes.py
Developed by: Haris Zafeiropoulos
Framework: PREGO - WP5
Usage: ./get_core_met_processes.py
'''

import os
import sys

path = os.getcwd()
input_file = path + "/more_500kos_ncbi_ids.tsv"
output_file_1 = path + "/rare_met_processes.tsv"
output_file_2 = path + "/core_met_processes.tsv"
# We are working with 21,139 genomes that have at least 500 KO terms 
data = open(input_file, "r")

# Example: {K01719: YYYY, ZZZZ, FFFFF}
process_ncbi_ids = {}

for line in data:

   line = line[:-1]
   ncbi_id = line.split("\t")[0]
   processes = line.split("\t")[1].split("|")

   for term in processes: 

      if term not in process_ncbi_ids.keys():

         process_ncbi_ids[term] = [ncbi_id]

      else: 

         process_ncbi_ids[term].append(ncbi_id)

# Example: [(KO1719, 324), (KO2342, 342),...]
ko_number_of_genomes = []

for key, value in process_ncbi_ids.items():
   ko_number_of_genomes.append( (key, len(value)) )

with open(output_file_1, "w") as f:
   for pair in ko_number_of_genomes:

      if pair[1] < 10:

         ncbi_ids = ""
         for id in process_ncbi_ids[pair[0]]:
            ncbi_ids += id + "|"

         ncbi_ids = ncbi_ids[:-1]

         f.write(pair[0] + "\t" + ncbi_ids + "\n")


with open(output_file_2, "w") as f:
   for pair in ko_number_of_genomes:

      if pair[1] > 18000:

         f.write(pair[0] + "\t" + str(pair[1]) + "\n")



