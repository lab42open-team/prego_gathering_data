#!/usr/bin/python3.5

'''
ABOUT
Script name: build_species_ko_ref_db.py
Path on oxygen: /data/databases/struo
Developed by: Haris Zafeiropoulos
Framework: PREGO - WP2

GOAL
Map the UniRef terms annotated in GTDB oriented genomes 
to KO terms if related to metabolsm

NOTES
The script uses the multiprocessing library to run faster.  
You can set the number of processes you want by changing 
the value of the number_processes variable. 

USAGE
Change the values of the directory_path and ref_species_ko_terms_db variables
and then run
./build_species_ko_ref_db.py
'''


import os, re, time, random, json
import multiprocessing


# Input and output files 
directory_path = "/data/databases/struo/genomes_13/"
ref_species_ko_terms_db = "/data/databases/struo/in_house_ref_ko_species_db_genomes_13.tsv"

# Files with the mapping terms and the KO terms that are metabolism realted
corresponding_terms_of_uniprot_to_ko = "/data/databases/struo/uniref_ko.json"
metabolism_related_ko_terms_file = "/data/databases/struo/only_metabolism_KO_terms.txt"

# Keep all metabolism related terms in a set
metabol_ko_terms_global = set(line.strip() for line in open(metabolism_related_ko_terms_file))


# Main function
def map_uniref_to_ko(filename):

	dictionary = corresponding_terms_of_uniprot_to_ko
	directory = directory_path
	metabol_mapping_file = metabol_ko_terms_global 
	output_file = ref_species_ko_terms_db

   # create the variable you will get the info you need
	my_string = ''
	
   # read the json file with the corresponding terms
	with open(dictionary) as json_file:
		corresponding_terms_dictionary = json.load(json_file)
	
      # each file corresponds to a ncbi taxonomy id
		taxon_id = filename.split("_")[1].split(".")[0]
		file = directory + filename
		read_file = open(file, "r")
		
		print("the taxon id under study is: " + taxon_id)
		
      # We need to find the list of KO terms regarding the species under study in each file 
		list_with_ko_terms_of_the_species = set()
		
      # Each line is a protein coding gene of the species under study; thus, every line has a uniref id. 
      # Our approach is to get this id and check if that maps to a KO correspondint term
		counter = 0
		for line in read_file:
	
         # Get the info regarding each every gene
			if line[0] == ">":

            # Keep the uniref id of the gene			
				uniref = line.split("|")[0].split("_")[1]
				counter += 1

	###################     DICTIONARY FOR UNIREF - KO  LIST          #####################

            # Iterate through the dictionary of uniref - ko terms				
				for pair in corresponding_terms_dictionary:

               # Chech if there is the uniref id that there is on the sequence title under study 
					if pair['uniref_id'] == uniref:						
						ko_term_found = pair['ko_term']

                  # If the term corresponds to a KO one, then, check if this term is related with metabolism
						if ko_term_found in metabol_mapping_file:
							print(uniref)
							print(ko_term_found)
							
                     # If yes, then keep the KO term as on related with the genome under study; 
                     # add it to a list which will be assigned to it in the end
							list_with_ko_terms_of_the_species.add(ko_term_found)
							break

                  # If it is not a metabolism reletad, then go to the next sequence of the genome
						else:
							break

               # If it is not a match, then keep searching for the corresponding term; if it exists
					else:
						continue
				
				modulo = int(counter%100)	
				if modulo == 0:
					print("for the genome: " + taxon_id + " " +  str(counter) + " sequences have been parsed" )

      # When all sequences of the genome have been parsed, then print that it's ok..
		print("the genome: " + taxon_id  + " has been done! Here you may see the KO terms found: " + "\n")

      # ... and add a new entry in the dictionary with the genome under study
      #  my_dictionary[taxon_id] = list_with_ko_terms_of_the_species

		with open(output_file, "a") as outfile:
			for ko in list_with_ko_terms_of_the_species:
				my_string += ko + "|"
			outfile.write(taxon_id + "\t" + my_string + "\n")
	
	json_file.close()
	return my_string


# Now that we have built everything we need, we can set a pool with processes
# we make a list with the items that our function will run for
list_with_files = []
for filename in  os.listdir(directory_path):
	list_with_files.append(filename)

# and we build a pool... 
if __name__ == '__main__':
   
   # Set the maximum number of processes that our script is going to use
	number_processes = 10
   
   # Start a pool
	pool = multiprocessing.Pool(number_processes)
   
   # Map the function that is about to run with the list of items for which it has to iterate
	pool.map(map_uniref_to_ko, list_with_files)
   
   # Can't submit new tasks to pool of worker processes. 
   # Once all tasks complet, the worker processes will exit.
	pool.close()
   
   # The code in __main__ must wait until all our tasks are complete before continuing
	pool.join()


'''
ATTENTION! In the first place we had a ThreadPoolExecutor. 
However, it turns out that parsing files is a CPU bound and not I/O as I thought.
That said, we need to replace our initial executor with a ProcessPoolExecutor
with concurrent.futures.ProcessPoolExecutor(max_workers = 4) as executor:	

Here is where the concurrent.futures part is taking place:
	future_to_url = {executor.submit(map_uniref_to_ko, filename): filename for filename in  os.listdir(directory_path)}
'''
