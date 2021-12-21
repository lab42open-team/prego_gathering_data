#!/usr/bin/python

import sys, os, getopt



def map_gtdb_to_ncbi_id():

   metadata_file = open("../metadata_1per-GTDB-Spec_gte50comp-lt5cont_wtaxID.tsv","r")

   gtdb_taxid_ncbi_id  = {}
   ncbi_id_gids = {}
   for line in metadata_file:
      line       = line.split("\t")
      ncbi_id    = line[73]
      gid        = line[55]
      gtdb_taxid = line[-1][:-1]
      if ncbi_id not in ncbi_id_gids:
         ncbi_id_gids[ncbi_id] = []
         ncbi_id_gids[ncbi_id].append(gid)
      else:
         ncbi_id_gids[ncbi_id].append(gid)

      gtdb_taxid_ncbi_id[gtdb_taxid] = ncbi_id

   return gtdb_taxid_ncbi_id, ncbi_id_gids


def parse_arguments():

   opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["inputDir=", "output="])
   current_directory = os.getcwd()

   for i in opts:

      if i[0] == "--inputDir" or i[0] == "-i":
         dir_path = i[1]

      elif i[0] == "--output" or i[0] == "-o":
         output_file = current_directory + "/" + i[1]


   file_list = []
   try:
      for filename in  os.listdir(dir_path):
         file_list.append(filename)
   except:
      print("ValueError: No directory was given.")
      raise 

   return dir_path, output_file, file_list


def parse_list_of_files(outputfile, path_to_files, files_to_parse, gtdb_to_ncbi_maps): 

   o = open(outputfile, "w")

   for filename in files_to_parse: 
      
      print(filename)
      
      genome_file = open(path_to_files + filename, "r")
      gtdb_id     = filename.split("_")[1].split(".")[0]
      ncbi_id     = gtdb_to_ncbi_maps[gtdb_id]
      for line in genome_file:

         if line[0] != ">": 
            continue
         
         else:
            uniref = line.split("|")[0].split("_")[1]
            
            o.write(ncbi_id + "\t" + uniref + "\n")


# and we build a pool... 
if __name__ == '__main__':
   
   # Run the general
   input_path, ofile, lfiles = parse_arguments()

   gtdb_taxid_ncbi_id, ncbi_id_gids = map_gtdb_to_ncbi_id()

   parse_list_of_files(ofile, input_path, lfiles, gtdb_taxid_ncbi_id)
   
