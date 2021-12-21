#!/usr/bin/python3

# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
# usage: ./build_associations_file.py

"""
This script parses the all_genomes_gtdb_ids.tsv file, a .tsv file 
with 2 cols where in the 1st one there is a GTDB id and on the 2nd
we have all the KOs found on the representative genome of that id 
in the GTDB v.89. By making use of the GTDB metadata file 
it makes a single set of KOs per NCBI id and links those with their
corresponding GTDB ids. 

e.g. -2	2162557	-90	KO:K01696	Struo	Genome annotation	0.5	GCA_003212335.1
     -2	43263	-90	KO:K01928	Struo	Genome annotation	0.5	GCA_000467105.1,GCA_900156545.1,GCA_003205495.1,GCA_000474255.1
"""


import sys

gtdb_file_input = open("all_genomes_gtdb_ids.tsv","r")
metadata_file = open("../metadata_1per-GTDB-Spec_gte50comp-lt5cont_wtaxID.tsv","r")
output_file   = open("struo_database_pairs.tsv", "w")

# gid: ncbi_genbank_assembly_accession
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

ncbi_id_kos = {}
for line in gtdb_file_input:
    line     = line.split("\t")
    gtdb_id  = str(line[0])
    kos      = line[1].split("|")
    ncbi_id  = gtdb_taxid_ncbi_id[gtdb_id]

    if ncbi_id not in ncbi_id_kos:
        ncbi_id_kos[ncbi_id] = set()
    for ko in kos:
        ncbi_id_kos[ncbi_id].add(ko)

for ncbi_id, kos in ncbi_id_kos.items():

    gids     = ncbi_id_gids[ncbi_id]

    if gids:
       
        ncbi_url = 'https://www.ncbi.nlm.nih.gov/assembly/?term='
        for gid in gids: 
            ncbi_url += gid + "+"
        ncbi_url = ncbi_url[:-1]
 
        for ko in kos:
            if ko[0] == "K":
                output_file.write("-2" + "\t" + ncbi_id + "\t" + "-90" + "\t" + "KO:" + ko + "\t" +\
                        "Struo"+ "\t" + "Genome annotation" + "\t" +\
                        "50" + "\t" + "TRUE" + "\t" + ncbi_url + "\n")
