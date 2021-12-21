#!/bin/bash

# ABOUT:
# Script name: map_gtdb_ids_to_ncbi.sh
# Developed by: Haris Zafeiropoulos
# Framework: PREGO - WP2

# GOAL:
# This script makes use of the GTDB metadata file to map the internal organisms ids to NCBI Taxonomy IDs

#USAGE
#./map_gtdb_ids_to_ncbi.sh <FILE>  <FILE_TO_SEARCH: metadata_1per-GTDB-Spec_gte50comp-lt5cont_wtaxID.tsv > > <NEW_FILE: head gtdb_ncbi_pair_ids.tsv>


FILE=$1
FILE_TO_SEARCH=$2

while read LINE; do

   INIT=$LINE
   LINE="^$LINE$"

   awk -F "\t" -v lvar="$LINE" '$NF~lvar {print $0}' $FILE_TO_SEARCH | awk -F "\t" -v init="$INIT" '{print init "\t" $78}'

done < $FILE




