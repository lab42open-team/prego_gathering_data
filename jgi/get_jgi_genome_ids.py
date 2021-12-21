#!/usr/bin/python3.5

########################################################################################
# script name: get_jgi_genome_ids.py
# path on oxygen: /data/databases/scripts/gathering_data/jgi
# developed by: Haris Zafeiropoulos
# framework: PREGO - WP2
########################################################################################
# GOAL
# Aim of this script is to get all the genome IDs that are available on JGI portal and
# coming from isolates or Single Amplified Genomes (SAGs) or Metagenome Assembled Genomes
# (MAGs).
########################################################################################
#
# usage: ./get_jgi_genome_ids.py
#
########################################################################################

from functions_jgi import the_jgi_genome_id_function
import datetime, sys


# Regarding the .log file 
start = datetime.datetime.now()
date =  datetime.date.today()
old_stdout = sys.stdout
log_file_path = "/data/databases/scripts/gathering_data/logfiles/" 
log_file_name = log_file_path + "jgi_getting_genome_ids"+ str(date) + ".log"
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

# Url bases that are going to be used for the drivers
base_url_bacteria_isolates = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonList&page=isolateList&domain=Bacteria%20isolates&seq_center=jgi#taxontable=results%3D100%26startIndex%3D" 
base_url_bacteria_sags = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonList&page=isolateList&domain=Bacteria%20SAGs&seq_center=jgi#taxontable=results%3D100%26startIndex%3D" 
base_url_bacteria_mags = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonList&page=isolateList&domain=Bacteria%20MAGs&seq_center=jgi#taxontable=results%3D100%26startIndex%3D" 
base_url_archaea_isolates = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonList&page=isolateList&domain=Archaea%20isolates&seq_center=jgi#taxontable=results%3D100%26startIndex%3D"
base_url_archaea_sags = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonList&page=isolateList&domain=Archaea%20SAGs&seq_center=jgi#taxontable=results%3D100%26startIndex%3D"
base_url_archaea_mags = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonList&page=isolateList&domain=Archaea%20MAGs&seq_center=jgi#taxontable=results%3D100%26startIndex%3D"

isolates_and_mags = ((base_url_bacteria_isolates, jgi_genome_ids_file_bacteria_isolates),
    (base_url_bacteria_mags, jgi_genome_ids_file_bacteria_mags),
    (base_url_archaea_isolates, jgi_genome_ids_file_archaea_isolates),
    (base_url_archaea_mags, jgi_genome_ids_file_archaea_mags)
    )

sags = ((base_url_archaea_sags, jgi_genome_ids_file_archaea_sags), (base_url_bacteria_sags, jgi_genome_ids_file_bacteria_sags))

#  Perform the above function for both Bacteria and Archaea for the case of isolates and MAGs.
for tuple_entry in isolates_and_mags:
    the_jgi_genome_id_function(tuple_entry[0], tuple_entry[1])
print("Done with the isolates and MAGs.")

# Now perform the same for the SAGs
for tuple_entry in sags:
    the_jgi_genome_id_function(tuple_entry[0], tuple_entry[1])
print("All JGI Genome IDs needed are now in the corresponding files.")



print("script execution stared at:", start)
end = datetime.datetime.now()
print("Script execution ended at:", end)
total_time = end - start
print("Script totally ran for :", total_time)

sys.stdout = old_stdout
log_file.close()

# ----------------------------------------------------------------------------------------------------------------------------

####     apo Vangeli - JavaScript sti consola
# for(i = 0;i < document.getElementsByTagName("td").length; i++)
# {
# if (document.getElementsByTagName("td")[i].headers == "yui-dt0-th-IMGGenomeID ") { console.log(i+":"+ document.getElementsByTagName("td")[i].innerText + "\t" + document.getElementsByTagName("td")[i-2].innerText + "\t" + document.getElementsByTagName("td")[i-3].innerText) }
# }


