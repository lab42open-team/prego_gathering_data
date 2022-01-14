# PREGO gathering data scripts

[![License](https://img.shields.io/badge/License-BSD_2--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)

In this repository are all the scripts for [PREGO](https://prego.hcmr.gr/Search) gathering data module.
Scripts are organized based on the resources, each resource has its' own collection of scripts.
Included are functions for API calls, FTP access urls and association extraction.

## Structure

Currently, 24/12/2021, there are five supported resousces, 
[BioProject](https://www.ncbi.nlm.nih.gov/bioproject/), [JGI IMG](https://img.jgi.doe.gov), [MG-RAST](https://www.mg-rast.org), 
[MGnify](https://www.ebi.ac.uk/metagenomics/) and [Struo pipeline](https://github.com/leylabmpi/Struo).


Each resource has scripts with the prefix `*get*`, `*extract*` and `*functions*`.  
`*get*` contains the API calls and transformations of data and 
`*extract*` scripts perform the Named Entity Recognition using the EXTRACT tagger, 
find the associations and then calculate the score for each association.
