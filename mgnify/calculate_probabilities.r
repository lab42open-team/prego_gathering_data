#!/usr/bin/Rscript

# Dependencies
library(Rfast)

######################################################
### Regarding the NCBI IDs and their probabilities
######################################################

# read the background file
x<- read.csv("/data/databases/scripts/gathering_data/mgnify/ncbi_background.tsv", header = FALSE, sep = "\t")

# keep the values of each id occurence
mydata_x <- x$V2

# build a model using the negative binomial distribution
mod_x <- negbin.mle(mydata_x)

# calculate the probability of each ncbi id to come with respect to the non binomial distribution we built
x$new_x <- dnbinom(x = mydata_x, size = mod_x$param[2], prob = mod_x$param[1])

# write an output file with the probabilities of each NCBI ID to come 
write.table(x, file='/data/databases/scripts/gathering_data/mgnify/probabilities_ncbi.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)


########################################################
### Regarding the metadata terms and their probabilities
########################################################

# read the background file
y <- read.csv("/data/databases/scripts/gathering_data/mgnify/metadata_background.tsv", header = FALSE, sep = "\t")

# keep the values of each id occurence
mydata_y <- y$V2

# build a model using the negative binomial distribution
mod_y <- negbin.mle(mydata_y)

# calculate the probability of each ncbi id to come with respect to the non binomial distribution we built
y$new_y <- dnbinom(x = mydata_y, size = mod_y$param[2], prob = mod_y$param[1])

# write an output file with the probabilities of each metadata ter to come 
write.table(y, file='/data/databases/scripts/gathering_data/mgnify/probabilities_metadata.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)
