all: all.tsv

jgiimg.tsv: jgiimg_metadata.tsv bioproject_publications.tsv
	./jgiimg.py > $@

all.tsv: jgiimg.tsv
	cat $^ > $@
