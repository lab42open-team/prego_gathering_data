#!/usr/bin/env python

import pg
import sys

publication_bioprojects = {}
for line in open('/data/databases/bioproject/bioproject_unicellular_pubmed_ncbi.tsv', 'r'):
    fields = line.strip().split('\t')
    if len(fields) > 1:
        bioproject = fields[0]
        for publication in fields[1].split(','):
            if publication not in publication_bioprojects:
                publication_bioprojects[publication] = set()
            publication_bioprojects[publication].add(bioproject)

db = pg.connect(host='localhost', port=5432, user='guest', dbname='textmining')
bioproject_entity_publications = {}
previous = []
#for match in db.query('SELECT * FROM matches WHERE document in (%s) AND type in (-21,-23,-25,-26,-27);' % ','.join(publication_bioprojects.keys())).getresult():
for match in db.query('SELECT * FROM matches WHERE document in (%s) AND type in (-21,-23,-25,-26,-27);' % '9738079').getresult():
    print match
    current = match[0:2]
    if current == previous:
        continue
    previous = current
    publication = str(match[0])
    entity = (match[3], match[4])
    for bioproject in publication_bioprojects[publication]:
        if bioproject not in bioproject_entity_publications:
            bioproject_entity_publications[bioproject] = {}
        if entity not in bioproject_entity_publications[bioproject]:
            bioproject_entity_publications[bioproject][entity] = set()
        bioproject_entity_publications[bioproject][entity].add(publication)
        print bioproject_entity_publications[bioproject][entity]
        sys.exit(0)

with open('bioproject_metadata.tsv', 'r') as bioproject: 
    next(bioproject)
    for line in img:
        fields = line.strip().split('\t')
        bioproject = fields[10]
        if bioproject != '' and bioproject != '-' and bioproject != '0':
            if bioproject in bioproject_entity_publications:
                for entity in bioproject_entity_publications[bioproject]:
                    publications = sorted(bioproject_entity_publications[bioproject][entity])
                    evidence = ' '.join(list(map(lambda x: 'PMID:'+str(x), publications)))
                    print '-2\t%d\t%d\t%s\tBioProject\t%s\t3\tTRUE\t' % (ncbi_taxon, entity[0], entity[1], evidence)
