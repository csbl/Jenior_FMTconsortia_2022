#!/usr/bin/python

import pickle


pathwayDict = pickle.load(open('/home/mjenior/Desktop/active_projects/CDI_FMT_project/fmt_metaG/gene_to_pathway.pickle', 'rb'))

outFile = open('compiled_mapping.pathways.tsv', 'w')
with open('compiled_mapping.tsv', 'r') as inFile:
	header = inFile.readline()
	header = header.strip() + '\tpathway\n'
	outFile.write(header)

	for line in inFile:

		try:
			pathway = pathwayDict[line.split()[0]]
		except KeyError:
			pathway = 'unknown'

		entry = line.strip() + '\t' + pathway + '\n'
		outFile.write(entry)


outFile.close()

