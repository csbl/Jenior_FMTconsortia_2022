#!/usr/bin/python

import pickle

title_dict = {}
with open('/media/mjenior/Jenior_HD/kegg/pathway/map_title.tab', 'r') as inFile:
	for line in inFile:
		line = line.split()
		entry = str(line[0])
		path = '_'.join(line[1:])
		try:
			title_dict[entry] += ':' + path
			print('found multi')
		except KeyError:
			title_dict[entry] = path

path_dict = {}
with open('/media/mjenior/Jenior_HD/kegg/genes/links/genes_pathway.list', 'r') as inFile:
	for line in inFile:
		line = line.split()
		gene = line[0]
		path = title_dict[str(line[1][-5:])]
		try:
			path_dict[entry] += '|' + path
		except KeyError:
			path_dict[gene] = path




pickle.dump(path_dict, open('gene_to_pathway.pickle', 'wb'))

