#!/usr/bin/python
'''USAGE: addTaxData.py file_list
Adds bacterial taxonomy information to bin annotation files
'''
import sys
import pickle

with open('/home/matt/Desktop/kegg_genome.pickle', 'rb') as f: 
	taxDict = pickle.load(f)

outFile = str(sys.argv[1]).rstrip('txt') + 'tax_gram.txt'
outFile = open(outFile,'w')
with open(sys.argv[1],'r') as annotation:

	# fmt.metaG.01044A.bin.149.KEGGprot.out | Blautia obeum | 3296 of 5051 (65.3%)
	for line in annotation:
		species = line.split(' | ')[1]
		species = species.replace(' ','_')
		line = line.split()
		filename = line[0]
		genus = line[2]
		try:
			taxonomy = taxDict[genus]
		except:
			continue
		taxonomy = taxDict[genus].split('_')
		entry = filename + '\t' + taxonomy[0] + '\t' + taxonomy[1] + '\t' + species + '\n'
		outFile.write(entry)

outFile.close()
