#!/usr/bin/python
'''USAGE: makeTax.py
Creates KEGG lineage and Gram stain information dictionary
'''
import pickle

gramDict = {'Proteobacteria':'negative', 'Tenericutes':'negative', 'Bacteroidetes':'negative','Fusobacteria':'negative','Verrucomicrobia':'negative',
'Actinobacteria':'positive', 'Spirochaetes':'negative', 'Chlamydiae':'negative', 'Thermotogae':'negative', 'Planctomycetes':'negative',
'Firmicutes':'positive', 'Cyanobacteria':'negative', 'Chloroflexi':'negative', 'Chlorobi':'negative', 'Deinococcus-Thermus':'positive'}

taxDict = {}
with open('/home/matt/Desktop/kegg_genome.txt', 'r', errors='ignore') as kegg:

	for line in kegg:
		line = line.split()
		
		if line[0] != 'LINEAGE':
			continue
		elif len(line) < 3:
			continue
		else:
			# LINEAGE   Bacteria; Proteobacteria; Gammaproteobacteria; Pasteurellales; Pasteurellaceae; Haemophilus
			if 'unclassified' == line[-2]:
				genus = line[-3]
				tax = ''.join(line[1:-2])
			else:
				genus = line[-1]
				tax = ''.join(line[1:])
			phylum = line[2].replace(';','')
			if not phylum in list(gramDict.keys()):
				continue
			else:
				entry = tax + '_' + gramDict[phylum]
				taxDict[genus] = entry

with open('/home/matt/Desktop/kegg_genome.pickle', 'wb') as outFile:
    pickle.dump(taxDict, outFile)
