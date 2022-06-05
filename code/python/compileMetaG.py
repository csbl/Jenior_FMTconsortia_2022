#!/usr/bin/python
'''
USAGE: python compileMetaG.py mapping_files.txt compiled_mapping.tsv
'''

import sys

outFile = open(sys.argv[2], 'w')
outFile.write('keggID\tgene')

inFiles = []
with open(sys.argv[1], 'r') as mapping:
	for line in mapping:
		outFile.write('\t')
		newFile = line.strip()
		outFile.write(newFile)
		inFiles.append(newFile)
outFile.write('\n')

mappingDict = {}
all_genes = set()
for inFile in inFiles:
	locals()[inFile + '.dict'] = {}
	with open(inFile, 'r') as current:
		header = current.readline()
		for line in current:
			if 'unknown:' in line:
				continue
			elif 'mgy:MGMSRv2' in line:
				gene = '__'.join(line.split()[0].split('|')[0:2]) + '|'
				gene += line.split()[0].split('|')[2]
				abund = int(line.split()[-1])
			else:
				gene = '|'.join(line.split()[0].split('|')[0:2])
				abund = int(line.split()[-1])
				locals()[inFile + '.dict'][gene] = abund
				all_genes |= set([gene])


for gene in all_genes:
	kegg = gene.split('|')[0]
	if len(kegg) > 50: continue
	gene_name = gene.split('|')[1]



	if len(gene_name) > 133:
		gene_name = gene_name[0:134]
		if len(gene_name) > 133:
			gene_name = gene_name.split('\\')[0]
			if len(gene_name) > 133:
				gene_name = gene_name.split('/')[0]


	outFile.write(kegg + '\t' + gene_name)

	for inFile in inFiles:
		outFile.write('\t')

		try:
			abund = locals()[inFile + '.dict'][gene]
		except KeyError:
			abund = 0
		outFile.write(str(abund))

	outFile.write('\n')


outFile.close()





