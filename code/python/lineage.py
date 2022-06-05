#!/usr/bin/python
# Rewrites kegg_genome information file to just lineages excluding Archaea 

outFile = open('/home/matt/Desktop/kegg_lineage.txt','w')
with open('/home/matt/Desktop/kegg_genome.txt','r') as kegg:

	for line in kegg:
		line = line.split()
		print(line)
		
		if line[0] != 'LINEAGE':
			continue
		elif line[1] == 'Archaea':
			continue
		else:
			tax = ' '.join(line[1:])
			outFile.write(tax)

outFile.close()
