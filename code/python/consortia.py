#!/usr/bin/python

# Dependencies
import cobra
import warnings
import symengine
from sys import argv
from random import shuffle
from copy import deepcopy
from cobra.util import solver
from cobra.manipulation.delete import *


# Sets reaction bounds based on the species abundance within the consortia
def constrain_by_abundance():
	pass


# Function to add reactions to consortia super-organism
def grow_genre(community, member):
	add_rxns = []
	community_rxns = set([x.id for x in community.reactions])
	for rxn in member.reactions:
		if rxn.id in community_rxns:
			continue
		elif rxn.id in c('biomass_GmPos','biomass_GmNeg'):
			continue
		else:
			add_rxns.append(deepcopy(rxn))

	community.add_reactions(add_rxns)

	return community


# Get list of GENRE files to include
genre_files = []
with open(argv[1], 'r') as genres:
	for line in genres:
		genre_files.append(line.strip())


# Read in GENREs and add to new growing consortia model
community_model = cobra.Model('community_model')
for genre in genre_files:
	genre = cobra.io.read_sbml_model(genre)
	community_model = grow_genre(community_model, genre)


# Add community objective
biomass_rxn = cobra.Reaction('biomass_bacterial')
biomass_rxn.name = 'Biomass Reaction'
biomass_rxn.lower_bound = 0.
biomass_rxn.upper_bound = 1000.
biomass_rxn.add_metabolites({
    cpd17041_c: -0.4,  # Protein
    cpd17043_c: -0.15, # RNA
    cpd17042_c: -0.05, # DNA
    cpd11852_c: -0.05, # Lipid
    cellwall_c: -0.2,
    cofactor_c: -0.2,
    cpd00001_c: -20.0,
    cpd00002_c: -20.0,
    cpd00008_c: 20.0,
    cpd00009_c: 20.0,
    cpd11416_c: 1.0 # Biomass
})
community_model.add_reactions([biomass_rxn])
community_model.objective = biomass_rxn
test_flux = round(community_model.slim_optimize(), 4)


# Check stats for new model
print('Consortia GENRE has', len(community_model.genes), 'genes,', len(community_model.reactions), 'reactions, and', len(community_model.metabolites), 'metabolites')
print('Biomass flux:', test_flux)


