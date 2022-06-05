#!/usr/bin/python

import cobra
from copy import deepcopy
from cobra.util import solver
from random import shuffle
from cobra.medium import minimal_medium
from math import ceil
from cobra.flux_analysis import flux_variability_analysis as fva


def _minimize_exchanges_coop(model, fraction, avail_lb=-1000.):

    newModel = deepcopy(model)
    newModel = _minimum_contraint_coop(newModel, fraction)
    for rxn in newModel.exchanges: rxn.bounds = (-1000., 1000.)

    minObj = float(ceil(newModel.slim_optimize() * fraction))
    minmedia = set(list(minimal_medium(newModel, min_objective_value=minObj, minimize_components=True).index))

    for rxn in newModel.exchanges:
        if rxn.id in minmedia:
            rxn.lower_bound = float(avail_lb)
        else:
            rxn.lower_bound = 0.

    mingrowth = newModel.slim_optimize()

    return newModel, minmedia, mingrowth

def _minimum_contraint_coop(model, fraction):

    mingrowth = model.slim_optimize() * fraction
    objConstr = model.problem.Constraint(model.objective.expression, lb=mingrowth, ub=1000.)
    model.add_cons_vars([objConstr])
    model.solver.update()

    return model

def _fva_produce(model):
    
    model_exch = [rxn.id for rxn in model.exchanges]
    exch_fva = fva(model, reaction_list=model_exch)
    
    consumed = []
    for exch, row in exch_fva.iterrows():
        if row['maximum'] > 1e-5: consumed.append(exch)
            
    return set(consumed)


# Iterative calculations of possible metabolite sharing between GENREs
def mutualism(raw_model1, raw_model2, iters=5, exch_dict=False, fraction=0.25, avail_lb=-1000., step=100):
    
    model1, minmedia1, mingrowth1 = _minimize_exchanges_coop(raw_model1, fraction, avail_lb)
    model2, minmedia2, mingrowth2 = _minimize_exchanges_coop(raw_model2, fraction, avail_lb)
    model1_exchs = set([rxn.id for rxn in model1.exchanges])
    model2_exchs = set([rxn.id for rxn in model2.exchanges])
    minmedia = minmedia1.union(minmedia2)

    newMedia = {}
    for exch in minmedia:
        newMedia[exch] = [0, 'base_medium']
        try: 
            model1.reactions.get_by_id(exch).lower_bound = float(avail_lb)
        except:
            pass
        try:
            model2.reactions.get_by_id(exch).lower_bound = float(avail_lb)
        except:
            pass

    one_to_two = 0
    two_to_one = 0
    for i in range(1, iters+1):
        current_media = set(list(newMedia.keys()))
        additions1 = _fva_produce(model1)
        additions2 = _fva_produce(model2)
        curr_additions = additions1.union(additions2)
        
        for rxn in additions1:
            newMedia[rxn] = [i]
            newMedia[rxn].append('1-->2')
            one_to_two += 1
            try:
                if model2.reactions.get_by_id(rxn).lower_bound > -1000.:
                    model2.reactions.get_by_id(rxn).lower_bound -= float(step)
            except:
                continue
        for rxn in additions2:
            newMedia[rxn] = [i]
            newMedia[rxn].append('2-->1')
            two_to_one += 1
            try:
                if model1.reactions.get_by_id(rxn).lower_bound > -1000.:
                    model1.reactions.get_by_id(rxn).lower_bound -= float(step)
            except:
                continue

        test = len(newMedia.keys()) - len(current_media)
        if test == 0:
            print('Failed to find any new metabolites in iteration', i)
            break
        else:
            print('Added', test,'metabolites in iteration', i)

    totalNew = len(newMedia.keys()) - len(minmedia)
    print('A total of', totalNew, 'new metabolites added to media')
    print('\t', one_to_two, 'metabolites passed from GENRE 1 to GENRE 2')
    print('\t', two_to_one, 'metabolites passed from GENRE 2 to GENRE 1')
    newgrowth1 = model1.slim_optimize()
    newgrowth1 = newgrowth1 - mingrowth1
    print('GENRE 1 objective flux increased by', round(newgrowth1,4))
    newgrowth2 = model2.slim_optimize()
    newgrowth2 = newgrowth2 - mingrowth2
    print('GENRE 2 objective flux increased by', round(newgrowth2,4))
    
    ids = []
    names = []
    iterations = []
    crossfed = []
    for exch in newMedia.keys():
        ids.append(exch)
        iterations.append(newMedia[exch][0])
        crossfed.append(newMedia[exch][1])
        try:
            names.append(list(raw_model1.reactions.get_by_id(exch).metabolites)[0].name)
        except:
            names.append(list(raw_model2.reactions.get_by_id(exch).metabolites)[0].name)
        
    cooperation = {'exchange':  ids, 
                  'metabolite': names, 
                  'iteration': iterations, 
                  'direction': crossfed}
    cooperation = pandas.DataFrame(cooperation, columns = ['exchange','metabolite','iteration','direction'])
        
    return cooperation


def _minimum_contraint_comp(model, fraction):
    
    newModel = deepcopy(model)
    for rxn in newModel.exchanges: rxn.bounds = (-1000., 1000.)
    mingrowth = newModel.slim_optimize() * fraction
    objConstr = model.problem.Constraint(newModel.objective.expression, lb=mingrowth, ub=1000.)
    newModel.add_cons_vars([objConstr])
    newModel.solver.update()

    return newModel

def _fva_consume(model):
    
    model_exch = [rxn.id for rxn in model.exchanges]
    exch_fva = fva(model, reaction_list=model_exch)
    
    consumed = []
    for exch, row in exch_fva.iterrows():
        if row['minimum'] < -1.: consumed.append(exch) # Needs to be larger, trace is not relevant
            
    return set(consumed)


# Iterative calculation of potential substrate competition between GENREs
def competition(model1, model2, fraction=0.25, avail_lb=-1000., iters=9, step=0.5):
    newModel1 = _minimum_contraint_comp(model1, fraction)
    newModel2 = _minimum_contraint_comp(model2, fraction)
    
    # Set availability for edges of competition
    shared_minMedia = set([x.id for x in newModel1.exchanges]).intersection(set([y.id for y in newModel2.exchanges]))
    for exch in newModel1.reactions:
        if exch.id in shared_minMedia: exch.lower_bound = float(avail_lb)
    for exch in newModel2.reactions:
        if exch.id in shared_minMedia: exch.lower_bound = float(avail_lb)
    
    pre_flux1 = round(newModel1.slim_optimize(), 4)
    pre_flux2 = round(newModel2.slim_optimize(), 4)
    print('Starting flux', pre_flux1, pre_flux2)
    competition = {}
    for i in range(1, iters+1):
        consumed1 = _fva_consume(newModel1)
        consumed2 = _fva_consume(newModel2)
        contested = consumed1.intersection(consumed2)  
        for x in contested:
            newModel1.reactions.get_by_id(x).lower_bound *= step
            newModel2.reactions.get_by_id(x).lower_bound *= step
            try:
                competition[x][1] += 1
            except:
                competition[x] = [i, 1]
            
        model1_current = round(newModel1.slim_optimize(), 4)
        model2_current = round(newModel2.slim_optimize(), 4)
        print('Iteration', i, '|', len(contested), 'contested metabolites |', model1_current, model2_current)
        if str(model1_current) == 'nan': 
            break
        elif str(model2_current) == 'nan':
            break
        elif pre_flux1 == model1_current and pre_flux2 == model2_current:
            break
        else:
            pre_flux1 = model1_current
            pre_flux2 = model2_current
    
    ids = []
    names = []
    iterations = []
    contested = []
    for exch in competition.keys():
        ids.append(exch)
        iterations.append(competition[exch][0])
        contested.append(competition[exch][1])
        try:
            names.append(list(model1.reactions.get_by_id(exch).metabolites)[0].name)
        except:
            names.append(list(model2.reactions.get_by_id(exch).metabolites)[0].name)
        
    competition = {'exchange':  ids, 
                  'metabolite': names, 
                  'iteration': iterations, 
                  'times_contested': contested}
    competition = pandas.DataFrame(competition, columns = ['exchange','metabolite','iteration','times_contested'])
            
    
    return competition

#---------------------------------------------------------------#

# Under construction
# Increase number by a defined fraction of intracellular metabolites with transports and extracellular exchanges
def expand_transporters(model, min_rxns=3, fraction = 0):    
    transporters = []
    exchs = set([x.id for x in model.exchanges])
    for rxn in model.reactions:
        if rxn.id in exchs: 
            continue
        else:
            compartments = set()
            for x in list(rxn.metabolites): compartments |= set([x.compartment])
            if len(compartments) > 1: transporters.append(rxn.id)
    transported = set()
    for x in transporters:
        transported |= set([y.id for y in model.reactions.get_by_id(x).metabolites])
    
    new_exchs = []
    new_rxns = []
    curr_new_rxn = 0
 
    for cyt_cpd in model.metabolites:
        if cyt_cpd.id in transported:
            continue
        elif len(cyt_cpd.reactions) < min_rxns:
            continue
        else:
            ext_id = cyt_cpd.id.split('_')[0] + '_e'
            try:
                ext_cpd = model.metabolites.get_by_id(ext_id)
            except:
                ext_cpd = deepcopy(cyt_cpd)
                ext_cpd.compartment = 'extracellular'
                ext_cpd.id = ext_id
                model.add_boundary(ext_cpd, type='exchange', reaction_id='EX_'+ext_cpd.id, lb=0., ub=1000.)
                new_exchs.append('EX_'+ext_cpd.id)
                
            # Add new transporters
            curr_new_rxn += 1
            new_rxn = cobra.Reaction('transporter_' + str(curr_new_rxn))
            new_rxn.name = ext_cpd.name + ' diffusion'
            new_rxn.lower_bound = -1000.
            new_rxn.upper_bound = 1000.
            new_rxn.add_metabolites({ext_cpd: -1.0, cyt_cpd: 1.0})
            new_rxns.append(new_rxn)
    
    # Open indicated fraction of new exchanges
    if fraction != 0:
        shuffle(new_exchs)
        perc_changed = round(len(new_exchs)*fraction)
        for x in new_exchs[0:perc_changed]: model.reactions.get_by_id(x).lower_bound = -1000.
    model.add_reactions(new_rxns)
    
    total_added = len(new_rxns) + len(new_exchs)
    print('Total added reactions:', total_added)
    if fraction != 0: print('Opened exchanges:', perc_changed)
    
    return model


