#!/usr/bin/python

import numpy
import pandas
import cobra
import symengine
from cobra.flux_analysis import flux_variability_analysis

# gapsplit flux sampler
# Keaty TC & Jensen PA (2019). gapsplit: Efficient random sampling for non-convex constraint-based models.
# bioRxiv 652917; doi: https://doi.org/10.1101/652917 

def gapsplit(model, depth=500):
    fva = flux_variability_analysis(model, model.reactions, fraction_of_optimum=0.001)

    # only split reactions with feasible range >= min_range
    idxs = (fva.maximum - fva.minimum >= 1e-5).to_numpy().nonzero()[0]
    weights = (1.0 / (fva.maximum - fva.minimum) ** 2).to_numpy()
    samples = numpy.zeros((depth, len(model.reactions)))
    k = 0
    for try_ in range(1000):
        relative, target, width = _maxgap(samples[0:k,idxs], fva.iloc[idxs,:])
        primary_var = numpy.argmax(relative)
        primary_target = target[primary_var]
        primary_lb = primary_target - 0.001*width[primary_var]
        primary_ub = primary_target + 0.001*width[primary_var]

        new_sample = _generate_sample(model, idxs[primary_var], primary_lb, primary_ub)
        if new_sample is not None:
            new_sample[new_sample > fva.maximum] = fva.maximum[new_sample > fva.maximum]
            new_sample[new_sample < fva.minimum] = fva.minimum[new_sample < fva.minimum]
            samples[k,:] = new_sample
            k += 1
        if k >= depth: break

    if k < depth:
        # max_tries reached; return fewer samples
        samples = samples[:k,:]

    return pandas.DataFrame(data=samples, columns=fva.maximum.index)

def _generate_sample(model, primary_var, primary_lb, primary_ub):
    
    # Formulate a [MI]QP to find a single solution
    with model:
        model.reactions[primary_var].lower_bound = primary_lb
        model.reactions[primary_var].upper_bound = primary_ub
        model.objective = symengine.RealDouble(0)
        solution = model.optimize()
        if solution.status != 'optimal':
            return None
        else:
            return solution.fluxes

def _maxgap(points, fva):
    points = points.copy()
    points = numpy.vstack((fva.minimum, points, fva.maximum))
    points.sort(0)

    gaps = points[1:,:] - points[0:-1,:]
    width = gaps.max(0)
    loc = gaps.argmax(0)
    left = numpy.zeros(width.size)
    for i in range(width.size):
        left[i] = points[loc[i],i]

    relative = width / (points[-1,:] - points[0,:])
    target = left + width / 2.0

    return relative, target, width

