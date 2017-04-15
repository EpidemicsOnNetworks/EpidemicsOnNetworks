r'''
EoN
===

EoN (Epidemics on Networks)

EoN is a Python package for the simulation of epidemics on networks 
and ODE models of disease spread.

The algorithms are based on the book
        
`Mathematics of epidemics on networks: from exact to approximate 
models`
by Kiss, Miller & Simon
        http://www.springer.com/book/9783319508047
        
Please cite the book if using these algorithms

For simulations, we assume that input networks are **NetworkX** 
graphs; see https://networkx.github.io/




Distributed under MIT license.  See license.txt for full details.

'''

print("warning - EoN is not yet tested in Python 3")

print("warning - EoN is currently under significant development.  Interface"
      +" may change with little if any warning until it version 1.0")


__author__ = "Joel C. Miller, Istvan Z. Kiss, and Peter Simon"
__version__ = "0.92"



class EoNError(Exception):
    r'''
    this will be the basic error type for EoN
    '''
    pass


import simulation
from simulation import *
import analytic
from analytic import *



def subsample(report_times, times, status1, status2=None, 
                status3 = None):
    r'''
    Given 
      S, I, and/or R as lists (or other iterable) of numbers of nodes of
      given status at given times
    returns them 
      subsampled at specific report_times.
      
    If more than one argument is given, does so as a list in order given, but 
    skipping whichever was not included (if any not included)

    If only one is given then returns just that.

    If report_times goes longer than times, then this simply assumes the 
    system freezes in the final state.
    
    This uses a recursive approach if multiple arguments are defined.

    :INPUTS:

    report_times : iterable (ordered)
                   times at which we want to know state of system
                   
    times : iterable (ordered)
            times at which we have the system state (assumed no change 
            between these times)
            
    statusX (X one of 1, 2 or 3) : iterable (order corresponds to times)
                          generally S, I, or R
                          number of nodes in given status.
    :RETURNS:

    If only status1 is defined
        report_status1 : scipy array gives status1 subsampled just at 
                         report_times.
                     
    If more are defined then it returns a list, either
        [report_status1, report_status2]
    or
        [report_status1, report_status2, report_status3]
        
    '''
    if report_times[0] < times[0]:
        raise EoNError("report_times[0]<times[0]")
        
    report_status1 = []
    next_report_index = 0
    next_observation_index = 0
    while next_report_index < len(report_times):
        while next_observation_index < len(times) and \
              times[next_observation_index]<= report_times[next_report_index]:
            candidate = status1[next_observation_index]
            next_observation_index += 1
        report_status1.append(candidate)
        next_report_index +=1
        
    report_status1= scipy.array(report_status1)
    
    if status2 is not None:
        if status3 is not None:
            report_status2, report_status3 = subsample(report_times, times, status2, status3)
            return report_status1, report_status2, report_status3
        else:
            report_status2 = subsample(report_times, times, status2)
            return report_status1, report_status2
    else:
        return report_status1



def get_time_shift(times, L, threshold):
    r'''
    Identifies the first time at which L crosses a threshold.  
    Useful for shifting times.
    
    :INPUTS:

    times : list or scipy array (ordered)
            the times we have observations
    L : a list or scipy array
        order of L corresponds to times
    threshold : number
        a threshold value

    :RETURNS:

    t : number
        the first time at which L reaches or exceeds a threshold.
    '''
    for index, t in enumerate(times):
        if L[index]>= threshold:
            break
    return t




'''
These are the systems I want to include based on their numbering in the 
book:

coded (3.7) SIS individual based
(3.30) SIR individual based
NOT coded (3.26) SIS pair based
(3.39) SIR pair based

chapter 4?

(5.13) SIS heterogeneous pairwise
(5.15) SIR heterogeneous pairwise
(5.18) SIS compact pairwise
(5.19) SIR compact pairwise
(5.20) SIS super compact pairwise
(5.22) SIR super compact pairwise
(5.36) SIS effective degree
(5.38) SIR effective degree
(5.43) SIR compact effective degree
(5.44) SIS compact effective degree = SIS compact pairwise

(6.2) Epidemic probability discrete time
(6.3) Epidemic probability continuous time
(6.5) Epidemic probability non-Markovian
(6.6) Epidemic size discrete time
(6.7) Epidemic size continuous time
(6.8) Epidemic size non-Markovian
(6.10) Epidemic size discrete time (large IC)
(6.11) Discrete-time EBCM model
(6.12) Continuous time EBCM model

(8.1) SIS pairwise contact conserving rewiring
(8.5) SIS eff. deg. contact conserving rewiring
(8.7) SIS pairwise random activation/deletion
(8.13) SIS eff. deg. random activation/deletion
(8.15) SIS pairwise link-status dependent act/del
(8.16) SIS link deactivation-activation on fixed networks.
(8.19) EBCM dynamic network

(9.5) SI^{K}R multistage pairwise for homogeneous
(9.27) SIR pairwise, constant infection duration.
(9.35) SIR homogeneous pairwise, general recovery
(9.36) SIR EBCM non-Markovian trans/recovery

add models that take in graph, measure degree distribution and run EBCM
similarly for EBCM with neighbor degrees (see barabasi_SIR.py)

consider explicitly defining toast graph etc.

SIR_pair_based --- which of 2 versions to keep, and comments need to 
                   explain it a bit better.
'''

