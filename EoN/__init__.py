r'''
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



EoN consists of two sets of algorithms.  

- The first deals with simulation of epidemics on networks.  The most significant of these are `fast_SIS` and `fast_SIR` which significantly outperform Gillespie algorithms (also included).  These algorithms are discussed in more detail in the appendix of the book.


- The second deals with solution of systems of equations derived in the book.  For these it is possible to either provide the degree distribution, or simply use a network and let the code determine the degree distribution.


- There are a few additional algorithms which are not described in the book, but which we believe will be useful. Most notably, the function `visualize` which creates a sequence of images of a network which are appropriate for creating a movie showing disease spread.

Distributed under MIT license.  See :download:`license.txt<../license.txt>` for full details.



'''

print("warning - EoN is currently under significant development.  Interface"
      +" may change with little if any warning until version 1.0")
print("Warning - testing in Python 3 is limited")

__author__ = "Joel C. Miller, Istvan Z. Kiss, and Peter Simon"
__version__ = "0.93"

#__all__ = 

class EoNError(Exception):
    r'''
    this will be the basic error type for EoN
    '''
    pass

def _get_rate_functions(G, tau, gamma, transmission_weight = None, 
                        recovery_weight=None):
    r'''
    Arguments:
        G : networkx Graph
            the graph disease spread on

        tau : number
            disease parameter giving edge transmission rate (subject to edge scaling)

        gamma : number (default None)
            disease parameter giving typical recovery rate, 
        
        transmission_weight : string (default None)
            `G.edge[u][v][transmission_weight]` scales up or down the recovery rate.

        recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                `gamma_i = G.node[i][recovery_weight]*gamma`
    Returns:
        : trans_rate_fxn, rec_rate_fxn
            Two functions such that 
            - `trans_rate_fxn(u,v)` is the transmission rate from u to v and
            - `rec_rate_fxn(u)` is the recovery rate of u.
'''
    if transmission_weight is None:
        trans_rate_fxn = lambda x, y: tau
    else:
        trans_rate_fxn = lambda x, y: tau*G.edge[x][y][transmission_weight]

    if recovery_weight is None:
        rec_rate_fxn = lambda x : gamma
    else:
        rec_rate_fxn = lambda x : gamma*G.node[x][recovery_weight]


    return trans_rate_fxn, rec_rate_fxn



import EoN.simulation
from EoN.simulation import *
import EoN.analytic
from EoN.analytic import *



def subsample(report_times, times, status1, status2=None, 
                status3 = None):
    r'''
    Given 
      S, I, and/or R as lists of numbers of nodes of the given status
      at given times

    returns them 
      subsampled at specific report_times.
    

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

    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        import scipy
        import matplotlib.pyplot as plt

        """ in this example we will run 100 stochastic simulations.
            Each simulation will produce output at a different set
            of times.  In order to calculate an average we will use
            subsample to find the epidemic sizes at a specific set
            of times given by report_times.
        """

        G = nx.fast_gnp_random_graph(10000,0.001)
        tau = 1.
        gamma = 1.
        report_times = scipy.linspace(0,5,101)
        Ssum = scipy.zeros(len(report_times))
        Isum = scipy.zeros(len(report_times))
        Rsum = scipy.zeros(len(report_times))
        iterations = 100
        for counter in range(iterations): 
            t, S, I, R = EoN.fast_SIS(G, tau, gamma, initial_infecteds = range(10))
            #t, S, I, and R have an entry for every single event.
            newS, newI, newR = EoN.subsample(report_times, t, S, I, R)
            plt.plot(t, newI, linewidth=1, alpha = 0.4)
            Ssum += newS
            Isum += newI
            Rsum += newR
        Save = Ssum / float(iterations)
        Iave = Isum / float(iterations)
        Rave = Rsum / float(iterations)
        plt.plot(report_times, Save, "--", linewidth = 3, label = "average")
        plt.plot(report_times, Iave, "--", linewidth = 3)
        plt.plot(report_times, Rave, "--", linewidth = 3)
        plt.legend(loc = "upper right")
        plt.savefig("tmp.pdf")

    If only one of the sample times is given then returns just that.

    If report_times goes longer than times, then this simply assumes the 
    system freezes in the final state.
    
    This uses a recursive approach if multiple arguments are defined.


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
    
    Arguments:
        times : list or scipy array (ordered)
            the times we have observations
        L : a list or scipy array
            order of L corresponds to times
        threshold : number
            a threshold value

    Returns:
        : 
        t : number
            the first time at which L reaches or exceeds a threshold.

    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN
        import matplotlib.pyplot as plt

        """ in this example we will run 20 stochastic simulations.
            We will produce one plot showing the unshifted
            curves and one with them shifted so that t=0 when 1%
            are in the I class.
        """
        N=1000000
        kave = 10.
        G = nx.fast_gnp_random_graph(100000,kave/(N-1.))
        tau = 1.
        gamma = 1.
        report_times = scipy.linspace(0,5,101)
        Ssum = scipy.zeros(len(report_times))
        Isum = scipy.zeros(len(report_times))
        Rsum = scipy.zeros(len(report_times))
        iterations = 20
        for counter in range(iterations):
            R=[0]
            while R[-1]<1000: #if an epidemic doesn't happen, repeat
                t, S, I, R = EoN.fast_SIS(G, tau, gamma, initial_infecteds = range(10))
            plt.plot(t, I, linewidth = 1, color = 'gray', alpha=0.4)
            tshift = EoN.get_time_shift(t, I, 0.01*kave)
            plt.plot(t-tshift, I, color = 'red', linewidth = 1, alpha = 0.4)
        plt.savefig("tmp.pdf")
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
'''

