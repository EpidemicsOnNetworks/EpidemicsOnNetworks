QuickStart Guide
================

The code here provides an example of creating a Barabasi-Albert network.  Then it performs several simulations of an SIR epidemic starting with a fraction rho randomly infected initially.  Finally it uses several analytic models to predict the spread of an epidemic in a random network with the given properties.

>>> import networkx as nx
>>> import matplotlib.pyplot as plt
>>> import EoN
>>> import random
>>> 
>>> N=10**5
>>> G=nx.barabasi_albert_graph(N, 5) #create a barabasi-albert graph
>>> 
>>> tmax = 20
>>> iterations = 5  #run 5 simulations
>>> tau = 0.1           #transmission rate
>>> gamma = 1.0    #recovery rate
>>> rho = 0.005      #random fraction initially infected
>>> 
>>> for counter in range(iterations): #run simulations
>>>     initial_infections = random.sample(G.nodes(), int(round(rho*N))) 
>>>     t, S, I, R = EoN.fast_SIR(G, tau, gamma, initial_infecteds=initial_infections, tmax = tmax)
>>>     plt.plot(t, I, color = 'k', alpha=0.3)
>>>         
>>> t, S, I, R = EoN.SIR_homogeneous_pairwise_from_graph(G, tau, gamma, rho=rho, tmax = tmax)
>>> plt.plot(t, I, '-.', label = 'Homogeneous pairwise', linewidth = 5)
>>> 
>>> t, S, I, R = EoN.SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, rho=rho, tmax=tmax)
>>> plt.plot(t, I, ':', label = 'Heterogeneous meanfield', linewidth = 5)
>>> 
>>> t, S, I, R = EoN.EBCM_from_graph(G, tau, gamma, rho=rho, tmax = tmax)
>>> plt.plot(t, I, '--', label = 'EBCM approximation', linewidth = 5)
>>> 
>>> plt.legend(loc = 'upper right')
>>> plt.savefig('SIR_BA_model_vs_sim.pdf')
>>> 
>>> plt.clf()
>>>
>>> #Now run for SIS.   Simulation is much slower so need smaller network
>>> N=10**4  
>>> G=nx.barabasi_albert_graph(N, 5) #create a barabasi-albert graph
>>> print 'got graph'
>>> for counter in range(iterations):
>>>     initial_infections = random.sample(G.nodes(), int(round(rho*N))) 
>>>     t, S, I = EoN.fast_SIS(G, tau, gamma, initial_infecteds=initial_infections, tmax = tmax)
>>>     plt.plot(t, I, color = 'k', alpha=0.3)
>>>         
>>> t, S, I = EoN.SIS_homogeneous_pairwise_from_graph(G, tau, gamma, rho=rho, tmax = tmax)
>>> plt.plot(t, I, '-.', label = 'Homogeneous pairwise', linewidth = 5)
>>> 
>>> t, S, I = EoN.SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, rho=rho, tmax=tmax)
>>> plt.plot(t, I, ':', label = 'Heterogeneous meanfield', linewidth = 5)
>>> 
>>> t, S, I = EoN.SIS_compact_pairwise_from_graph(G, tau, gamma, rho=rho, tmax=tmax)
>>> plt.plot(t, I, '--', label = 'Compact pairwise', linewidth = 5)
>>> 
>>> plt.legend(loc = 'lower right')
>>> plt.savefig('SIS_BA_model_vs_sim.pdf')

