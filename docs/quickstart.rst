QuickStart Guide
================


>>>import networkx as nx
>>>import matplotlib.pyplot as plt
>>>import EoN
>>>
>>>G=nx.barabasi_albert_graph(10**6, 5) #create a barabasi-albert graph
>>>
>>>iterations = 10  #run 10 simulations
>>>tau = 0.1           #transmission rate
>>>gamma = 1.0    #recovery rate
>>>rho = 0.005      #random fraction initially infected
>>>
>>>for counter in range(iterations):
>>>    t, S, I, R = EoN.fast_SIR(G, tau, gamma, rho=rho, tmax = 20)
>>>    plt.plot(t, I)
>>>
>>>t, S, I, R = EoN.EBCM_from_graph(G, tau, gamma, rho=rho)
>>>plt.plot(t, I, ':', label = 'EBCM approximation')
>>>
>>>plt.legend(loc = 'upper right')
>>>plt.savefig('BA_EBCM_vs_sim.pdf')
