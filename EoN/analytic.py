from scipy import integrate
from scipy.ndimage.interpolation import shift
import scipy
import networkx as nx

#######################
#                     #
#   Auxiliary stuff   #
#                     #
#######################



def my_odeint(dfunc, V0, times, args=()):
    r'''For some of the systems odeint will switch to the BDF solver.
    In large enough systems, it then gets stuck trying to estimate the 
    Jacobian.

    This routine has identical inputs to integrate.odeint, but relies on 
    integrate.ode.  It avoids BDF.

    In particular, this seems to be important for SIS heterogeneous 
    pairwise where the number of equations is very large.  I have found 
    that near equilibrium, this often is interpreted as being a stiff 
    system and it switches to bdf, which requires calculating a 
    Jacobian.  In some systems this is impractically large.
    
    See this question: 
        http://stackoverflow.com/q/40317096/2966723,
    with the answer by 
        Phillip: http://stackoverflow.com/users/1881610/phillip
        
    INPUT and OUTPUT are as integrate.odeint
    '''

    r = integrate.ode(lambda t, X: dfunc(X, t, *args))
    r.set_integrator('vode', method='adams')
    r.set_initial_value(V0,times[0])
    V=[V0]
    for time in times[1:]:
        V.append(r.integrate(time))
    V = scipy.array(V)
    return V


def get_Nk_and_IC_as_arrays(G, rho, SIR=True):
    r'''
    Given the graph and initial proportion infected this finds the 
    initial conditions and number of nodes of each degree as needed
    by many of the differential equations models.
    
    :INPUTS:

    G : networkx graph

    rho : number between 0 and 1
          fraction of nodes to infect at time 0.

    SIR : boolean
          says whether the system will be SIR or SIS.

    :RETURNS:

    Nk : scipy array
         NUMBER (not proportion) of nodes of each degree.

    Sk0 : scipy array
          NUMBER of susceptible nodes of each degree at t=0, 
          = (1-rho)Nk

    Ik0 : scipy array    
          NUMBER of infected nodes of each degree at t=0,   
          = rho Nk

    if SIR, also returns
    Rk0 : scipy array
          NUMBER of recovered nodes of each degree at t=0,    
          = 0 Nk
    '''
    
    degree_count = Counter(G.degree().values())
    maxk = max(degree_count.keys())
    
    Nk = scipy.array([degree_count[k] for k in range(maxk+1)])
    Sk0 = (1-rho)*Nk
    Ik0 = rho*Nk
    Rk0 = 0*Nk
    
    if SIR:
        return Nk, Sk0, Ik0, Rk0
    else:
        return Nk, Sk0, Ik0

def get_NkNl_and_IC_as_arrays(G, rho, withKs = False, SIR=True):
    r'''
    In some of the differential equations models, we need to know how
    many edges exist between nodes of given degrees in a graph.  
    
    This finds that and the initial conditions for numbers of the
    various edges assuming a fraction rho is initially infected.
    
    INPUT
    -----
    G : networkx graph
    rho : number between 0 and 1
          fraction of nodes to infect at time 0.
    withKs : boolean
             flag to say whether we are restricting our attention to 
             just those degrees observed in the network or to all 
             degrees.
             If True, 
                 then we only consider those degrees that are observed.
             If False, 
                 then we treat it as if all degrees from 0 to kmax are 
                 observed.
    SIR : boolean
          says whether the system will be SIR or SIS.

    :RETURNS:

    NkNl : 2D scipy array
         NUMBER (not proportion) of edges between each pair of degrees.
    SkSl0 : 2D scipy array
          initial NUMBER of edges between pair of susceptibel nodes of 
          each degree type.
          = (1-rho)^2 NkNl
    SkIl0 : 2D scipy array    
          initial NUMBER of edges from a susceptible to an infected node 
          of the given degrees.
          = rho(1-rho) NkNl

    if not SIR, also returns
    IkIl0 : 2D scipy array
          initial NUMBER of edges between 2 infected nodes.  This is not 
          needed for SIR model.
          = rho^2*NkNl

    if withKs, also returns
    Ks : scipy array
         The observed degrees in the population.
    ''' 
    if withKs:
        Ks = sorted(list(set(G.degree().values())))
        klength = len(Ks)
    else:
        klength = max(G.degree().values())+1
    NkNl = scipy.zeros(shape=(klength,klength))
    NkNl = scipy.zeros(shape=(klength,klength))
    NkNl = scipy.zeros(shape=(klength,klength))
    for u,v in G.edges():
        k = G.degree(u)
        l = G.degree(v)
        NkNl[Ks.index(k)][Ks.index(l)] += 1
        NkNl[Ks.index(l)][Ks.index(k)] += 1
    SkSl0 = (1-rho)*(1-rho)*NkNl
    SkIl0 = (1-rho)*rho*NkNl
    IkIl0 = rho*rho*NkNl
    if withKs:
        if SIR:
            return NkNl, SkSl0, SkIl0, scipy.array(Ks)
        else:
            return NkNl, SkSl0, SkIl0, IkIl0, scipy.array(Ks)
    else:
        if SIR:
            return NkNl, SkSl0, SkIl0
        else:
            return NkNl, SkSl0, SkIl0, IkIl0

def get_Pk(G):
    r'''
    Used in several places so that we can input a graph and then we 
    can call the methods that depend on the degree distribution

    :INPUTS:

    G : networkx Graph

    :RETURNS:

    Pk : dict
         Pk[k] is the proportion of nodes with degree k.
    '''

    degree_count = Counter(G.degree().values())
    Pk = {x:degree_count[x]/float(G.order()) for x in degree_count.keys()}
    return Pk

def get_Psi(Pk):
    r'''
    Given a degree distribution (as a dict), returns the function psi
    
    :INPUTS:

    Pk : dict
         Pk[k] is the proportion of nodes with degree k.

    :RETURNS:

    psi : function.
          psi(x) = \sum_k Pk[k] x^k
    '''
    maxk = max(Pk.keys())
    Pkarray = scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    return lambda x: Pkarray.dot(x**k)

def get_PsiPrime(Pk):
    r'''
    Given a degree distribution (as a dict) returns the function
    dPsi(x)/dx
    
    :INPUTS:

    Pk : dict
         Pk[k] is the proportion of nodes with degree k.

    :RETURNS:

    psiPrime : function.
          \sum_k k Pk[k] x^{k-1}
    '''
    maxk = max(Pk.keys())
    Pkarray = scipy.array([Pk.get(k,0) for k in range(maxk+1)])

    return lambda x: Pkarray[k].dot(k*x**(k-1))

def get_PsiDPrime(Pk):
    r'''
    Given a degree distribution (as a dict) returns the function 
    
    d^2Psi(x)/dx^2
    
    :INPUTS:

    Pk : dict
         Pk[k] is the proportion of nodes with degree k.

    :RETURNS:

    psiDPrime : function.
          \sum_k k(k-1)Pk[k] x^{k-2}
    '''
    maxk = max(Pk.keys())
    Pkarray = scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    return lambda x: Pkarray[k].dot(k*(k-1)*x**(k-2))



##################
#                #
#    ODE CODE    #
#                #
##################


########      INDIVIDUAL BASED code -
########  given node and who its neighbors are, we track probability of
########  having given status based on probabilities of neighbors.  Assumes
########  independence.
    
def _dSIS_individual_based_(Y, t, G, nodelist, trans_rate_fxn, rec_rate_fxn):
    N = len(nodelist)
    dY = scipy.zeros(N)
    for index, (node, Yi) in enumerate(zip(nodelist,Y)):
        #This would probably be faster if it were done as a
        #matrix multiplication, but then I'd need to have
        #the matrix Gij explicitly included.  Perhaps that
        #would be better.  Networkx has something to create
        #numpy sparse matrices.  Perhaps that works?
        #No plan to do premature optimization.  Let's get it
        #working and then see if it's slow.
        dY[index] = sum(trans_rate_fxn(node,nbr)*(1-Y[node])*Y[nbr] 
                            for nbr in G.neighbors(node)) - rec_rate_fxn(node)*Yi
    return dY

def _dSIR_individual_based_(V, t, G, nodelist, index_of_node, trans_rate_fxn, 
                            rec_rate_fxn):
    '''    <\dot{X}_i> = - tau sum_j g_{ij} <Xi><Yj>
    <\dot{Y}_i> = tau sum_j g_{ij} <Xi><Yj> - gamma_i <Y_i>
    Z_i = 1-X_i-Y_i
    '''
    N = len(nodelist)
    X = V[:N]
    Y = V[N:]
    dX = scipy.zeros(N)
    dY = scipy.zeros(N)
    for index, (node, Xi, Yi) in enumerate(zip(nodelist,X, Y)):
        #This would probably be faster if it were done as a
        #matrix multiplication, but then I'd need to have
        #the matrix Gij explicitly included.  Perhaps that
        #would be better.  Networkx has something to create
        #numpy sparse matrices.  Perhaps that works?
        #No plan to do premature optimization.  Let's get it
        #working and then see if it's slow.
        
        dX[index] = -Xi*sum(trans_rate_fxn(node,nbr)*Y[index_of_node[nbr]] 
                                for nbr in G.neighbors(node))
        dY[index] =  -dX[index] - rec_rate_fxn(node)*Yi
    dV = scipy.concatenate((dX,dY), axis=0)
    return scipy.array(dV)

def SIS_individual_based(G, nodelist, Y0, tau, gamma, tmin = 0, 
                            tmax = 100, tcount = 1001, transmission_weight=None, 
                            recovery_weight=None, return_full_data = False):
    #tested in test_SIS_individual_based
    '''Encodes System (3.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention 
        strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology

    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    :INPUTS:

    G : Networkx graph
    
    Y0 : scipy array
         the array of initial infection probabilities

    nodelist : list
         list of nodes in G in the same order as in Y0

    tau : number
          transmission rate of disease

    gamma : number 
            global recovery rate 

    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    transmission_weight : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    return_full_data       (default False)
            If True, returns times, Ss, Is
            if False, returns times, S, I

    :RETURNS:

    if return_full_data is True:
        returns times, Ss, Is
           where times is a scipy array of times, Ss is a 2D scipy array
           Ss[i,j] gives probability nodelist[i] is susceptible at time
           times[j].
           Similarly for Is.
    if return_full data is False:
        returns times, S, I
             all are scipy arrays.  gives times, and expected number 
             susceptible and expected number infected.
             
    :SAMPLE USE:

    ::

        import networkx as nx
        import EoN as EoN
        import scipy
        
        G = nx.configuration_model([3,10]*1000)
        nodelist = G.nodes()
        N = G.order()
        rho = 1./N
        Y = rho*scipy.ones(N)
        t, S, I = EoN.SIS_individual_based(G, nodelist, Y, 0.3, gamma=1, 
                    tmax = 20)
    '''
    trans_rate_fxn, rec_rate_fxn = _get_rate_functions(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    times = scipy.linspace(tmin, tmax, tcount)
    Y = integrate.odeint(_dSIS_individual_based_, Y0, times,  
                                    args =(G, nodelist, trans_rate_fxn, rec_rate_fxn))
    Is = Y.T
    Ss = scipy.ones(len(Is))[:,None]-Is 
    
    if return_full_data:
        return times, Ss, Is
    else:
        return times, sum(Ss), sum(Is)


def SIR_individual_based(G, nodelist, X0, Y0, tau, gamma, tmin = 0, 
                            tmax = 100, tcount = 1001, transmission_weight=None, 
                            recovery_weight=None, return_full_data = False):
    '''
    Encodes System (3.30) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    See also:

    :INPUTS:

    G : Networkx graph

    X0 : scipy array
         the array of initial susceptibility probabilities
    Y0 : scipy array
         the array of initial infection probabilities

    nodelist : list
         list of nodes in G in the same order as in X0 and Y0

    tau : number
          transmission rate of disease

    gamma : number      (default None)
            global recovery rate  

    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    transmission_weight : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    return_full_data       (default False)
            If True, returns times, S, I, R, Ss, Is, Rs
            if False, returns times, S, I, R

    :RETURNS:

    if return_full_data is True:
        returns times, Ss, Is, Rs
           where times is a scipy array of times, Ss is a 2D scipy array
           Ss[i,j] gives probability nodelist[i] is susceptible at time
           times[j].
           Similarly for Is ans Rs
    if return_full data is False:
        returns times, S, I, R
             all are scipy arrays.  gives times, and expected number 
             susceptible, expected number infected, and expected number
             recovered
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        import scipy
        import matplotlib.pyplot as plt
        
        G = nx.configuration_model([3,10]*10000)
        tau = 0.3
        gamma = 1
        N = G.order()
        rho = 1./N

    nodelist = G.nodes()
    X = (1-rho)*scipy.ones(N)
    Y = rho*scipy.ones(N)
    
    t, S, I, R = EoN.SIR_individual_based(G, nodelist, X, Y, tau, gamma=gamma, tmax = 20)
    plt.plot(t,I)
    '''

    trans_rate_fxn, rec_rate_fxn = _get_rate_functions(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    index_of_node = {}
    for i, node in enumerate(nodelist):
        index_of_node[node] = i

    N = len(X0)
    times = scipy.linspace(tmin, tmax, tcount)
    V0 = scipy.concatenate((X0,Y0), axis=0)
    V = integrate.odeint(_dSIR_individual_based_, V0, times, 
                            args = (G, nodelist, index_of_node, 
                                    trans_rate_fxn, rec_rate_fxn))
    Ss = V.T[:N]
    S = Ss.sum(axis=0)
    Is = V.T[N:]
    I = Is.sum(axis=0)
    Rs = scipy.ones(N)[:, None] - Ss - Is
    R = Rs.sum(axis=0)
    if return_full_data:
        return times, S, I, R, Ss, Is, Rs
    else:
        return times, S, I, R


def SIS_individual_based_pure_IC(G, index_nodes, nodelist, tau, gamma, 
                                    tmin = 0, tmax = 100, tcount = 1001, 
                                    transmission_weight=None, 
                                    recovery_weight=None, 
                                    return_full_data = False):
    '''Encodes System (3.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The difference between this and SIS_individual_based is that this 
    one assumes a "pure initial condition", that is, we know exactly 
    what the statuses of the nodes are at the initial time.  
    
    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    :INPUTS:

    G : Networkx graph
    index_nodes : list or set
      the set of nodes initially infected
    nodelist : list
         list of nodes in G in the same order as in Y0
    tau : number
          transmission rate of disease
    gamma : number      (default None)
            global recovery rate  

    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    transmission_weight : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    return_full_data : boolean      (default False)


    RETURNS:
    --------
            if return_full_data is True,
                returns times, Ss, Is
            if return_full_data is False,
                returns times, S, I
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        import scipy
        import matplotlib.pyplot as plt
        
        G = nx.configuration_model([3,10]*1000)
        nodelist = G.nodes()
        index_nodes = range(100)
        t, S, I = EoN.SIS_individual_based(G, index_nodes, nodelist, 0.3,  
                    gamma=1, tmax = 20)
        plt.plot(t,I)

    '''
    #make Y0[u] be 1 if infected 0 if not
    Y0 = scipy.array([1 if u in index_nodes else 0 for u in nodelist])

    return SIS_individual_based(G, nodelist, Y0, tau, gamma, tmin, tmax, tcount,
                                transmission_weight, recovery_weight, return_full_data)
        



def SIR_individual_based_pure_IC(G, index_nodes, nodelist, tau, gamma, 
                                    initial_susceptible=None, tmin = 0, 
                                    tmax = 100, tcount = 1001, 
                                    transmission_weight=None, 
                                    recovery_weight=None, 
                                    return_full_data = False):
    '''Encodes System (3.30) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The difference between this and SIR_individual_based is that this 
    one assumes a "pure initial condition", that is, we know exactly 
    what the statuses of the nodes are at the initial time.
    
    <\dot{Y}_i> = tau \sum_j g_{ij} (1-<Y_i>)<Y_j>  -  gamma_i <Y_i>

    :INPUTS:

    G : Networkx graph
    
    index_nodes : list or set
      the set of nodes initially infected
      
    nodelist : list
         list of nodes in G in the same order as in Y0

    tau : number
          transmission rate of disease

    gamma : number      (default None)
            global recovery rate  

    initial_susceptible : list or set  (default None)
      initially susceptible nodes
      if equal to None, then all non-index nodes are initially 
      susceptible.

    tmin : number       (default 0)
           minimum report time

    tmax : number       (default 100)
           maximum report time

    tcount : integer       (default 1001)
             number of reports

    transmission_weight : string       (default None)
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    return_full_data : boolean      (default False)

    RETURNS:
    -------
    if return_full_data is True,
        returns times, S, I, R, Ss, Is, Rs
    if return_full_data is False,
        returns times, S, I, R
    
    '''
    nodelist = G.nodes()
    N = len(nodelist)
    #make Y0[u] be 1 if infected 0 if not
    Y0 = scipy.array([1 if u in index_nodes else 0 for u in nodelist])
    if initially_susceptible is None:
        X0 = scipy.ones - Y0
    else:
        X0 = scipy.array([1 if u in initially_susceptible else 0 
                            for u in nodelist])
    
    return SIR_individual_based(G, nodelist, X0, Y0, tau, gamma, tmin, 
                                tmax, tcount, transmission_weight, 
                                recovery_weight, return_full_data)

########   PAIR BASED

def _dSIS_pair_based_(V, t, G, nodelist, index_of_node, trans_rate_fxn, rec_rate_fxn):
    '''
    <\dot{Y}_i> = tau \sum_j g_{ij} <XiYj>  -  gamma_i <Yi>
    <\dot{XY}_ij> = tau sum_{k \neq i} g_{jk} <XiXj><XjYk>/<Xj>
                   - tau sum_{k neq j} g_{ik} <YkXi><XiYj>/<Xi>
                   - tau g_{ij}<XiYj> - gamma_j <XiYj>  
                   + ***gamma_i <YiYj>***
    <\dot{XX}_ij> = - tau sum_{k\neq i} g_{jk} <XiXj><XjYk>/<Xj>
                    - tau sum_{k neq j} g_{ik} <YkXi><XiXj>/<Xi>
                    + **** \gamma_i <YiXj> + gamma_j <XiYj>****

    <Xi>=1-<Yi>
    <YiYj> = 1 - <XiXj> - <XiYj> - <XjYi>
    <YiXj> = <XjYi>

    (Starred terms differ from SIR)
    
    The equations as coded involve all pairs rather than just the
    pairs that are in edges.  Those that are not part of an edge are
    set to zero and their derivatives are zero.  So the code could run
    faster if we took these out of the calculation.  That is a
    potential future improvement.

    I think for most cases this is a small contribution.
    All derivatives are initialized to 0, and then the loop only makes 
    changes for those terms where an edge exists.
    '''
    N=G.order()
    Y = V[0:N] #infecteds
    X = 1-Y    #susceptibles
    Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) 
            #there are places where we divide by X[i] which may = 0.
            #In those cases the numerator is (very) 0, so it's easier
            #to set this up as mult by inverse with a dummy value when
            #it is 1/0.
    Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[N: N+N**2]
    XX = V[N+N**2:]

    XY.shape = (N,N)
    XX.shape = (N,N)
    

    YX = XY.T   #not really needed, 
                #but helps keep consistent with equations as written.

    YY = 1 - XY-XX-YX

    dY = scipy.zeros(N)
    dXY = scipy.zeros((N,N))
    dXX = scipy.zeros((N,N))

    
    #I could make the below more efficient, but I think this sequence of for 
    # loops is easier to read, or at least understand.
    #I expect this isn't the bottleneck.  
    #Will avoid (premature) optimization for now.
    for u in nodelist:
        i = index_of_node[u]
        dY[i] += -rec_rate_fxn(u)*Y[i] 
        for v in G.neighbors(u):
            j = index_of_node[v]
            dY[i] += trans_rate_fxn(u,v)*XY[i,j]
            
            dXY[i,j] +=  - (trans_rate_fxn(u,v)+rec_rate_fxn(v))*XY[i,j] \
                            + rec_rate_fxn(u)*YY[i,j]
            dXX[i,j] +=  rec_rate_fxn(u)*YX[i,j] + rec_rate_fxn(v)*XY[i,j]
            #all the pure pairs are dealt with.  Now the triples
            for w in G.neighbors(u):
                if w == v: #skip these
                    continue
                #so w != v. 
                k= index_of_node[w]

                dXY[i,j] += trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j]  \
                            -  trans_rate_fxn(u,w) * YX[k,i] * XY[i,j]*Xinv[i]
                dXX[i,j] += -trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j] \
                            -  trans_rate_fxn(u,w) * YX[k,i] * XX[i,j]*Xinv[i]

    dXY.shape = (N**2,1)
    dXX.shape = (N**2,1)

    dV = scipy.concatenate((dY[:,None], dXY, dXX), axis=0).T[0]
    return dV

def _dSIR_pair_based_(V, t, G, nodelist, index_of_node, trans_rate_fxn, rec_rate_fxn):
    '''
    <\dot{X}_i> = -tau sum_j g_{ij} <XiYj>
    <\dot{Y}_i> = tau \sum_j g_{ij} <XiYj>  -  gamma_i <Y_i>
    <\dot{XY}_ij> = tau sum_{k \neq i} g_{jk} <XiXj><XjYk>/<Xj>
                   - tau sum_{k neq j} g_{ik} <YkXi><XiYj>/<Xi>
                   - tau g_{ij}<XiYj> - gamma_j <XiYj> 
    <\dot{XX}_ij> = -tau sum_{k\neq j} gik <YkXi><XiXj>/<Xi>
                    -tau sum_{k neq i} gjk <XiXj><XjYk>/<Xj>
    <>
    The equations as coded involve all pairs rather than just the
    pairs that are in edges.  Those that are not part of an edge are
    set to zero and their derivatives are zero.  So the code could run
    faster if we took these out of the calculation. I think for most
    cases this is a small contribution.  Before I forced the initial
    conditions for these nonedges to be 0, they caused quite a bit of
    numerical headaches.
    '''
    N=G.order()
    X = V[0:N] #susceptibles
    Y = V[N:2*N] #infecteds
    Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) 
            #there are places where we divide by X[i] which may = 0.
            #In those cases the numerator is (very) 0, so it's easier
            #to set this up as mult by inverse with a dummy value when
            #it is 1/0.
    Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[2*N: 2*N+N**2]
    XX = V[2*N+N**2:]

    #print X.shape, Y.shape, XY.shape, XX.shape, N
    XY.shape = (N,N)
    XX.shape = (N,N)
    

    YX = XY.T   #not really needed, 
                #but helps keep consistent with equations as written.

    dX = scipy.zeros(N)
    dY = scipy.zeros(N)
    dXY = scipy.zeros((N,N))
    dXX = scipy.zeros((N,N))

    
    #I could make the below more efficient, 
    #but I think this sequence of for loops is easier to read, 
    #or at least understand.
    #I expect it to run quickly regardless.  Will avoid (premature) 
    #optimization for now.
    for u in nodelist:
        i = index_of_node[u]
        dY[i] += -rec_rate_fxn(u)*Y[i] 
        for v in G.neighbors(u):
            j = index_of_node[v]
            dX[i] += -trans_rate_fxn(u,v)*XY[i,j]
            dY[i] += trans_rate_fxn(u,v)*XY[i,j]
            
            dXY[i,j] +=  - (trans_rate_fxn(u,v)+rec_rate_fxn(v))*XY[i,j] 

            #all the pure pairs are dealt with.  Now the triples
            for w in G.neighbors(v):
                if w == u: #skip these
                    continue
                #so w != u.  
                k= index_of_node[w]
                #i corresponds to u, j to v and k to w.
                dXY[i,j] += trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j]  
                dXX[i,j] += -trans_rate_fxn(v,w) * XX[i,j] * XY[j,k]*Xinv[j] 
            for w in G.neighbors(u):
                if w == v:
                    continue #skip these
                k = index_of_node[w]
                dXY[i,j] += -  trans_rate_fxn(u,w) * YX[k,i] * XY[i,j]*Xinv[i]
                dXX[i,j] += -  trans_rate_fxn(u,w) * YX[k,i] * XX[i,j]*Xinv[i]

    dXY.shape = (N**2,1)
    dXX.shape = (N**2,1)

    dV = scipy.concatenate((dX[:, None], dY[:,None], dXY, dXX), axis=0).T[0]
    return dV


def SIS_pair_based(G, nodelist, Y0, tau, gamma, XY0=None, XX0 = None, 
                    tmin = 0, tmax = 100, tcount = 1001, transmission_weight=None, 
                    recovery_weight=None, return_full_data = False):
    r'''
    Encodes System (3.26) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    This system solves equations for an SIS disease model spreading on a 
    given graph.  
    
    It captures the dependence with pairs, but not triples.

    It does not include corrections for triangles (or any other cycles).  
    
    The corrections for triangles are provided in the text, but not 
    implemented here.

    There are some inefficiencies in the implementation:
        we track all pairs, rather than just those pairs in edges, but 
        this is unlikely to significantly affect the calculation time.  
        
        This makes it much easier to vectorize things.
        
        We track pairs in both directions: e.g., XX[1,2] and XX[2,1].


    INPUT:
    ------
    G : Networkx graph

    nodelist : list
         list of nodes in G in the same order as in Y0

    Y0 : scipy array
         the array of initial infection probabilities for each node in 
         order as in nodelist

    tau : number
          transmission rate of disease

    gamma : number (default None)
            global recovery rate  

    XY0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XY0[i,j] is probability node i is susceptible and j is 
            infected.
            if None, then assumes that infections are introduced 
            randomly according to Y0.

    XX0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XX0[i,j] is probability nodes i and j are susceptible.
            if None, then assumes that infections are introduced 
            randomly according to Y0.

    tmin : number (default 0)
           minimum report time

    tmax : number (default 100)
           maximum report time 

    tcount : integer (default 1001)
             number of reports

    transmission_weight : string
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    return_full_data : boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R

    :RETURNS:

    if return_full_data is True:
        returns times, S, I, Xs, Ys, XY, XX
    if False:
        returns times, S, I
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        
        G = nx.fast_gnp_random_graph(1000,0.004)
        nodelist = G.nodes()
        Y0 = scipy.array([1 if node<10 else 0 for node in nodelist]) #infect first 10
        t, S, I = EoN.SIS_pair_based(G, nodelist, Y0, 2, 0.5, tmax = 4, tcount = 101)
        plt.plot(t,I)
        
'''
    trans_rate_fxn, rec_rate_fxn = _get_rate_functions(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    times = scipy.linspace(tmin,tmax,tcount)


    N = len(Y0)
    X0=1-Y0

    if XY0 is None:
        XY0 = X0[:,None]*Y0[None,:]
    else:
        if XY0.shape != (N,N):
            raise EoNError("incompatible lengths for XY0 and Y0")

    if XX0 is None:
        XX0 = X0[:,None]*X0[None,:]
    else:
        if XX0.shape != (N,N):
            raise EoNError("incompatible lengths for XX0 and Y0")
    A = nx.adjacency_matrix(G).toarray()
    XY0 = XY0*A  #in principle the equations should still work for pairs not
    XX0 = XX0*A  #in an edge, but this led to the error with odeint.  
                 #Multiplying by A restricts attention just to present edges.
    
    XY0.shape=(N**2,1)
    XX0.shape=(N**2,1)
    
    V0 = scipy.concatenate((Y0[:,None], XY0, XX0), axis=0).T[0]
    index_of_node = {node:i for i, node in enumerate(nodelist)}

    V = integrate.odeint(_dSIS_pair_based_, V0, times, 
                            args = (G, nodelist, index_of_node, trans_rate_fxn, 
                                    rec_rate_fxn))
    Ys = V.T[0:N]
    I = Ys.sum(axis=0)
    Xs = scipy.ones(N)[:,None]-Ys
    S = Xs.sum(axis=0)
    if return_full_data:
        XY = V.T[N: N+N**2]
        XX = V.T[N+N**2:]
        XY.shape = (N,N,tcount)
        XX.shape = (N,N,tcount)
        return times, S, I, Xs, Ys, XY, XX
    else:
        return times, S, I






def SIR_pair_based(G, nodelist, Y0, tau, gamma, X0 = None, XY0=None, 
                    XX0 = None, tmin = 0, tmax = 100, tcount = 1001, 
                    transmission_weight=None, recovery_weight=None, 
                    return_full_data = False):
    '''
    Encodes System (3.39) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    This system solves equations for an SIR disease model spreading on a
    given graph.  It captures the dependence with pairs, but not 
    triples.
    
    It will be exact for a tree.

    There are NO CORRECTIONS for the existence of TRIANGLES or any other
    CYCLES.
    
    Some corrections for triangles are provided in the text, but not 
    implemented here.
    
    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention 
        strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology
    

    INPUT:
    ------
    G : Networkx graph
    
    nodelist : list
         list of nodes in G in the same order as in Y0

    Y0 : scipy array
         the array of initial infection probabilities for each node in 
         order as in nodelist

    tau : number
          transmission rate of disease

    gamma : number (default None)
            global recovery rate  

    X0 : scipy array (default None)
            probability a random node is initially susceptible.
            the probability of initially recovered will be 1-X0-Y0.  By 
            default we assume no initial recoveries, so X0=1-Y0 will be 
            assumed in this case.

    XY0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XY0[i,j] is probability node i is susceptible and j is 
                infected.
            if None, then assumes that infections are introduced 
                randomly according to Y0.

    XX0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XX0[i,j] is probability nodes i and j are susceptible.
            if None, then assumes that infections are introduced 
                randomly according to Y0.

    tmin : number (default 0)
           minimum report time

    tmax : number (default 100)
           maximum report time 

    tcount : integer (default 1001)
             number of reports

    transmission_weight : string
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    return_full_data : boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R

    :RETURNS:

    if return_full_data is True:
        returns times, S, I, R, Xs, Ys, Zs, XY, XX
    if False:
        returns times, S, I, R

    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        
        G = nx.fast_gnp_random_graph(1000,0.004)
        nodelist = G.nodes()
        Y0 = scipy.array([1 if node<10 else 0 for node in nodelist]) #infect first 10
        t, S, I, R = EoN.SIR_pair_based(G, nodelist, Y0, 2, 0.5, tmax = 4, tcount = 101)
        plt.plot(t,I)
    '''


    trans_rate_fxn, rec_rate_fxn = _get_rate_functions(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    times = scipy.linspace(tmin,tmax,tcount)
    if X0 is None:
        X0 = 1-Y0
    N = len(Y0)

    if XY0 is None:
        XY0 = X0[:,None]*Y0[None,:]
    else:
        if XY0.shape != (N,N):
            raise EoNError("incompatible lengths for XY0 and Y0")
    if XX0 is None:
        XX0 = X0[:,None]*X0[None,:]
    else:
        if XX0.shape != (N,N):
            raise EoNError("incompatible lengths for XX0 and Y0")
    A = nx.adjacency_matrix(G).toarray()
    XY0 = XY0*A  #in principle the equations should still work for pairs not
    XX0 = XX0*A  #in an edge, but this led to the error with odeint.  
                 #Multiplying by A restricts attention just to present edges.
    
    XY0.shape=(N**2,1)
    XX0.shape=(N**2,1)
    
    V0 = scipy.concatenate((X0[:,None], Y0[:,None], XY0, XX0), axis=0).T[0]
    #print V0.shape
    index_of_node = {node:i for i, node in enumerate(nodelist)}
    #    index_of_node = {}
    #    for i, node in enumerate(nodelist):
    #        index_of_node[node] = i


    V = integrate.odeint(_dSIR_pair_based_, V0, times, 
                            args = (G, nodelist, index_of_node, trans_rate_fxn, 
                                    rec_rate_fxn))
    Xs = V.T[0:N]
    S = Xs.sum(axis=0)
    Ys = V.T[N:2*N]
    I = Ys.sum(axis=0)
    Zs = scipy.ones(N)[:,None]-Xs-Ys
    R = Zs.sum(axis=0)
    if return_full_data:
        XY = V.T[2*N: 2*N+N**2]
        XX = V.T[2*N+N**2:]
        #print len(XX)
        XY.shape = (N,N,tcount)
        XX.shape = (N,N,tcount)
        return times, S, I, R, Xs, Ys, Zs, XY, XX
    else:
        return times, S, I, R


def _dSIR_pair_based2_(V, t, G, nodelist, index_of_node, edgelist, 
                        index_of_edge, trans_rate_fxn, rec_rate_fxn):
    #X, Y, XY, YX, XX
    N=len(nodelist)
    E=len(edgelist)
    
    X = V[0:N]
    Y = V[N:2*N]
    def Xinv(u):
        '''
        There are places where we divide by X[i] which may = 0.  
        
        In those cases the numerator is (very) 0, so it's easier
        to set this up as mult by inverse with a dummy value when it 
        is 1/0.
        '''
        index = index_of_node[u]
        if X[index]==0:
            return 0
        else:
            return 1/X[index]
    def Yinv(u):
        index = index_of_node[u]
        if Y[index]==0:
            return 0
        else:
            return 1/Y[index]
    #Xinv = scipy.array([1/v if v!=0 else 0 for v in X]) 
    #Yinv = scipy.array([1/v if v!=0 else 0 for v in Y])
    
    XY = V[2*N: 2*N+E]
    YX = V[2*N+E:2*N+2*E]
    XX = V[2*N+2*E:]

    #print X.shape, Y.shape, XY.shape, XX.shape, N

    dX = scipy.zeros(N)
    dY = scipy.zeros(N)
    dXY = scipy.zeros(E)
    dYX = scipy.zeros(E)
    dXX = scipy.zeros(E)

    def my_XY(u,v): #either XY[i,j] or YX[i,j], depending on which is defined.  
                    #Might be cleaner with Try, Except
        if index_of_edge.has_key((u,v)):
            index = index_of_edge[(u,v)]
            return XY[index]
        else:
            index = index_of_edge[(v,u)]
            return YX[index]
    def my_YX(u,v):
        return my_XY(v,u)
    def my_XX(u,v):
        if index_of_edge.has_key((u,v)):
            index = index_of_edge[(u,v)]
        else:
            index = index_of_edge[(v,u)]
        return XX[index]
    
    #I could make the below more efficient, but I think this sequence 
    #of for loops is easier to read, or at least understand.
    #I expect it to run quickly regardless.  
    #Will avoid (premature) optimization for now.
    for i, u in enumerate(nodelist):
        dY[i] += -rec_rate_fxn(u)*Y[i] 
    for edgeindex, (u,v) in enumerate(edgelist):
        i = index_of_node[u]
        j = index_of_node[v]

        dX[i] += -trans_rate_fxn(u,v)*XY[edgeindex] 
        dY[i] += trans_rate_fxn(u,v)*XY[edgeindex]  

        dX[j] += -trans_rate_fxn(u,v)*YX[edgeindex]
        dY[j] += trans_rate_fxn(u,v)*YX[edgeindex]

        dXY[edgeindex] += - (trans_rate_fxn(u,v)+rec_rate_fxn(v))*XY[edgeindex]
        dYX[edgeindex] += - (trans_rate_fxn(u,v)+rec_rate_fxn(v))*YX[edgeindex]

        for w in G.neighbors(u):
            if w==v:  #skip this case
                continue
            k=index_of_node[w]
            dXY[edgeindex] += -trans_rate_fxn(w,u)*my_YX(w,u) * my_XY(u,v)*Xinv(u)
            dYX[edgeindex] += trans_rate_fxn(w,u)*my_XX(v,u)*my_XY(u,w)*Xinv(u)
            dXX[edgeindex] += -trans_rate_fxn(w,u) * my_XX(u,v) \
                                * my_XY(u,w)*Xinv(u) 
        for w in G.neighbors(v):
            if w==u:
                continue
            #w transmits to v.
            dXY[edgeindex] += trans_rate_fxn(w,v)*my_XX(u,v) *my_XY(v,w)*Xinv(v)
            dYX[edgeindex] += -trans_rate_fxn(w,v)*my_XY(v,w)*my_YX(u,v)*Xinv(v)
            dXX[edgeindex] += -trans_rate_fxn(w,v)*my_XX(u,v)*my_XY(v,w)*Xinv(v)

    dV = scipy.concatenate((dX, dY, dXY, dYX, dXX), axis=0)
    return dV


def _SIR_pair_based_initialize_node_data(G, rho, nodelist, X0, Y0):
    #inputs must define either rho or Y0.  In first case nodelist is optional.
    if (rho and Y0 is not None) or (not rho and Y0 is None):
        raise EoNError("need rho or Y0 defined for initial condition, \
                        but not both")
    if Y0 is not None and nodelist is None:
        raise EoNError("order in Y0 is ambiguous if nodelist is not given.")

    if not nodelist: #then rho is defined but not Y0
        nodelist = list(G.nodes())
    if rho: #Y0 not defined
        Y0 = rho*scipy.ones(len(nodelist))
        X0 = 1-Y0 #assume X0=0
    else:  #Y0 is defined
        if not X0:
            X0=1-Y0  #assume Z0=0
        #otherwise X0 is given and Z0=1-X0-Y0
    return nodelist, X0, Y0

def _SIR_pair_based_initialize_edge_data(G, edgelist, nodelist, XY0, YX0, 
                                            XX0, X0, Y0, index_of_node):
    if (not XY0 is None or YX0 is None or XX0 is None) \
            and (XY0 is not None or YX0 !=None  or XX0 is not None):  
            #at least one defined and one not defined
        raise EoNError("must define all of XY0, YX0, and XX0 or none of them")
    if not edgelist:
        if XY0:
            raise EoNError("order in XY0, YX0, and XX0 is ambiguous if \
                            edgelist is not given.")
        else:
            edgelist = list(G.edges())

    if XY0:
        #test that XY0 <= X0Y0, same for  YX0 and XX0
        for index,(u,v) in enumerate(edgelist):
           i_u = index_of_node[u]
           i_v = index_of_node[v]
           if XY0[index] >X0[i_u]*Y0[i_v] or YX0[index]>Y0[i_u]*X0[I_v] \
                                or XX0[index]>X0[i_u]*X0[i_v]:
               raise EoNError("edge probabilities inconsistent with node \
                                probabilities")
    else:
        XY0 = scipy.array([X0[index_of_node[u]]*Y0[index_of_node[v]] 
                            for u,v in edgelist])
        YX0 = scipy.array([Y0[index_of_node[u]]*X0[index_of_node[v]] 
                            for u,v in edgelist])
        XX0 = scipy.array([X0[index_of_node[u]]*X0[index_of_node[v]] 
                            for u,v in edgelist])
    return edgelist, XY0, YX0, XX0

    
def SIR_pair_based2(G, tau, gamma, rho = None, nodelist=None, X0=None, 
                    Y0=None, edgelist = None, XY0=None, YX0 = None, 
                    XX0 = None, tmin = 0, tmax = 100, tcount = 1001, 
                    transmission_weight=None, recovery_weight=None, 
                    return_full_data = False):
    '''
    Encodes System (3.39) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    This system solves equations for an SIR disease model spreading on a 
    given graph.  It captures the dependence with pairs, but not 
    triples.

    It will be exact for a tree.

    There are NO CORRECTIONS for the existence of TRIANGLES or any other 
    CYCLES.
    
    Some corrections for triangles are provided in the text, but not 
    implemented here.

    See also:
    Hadjichrysanthou and Sharkey
    Epidemic control analysis: Desigining targeted intervention 
    strategies against epidemics propagated on contact networks,
    Journal of Theoretical Biology

    <\dot{Y}_i> = tau \sum_j g_{ij} <XY>>  -  gamma_i <Y_i>
    <\dot{XY}> = tau sum_{k \neq i} g_{jk} <XX><XjYk>/<Xj>
                   - tau sum_{k neq j} g_{ik} <YkXi><XY>/<Xi>
                   - tau g_{ij}<XY> - gamma_j <XIYj> 
    <\dot{XX}> = 
    <>
    The equations as coded involve all pairs rather than just the pairs 
    that are in edges.  
    
    Those that are not part of an edge are set to zero and their 
    derivatives are zero.  
    
    So the code could run faster, but I think for most cases this is a 
    small contribution.  
    
    Before I forced the initial conditions for these nonedges to be 0, 
    they caused quite a bit of numerical headaches.

    :INPUTS:
        
    G : Networkx graph

    nodelist : list
         list of nodes in G in the same order as in Y0

    Y0 : scipy array
         the array of initial infection probabilities for each node in 
         order as in nodelist

    tau : number
          transmission rate of disease

    gamma : number (default None)
            global recovery rate  

    XY0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XY0[i,j] is probability node i is susceptible and j is 
            infected.
            if None, then assumes that infections are introduced 
                randomly according to Y0.

    XX0 : 2D scipy array (default None)
            (each dimension has length number of nodes of G)
            XX0[i,j] is probability nodes i and j are susceptible.
            if None, then assumes that infections are introduced 
                randomly according to Y0.

    tmin : number (default 0)
           minimum report time

    tmax : number (default 100)
           maximum report time 

    tcount : integer (default 1001)
             number of reports

    transmission_weight : string
            the label for a weight given to the edges.
            G.edge[i][j][transmission_weight] = g_{ij}

    recovery_weight : string       (default None)
            a label for a weight given to the nodes to scale their 
            recovery rates
                gamma_i = G.node[i][recovery_weight]*gamma

    return_full_data : boolean      (default False)
            if True:
                returns times, S, I, R, Xs, Ys, Zs, XY, XX
            if False:
                returns times, S, I, R
    '''

    #note: we do not test whether the nodelist and edgelist are in 
    #fact lists of nodes or edges
    nodelist, X0, Y0 = _SIR_pair_based_initialize_node_data(G, rho, nodelist, 
                                                                X0, Y0)

    index_of_node = {node:i for i, node in enumerate(nodelist)}

    edgelist, XY0, YX0, XX0 = _SIR_pair_based_initialize_edge_data(G, 
                                                edgelist, nodelist, XY0, YX0, 
                                                XX0, X0, Y0, index_of_node)

    index_of_edge = {edge:i for i, edge in enumerate(edgelist)}
    N = len(nodelist)
    E = len(edgelist)
    #now we define functions which give the transmission rate of edges 
    #and recovery rate of nodes.  
    trans_rate_fxn, rec_rate_fxn = _get_rate_functions(G, tau, gamma, 
                                                transmission_weight,
                                                recovery_weight)

    times = scipy.linspace(tmin,tmax,tcount)

    
    V0 = scipy.concatenate((X0, Y0, XY0, YX0, XX0),axis=0)
    V = integrate.odeint(_dSIR_pair_based2_, V0, times, 
                            args = (G, nodelist, index_of_node, edgelist, 
                                    index_of_edge, trans_rate_fxn, rec_rate_fxn))
                                    
    #dX, dY, dXY, dYX, dXX
    Xs = V.T[0:N]
    S = Xs.sum(axis=0)
    Ys = V.T[N:2*N]
    I = Ys.sum(axis=0)
    Zs = scipy.ones(N)[:,None]-Xs-Ys
    R = Zs.sum(axis=0)
    if return_full_data:
        XY = V.T[2*N: 2*N+E]
        YX = V.T[2*N+E:2*N+2*E]
        XX = V.T[2*N+2*E:]
        YY = 1 - XY - YX-XX
        return times, S, I, R, Xs, Ys, Zs, XY, YX, XX, YY,edgelist, nodelist
    else:
        return times, S, I, R




######    HOMOGENEOUS MEANFIELD

def _dSIS_homogeneous_meanfield_(X, t, n_over_N, tau, gamma):
    S, I = X
    dSdt = gamma * I - tau*n_over_N*S*I
    dIdt = -dSdt
    dX = scipy.array([dSdt,dIdt])
    return dX

def _dSIR_homogeneous_meanfield_(X, t, n_over_N, tau, gamma):
    S, I = X
    dSdt = -tau*n_over_N*S*I
    dIdt = tau*n_over_N*S*I - gamma*I
    dX = scipy.array([dSdt, dIdt])
    return dX

           
def SIS_homogeneous_meanfield(S0, I0, n, tau, gamma, tmin=0, tmax=100, 
                                tcount=1001):
    '''Encodes System (4.8) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    In the text this is often referred to as the 
    "mean-field model closed at the level of pairs"

    [\dot{S}] = \gamma [I] - tau n[S][I]/N
    [\dot{I}] = \tau n[S][I]/N - \gamma [I]


    :INPUTS:

    S0 : number
         initial number susceptible
    I0 : number
         initial number infected
    n : integer
        degree of all nodes.
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports

    :RETURNS:

    times, S, I, all scipy arrays
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        S0 = 999
        I0 = 1
        n = 4 #degree
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_homogeneous_meanfield(S0, I0, n, tau, gamma)
    '''

    N=S0+I0
    X0=scipy.array([S0,I0])
    times=scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_homogeneous_meanfield_, X0,times, 
                            args=(float(n)/N, tau, gamma))
    S, I= X.T
    return times, S, I

def SIR_homogeneous_meanfield(S0, I0, R0, n, tau, gamma, tmin=0, tmax=100, 
                                tcount=1001):
    '''Encodes System (4.9) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "mean-field model closed at the level of pairs"

    [\dot{S}] = - tau n[S][I]/N
    [\dot{I}] = \tau n[S][I]/N - \gamma [I]
    [\dot{R}] = \gamma [I]


    :INPUTS:

    S0 : number
         initial number susceptible
    I0 : number
         initial number infected
    R0 : number
         initial number recovered
    n : integer
        degree of each node
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    
    :RETURNS:

    times, S, I, R : all scipy arrays
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        S0 = 999
        I0 = 1
        n = 4 #degree
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_homogeneous_meanfield(S0, I0, 0, n, tau, gamma)
        
    '''

    N=S0+I0+R0
    X0= scipy.array([S0,I0])
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIR_homogeneous_meanfield_, X0, times, 
                            args=(float(n)/N, tau, gamma))
    S, I= X.T
    R = N-S-I
    return times, S, I, R


####   HOMOGENEOUS PAIRWISE

def _dSIS_homogeneous_pairwise_(X, t, N, n, tau, gamma):
    r'''
    [\dot{S}] = gamma [I] - tau [SI]
    [\dot{I}] = \tau [SI] - \gamma [I] = -[\dot{S}]
    [\dot{SI}] = \gamma([II]-[SI])+ \tau ((n-1)/n) [SI]([SS]-[SI])/[S] 
                 - \tau [SI]
    [\dot{SS}] = 2\gamma[SI] - 2\tau ((n-1)/n) [SI][SS]/[S]
    [\dot{II}] = -2\gamma[II] + 2\tau((n-1)/n) [SI]^2/[S] + 2\tau[SI]

    conserved quantities: [S]+[I]
                          [SS]+2[SI]+[II]

    n([S]+[I]) should equal [SS]+2[SI]+[II], so II will be calculated 
    based on this.
    '''
    S, SI, SS = X 
    I = N-S
    II = N*n - SS - 2*SI
    nm1_over_n = (n-1.)/n
    dSdt = gamma * I - tau*SI
    #dIdt = -dSdt
    dSIdt = gamma*(II-SI) + tau*nm1_over_n*SI*(SS-SI)/S - tau*SI
    dSSdt = 2*gamma*SI - 2*tau*nm1_over_n*SI*SS/S
    #dIIdt = -2*gamma*II + 2*tau*nm1_over_n*SI**2/S + 2*tau*SI
    dX = scipy.array([dSdt,dSIdt,dSSdt])
    return dX 


def _dSIR_homogeneous_pairwise_(X, t, n, tau, gamma):
    S,I,SI,SS = X
    nm1_over_n = (n-1.)/n
    dSdt = -tau*SI
    dIdt = tau*SI - gamma*I
    dSIdt = -gamma*SI + tau*nm1_over_n * SI*(SS-SI)/S - tau*SI
    dSSdt = -2*tau*nm1_over_n*SI*SS/S
    dX =  scipy.array([dSdt,dIdt,dSIdt,dSSdt])
    return dX

def SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, tmin = 0, 
                                tmax=100, tcount=1001, 
                                return_full_data=False):
    r'''Encodes System (4.10) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "mean-field model closed at the level of triples"

    :INPUTS:

    S0 : number
         initial number susceptible
    I0 : number
         initial number infected
    SI0 : number
          initial number of SI edges
    SS0 : number
          initial number of SS edges
    n : integer
          (common) degree of nodes.
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I or 
                       all calculated data.
    
    :RETURNS:

    if return_full_data is True:
        t, S, I, SI, SS, II
    if return_full_data is False:
        t, S, I
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        S0 = 990
        I0 = 10
        SI0 = 50
        SS0 = 4900
        n = 5
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, 
                                                tmax = 20)
    '''
    N = S0+I0

    if SS0 + SI0*2>n*N:
        raise EoNError('Initial condition has more SS, SI, and IS edges than allowed')

    X0 = scipy.array([S0, SI0, SS0])
    
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_homogeneous_pairwise_, X0, times, 
                            args=(N, n, tau, gamma))
    S, SI, SS= X.T
    I = N-S
    
    if return_full_data:
	II = N*n - SS-2*SI
        return times, S, I, SI, SS, II
    else:
        return times, S, I
    
    
def SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, tmin = 0, 
                                tmax=100, tcount=1001, 
                                return_full_data=False):
    '''Encodes System (4.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "mean-field model closed at the level of triples"

    [\dot{S}] = - tau [SI]
    [\dot{I}] = \tau [SI] - \gamma [I]
    [\dot{R}] = \gamma [I]    ;    [R] = N-[S]-[I]
    [\dot{SI}] = -\gamma [SI]+ \tau ((n-1)/n) [SI]([SS]-[SI])/[S] 
                 - \tau [SI]
    [\dot{SS}] = - 2\tau ((n-1)/n) [SI][SS]/[S]

    conserved quantities: [S]+[I]+[R]  also 
                          [SS]+[II]+[RR] + 2([SI] + [SR] + [IR])

    :INPUTS:

    S0 : Initial number suusceptible
    I0 : Initial number infected
    R0 : Initial number recovered
    SI0 : Initial number of SI edges
    SS0 : Initial number of SS edges
    n : Degree of nodes
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       tells whether to just return times, S, I, R or 
                       all calculated data.
                       if True, then returns times, S, I, R, SI, SS
    :RETURNS:

    if return_full_data is True:
        times, S, I, R, SI, SS
    if return_full_data is False:
        times, S, I, R 

    :SAMPLE USE:


    ::

        import networkx as nx
        import EoN
        S0 = 990
        I0 = 10
        R0 = 1
        SI0 = 45
        SS0 = 4900
        n = 5
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, 
                                                    tmax = 20)

    '''
    N = S0+I0+R0
    if SS0 + 2*SI0 > n*N:
        raise EoNError('Initial condition has more SS, SI, and IS edges than allowed')
    X0 = scipy.array([S0, I0, SI0, SS0])
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIR_homogeneous_pairwise_, X0, times, 
                            args=(n, tau, gamma))
    S, I, SI, SS = X.T
    R = N-S-I
    if return_full_data:
        return times, S, I, R, SI, SS
    else:
        return times, S, I, R


def SIS_homogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, 
                                        tmax=100, tcount=1001, 
                                        return_full_data=False):
    r'''
    Calls SIS_homogeneous_pairwise with a graph, disease parameters, and
    a random fraction rho initially infected.

    :INPUTS:

    G : networkx Graph
    tau : number
          transmission rate
    gamma : number
            recovery rate
    rho : number (default 1/N)
          initial fraction infected 
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       tells whether to just return times, S, I, or 
                       all calculated data.
                       if True, then returns times, S, I, SI, SS
    :RETURNS:

    if return_full_data is True:
        t, S, I, SI, SS, II
    if return_full_data is False:
        t, S, I
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        G = nx.fast_gnp_random_graph(10000,0.0005)
        tau = 1
        gamma = 3
        rho = 0.02
        t, S, I = EoN.SIS_homogeneous_pairwise_from_graph(G, tau, gamma, rho, tmax = 20)
    '''
    
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    n = sum(k*Pk[k] for k in Pk.keys())
    N=G.order()
    S0 = (1-rho)*N
    I0 = rho*N
    SI0 = (1-rho)*N*n*rho
    SS0 = (1-rho)*N*n*(1-rho)
    return SIS_homogeneous_pairwise(S0, I0, SI0, SS0, n, tau, gamma, tmin, 
                                        tmax, tcount, return_full_data)

def SIR_homogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, 
                                        tmax=100, tcount=1001, 
                                        return_full_data=False):
    r'''
    Calls SIS_homogeneous_pairwise with a graph, disease parameters, and
    a random fraction rho initially infected.

    :INPUTS:

    G : networkx Graph
    tau : number
          transmission rate
    gamma : number
            recovery rate
    rho : number (default 1/N)
          initial fraction infected 
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       tells whether to just return times, S, I, R or 
                       all calculated data.
                       if True, then returns times, S, I, R, SI, SS
    :RETURNS:

    if return_full_data is True:
        t, S, I, SI, SS, II
    if return_full_data is False:
        t, S, I
        
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        G = nx.fast_gnp_random_graph(10000,0.0005)
        tau = 1
        gamma = 3
        rho = 0.02
        t, S, I, R = EoN.SIR_homogeneous_pairwise_from_graph(G, tau, gamma, rho, 
                                                                tmax = 20)
    '''
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    n = sum(k*Pk[k] for k in Pk.keys())
    N=G.order()
    S0 = (1-rho)*N
    I0 = rho*N
    SI0 = (1-rho)*N*n*rho
    SS0 = (1-rho)*N*n*(1-rho)
    R0=0
    return SIR_homogeneous_pairwise(S0, I0, R0, SI0, SS0, n, tau, gamma, tmin,
                                    tmax, tcount, return_full_data)




#######     HETEROGENEOUS MEAN FIELD

def _dSIS_heterogeneous_meanfield_(X, t, kcount, tau, gamma):
    ks = scipy.arange(kcount)
    S = scipy.array(X[:kcount])
    I = scipy.array(X[kcount:])
    pi_I = ks.dot(I)/ ks.dot(I+S)

    Sdot = gamma*I - tau * ks*S*pi_I
    Idot = tau*ks*S*pi_I - gamma*I
    dX = scipy.concatenate((Sdot,Idot), axis=0)
    return dX
    
def _dSIR_heterogeneous_meanfield_(X, t, S0, Nk, tau, gamma):
    theta = X[0]
    Rk = scipy.array(X[1:])
    ks = scipy.arange(len(Rk))
    Sk = S0*(theta**ks)
    Ik = Nk - Sk - Rk
    pi_I = ks.dot(Ik)/ks.dot(Nk)
    dRkdt = gamma*Ik
    dThetadt = - tau *pi_I * theta

    dX = scipy.concatenate(([dThetadt],dRkdt), axis=0)
    return dX




def SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmin = 0, tmax=100, 
                                tcount=1001, return_full_data=False):
    '''Encodes System (5.10) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "heterogeneous mean-field model closed at the level of pairs"

    a few notes on the inputs:
    Sk0 is an array (or a list). 
    
    It is not a dict.  
    
    Sk0[k] is the *number* of nodes that are susceptible and have degree 
    k (even if some degrees missing).  
    
    A dict like this can be converted into an array by
    Sk0 = scipy.array([Sk0dict.get(k,0) 
                       for k in xrange(max(Sk0dict.keys())+1)])

    Ik0 is similar to Sk0.


    [\dot{S}_k] = \gamma [I_k] - \tau k [S_k] \pi_I
    [\dot{I}_k] = -(above)
    \pi_I = \sum_k k [I_k] / \sum_k  k [N_k]



    :INPUTS:

    Sk0 : scipy array
          number susceptible for each k
    Ik0 : scipy array
          number infected for each k
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       tells whether to just return times, S, I or all 
                       calculated data.
                       if True, returns t, S, I, Sk, Ik

    :RETURNS:

    
    if return_full_data is True:
        times, S, I, Sk  (Sk is scipy 2D arrays)
    if return_full_data is False:
        times, S, I      (all scipy arrays)
        
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        Sk0 = [995, 995, 995, 995, 995]
        Ik0 = [5, 5, 5, 5, 5]
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmax = 10)
    '''
    if len(Sk0) != len(Ik0):
        raise EoNError('length of Sk0 not equal to length of Ik0')
    Sk0 = scipy.array(Sk0)
    Ik0 = scipy.array(Ik0)

    kcount = len(Sk0)
    
    X0 = scipy.concatenate((Sk0,Ik0), axis=0)
    Nk = Sk0+Ik0
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIS_heterogeneous_meanfield_, X0, times, 
                            args=(kcount,tau,gamma))
    Sk = scipy.array(X.T[:kcount])
    Ik = scipy.array(X.T[kcount:])
    S = Sk.sum(axis=0)
    I = Ik.sum(axis=0)
    if return_full_data:
	return times, S, I, Sk, Ik
    else:
        return times, S, I
	
    
def SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, tmin = 0, tmax=100, 
                                    tcount=1001, return_full_data=False):
    '''
    Encodes System (5.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "heterogeneous mean-field model closed at the level of pairs"

    Ik0 and Rk0 are similar to Sk0.

    [S_k] = [S_k](0) theta^k
    [I_k] = [N_k] - [S_k] - [R_k]
    [\dot{R}_k] = \gamma [I_k]
    pi_I = \sum_k k[I_k]


    :INPUTS:

    Sk0 : array
          Sk0[k] is the number of
          nodes that are susceptible and have degree k (even if some degrees 
          missing).
    Ik0 : array
    Rk0 : array
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or 
                       all calculated data.

    :RETURNS:

    
    if return_full_data is True:
        times, S, I, R, Sk, Ik, Rk (the Xk are scipy 2D arrays)
    if return_full_data is False:
        times, S, I, R          (all scipy arrays)
        
    :SAMPLE USE:


    ::

        import networkx as nx
        import EoN
        Sk0 = [995, 995, 995, 995, 995]
        Ik0 = [5, 5, 5, 5, 5]
        Rk0 = [0,0,0,0,0]
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, 
                                                        tmax = 10)
    '''
    if len(Sk0) != len(Ik0) or len(Sk0) != len(Rk0):
        raise EoNError('length of Sk0, Ik0, and Rk0 must be the same')

    theta0=1
    Sk0 = scipy.array(Sk0)
    Ik0 = scipy.array(Ik0)
    Rk0 = scipy.array(Rk0)
    Nk = Sk0+Ik0 +Rk0
    X0 = scipy.concatenate(([theta0],Rk0), axis=0)
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIR_heterogeneous_meanfield_, X0, times, 
                            args = (Sk0, Nk, tau, gamma))
    
    theta = X.T[0]
    Rk = X.T[1:]
    ks = scipy.arange(len(Rk))
    L=(theta[None,:]**ks[:,None])
    Sk=Sk0[:,None]*L
    Ik = Nk[:,None]-Sk-Rk
    if not return_full_data:
	return times, sum(Sk), sum(Ik), sum(Rk)
    else:
	return times, Sk, Ik, Rk
    
def SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, rho = None, 
                                            tmin = 0, tmax=100, tcount=1001, 
                                            return_full_data=False):
    r'''
    Takes a graph and an initial proportion infected rho.  
    Calculates Sk0 and Ik0 and calls the heterogeneous meanfield model

    :INPUTS:

    G : networkx Graph
    tau : number
          transmission rate
    gamma : number
            recovery rate
    rho : number between 0 and 1  (default None)
          the fraction to be randomly infected at time 0
          If None, then rho=1/N is used where N = G.order()
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or 
                       all calculated data.

    :RETURNS:

    
    if return_full_data is True:
        times, S, I, Sk, Ik, (the Xk are scipy 2D arrays)
    if return_full_data is False:
        times, S, I         (all scipy arrays)
        
    :SAMPLE USE:


    ::

        import networkx as nx
        import EoN
        G = nx.configuration_model([1,2,3,4]*1000)
        tau = 1
        gamma = 2
        t, S, I = EoN.SIS_heterogeneous_meanfield_from_graph(G, tau, gamma, 
                                                                tmax = 15)

    '''
    if rho is None:
        rho = 1./G.order()
    #Pk = get_Pk(G)
    #maxk = max(Pk.keys())
    Nk, Sk0, Ik0 = get_Nk_and_IC_as_arrays(G, rho, SIR=False)#G.order()*scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    return SIS_heterogeneous_meanfield(Sk0, Ik0, tau, gamma, tmin, tmax, tcount, return_full_data)

def SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, rho = None, 
                                            tmin = 0, tmax=100, tcount=1001, 
                                            return_full_data=False):
    r'''
    Takes a graph and an initial proportion infected rho.  
    Calculates Sk0 and Ik0 and calls the heterogeneous meanfield model

    :INPUTS:

    G : networkx Graph
    tau : number
          transmission rate
    gamma : number
            recovery rate
    rho : number between 0 and 1  (default None)
          the fraction to be randomly infected at time 0
          If None, then rho=1/N is used where N = G.order()
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or 
                       all calculated data.


    :RETURNS:

    if return_full_data is True
        times, Sk, Ik, Rk (the Xk are scipy 2D arrays)
    if False,
        times, S, I, R (all scipy arrays)
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        G = nx.configuration_model([1,2,3,4]*1000)
        tau = 1
        gamma = 2
        t, S, I, R = EoN.SIR_heterogeneous_meanfield_from_graph(G, tau, gamma, 
                                                                tmax = 10)
    
    '''
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    
    return SIR_heterogeneous_meanfield(Sk0, Ik0, Rk0, tau, gamma, tmin, tmax, 
                                        tcount, return_full_data=False)


#######      HETEROGENEOUS PAIRWISE
def _dSIS_heterogeneous_pairwise_(X, t, Nk, NkNl, tau, gamma, Ks):
    '''
    Gives the derivatives


    [\dot{Sk}] = gamma [Ik] - tau [SkI]
    [\dot{Ik}] = tau [SkI] - gamma [Ik]   = - [\dot{Sk}]
    [\dot{SkIl}] = gamma([IkIl] - [SkIl]) + tau([SkSlI] - [ISkIl] 
                                                 - [SkIl])
    [\dot{SkSl}] = gamma([SkIl]+[IkSl]) - tau ([SkSlI] + [ISkSl])
    [\dot{IkIl}] = -2 gamma[IkIl] + tau ([SkIl]+[IkSl]+[ISkIl]+[IkSlI])
    [AlSkI] = ((k-1)/k) [AlSk][SkI]/[Sk]
    [ISkAl] = ((k-1)/k) [ISk][SkAl]/[Sk]

    So: [SkSlI] = ((l-1)/l) [SkSl][SlI]/[Sl] 
        [ISkIl] = ((k-1)/k) [ISk] [SkIl]/[Sk]
        [ISkSl] = ((k-1)/k) [ISk] [SkSl]/[Sk]
        [IkSlI] = ((l-1)/l) [IkSl][SlI]/[Sl]

    conserved quantities : [Sk]+[Ik],  [SkSl] + [SkIl] + [IkIl] + [IkSl]

    identities: [IkSl] = [SlIk], so where IkSl needed, use SlIk.T


    :INPUTS:

    X : current values of variables
    t : current time
    Nk : The number of nodes of each degree; Nk[i] is the number of 
         nodes with degree Ks[i]
    NkNl : number of edges of various types.  NkNl[i,j] corresponds to 
           Ks[i] and Ks[j].
    tau : transmission rate
    gamma : recovery rate
    Ks : a scipy array --- gives the observed degrees in increasing 
         order.
    '''
    kcount = len(Ks)
    Sk = X[:kcount]
    SkSl = X[kcount: kcount + kcount**2]
    SkIl = X[kcount+kcount**2: kcount + 2*kcount**2]

    SkSl.shape=(kcount,kcount)
    SkIl.shape=(kcount,kcount)

    Ik = Nk - Sk
    SkI = SkIl.sum(1)#sum(SkIl.T)
    IkSl = SkIl.T
    IkIl = NkNl - SkSl - SkIl - IkSl

    Ls = Ks #easier to keep as l for readability of eqns

    SlI = SkI #easier in this form for readability
    ISk = SkI

    tmpSk = 1.*Sk #the 1.* is to make it a copy, not the same object
    tmpSk[tmpSk==0] = 1

    SkSlI = SkSl * ((Ls-1))*SlI/ (Ls*tmpSk)
    ISkIl = (ISk * ((Ks-1))*SkIl.T/ (Ks*tmpSk)).T
    ISkSl = SkSlI.T 
    IkSlI = ISkIl.T 
    
    dSk = gamma*Ik - tau *SkI
    dSkSl = gamma*(SkIl + IkSl) - tau*(SkSlI + ISkSl)
    dSkIl = gamma*(IkIl-SkIl) + tau*(SkSlI - ISkIl - SkIl)

    dSkIl.shape=(kcount**2,1)
    dSkSl.shape=(kcount**2,1)
    dX = scipy.concatenate((dSk[:,None], dSkSl, dSkIl),axis=0).T[0]
    return dX

def _dSIR_heterogeneous_pairwise_(X, t, tau, gamma, Nk, Ks):
    kcount = len(Ks)
    Sk = X[:kcount]
    Ik = X[kcount:2*kcount]
    SkSl = X[2*kcount:2*kcount+kcount**2]
    SkIl = X[2*kcount+kcount**2:2*kcount+2*kcount**2]

    SkSl.shape=(kcount,kcount)
    SkIl.shape=(kcount,kcount)

    SkI = SkIl.sum(1)
    IkSl = SkIl.T

    Ls = Ks #easier to keep as l for readability of eqns
    SlI = SkI
    ISk = SkI

    tmpSk = 1.*Sk #the 1.* is to make it a copy, not the same object
    tmpSk[tmpSk==0] = 1
    tmpSl = tmpSk

    SkSlI = SkSl * ((Ls-1))*SlI/ (Ls*tmpSl)
    ISkIl = (ISk * ((Ks-1))*SkIl.T/ (Ks*tmpSk)).T
    ISkSl = SkSlI.T 
    IkSlI = ISkIl.T 

    dSk = -tau*SkI
    dIk = tau*SkI - gamma*Ik
    dSkIl = -gamma*SkIl + tau*(SkSlI - ISkIl - SkIl)
    dSkSl = -tau*(SkSlI + ISkSl)

    dSkIl.shape=(kcount**2,1)
    dSkSl.shape = (kcount**2,1)

    dX = scipy.concatenate((dSk[:,None], dIk[:,None], 
                                dSkSl, dSkIl),axis=0).T[0]
    return dX



def SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, gamma, 
                                tmin = 0, tmax=100, tcount=1001, 
                                return_full_data = False, Ks = None):
    '''Encodes System (5.13) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "heterogeneous mean-field model closed at the level of triples"

    :INPUTS:

    Sk0 : array.  Sk0[k] is the number of
                 nodes that are susceptible and have degree k.  If one 
                 is empty, it becomes 0.
                 (if Ks is defined, the definition changes slightly, see
                  below)

    Ik0 : array
         similar to Sk0, but for infected.
        (if Ks is defined, the definition changes slightly, see below)

    SkSl0, SkIl0, and IkIl0 : 2D arrays
         SkIl0[k][l] is [S_kI_l]
             see below for constraints these should satisfy related to 
             Sk0 and Ik0.  
             The code does not enforce these constraints.
             (if Ks is defined, the definition changes slightly, 
              see below)

    tau : number
          transmission rate

    gamma : number
            recovery rate

    tmin : number (default 0)
           minimum report time

    tmax : number (default 100)
           maximum report time 

    tcount : integer (default 1001)
             number of reports

    return_full_data : boolean (default False)
                       If True, return times, Sk, Ik, SkIl, SkSl, IkIl
                       If False, return times, S, I

    Ks : scipy array. (default None)
         (helps prevent memory errors) if some degrees are not
         observed, then the corresponding entries of these arrays are
         zero.  This can lead to memory errors in the case of a
         network with many missing degrees.  So Ks is an (assumed)
         ordered vector stating which Ks are actually observed.  Then
         the Sk0[i] is the number of nodes that are susceptible and
         have degree Ks[i].  Similarly for Ik0 and SkIl0 etc.

    In principle, there are constraints relating Sk with SkSl and SkIl 
    and similarly relating Ik with IkIl and SkIl.T.  
    
    No attempt is made to enforce these.  
    
    It is assumed the user will ensure acceptible inputs.

    We could also remove Sk0 and Ik0 as inputs and infer them from the 
    others, but for consistency with elsewhere, this is not done here.
    
    
    :RETURNS:

    if return_full_data is True:
        returns times, S, I, Sk, Ik, SkIl, SkSl, IkIl
    if return_full_data is False:
        returns times, S, I
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        import scipy
        Sk0 = 100 * scipy.ones(4)
        Ik0 = scipy.zeros(4)
        Ik0[3]=1
        SkSl0 = scipy.matrix([[0, 0,0,0],[0,100,0,0],[0,0,200,0],[0,0,0,294]])
        #only interact within a degree class, so the deg 1 and 2 are safe.
        SkIl0 = scipy.zeros((4,4))
        SkIl0[3,3] = 3
        IkIl0 = scipy.zeros((4,4))
        tau = 1
        gamma = 1
        
        t, S, I = EoN.SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, 
                                                    gamma)
    
    '''

    if Ks is None:
        Ks = scipy.array(range(len(Sk0)))
    times = scipy.linspace(tmin,tmax,tcount)
    Nk = Sk0+Ik0
    kcount = len(Nk)
    NkNl = SkSl0 + SkIl0 + IkIl0 + SkIl0.T
    SkSl0 = SkSl0.copy()
    SkIl0 = SkIl0.copy()
    SkSl0.shape = (kcount**2,1)
    SkIl0.shape = (kcount**2,1)
    X0 = scipy.concatenate((Sk0[:,None], SkSl0, SkIl0), axis=0).T[0]

    X = my_odeint(_dSIS_heterogeneous_pairwise_, X0, times, 
                    args = (Nk, NkNl, tau, gamma, Ks))

    kcount = len(Nk)
    Sk = X.T[:kcount]
    print Sk.size
    S = Sk.sum(axis=0)
    print S.size
    Ik = Nk[:,None] - Sk
    I = Ik.sum(axis=0)
    if return_full_data:
	SkSl = X.T[kcount:kcount+kcaount**2]
	SkIl = X.T[kcount+kcount**2:]
        SkSl.shape = (kcount,kcount,tcount)
        SkIl.shape = (kcount,kcount,tcount)
	IkIl = NkNl - SkSl - SkIl - SkIl.T
	return times, S, I, Sk, Ik, SkIl, SkSl, IkIl
    else:
	return times, S, I

def SIR_heterogeneous_pairwise(Sk0, Ik0, Rk0, SkSl0, SkIl0, tau, gamma, 
                                tmin = 0, tmax=100, tcount=1001, 
                                return_full_data=False, Ks = None):
    '''Encodes System (5.15) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    In the text this is often referred to as the 
    "heterogeneous mean-field model closed at the level of triples"

    [\dot{S}_k] = -tau [S_k I]
    [\dot{I}_k] = tau [S_k I] - gamma [I_k]
    [\dot{R}_k] = gamma [I_k]  (but using Rk=Nk-Sk-Ik for this equation)
    [\dot{S_kI_l}] = -gamma[S_k I_l] + tau([S_k S_l I] - [I S_k I_l] 
                                           - [S_k I_l])
    [\dot{S_kS_l}] = -tau([S_k S_l I] + [I S_k S_l])

    [A_l S_k I] = ((k-1)/k) [A_l S_k] [S_k I]/ [S_k]
    [I S_k A_l] = ((k-1)/k) [I S_k] [S_k A_l]/ [S_k]

    :INPUTS:

    Sk0 : scipy array, Sk0[k] is number of degree k susceptible at 
          time 0.
    Ik0 : scipy array
    Rk0 : scipy array
    SkIl0 : scipy 2D array
    SkSl0 : scipy 2D array

    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       If True, return times, Sk, Ik, Rk, SkIl, SkSl
                       If False, return times, S, I, R


    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or 
                       all calculated data.
    Ks : scipy array. (default None)
         (helps prevent memory errors) if some degrees are not
         observed, then the corresponding entries of these arrays are
         zero.  This can lead to memory errors in the case of a
         network with many missing degrees.  So Ks is an (assumed)
         ordered vector stating which Ks are actually observed.  Then
         the Sk0[i] is the number of nodes that are susceptible and
         have degree Ks[i].  Similarly for Ik0 and SkIl0 etc.        

    :RETURNS:

    if return_full_data is True
        returns times, S, I, R, Sk, Ik, Rk, SkIl, SkSl
    if return_full_data is False
        return times, S, I, R
    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        
'''
    
    if Ks is None:
        Ks = scipy.array(range(len(Sk0)))
#    print "Ks is ", Ks
    times = scipy.linspace(tmin,tmax,tcount)

    Nk = Sk0+Ik0+Rk0
    kcount = len(Ks)
    SkSl0.shape = (kcount**2,1)
    SkIl0.shape = (kcount**2,1)

    X0 = scipy.concatenate((Sk0[:,None], Ik0[:,None], SkSl0, SkIl0), 
                                axis=0).T[0]
    X = integrate.odeint(_dSIR_heterogeneous_pairwise_, X0, times, 
                            args = (tau, gamma, Nk, Ks))

    Sk = X.T[:kcount]
    S = Sk.sum(axis=0)
    Ik = X.T[kcount:2*kcount]
    I = Ik.sum(axis=0)
    SkIl = X.T[2*kcount:2*kcount+kcount**2]
    SkSl = X.T[2*kcount+kcount**2: 2*kcount+2*kcount**2]

    Rk = Nk[:,None] - Sk - Ik
    R = Rk.sum(axis=0)
    
    if return_full_data:
        SkIl.shape = (kcount,kcount,tcount)
        SkSl.shape = (kcount,kcount,tcount)
        return times, S, I, R, Sk, Ik, Rk, SkIl, SkSl
    else:
        return times, S, I, R


def SIS_heterogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, 
                                            tmax=100, tcount=1001, 
                                            return_full_data = False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0 = get_Nk_and_IC_as_arrays(G, rho, SIR=False)

    NkNl, SkSl0, SkIl0, IkIl0, Ks = \
                get_NkNl_and_IC_as_arrays(G, rho, withKs = True, SIR=False)

    Sk0 = scipy.array([Sk0[k] for k in Ks])
    Ik0 = scipy.array([Ik0[k] for k in Ks])
    return SIS_heterogeneous_pairwise(Sk0, Ik0, SkSl0, SkIl0, IkIl0, tau, 
                                        gamma, tmin, tmax, tcount, 
                                        return_full_data, 
                                        Ks = scipy.array(Ks))

    
def SIR_heterogeneous_pairwise_from_graph(G, tau, gamma, rho = None, tmin=0, 
                                            tmax=100, tcount = 1001, 
                                            return_full_data=False):
    if rho==None:
        rho = 1./G.order()    
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    NkNl, SkSl0, SkIl0, Ks = get_NkNl_and_IC_as_arrays(G, rho, withKs = True, 
                                                        SIR = True)

    Sk0 = scipy.array([Sk0[k] for k in Ks])
    Ik0 = scipy.array([Ik0[k] for k in Ks])
    Rk0 = scipy.array([Rk0[k] for k in Ks])
    return SIR_heterogeneous_pairwise(Sk0, Ik0, Rk0, SkSl0, SkIl0, tau, gamma,
                                        tmin = tmin, tmax=tmax, tcount=tcount,
                                        return_full_data=return_full_data,
                                        Ks = Ks)




#######    COMPACT PAIRWISE


def _dSIS_compact_pairwise_(X, t, Nk, twoM, tau, gamma):
    SI, SS= X[-2:]
    Sk = X[:-2]

    Ik = Nk - Sk
    II = twoM - SS - 2*SI

    ks = scipy.arange(len(Sk))
    SX = float(SI + SS)
    Q = (1./SX**2) * (ks*(ks-1)).dot(Sk)

    dSk = gamma *Ik - tau * ks *Sk *SI/SX
    dSI = gamma*(II-SI) + tau*(SS-SI)*(SI)*Q - tau*SI
    dSS = 2*gamma*SI - 2*tau*SS*SI*Q

    dX = scipy.concatenate((dSk, [dSI, dSS]), axis=0)
    return dX

def _dSIR_compact_pairwise_(X, t, N, tau, gamma):
    SS, SI, R = X[-3:]
    Sk = X[:-3]
    ks = scipy.arange(len(Sk))
    SX = ks.dot(Sk)
    Q = (ks*(ks-1)).dot(Sk)/SX**2
    I = N - sum(Sk) - R
    
    dSk = -tau * ks * Sk *SI/SX
    dSS = -2*tau*SS*SI*Q
    dSI = -gamma*SI + tau*(SS-SI)*SI*Q-tau*SI
    dR = gamma*I

    dX = scipy.concatenate((dSk, [dSS,dSI,dR]),axis=0)
    return dX 


def SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin = 0, 
                            tmax=100, tcount=1001, return_full_data=False):
    '''Encodes system (5.18) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [dot{S}_k] = gamma[I_k] - tau k [S_k] [SI]/[SX]
    [dot{I}_k] = tau * k *[S_k] [SI]/SX - gamma [I_k] = -[dot{S}_k]
    [dot{SI}] = gamma([II]-[SI]) + tau([SS]-[SI])[SI]Q - tau[SI]
    [dot{SS}] = 2 gamma[SI] - 2 tau [SS] [SI] Q
    [dot{II}] = 2 tau[SI] = 2 gamma[II] + 2 tau[SI]^2Q
    [SX] = sum_k k [S_k]
    Q = (1/[SX]^2) sum_k (k-1)k[S_k]

    INPUT
    -----
    Sk0 : scipy array
          number susceptible for each k
    Ik0 : scipy array
          number infected for each k
    SI0 : number
          number of SI edges
    SS0 : number
          number of SS edges
    II0 : number
          number of II edges
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       if True, return times, S, I, Sk, Ik, SI, SS, II
		       if False,  return times, S, I




    conserved quantities:  [Sk]+[Ik]     ;     SS + II + 2SI
    
    
    :RETURNS:

    
    :SAMPLE USE:

    ::


        import networkx as nx
        import EoN
        
    '''
    Nk = Sk0+Ik0
    twoM = SS0+II0+2*SI0
    times = scipy.linspace(tmin,tmax,tcount)
    X0 = scipy.concatenate((Sk0, [SI0,SS0]), axis=0)
    X = integrate.odeint(_dSIS_compact_pairwise_, X0, times, 
                            args = (Nk, twoM, tau, gamma))

    SI, SS = X.T[-2:]
    Sk = X.T[:-2]
    S = Sk.sum(axis=0)

    Ik = Nk[:,None] - Sk
    I = Ik.sum(axis=0)

    II = twoM - SS - 2*SI

    if return_full_data:
	return times, S, I, Sk, Ik, SI, SS, II
    else:
	return times, S, I

def SIR_compact_pairwise(Sk0, I0, R0, SS0, SI0, tau, gamma, tmin=0, tmax=100,
                            tcount=1001, return_full_data=False):
    '''Encodes system (5.19) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    [dot{S}_k] = -tau k [S_k] [SI]/[SX]
    [dot{SS}] = -2 tau [SS] [SI] Q
    [dot{SI}] = -gamma [SI] + tau([SS]-[SI])[SI]Q - tau [SI]
    [dot{R} = gamma [I]
    [SX] = sum_k k[S_k]
    Q = (1/[SX]^2) sum_k (k-1) k [S_k]
    [S] = sum [S_k]
    I = N-[S]-R

    :INPUTS:

    Sk0 : scipy array
          initial number of suscetibles of each degree k
    I0 : number
         initial number infected
    R0 : number
         initial number recovered
    SS0 : number
          initial number of SS edges
    SI0 : number
          initial number of SI edges
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or 
                       all calculated data.

    :RETURNS:


    :SAMPLE USE:

    ::

       
        import networkx as nx
        import EoN
    '''
    
    times = scipy.linspace(tmin,tmax,tcount)
    N = I0+R0+sum(Sk0)
    X0 = scipy.concatenate((Sk0, [SS0, SI0, R0]), axis=0)
    X = integrate.odeint(_dSIR_compact_pairwise_, X0, times, 
                            args = (N, tau, gamma))
    SI, SS, R = X.T[-3:]
    Sk = X.T[:-3]
    S = Sk.sum(axis=0)
    I = N - R - S
    if return_full_data:
	return times, Sk, I, R, SS, SI
    else:
	return times, S, I, R


def SIS_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    N = G.order()
    maxk = max(Pk.keys())
    Nk = G.order()*scipy.array([Pk.get(k,0) for k in range(maxk+1)])
    Sk0 = (1-rho)*Nk
    Ik0 = rho*Nk
    SI0 = sum(Nk[k]*k*(1-rho)*rho for k in range(maxk+1))
    SS0 = sum(Nk[k]*k*(1-rho)*(1-rho) for k in range(maxk+1))
    II0 = sum(Nk[k]*k*rho*rho for k in range(maxk+1))
    return SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin, 
                                    tmax, tcount, return_full_data)

    
def SIR_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]
    
    I0 = sum(Ik0)
    R0 = 0

    SX0 = scipy.dot(Sk0,ks)
    SS0 = (1-rho)*SX0
    SI0 = rho*SX0
    
    return SIR_compact_pairwise(Sk0, I0, R0, SS0, SI0, tau, gamma, tmin=tmin,
                                tmax=tmax, tcount=tcount, 
                                return_full_data=return_full_data)




#######SUPER COMPACT PAIRWISE

def _dSIS_super_compact_pairwise_(X, t, tau, gamma, N, k_ave, ksquare_ave, 
                                    kcube_ave):
    '''    
    [\dot{I}] = tau [SI] - gamma [I]
    [\dot{SS}] = 2 gamma [SI] - 2 tau [SI] [SS] Q
    [\dot{SI}] = gamma ([II]-[SI]) + tau [SI] ([SS]-[SI])Q - tau [SI]
    [\dot{II}] = -2 gamma [II] + 2 tau [SI]^2 Q + 2 tau [SI]
    Q = ((<K^2>(<K^2>-n_S<K>) + <K^3>(n_S-<K>))
        /
        (n_S(<K^2>-<K>^2)) - 1)/n_S[S]
    n_S = ([SI] + [SS])/(N-[I])
    '''

    I, SS, SI, II = X
    S = N-I
    
    n_S = (SS+SI)/(S)
    
    Q = ((ksquare_ave*(ksquare_ave-n_S*k_ave) \
            + kcube_ave*(n_S-k_ave)) *(n_S*(ksquare_ave-k_ave**2))-1)/(S*n_S)
    dIdt = tau*SI - gamma*I
    dSSdt = 2*gamma*SI - 2*tau*SI*SS*Q
    dSIdt = gamma*(II - SI) + tau*SI*(SS-SI)*Q - tau*SI
    dIIdt = -2*gamma*II + 2*tau*SI**2*Q + 2*tau*SI

    dX = scipy.array([dIdt, dSSdt, dSIdt, dIIdt])
    return dX

def _dSIR_super_compact_pairwise_(X, t, tau, gamma, psihat, psihatPrime, 
                                    psihatDPrime, N):
    theta = X[0]
    SS = X[1]
    SI = X[2]
    R = X[3]
    S = N * psihat(theta)
    I = N - S - R

    Q = psihatDPrime(theta)/(N*psihatPrime(theta)**2)
    dThetadt = - tau*SI/(N*psihatPrime(theta))
    dSSdt = -2*tau*SS*SI*Q
    dSIdt = - gamma*SI + tau*(SS-SI)*SI*Q - tau*SI
    dRdt = gamma*I

    dX = scipy.array([dThetadt, dSSdt, dSIdt, dRdt])
    return dX


def SIS_super_compact_pairwise(S0, I0, SS0, SI0, II0, tau, gamma, k_ave, 
                                ksquare_ave, kcube_ave, tmin = 0, tmax=100, 
                                tcount=1001, return_full_data=False):
    '''
    Encodes system (5.20) of Kiss, Miller, & Simon.  Please cite the 
    book if using this algorithm.

    :INPUTS:

    S0 : 
    I0 :
    SS0 :
    SI0 :  
    II0 : 
    tau : number
          transmission rate
    gamma : number
            recovery rate
    k_ave : number
            average value of k
    ksquare_ave : number
            average value of k**2
    kcube_ave : number
            average value of k**3
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
             tells whether to just return times, S, I, R or all 
             calculated data.

    :RETURNS:

    if return_full_data is True
        returns times, S, I, SS, SI, II
    if return_full_data is False
        returns times, S, I
        
    :SAMPLE USE:

    import networkx as nx
    import EoN
    '''
    X0 = [I0, SS0, SI0, II0]
    N = S0+I0
    times = scipy.linspace(tmin,tmax,tcount)
    X = integrate.odeint(_dSIS_super_compact_pairwise_, X0, times, 
                            args = (tau, gamma, N, k_ave, ksquare_ave, 
                                    kcube_ave))
    I, SS, SI, II = X.T
    S = N-I
    if return_full_data:
        return times, S, I, SS, SI, II
    else:
        return times, S, I



def SIR_super_compact_pairwise(SS0, SI0, R0, N, tau, gamma, psihat, 
                                psihatPrime, psihatDPrime, tmin = 0, 
                                tmax = 100, tcount = 1001, 
                                return_full_data = False):
    '''
    Encodes system (5.22) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{theta} = -tau [SI]/N*psihat(theta)
    [dot{SS}] = -2 tau [SS] [SI] Q
    [dot{SI}] = -gamma[SI] + tau([SS]-[SI])[SI]Q - tau*[SI]
    [dot{R}] = gamma*[I]
    [S] = N psihat(theta)
    [I] = N-[S]-[R]
    Q = psihat_xx(theta)/N(psihat_x(theta))^2

    
    :INPUTS:

    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
            tells whether to just return times, S, I, R or all 
            calculated data.
'''
    times = scipy.linspace(tmin,tmax,tcount)
    X0 = scipy.array([1., SS0, SI0, R0])
    X = integrate.odeint(_dSIR_super_compact_pairwise_, X0, times, 
                            args = (tau, gamma, psihat, psihatPrime, 
                                    psihatDPrime, N))
    theta, SS, SI, R = X.T
    S = N*psihat(theta)
    I = N-S-R
    if return_full_data:
	return times, S, I, R, SS, SI
    else:
	return times, S, I, R

def SIS_super_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0,
                                            tmax=100, tcount=1001, 
                                            return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0 = get_Nk_and_IC_as_arrays(G, rho, SIR=False)
    ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]

    S0 = sum(Sk0)
    I0 = sum(Ik0)
    R0 = 0

    SX0 = scipy.dot(Sk0,ks)
    SS0 = (1-rho)*SX0
    SI0 = rho*SX0
    II0 = scipy.dot(Nk,ks)-SX0
    Pk = get_Pk(G)
    Pks = scipy.array([Pk.get(k,0) for k in ks])
    k_ave = scipy.dot(Pks, ks)
    ksquare_ave = scipy.dot(Pks, ks*ks)
    kcube_ave = scipy.dot(Pks, ks*ks*ks)
    
    return SIS_super_compact_pairwise(S0, I0, SS0, SI0, II0, tau, gamma, 
                                        k_ave, ksquare_ave, kcube_ave, 
                                        tmin=tmin, tmax=tmax, tcount=tcount, 
                                        return_full_data=return_full_data)

def SIR_super_compact_pairwise_from_graph(G, tau, gamma, rho = None, tmin = 0,
                                            tmax=100, tcount=1001, 
                                            return_full_data=False):
    if rho is None:
        rho = 1./G.order()
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    ks = scipy.array(range(len(Nk))) #[0,1,2,3,...,k]

    N = G.order()
    R0 = 0

    SX0 = Sk0.dot(ks)
    SS0 = (1-rho)*SX0
    SI0 = rho*SX0

    Pk = get_Pk(G)
    def psihat(x): #probably faster if vectorized, 
                   #but need to be careful with broadcasting...
        return (1-rho)*sum(Pk[k]*x**k for k in Pk)
    def psihatPrime(x):
        return (1-rho)*sum(k*Pk[k]*x**(k-1) for k in Pk)
    def psihatDPrime(x):
        return (1-rho)*sum(k*(k-1)*Pk[k]*x**(k-2) for k in Pk)

    return  SIR_super_compact_pairwise(SS0, SI0, R0, N, tau, gamma, psihat, 
                                        psihatPrime, psihatDPrime, 
                                        tmin = tmin, tmax = tmax, 
                                        tcount = tcount, 
                                        return_full_data = return_full_data)






########     EFFECTIVE DEGREE


def _dSIS_effective_degree_(X, t, original_shape, tau, gamma):
    '''
    \dot{S}_{s,i} = - tau i S_{s,i} + gamma*I_{s,i}
                    + gamma((i+1)S_{s-1,i+1}-iS_{s,i})
                    + tau[ISS]((s+1)S_{s+1,i-1} - sS_{s,i})/[SS]

    \dot{I}_{s,i} = tau i S_{s,i} - gamma I_{s,i}
                    + gamma((i+1)I_{s-1,i+1} - iI_{s,i})
                    + tau([ISI]/[SI] + 1)((s+1)I_{s+1,i-1} - sI_{s,i})
    S = sum S_{s,i}
    I = sum I_{s,i}
    '''
    #note that dSIR_effective_degree has some commands for vectorizing 
    #tentatively coded (and commented out)
    ksq = original_shape[0]*original_shape[1]
    Ssi = X[:ksq]
    Isi = X[ksq:]
    Ssi.shape = original_shape
    Isi.shape = original_shape

    ISS = sum([sum([i*s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    SS = sum([sum([s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    ISI = sum([sum([i*(i-1)*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    SI = sum([sum([i*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])

    g1 = scipy.zeros(original_shape)
    g2 = scipy.zeros(original_shape)
    t1 = scipy.zeros(original_shape)
    t2 = scipy.zeros(original_shape)
    
    dSsi = scipy.zeros(original_shape)
    dIsi = scipy.zeros(original_shape)
    for s in xrange(original_shape[0]):
        for i in xrange(original_shape[1]):
            if s==0 or i+1 == original_shape[1]:
                Ssm1ip1 = 0
                Ism1ip1 = 0
            else:
                Ssm1ip1 = Ssi[s-1,i+1]
                Ism1ip1 = Isi[s-1,i+1]
            if i==0 or s+1 == original_shape[0]:
                Ssp1im1 = 0
                Isp1im1 = 0
            else:
                Ssp1im1 = Ssi[s+1,i-1]
                Isp1im1 = Isi[s+1,i-1]
            
            dSsi[s,i] = -tau*i*Ssi[s,i] + gamma*Isi[s,i] \
                        + gamma*((i+1)*Ssm1ip1 - i*Ssi[s,i]) \
                        + tau*ISS*((s+1)*Ssp1im1 - s*Ssi[s,i])/SS
            dIsi[s,i] =  tau*i*Ssi[s,i] - gamma*Isi[s,i] \
                        + gamma*((i+1)*Ism1ip1 - i*Isi[s,i]) \
                        + tau*(ISI/SI + 1)*((s+1)*Isp1im1 - s*Isi[s,i])# 

    dSsi.shape = (original_shape[0]*original_shape[1])
    dIsi.shape = (original_shape[0]*original_shape[1])

    dX = scipy.concatenate((dSsi, dIsi), axis=0)
    return dX


def _dSIR_effective_degree_(X, t, N, original_shape, tau, gamma):
    R = X[-1]
    Ssi = X[:-1]
    Ssi.shape=original_shape
    ISS = sum([sum([i*s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    SS = sum([sum([s*Ssi[s,i] for i in range(original_shape[1])]) 
                for s in range(original_shape[0])])
    
    #commenting out commands for vectorizing this.  
    #I should do this eventually, but not now.  Apply to SIS version as well.
    #ultimately I think right option is to pad with zeros and then 
    #cut out appropriate section.
    #ivec = scipy.array(range(original_shape[1]))
    #ivec.shape = (1,original_shape[1])
    #svec = scipy.array(range(original_shape[0]))
    #svec.shape = (original_shape[0],1)
    #
    #Ssip1 = Ssi[:,1:]
    #scipy.pad(Ssip1, pad_width=((0,0),(0,1)), mode = 'constant', 
    #          constant_values=0)
    #Ssp1im1 = Ssi[1:,:-1]
    #scipy.pad(Ssp1im1, pad_width=((0,1),(1,0)), mode = 'constant', 
    #           constant_values=0)
    #dSsi = - tau* ivec*Ssi + gamma*((ivec+1)*Ssip1 - i*Ssi) 
    #       + tau *ISS*((svec+1)*Sp1im1 - svec*Ssi)/SS

    dSsi = scipy.zeros(original_shape)
    for s in xrange(original_shape[0]):
        for i in xrange(original_shape[1]):
            if i+1 == original_shape[1]:
                Ssip1 = 0
            else:
                Ssip1 = Ssi[s,i+1]
            if s+1 == original_shape[0] or i == 0:
                Ssp1im1=0
            else:
                Ssp1im1 = Ssi[s+1,i-1]            
            dSsi[s,i] = -tau*i*Ssi[s,i] + gamma*((i+1)*Ssip1 - i*Ssi[s,i]) \
                        + tau*ISS*((s+1)*Ssp1im1 - s*Ssi[s,i])/SS
    S = Ssi.sum() 
    I = N-S-R
    dR = gamma*I

    dSsi.shape = (original_shape[0]*original_shape[1])
    dX = scipy.concatenate((dSsi, [dR]), axis=0)
    return dX



def SIS_effective_degree(Ssi0, Isi0, tau, gamma, tmin = 0, tmax=100, 
                            tcount=1001, return_full_data=False):
    '''Encodes system (5.36) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.


    :INPUTS:

    Ssi0 and Isi0 : (square) numpy 2D arrays of same shape.
                      Entries are initial number susceptible or infected 
                      with given initial number of susceptible/infected 
                      neighbors.
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean (default False)
                       tells whether to just return times, S, I, R or 
                       all calculated data.
                       if True, 
                           return times, S, I, Ssi, Isi
                       if False, 
                           return times, S, I
                       
    '''
    times = scipy.linspace(tmin,tmax,tcount) 
    original_shape = Ssi0.shape
    ksq = original_shape[0]*original_shape[1]
    Ssi0.shape = (1,ksq)
    Isi0.shape = (1,ksq)
    
    X0= scipy.concatenate((Ssi0[0], Isi0[0]), axis=0)
    X = integrate.odeint(_dSIS_effective_degree_, X0, times, 
                            args = (original_shape, tau, gamma))
    Ssi = X.T[0:ksq]
    Isi = X.T[ksq:]
    S = Ssi.sum(axis=0)
    I = Isi.sum(axis=0)
    if return_full_data:
        Ssi.shape = (original_shape[0],original_shape[1],tcount)
        Isi.shape = (original_shape[0],original_shape[1],tcount)
        return times, S, I, Ssi, Isi
    else:
        return times, S, I

def SIR_effective_degree(S_si0, I0, R0, tau, gamma, tmin=0, tmax=100, 
                            tcount=1001, return_full_data=False):
    '''Encodes system (5.38) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{S}_{s,i} = - tau i S_{s,i}  + gamma((i+1)S_{s,i+1} - i S_{s,i})
                    + tau [ISS]((s+1)S_{s+1,i-1} - sS_{s,i})/[SS]
    \dot{R} = gamma I
    S = \sum_{s,i} S_{s,i}
    I = N-S-R

    :INPUTS:

    S_si0 : (square) numpy 2-D array
            S_{s,i} at time 0
    I0 : number
         number of infected individuals at time 0
    R0 : number
         number of recovered individuals at time 0
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or 
                       all calculated data.
    '''
    times = scipy.linspace(tmin,tmax, tcount)
    N = S_si0.sum()+I0+R0
    original_shape = S_si0.shape
    S_si0.shape = (original_shape[0]*original_shape[1]) 
    #note this makes it array([[values...]])
    R0=scipy.array([R0])
    R0.shape=(1)
    X0 = scipy.concatenate((S_si0, R0), axis=0)
    X = integrate.odeint(_dSIR_effective_degree_, X0, times, 
                            args = (N, original_shape, tau, gamma))

    R = X.T[-1]
    S_si = X.T[:-1]
    S = S_si.sum(axis=0)
    I = N - R - S
    if return_full_data:
        S_si.shape = (original_shape[0], original_shape[1], tcount)
	return  times, S, I, R, S_si
    else:
        return times, S, I, R


def SIS_effective_degree_from_graph(G, tau, gamma, rho = None, tmin = 0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    S_si0 = scipy.zeros((len(Sk0),len(Sk0)))
    I_si0 = scipy.zeros((len(Sk0),len(Sk0)))

    for s in range(len(Sk0)):
        for i in range(len(Sk0)-s):
            S_si0[s,i] = Sk0[s+i] * scipy.special.binom(s+i,i) \
                            * (rho**i) * (1-rho)**s
            I_si0[s,i] = Ik0[s+i] * scipy.special.binom(s+i,i) \
                            * (rho**i) * (1-rho)**s

    return SIS_effective_degree(S_si0, I_si0, tau, gamma, tmin = tmin, 
                                tmax=tmax, tcount=tcount, 
                                return_full_data=return_full_data)


def SIR_effective_degree_from_graph(G, tau, gamma, rho = None, tmin = 0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    S_si0= scipy.zeros((len(Sk0),len(Sk0)))
    for s in range(len(Sk0)):
        for i in range(len(Sk0)-s):
            S_si0[s,i] = Sk0[s+i] * scipy.special.binom(s+i,i) \
                            * (rho**i) * (1-rho)**s
    I0 = sum(Ik0)
    R0 = sum(Rk0)
    return SIR_effective_degree(S_si0, I0, R0, tau, gamma, tmin=tmin, 
                                tmax=tmax, tcount=tcount, 
                                return_full_data=return_full_data)




#######     COMPACT EFFECTIVE DEGREE


def SIS_compact_effective_degree(Sk0, Ik0, SI0, SS0, II0, tau, gamma, 
                                    tmin = 0, tmax=100, tcount=1001, 
                                    return_full_data=False):
    r'''
    Encodes system (5.44) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    This model is identical to the SIS compact pairwise model, so it 
    simply calls SIS_compact_pairwise()'''

    return SIS_compact_pairwise(Sk0, Ik0, SI0, SS0, II0, tau, gamma, tmin, 
                                    tmax, tcount, return_full_data)

def SIS_compact_effective_degree_from_graph(G, tau, gamma, rho = None, 
                                            tmin = 0, tmax=100, tcount=1001, 
                                            return_full_data=False):
                                            
    return SIS_compact_pairwise_from_graph(G, tau, gamma, rho = rho, 
                                            tmin = tmin, tmax=tmax, 
                                            tcount=tcount, 
                                            return_full_data=return_full_data)


def _dSIR_compact_effective_degree_(X, t, N, tau, gamma):
    Skappa = X[:-2]
    R, SI = X[-2:]
    I = N- R- Skappa.sum()
    kappas = scipy.arange(len(Skappa))
    effectiveI = float(SI) /Skappa.dot(kappas)
    dSkappa = effectiveI*(-(tau+gamma)*kappas*Skappa \
                + gamma*shift(kappas*Skappa,-1))
    dSI = -(tau+gamma)*SI \
            + tau*(effectiveI-2*effectiveI**2)*sum(kappas*(kappas-1)*Skappa)

    dR = gamma*I
    dX = scipy.concatenate((dSkappa, [dR, dSI]), axis=0) 
    return dX
    
def SIR_compact_effective_degree(Skappa0, I0, R0, SI0, tau, gamma, tmin=0, 
                                    tmax=100, tcount=1001, 
                                    return_full_data=False):
    '''
    Encodes system (5.43) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    \dot{S}_kappa = <I> [-(tau+gamma) kappa S_kappa 
                         + gamma(kappa+1)S_{kappa+1}
    [\dot{SI}] = -(tau+gamma)[SI] 
                + tau(<I> - 2 <I>^2) sum_{kappa} kappa(kappa-1) S_kappa
    \dot{R} = gamma I
    <I> = [SI]/sum_kappa kappa S_kappa
    S = sum_kappa S_kappa
    I = N - S - R

    :INPUTS:

    Skappa0 : scipy array
              from S_0(0) up to S_kappamax(0) of number susceptible with 
                  each effective degree
              Skappa = number of nodes that are susceptible and have 
                       kappa non-recovered neighbors
    I0 : number
         number of infected individuals at time 0
    R0 : number
         initial number recovered
    SI0 : number
          initial number of SI edges
    tau : number
          transmission rate
    gamma : number
            recovery rate
    tmin : number (default 0)
           minimum report time
    tmax : number (default 100)
           maximum report time 
    tcount : integer (default 1001)
             number of reports
    return_full_data : boolean
                       tells whether to just return times, S, I, R or 
                       all calculated data.
    
    RETURNS (if return_full_data==False):
    --------
    times : scipy.array of times
    S : scipy.array of number susceptible
    I : scipy.array of number infected
    R : scipy.array of number recovered

    if return_full_data==True
    --------------------------
    times : as before
    Skappa : array of arrays, each subarray gives particular S_kappa
    I : number infected
    R : number recovered
    SI : number of SI edges
    '''

    N = Skappa0.sum()+I0+R0
    X0= scipy.concatenate((Skappa0,[R0,SI0]), axis=0)
    times = scipy.linspace(tmin, tmax, tcount)
    X = integrate.odeint(_dSIR_compact_effective_degree_, X0, times, 
                            args = (N, tau, gamma))
    Skappa = X.T[:-2]
    S = Skappa.sum(axis=0)
    R, SI = X.T[-2:]
    I = N - S -R
    if return_full_data:
        return times, S, I, R, Skappa, SI
    else:
	return times, S, I, R

def SIR_compact_effective_degree_from_graph(G, tau, gamma, rho = None, 
                                            tmin = 0, tmax=100, tcount=1001, 
                                            return_full_data=False):
    Nk, Sk0, Ik0, Rk0 = get_Nk_and_IC_as_arrays(G, rho, SIR=True)
    Skappa0 = Sk0
    I0 = sum(Nk)*rho
    R0 = 0
    SI0 = sum([k*Skappa0[k]*rho for k in range(len(Sk0))])
    return SIR_compact_effective_degree(Skappa0, I0, R0, SI0, tau, gamma, 
                                        tmin=tmin, tmax=tmax, tcount=tcount, 
                                        return_full_data=return_full_data)






#######################################
#    EBCM and other related results   #
#######################################

def Epi_Prob_discrete(Pk, p, number_its = 100):
    '''Encodes System (6.2) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    :INPUTS:

    Pk : scipy array Pk[k] is probability a node has degree k.

    p : transmission probability

    number_its : number of iterations before assumed converged.
                 default value is 100

    :RETURNS:

    Calculated Epidemic probability (assuming configuration model)
    '''
    ks = scipy.arange(len(Pk))
    def psi(x):
        return Pk.dot(x**ks)
    def psiPrime(x):
        return (ks*Pk).dot(x**(ks-1))
    
    alpha = 1-p
    k_ave = psiPrime(1)
    for counter in range(number_its):
        alpha = 1-p +p *psiPrime(alpha)/k_ave
    return 1- psi(alpha)


def Epi_Prob_cts_time(Pk, tau, gamma, umin=0, umax = 10, ucount = 1001, 
                        number_its = 100):
    r'''Encodes System (6.3) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    The equations are rescaled by setting $u=\gamma T$.  Then it becomes

    P = 1- \int_0^\infty \psi(\alpha(u/\gamma)) e^{-u} du
    alpha_d(u/\gamma) = 1- p(u/\gamma)
                        + p(u/\gamma)
                          \int_0^\infty 
                                (\psiPrime(\alpha(\hat{u}/\gamma))/<K>)
                                e^{-u}du

    where p(u/\gamma) = 1 - e^{-\tau u/\gamma}

    Define 
        \hat{p}(u) = p(u/\gamma), and \hat{\alpha}(u) = \alpha(u/\gamma)
    and then drop hats to get

    P = 1-\int_0^\infty \psi(\alpha(u)) e^{-u} du
    \alpha(u) = 1-p(u) + p(u) 
                         \int_0^\infty 
                              (\psiPrime(\alpha(u))/<K>)e^{-u} du

    with initial guess 
        \alpha_1(u) = e^{-\tau u/\gamma} 
    and 
        p(u) = 1-e^{-\tau u/\gamma}
    
    :INPUTS:

    Pk : scipy array Pk[k] is probability a node has degree k.

    tau : transmission rate

    gamma : recovery rate

    umin : minimal value of \gamma T used in calculation
    umax : maximum value of \gamma T used in calculation
    ucount : number of points taken for integral.
             So this integrates from umin to umax using simple Riemann 
             sum.

    number_its : number of iterations before assumed converged.
                 default value is 100

    :RETURNS:

    Calculated Epidemic probability (assuming configuration model)
    '''
    ks = scipy.arange(len(Pk))
    def psi(x):
        return Pk.dot(x**ks)
    def psiPrime(x):
        ks*Pk
        x**(ks-1)
        return (ks*Pk).dot(x**(ks-1))
    us = scipy.linspace(umin, umax, ucount) 
    alpha = scipy.e**(-tau*us/gamma)  #initial guess for alpha(u)
    p = 1- scipy.e**(-tau*us/gamma)    #initial guess for p(u)
    exp_neg_u = scipy.e**(-us)        #e^{-u}
    for counter in xrange(number_its):
        alpha = 1 - p + p* (psiPrime(alpha)/kave).dot(exp_neg_u)/(ucount-1.)
    return 1 - psi(alpha).dot(exp_neg_u)/(ucount - 1)


def Epi_Prob_non_Markovian(Pk, tau, gamma, Pxidxi, po, umin=0, umax = 10, 
                            ucount = 1001, number_its = 100):
    r'''Encodes system (6.5) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    
    :INPUTS:

    Pk : scipy array Pk[k] is probability a node has degree k.

    tau : transmission rate

    gamma : recovery rate

    Pxidxi : a dict.  Returns P(xi)dxi for user-selected xi.  The 
             algorithm will replace the integral with
             \sum_{xi \in Pxidxi.keys()} \psi(\alpha(xi)) Pxidxi(xi)

    po : a function.
         returns p_o(xi), the probability a node will transmit to a 
         random neighbor given xi.
    umin : minimal value of \gamma T used in calculation
    umax : maximum value of \gamma T used in calculation
    ucount : number of points taken for integral.
             So this integrates from umin to umax using simple Riemann 
             sum.

    number_its : number of iterations before assumed converged.
                 default value is 100

    :RETURNS:

    Calculated Epidemic probability (assuming configuration model)
    '''
    ks = scipy.arange(len(Pk))
    def psi(x):
        return Pk.dot(x**ks)
    def psiPrime(x):
        return (ks*Pk).dot(x**(ks-1))
    xis = Pxidxi.keys()
    alpha = {xi: 1-po(xi) for xi in xis}
    for counter in xrange(number_its):
        newalpha = {}
        for xi in xis:
            newalpha[xi] = 1 - po(xi)  \
                            + po(xi)*sum(psiPrime(alpha[xihat])*Pxidxi(xihat) 
                                                for xihat in xis)/kave
        alpha = newalpha
    return 1 - sum(psi(alpha[xi])*Pxidxi[xi] for xi in xis)

def Attack_rate_discrete(Pk, p, number_its=100, rho = None, Sk0=None, 
                            phiS0=None, phiR0=0):
    r'''
    Encodes systems (6.6) and (6.10) of Kiss, Miller, & Simon.  Please 
    cite the book if using this algorithm.

    To use system (6.6), leave rho and Sk0 as None.

    :INPUTS:

    Pk : dict
        Pk[k] is the probability a randomly selected node has degree k.
    tau : number
        per-edge transmission rate.
    gamma : number
        per-node recovery rate
    number_its : The solution is found iteratively, so this determines 
                 the number of iterations.
    rho : Number (default None)
    Sk0 : dict (default None)
          only one of rho and Sk0 can be defined.  
          The other (or both) should remain None.
          rho gives the fraction of nodes randomly infected.
          Sk0 is a dict such that Sk0[k] is the probability that a 
              degree k node is susceptible at start.
    phiS0 : number (default None)
          Should only be used if Sk0 is not None.  
          If it is None, then assumes that initial introduction is 
          randomly introduced
    phiR0 : number (default 0)
          As with phiS0, only used if Sk0 is not None.
    :RETURNS:

    A : number
        the predicted fraction infected.
    '''

    if rho is not None and Sk0 is not None:
        raise EoNError("at most one of rho and Sk0 can be defined")
    if Sk0 is None:
        if rho is None or rho == 0:
            return Epi_Prob_discrete(Pk, p, number_its)
        else:
            Sk0 = {k: (1-rho) for k in Pk.keys()}
    def psihat(x):
        return sum(Pk[k]*Sk0[k]*x**k for k in Pk.keys())
    def psihatPrime(x):
        return sum(k*Pk[k]*Sk0[k]*x**(k-1) for k in Pk.keys())

    if phiS0 == None:
        phiS0 = psihatPrime(1)/sum(k*Pk[k] for k in Pk.keys())
    if phiR0 == None:
        phiR0 = 0

    theta  = 1
    for counter in xrange(number_its):
        theta = 1-p + p*(phiR0 +  phiS0*psihatPrime(theta)/psihatPrime(1))
    return 1 - psihat(theta)

#def Attack_rate_discrete_from_graph(G, p, number_its = 100, rho = None, 
#                                        Sk0 = None):
#    Pk = get_Pks(G)
#    return Attack_rate_discrete(Pk, p, number_its = number_its, 
#                                            rho = rho, Sk0 = Sk0)

def Attack_rate_cts_time(Pk, tau, gamma, number_its =100, rho = None, 
                            Sk0 = None, phiS0=None, phiR0=0):
    #tested in test_SIR_final_sizes
    '''Encodes system (6.7) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    
    This system predicts the fraction of nodes infected if an epidemic 
    occurs in a Configuration Model network assuming a continuous-time 
    Markovian SIR disease.  
    
    This gives the limit of the attack rate of epidemics as the initial 
    fraction infected approaches 0.

    If we look for the limit of a nonzero initial fraction infected, we 
    introduce rho or Sk0

    :INPUTS:
        
    Pk : dict
    
    the probability a randomly selected node has degree k.

    tau : number
    per-edge transmission rate.

    gamma : number
    per-node recovery rate

    number_its : int
    The solution is found iteratively, so this determines 
    the number of iterations.

    rho : number, optional
    The initial proportion infected (defaults to None)

    Sk0 : dict (default None)
    only one of rho and Sk0 can be defined.  
    The other (or both) should remain None.
    rho gives the fraction of nodes randomly infected.
    Sk0 is a dict such that Sk0[k] is the probability that a 
    degree k node is susceptible at start.

    :RETURNS:

    A : number
        the predicted fraction infected.
    '''

    if rho is not None and Sk0 is not None:
        raise EoNError("at most one of rho and Sk0 can be defined")
    if Sk0 is None:
        if rho is None:
            rho = 0
        Sk0 = {k: (1-rho) for k in Pk.keys()}
    def psihat(x):
        return sum(Pk[k]*Sk0[k]*x**k for k in Pk.keys())
    def psihatPrime(x):
        return sum(k*Pk[k]*Sk0[k]*x**(k-1) for k in Pk.keys())

    if phiS0 == None:
        phiS0 = psihatPrime(1)/sum(k*Pk[k] for k in Pk.keys())
    if phiR0 == None:
        phiR0 = 0

    kave = sum(Pk[k]*k for k in Pk.keys())
    omega = gamma/(gamma+tau)
    for counter in range(number_its):
        print omega
        omega = gamma/(gamma+tau) \
                + tau*phiS0*psihatPrime(omega)/(psihatPrime(1)*(gamma+tau)) \
                + tau*phiR0/(gamma+tau)
    return 1 - psihat(omega)

def Attack_rate_cts_time_from_graph(G,  tau, gamma, number_its =100, rho=None,
                                    Sk0 = None):
    r'''
    Given a graph, predicts the attack rate for Configuration Model 
    networks with the given degree distribution.  This does not account 
    for any structure in G beyond degree distribution.
    
    First calculates the degree distribution and then calls 
    Attack_rate_cts_time.
    
    SEE ALSO
    estimate_SIR_prob_size(G, p) - accounts for entire structure of G
    '''
    Pk = get_Pk(G)
    return Attack_rate_cts_time(Pk, tau, gamma, number_its = number_its, 
                                rho=rho, Sk0=Sk0)
    
def Attack_rate_non_Markovian():
    r'''
    Encodes system (6.8) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.
    '''
    pass

def EBCM_discrete(N, psihat, psihatPrime, p, phiS0, phiR0=0, R0=0, tmax = 100,
                    return_full_data = False):
    #tested in test_basic_discrete_SIR_epidemic
    r'''
    Encodes system (6.11) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    theta(t) = (1-p) + p(phi_R(0) 
               + phi_S(0) psihatPrime(theta(t-1))/psihatPrime(1))
    R(t) = R(t-1) + I(t-1)
    S(t) = N psihat(theta(t))
    I(t) = N-S-R

    :INPUTS:

    N : number
        size of population
    psihat : function
        psihat(x) = \sum_k S(k,0) x^k
    psihatPrime : function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    p : number
        per edge transmission probability
    phiS0 : number
            initial proportion of edges (of susceptible nodes) 
            connecting to susceptible nodes
    phiR0 : number
            initial proportion of edges (of susceptible nodes) 
            connecting to recovered nodes
    R0 : number
        number of recovered nodes at time 0
    tmax : number
        maximum time
    return_full_data : boolean
                       if False, 
                           return t, S, I, R
                       if True 
                           return t, S, I, R, and theta

    :RETURNS:

    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
    '''
    times = [0]
    theta = [1]
    R = [R0]
    S = [N*psihat(1)]
    I = [N-S[-1]-R[-1]]

    for time in range(1,tmax+1):
        times.append(time)
        newtheta = (1-p) + p *(phiR0 
                                + phiS0*psihatPrime(theta[-1])/psihatPrime(1))
        newR = R[-1]+I[-1]
        newS = N*psihat(newtheta)
        newI = N-newR-newS
        theta.append(newtheta)
        R.append(newR)
        S.append(newS)
        I.append(newI)
    if not return_full_data:
        return scipy.array(times), scipy.array(S), scipy.array(I), \
                    scipy.array(R)
    else:
        return scipy.array(times), scipy.array(S), scipy.array(I), \
                    scipy.array(R), scipy.array(theta)

def EBCM_discrete_uniform_introduction(N, psi, psiPrime, p, rho, tmax=100, 
                                        return_full_data=False):
    #tested in test_basic_discrete_SIR_epidemic
    '''
    Handles the case that the disease is introduced uniformly as opposed
    to depending on degree.

    :INPUTS:

    N : Number
        number of nodes
    psi : function
        psi(x) = \sum P(k) x^k
    psiPrime : function
        psiPrime(x)=d psi(x)/dx = \sum kP(k) x^{k-1}
    p : number
        per edge transmission probability
    rho : number
        initial proportion of infected nodes
    tmax : number
        maximum time
    return_full-data : boolean
        if False, 
            return t, S, I, R
        if True 
            return t, S, I, R, and theta

    :RETURNS:

    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
    '''
    def psihat(x):
        return (1-rho)*psi(x)
    def psihatPrime(x):
        return (1-rho)*psiPrime(x)
    
    return EBCM_discrete(N, psihat, psihatPrime, p, 1-rho, tmax=tmax, 
                            return_full_data=return_full_data)


def EBCM_discrete_from_graph(G, p, rho = None, tmin = 0, tmax=100, 
                                tcount=1001, return_full_data=False):
    #tested in test_basic_discrete_SIR_epidemic
    '''
    Takes a given graph, finds the degree distribution 
    (from which it gets psi),
    assumes a constant proportion of the population is infected at time 
    0, 
    and then uses the discrete EBCM model.
    
    :INPUTS:

    G : Networkx Graph
    p : number
        per edge transmission probability
    rho : number
        initial proportion of infected nodes
    tmax : number
        maximum time
    return_full-data : boolean
        if False, 
            return t, S, I, R and if True return t, S, I, R, and theta

    :RETURNS:

    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
 '''
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    def psi(x):
        return sum(Pk[k]*x**k for k in Pk)
    def psiPrime(x):
        return sum(k*Pk[k]*x**(k-1) for k in Pk)
    return EBCM_discrete_uniform_introduction(G.order(), psi, psiPrime, p, 
                                                rho, tmax=tmax, 
                                                return_full_data=False)


def _dEBCM_(X, t, N, tau, gamma, psihat, psihatPrime, phiS0, phiR0):
    theta = X[0]
    R = X[1]
    
    dtheta = -tau*theta + tau*phiS0*psihatPrime(theta)/psihatPrime(1) \
                + gamma*(1-theta) + tau*phiR0

    S = N*psihat(theta)
    I = N-S-R
    dR = gamma*I
    return scipy.array([dtheta, dR])
    
def EBCM(N, psihat, psihatPrime, tau, gamma, phiS0, phiR0=0, R0=0, tmin=0, 
            tmax=100, tcount=1001, return_full_data=False):
    '''
    Encodes system (6.12) of Kiss, Miller, & Simon.  Please cite the
    book if using this algorithm.

    note : R0 is R(0), not the reproductive number


    :INPUTS:

    N : number
        size of population
    psihat : function
             psihat(x) = \sum_k S(k,0) x^k
    psihatPrime : function
               psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    tau : number
        per edge transmission rate
    gamma : number
        per node recovery rate
    phiS0 : number
            initial proportion of edges (of susceptible nodes) 
            connecting to susceptible nodes
    phiR0 : number
            initial proportion of edges (of susceptible nodes) 
            connecting to recovered nodes
    R0 : number
        number of recovered nodes at time 0
    tmin : number
        start time
    tmax : number
        stop time
    tcount : integer
        number of distinct times to calculate
    return_full_data : boolean
        if False, 
            return t, S, I, R
        if True 
            return t, S, I, R, and theta

    :RETURNS:

    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
    '''
    times = scipy.linspace(tmin, tmax, tcount)
    X0 = scipy.array([1, R0])
    X = integrate.odeint(_dEBCM_, X0, times, 
                            args = (N, tau, gamma, psihat, psihatPrime, phiS0, 
                                        phiR0))
    theta = X.T[0]
    R = X.T[1]
    S = N*psihat(theta)
    I = N-S-R
    if not return_full_data:
        return times, S, I, R
    else:
        return times, S, I, R, theta

def EBCM_uniform_introduction(N, psi, psiPrime, tau, gamma, rho, tmin=0, 
                                tmax=100, tcount=1001, 
                                return_full_data=False):
    r'''
    Handles the case that the disease is introduced uniformly as opposed
    to depending on degree.

    :INPUTS:

    N : number
        size of population
    psi : function
        psihat(x) = \sum_k S(k,0) x^k
    psiPrime : function
        psihatPrime(x) = d psihat(x)/dx = sum_k k S(k,0) x^{k-1}
    tau : number
        per edge transmission rate
    gamma : number
        per node recovery rate
    rho : number
        initial proportion infected
    tmin : number
        start time
    tmax : number
        stop time
    tcount : integer
        number of distinct times to calculate
    return_full_data : boolean
        if False,
            return t, S, I, R 
        if True 
            return t, S, I, R, and theta

    :RETURNS:

    if return_full_data == False:
        returns t, S, I, R, all scipy arrays
    if ...== True
        returns t, S, I, R and theta, all scipy arrays
'''
    def psihat(x):
        return (1-rho)*psi(x)
    def psihatPrime(x):
        return (1-rho)*psiPrime(x)
    
    return EBCM(N, psihat, psihatPrime, tau, gamma, 1-rho, tmin=tmin, 
                tmax=tmax, tcount=tcount, return_full_data=return_full_data)

def EBCM_from_graph(G, tau, gamma, rho = None, tmin = 0, tmax=100, 
                        tcount=1001, return_full_data=False):
    #tested in test_SIR_dynamics
    if rho is None:
        rho = 1./G.order()
    Pk = get_Pk(G)
    def psi(x):
        return sum(Pk[k]*x**k for k in Pk)
    def psiPrime(x):
        return sum(k*Pk[k]*x**(k-1) for k in Pk)
    return EBCM_uniform_introduction(G.order(), psi, psiPrime, tau, gamma, 
                                        rho, tmax=tmax, 
                                        return_full_data=return_full_data)

    
def EBCM_deg_corr():
    r'''Nothing to see here 
    - just a place holder until this code is written'''
    pass
    
def EBCM_deg_corr_discrete():
    r'''Nothing to see here 
    - just a place holder until this code is written'''
    pass
    
def EBCM_deg_corr_from_graph(G):
    r'''Nothing to see here 
    - just a place holder until this code is written'''
    pass

def EBCM_deg_corr_discrete_from_graph():
    r'''Nothing to see here 
    - just a place holder until this code is written'''
    pass
