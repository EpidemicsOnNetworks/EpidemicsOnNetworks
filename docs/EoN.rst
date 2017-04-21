EoN module
==========

EoN
---
.. automodule:: EoN
    :show-inheritance:
    
Simulation Toolkit
------------------
This submodule deals with epidemic simulation.  The most notable of the algorithms are:

- Epidemic dynamics
    - fast_SIR  (a very efficient event-based algorithm)
    - fast_SIS  (a very efficient event-based algorithm)
    - basic_discrete_SIR_epidemic
    - 
    - Gillespie_SIR
    - Gillespie_SIS
- Percolation-based approaches (for SIR disease)
    - estimate_SIR_prob_size
    - estimate_directed_SIR_prob_size
- visualize

.. automodule:: EoN.simulation
    :members:
    :undoc-members:
    :show-inheritance:

Analytic Toolkit
----------------
This submodule deals with solution to systems of equations appearing in the book.
The majority of these also have a version that take a graph G.  There are 
additional functions that calculate properties which these need.

- Chapter 3
    - 3.7 SIS_individual_based and SIS_individual_based_pure_IC
    - 3.26 SIS_pair_based
    - 3.30 SIR_individual_based and SIR_individual_based_pure_IC
    - 3.39 SIR_pair_based
- Chapter 4
    - 4.8 SIS_homogeneous_meanfield
    - 4.9 SIR_homogeneous_meanfield
    - 4.10 SIS_homogeneous_pairwise
    - 4.11 SIR_homogeneous_pairwise
- Chapter 5
    - 5.10 SIS_heterogeneous_meanfield
    - 5.11 SIR_heterogeneous_meanfield
    - 5.13 SIS_heterogeneous_pairwise
    - 5.15 SIR_heterogeneous_pairwise
    - 5.18 SIS_compact_pairwise
    - 5.19 SIR_compact_pairwise
    - 5.20 SIS_super_compact_pairwise
    - 5.22 SIR_super_compact_pairwise
    - 5.36 SIS_effective_degree
    - 5.38 SIR_effective_degree
    - 5.43 SIR_compact_effective_degree
    - 5.44 SIS_compact_effective_degree
- Chapter 6
    - 6.2 Epi_Prob_discrete
    - 6.3 Epi_Prob_cts_time
    - 6.5 Epi_Prob_non_Markovian
    - 6.6 Attack_rate_discrete
    - 6.7 Attack_rate_cts_time
    - 6.8 Attack_rate_non_Markovian
    - 6.10 Attack_rate_discrete
    - 6.11 EBCM_discrete
    - 6.12 estimate_SIR_prob_size

- exercise 6.21  EBCM_pref_mix and EBCM_pref_mix_discrete

.. automodule:: EoN.analytic
    :members:
    :undoc-members:
    :show-inheritance:
