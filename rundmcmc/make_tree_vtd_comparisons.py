# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 14:14:52 2018

@author: MGGG
"""

import MCMC_partitions
import tree_tools
import networkx as nx

G = nx.grid_2d_graph(3,3)
for sample_size in [1000]:
    M = MCMC_partitions.maps(G, num_districts, m)
    M.district_list = M.district_maker() 
    T = tree_tools.random_spanning_tree(G)
    e = tree_tools.best_edge_for_equipartition(G,T)[0]
    partition = tree_tools.R(G,T,e)
    M.district_list = partition
    M.districts = M.district_list
    M.set_node_district_flags()
    M.initialize_boundary()
    M.set_boundary()
    found_districts = []
    
    rejections = 0
    for i in range(sample_size):
        current = M.district_list
        M.metropolis_hastings()
        found_districts.append(M.district_list)