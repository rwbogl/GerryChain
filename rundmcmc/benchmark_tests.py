# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 14:19:58 2018

@author: MGGG
"""

from rundmcmc.grid import Grid

from rundmcmc.chain import MarkovChain
from rundmcmc.proposals import propose_random_flip
from rundmcmc.validity import Validator, contiguous, no_vanishing_districts
from rundmcmc.accept import always_accept
import numpy as np
import networkx as nx


from benchmark_running import read_chain
'''
making historgram from a list of partitions 
'''
def count(x, visited_partitions):

    x_lens = np.sort([len(k) for k in x])
    count = 0
    for sample_nodes in visited_partitions:
        sample_lens = np.sort([len(k) for k in sample_nodes])
        #if (x_lens == sample_lens).all():
        if np.array_equal(x_lens , sample_lens):
            if x == sample_nodes:
                count += 1
    return count

def make_histogram(A, visited_partitions):
    #A is the list of all partitions,
    A_node_lists = [ set([frozenset(g.nodes()) for g in x]) for x in A]
    dictionary = {}
    for x in A_node_lists:
        dictionary[str(x)] = count(x,visited_partitions) / len(visited_partitions)
    return dictionary

####################
'''
These are functions that translate the between the different kinds of outputs
'''
def dictionary_to_nodes(dictionary):
    #Takes dictionary stored partition and returns the list of nodes format
    node_lists = []
    keylist = dictionary.keys()
    for i in set(dictionary.values()):
        node_list = []
        for x in keylist:
            if dictionary[x] == i:
                node_list.append(x)
        node_lists.append(frozenset(node_list))
    return set(node_lists)

def dictionary_list_to_node_set(dictionary_list):
    list_of_partitions = []
    for x in dictionary_list:
        list_of_partitions.append(dictionary_to_nodes(x))
    return list_of_partitions

#####################

def chain_test(grid_size, k_part, steps = 100, equi = False):
    from naive_graph_partitions import k_connected_graph_partitions
    k_part = 4
    grid = Grid(grid_size)
    G = grid.graph
    A = list(k_connected_graph_partitions(G, k_part))
    partitions = read_chain(grid, steps, equi)
    
    histogram = make_histogram(A, dictionary_list_to_node_set(partitions))
    
    total_variation = 0
    for k in histogram.keys():
        total_variation += np.abs( histogram[k] - 1 / len(A))
    print("total variation", total_variation)
    return [histogram, A, partitions]

def chain_walk(grid_size, k_part, steps = 100, equi = False):
    k_part = 4
    grid = Grid(grid_size)
    partitions = read_chain(grid, steps, equi)
    return partitions
#new = chain_test(3,3)
    
#chain_test((3,3),4)
