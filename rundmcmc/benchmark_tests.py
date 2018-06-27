# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 14:19:58 2018

@author: MGGG
"""

import benchmarktools as bmt
from rundmcmc.grid import Grid

from rundmcmc.chain import MarkovChain
from rundmcmc.proposals import propose_random_flip
from rundmcmc.validity import Validator, contiguous
from rundmcmc.accept import always_accept
import numpy as np

import benchmarktools

def count(x, visited_partitions):

    x_lens = np.sort([len(k) for k in x])
    count = 0
    for i in visited_partitions:
        sample_nodes = set([frozenset(g.nodes()) for g in i])
        sample_lens = np.sort([len(k) for k in sample_nodes])
        #if (x_lens == sample_lens).all():
        if np.array_equal(x_lens , sample_lens):
            if x == sample_nodes:
                count += 1
    return count

def make_histogram(A, visited_partitions):
    A_node_lists = [ set([frozenset(g.nodes()) for g in x]) for x in A]
    dictionary = {}
    for x in A_node_lists:
        dictionary[str(x)] = count(x,visited_partitions) / len(visited_partitions)
    return dictionary

def chain_test(grid_size, k_part, steps = 100):
    from naive_graph_partitions import k_connected_graph_partitions
    #k_part = 3
    G = nx.grid_graph(grid_size)
    A = list(k_connected_graph_partitions(G, k_part))
   
    is_valid = Validator([contiguous])
    chain = MarkovChain(propose_random_flip, is_valid, always_accept, graph, total_steps = iterations)
    partitions = []
  for step in chain: 
    #print('parent = ')
    #print(step.parent.assignment)
    #print('current assignment')
    #print(step.assignment)
    #if step.flips:
      #if not (list(step.flips.keys())[0][0] == list(step.flips.keys())[0][1]):
    partitions.append(step.assignment)
    histogram = make_histogram(A, visited_partitions)
    total_variation = 0
    for k in histogram.keys():
        total_variation += np.abs( histogram[k] - 1 / len(A))
    print("total variation", total_variation)
    return [histogram, A]
     