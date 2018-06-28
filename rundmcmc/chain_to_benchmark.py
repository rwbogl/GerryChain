# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 10:57:56 2018

@author: MGGG
"""

###All the tools to take the output of chain and make it useful:

from grid import Grid

from rundmcmc.chain import MarkovChain
from rundmcmc.proposals import propose_random_flip
from rundmcmc.validity import Validator, contiguous, no_vanishing_districts
from rundmcmc.accept import always_accept
import numpy as np
#import rundmcmc.benchmark_distances
import metagraph_distances as md

def dict_invert(dictionary):
    dict = {val: [key for key in dictionary.keys() if dictionary[key] == val] for val in dictionary.values()}
    return dict

def dictionary_to_sets(dictionary):
    #print(dictionary[1])
    return set([frozenset(list) for list in dict_invert(dictionary).values()])

def dictionary_list_to_sets_list(dict_list):
    return [dictionary_to_sets(x) for x in dict_list]

def chain_to_partitions(graph, iterations, skip_standing = False):
    is_valid = Validator([contiguous,no_vanishing_districts])
    chain = MarkovChain(propose_random_flip, is_valid, always_accept, graph, total_steps = iterations)
    partitions = []
  
    for step in chain:
        if step.flips:
            if skip_standing:
                if not (list(step.flips.keys())[0][0] == list(step.flips.keys())[0][1]):
                    partitions.append(step.assignment)
            else:
                partitions.append(step.assignment)
    return partitions

def test_distance(plan1, plan2):
    return 1

def build_distance_matrix(partitions, function):
    n = len(partitions)
    M = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            M[i][j] = function(partitions[i], partitions[j])
    return M

def partitions_to_distance(partition_list, distance_function):
    #Partition_list is a list of set partitions
    if type(partition_list[0]) == set:
        plan = partition_list
    if type(partition_list[0]) == frozenset:
        plan = partition_list
    if type(partition_list[0]) == dict:
        plan = dictionary_list_to_sets_list(partition_list)

    return build_distance_matrix(plan,distance_function)

def main():
    grid = Grid((10,10))
    print(grid)
    partitions = chain_to_partitions(grid, 10)
    D = partitions_to_distance(partitions, md.stupid_hamming_distance)
    print(D)

if __name__ == "__main__":
    main()
