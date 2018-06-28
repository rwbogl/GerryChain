# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 10:41:37 2018

@author: Assaf, Lorenzo, Anthony
"""
from itertools import product
import numpy as np
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

#from sklearn import manifold
#from sklearn.metrics import euclidean_distances
#from sklearn.decomposition import PCA


#from scipy.stats import gaussian_kde

def common_refinement(partition1, partition2):
    refinement = {part1.intersection(part2) for part1 in partition1 for part2 in partition2}
    refinement.remove(frozenset())
    return refinement

def partition_entropy(partition):
    """Returns the entropy of a partition"""
    number_of_units = len(set().union(*partition))
    area_ratios = {len(part)/number_of_units for part in partition}
    return sum([- prob * np.log(prob) for prob in area_ratios])

def shared_information_distance(partition1, partition2):
  return(2*partition_entropy(common_refinement(partition1,partition2)) - partition_entropy(partition1) - partition_entropy(partition2))
  
def stupid_hamming_distance(partition1, partition2):
    for part1 in partition1:
        stupid_hamming += min([part2.symmetric_difference(part1) for part2 in partition2])
    for part2 in partition2:
        stupid_hamming += min([part1.symmetric_difference(part2) for part1 in partition1])
  return(stupid_hamming)

def subgraph_list_to_dictionary(subgraph):
    m = len(subgraph)
    node_lists = [g.nodes() for g in subgraph]
    dictionary = {}
    for i in range(m):
        for x in node_lists[i]:
            dictionary[x] = i
    return dictionary

def partition_list_to_dictionary_list(partitions):
    dictionary_list = []
    for x in partitions:
        dictionary_list.append(subgraph_list_to_dictionary(x))
    return dictionary_list

def lp_distance(p1, p2, p):
    if set(p1.keys()) != set(p2.keys()):
        return "Keys do not match!"
    
    keys = p1.keys()
    difference_vector = []
    for k in keys:
            difference_vector.append( p1[k] - p2[k])
    return np.linalg.norm(difference_vector, p)

def build_distance_matrix(partitions):
    #Partition should be stored as dictionaries 
    n = len(partitions)
    M = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            M[i][j] = shared_information_distance(partitions[i], partitions[j])
            #M[i][j] = stupid_hamming_distance(partitions[i], partitions[j])
    return M
