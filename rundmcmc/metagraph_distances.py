# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 10:41:37 2018

@author: Assaf, Lorenzo, Anthony
"""
import numpy as np
partition_bad = 1
partition_bad = 1

def common_refinement(partition1, partition2):
    global partition_bad1, partition_bad2
    refinement = {part1.intersection(part2) for part1 in partition1 for part2 in partition2}
    if frozenset() in refinement:
        refinement.remove(frozenset())
        #Note that if each part in partitoin1 intersect each part in partition2
        #then there is never the emptyset in the intersection.
    #print(len(set().union(*refinement)))
    #print([len(set) for set in refinement])
    return refinement

def partition_entropy(partition):
    """Returns the entropy of a partition"""
    number_of_units = len(set().union(*partition))
    #print(number_of_units)
    area_ratios = [len(part)/number_of_units for part in partition]
    #print(area_ratios)
    #return sum([-(prob * np.log(prob)) for prob in area_ratios])
    return -np.dot(area_ratios,np.log(area_ratios))

def partition_sqent(partition):
    """Returns the entropy of a partition"""
    number_of_units = len(set().union(*partition))
    #print(number_of_units)
    area_ratios = [len(part)/number_of_units for part in partition]
    #print(area_ratios)
    #return sum([-(prob * np.log(prob)) for prob in area_ratios])
    return -np.dot(area_ratios,np.sqrt(area_ratios))

def shared_information_distance(partition1, partition2):
  return(2*partition_entropy(common_refinement(partition1,partition2)) - partition_entropy(partition1) - partition_entropy(partition2))
  
def shared_sqent_distance(partition1, partition2):
  return(2*partition_sqent(common_refinement(partition1,partition2)) - partition_sqent(partition1) - partition_sqent(partition2))

def stupid_hamming_distance(partition1, partition2):
    stupid_hamming=0
    for part1 in partition1:
        stupid_hamming += min([len(part2.symmetric_difference(part1)) for part2 in partition2])
    for part2 in partition2:
        stupid_hamming += min([len(part1.symmetric_difference(part2)) for part1 in partition1])
    return(stupid_hamming)
