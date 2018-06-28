# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 11:15:33 2018

@author: MGGG
"""

    
######Converting data types:
    
#These s

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
