# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 12:01:21 2018

@author: MGGG
"""

import networkx as nx
import random
import itertools
import numpy as np
import scipy.linalg
from scipy.sparse import csc_matrix
import scipy
from scipy import array, linalg, dot
import sys
#from naive_graph_partitions import k_connected_graph_partitions

######Tree counting

def log_number_trees(G, weight = False):
    #Kirkoffs is the determinant of the minor..
    #at some point this should be replaced with a Cholesky decomposition based algorithm, which is supposedly faster. 
    if weight == False:
        m = nx.laplacian_matrix(G)[1:,1:]
    if weight == True:
        m = nx.laplacian_matrix(G, weight = "weight")[1:,1:]
    m = csc_matrix(m)
    splumatrix = scipy.sparse.linalg.splu(m)
    diag_L = np.diag(splumatrix.L.A)
    diag_U = np.diag(splumatrix.U.A)
    S_log_L = [np.log(np.abs(s)) for s in diag_L]
    S_log_U = [np.log(np.abs(s)) for s in diag_U]
    LU_prod = np.sum(S_log_U) + np.sum(S_log_L)
    return  LU_prod


def score_tree_edges_pair(G,T,e):
    partition = R(G,T,e)
    tree_term = np.sum([log_number_trees(g) for g in partition])

    #Use the weighted version of Kirkoffs theorem-- weight each edge by the number of cut edges
    #I.e. build the adjacency graph of the partition, weight by cut size
    #This counts the number of ways to extend spanning trees on the subgraphs of the parittion
    #individually, to a spanning tree on all of G
    connector_graph = nx.Graph()
    connector_graph.add_nodes_from(partition)
    for x in partition:
        for y in partition:
            if x != y:
                cutedges = cut_edges(G, x,y)
                if cutedges != []:
                    connector_graph.add_edge(x,y, weight = len(cutedges))
    cut_weight = log_number_trees(connector_graph, True)
    return -1 * (tree_term + cut_weight)
    #YOU NEED THIS -1 -- the score is the inverse! Don't take it away!
    
#####For creating a spanning tree

def srw(G,a):
    wet = set([a])
    trip = [a]
    while len(wet) < len(G.nodes()):
        b = random.choice(list(G.neighbors(a)))
        wet.add(b)
        trip.append(b)
        a = b
    return trip

def forward_tree(G,a):
    walk = srw(G,a)
    edges = []
    for x in G.nodes():
        if (x != walk[0]):
            t = walk.index(x)
            edges.append( [walk[t], walk[t-1]])
    return edges

def random_spanning_tree(G):
    #It's going to be faster to use the David Wilson algorithm here instead.
    root = random.choice(list(G.nodes()))
    T_edges = forward_tree(G, root)
    T = nx.DiGraph()
    T.add_nodes_from(list(G.nodes()))
    T.add_edges_from(T_edges)

    return T

def random_spanning_tree_wilson(G):
    #The David Wilson random spanning tree algorithm
    
    return T
#####For lifting:

def cut_edges(G, G_A,G_B):
    #Finds the edges in G from G_A to G_B
    edges_of_G = list(G.edges())

    list_of_cut_edges = []
    for e in edges_of_G:
        if e[0] in G_A and e[1] in G_B:
            list_of_cut_edges.append(e)
        if e[0] in G_B and e[1] in G_A:
            list_of_cut_edges.append(e)
    return list_of_cut_edges

def random_lift(G, subgraphs):
    number_of_parts = len(subgraphs)
    trees = [random_spanning_tree(g) for g in subgraphs]
    
    #This builds a graph with nodes the subgraph, and they are connected
    #if there is an edge connecting the two subgraphs
    #and each edge gets 'choices' = to all the edges in G that connect the two subgraphs
    connector_graph = nx.Graph()
    connector_graph.add_nodes_from(subgraphs)
    for x in subgraphs:
        for y in subgraphs:
            if x != y:
                cutedges = cut_edges(G, x,y)
                if cutedges != []:
                    connector_graph.add_edge(x,y, choices = cutedges)
                    #need to worry about directendess!!???
                    
                    
    connector_meta_tree = random_spanning_tree(connector_graph)
    connector_tree = nx.Graph()
    for e in connector_meta_tree.edges():
        w = random.choice(connector_graph[e[0]][e[1]]['choices'])
        connector_tree.add_edges_from([w])
        
    
    T = nx.Graph()
    for x in trees:
        T.add_edges_from(x.edges())
    T.add_edges_from(connector_tree.edges())
    e = random.sample(list(T.edges()),number_of_parts - 1)
    return [T,e]

######Projection tools:
    
def R(G,T,edge_list):
    T.remove_edges_from(edge_list)
    components = list(nx.weakly_connected_components(T))
    T.add_edges_from(edge_list)
    subgraphs = [nx.induced_subgraph(G, A) for A in components]
    return subgraphs

def R_all(G,T,n):
    T_edges = set(T.edges())
    partitions = []

    for e in itertools.combinations(T_edges, n):
        partitions.append(R(G,T,list(e)))
    return partitions

def R_sample(G,T,n,m):
    T_edges = set(T.edges())
    partitions = []
    
    iteration = random.sample(list(itertools.combinations(T_edges, n)), m)
    
    for e in iteration:
        partitions.append(R(G,T,list(e)))
    return partitions

    

def best_edge_for_equipartition(G,T):
    list_of_edges = list(T.edges())
    best = 0
    candidate = 0
    for e in list_of_edges:
        score = equi_score_tree_edge_pair(G,T,e)
        if score > best:
            best = score
            candidate = e
    return [candidate, best]

def equi_score_tree_edge_pair(G,T,e):
    T.remove_edges_from([e])
    components = list(nx.connected_components(T))
    T.add_edges_from([e])
    A = len(components[0])
    B = len(components[1])
    x =  np.min([A / (A + B), B / (A + B)])
    return x

###Metropolis-Hastings tools
    
def propose_step(G,T):
    T_edges = list(T.edges())
    T_edges_t = [ tuple((e[1], e[0])) for e in T_edges]
    #Because of stupid stuff in networkx
    A = [e for e in G.edges() if e not in T_edges and e not in T_edges_t]
    e = random.choice(A)
    T.add_edges_from([e])
    C = nx.find_cycle(T, orientation = 'ignore')
    w = random.choice(C)
    U = nx.Graph()
    U.add_edges_from(list(T.edges()))
    U.remove_edges_from([w])
    T.remove_edges_from([e])
#    print(len(U.edges()))
#    print(U.edges())
    return U
    

def MH_step(G, T,e, equi = False, MH = True):
    n = len(e)
    U = propose_step(G,T)
    if equi == False:
        e2 = random.sample(list(U.edges()),n)
    if equi == True:
        e2 = best_edge_for_equipartition(G,U)[0]
    if MH == True:
        current_score = score_tree_edges_pair(G,T,e)
        new_score = score_tree_edges_pair(G, U, e2)
        if new_score > current_score:
            return [U,e2]
        else:
           p = np.exp(new_score - current_score)
           a = np.random.uniform(0,1)
           if a < p:
               return [U,e2]
           else:
               return [T,e]
    if MH == False:
        return [U,e2]

        
########Validation -- 
            
       
###Histogram creation tools

def count(x, visited_partitions):

    x_lens = np.sort([len(k) for k in x])
    count = 0
    for sample_nodes in visited_partitions:
        #sample_nodes = set([frozenset(g.nodes()) for g in i])
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
        
def subgraph_to_node(visited_partitions):
    partition_list  = []
    for partitions in visited_partitions:
        partition_list.append(set([frozenset(g.nodes()) for g in partitions]))
        
    return partition_list
##########################
    
def test(grid_size, k_part, steps = 100, equi = False, MH = True):
    from naive_graph_partitions import k_connected_graph_partitions
    #k_part = 3 nnum partitions
    G = nx.grid_graph(grid_size)
    A = list(k_connected_graph_partitions(G, k_part))
    
    ##Tree walks:
    T = random_spanning_tree(G)
    e = list(T.edges())[0:k_part - 1]
    visited_partitions = []
    for i in range(steps):
        new = MH_step(G, T, e, equi, MH)
        #This is the step that takes in the graph G, the spanning tree T, 
        #and the list of edges e that we are currently planning to delete.
        T = new[0]
        e = new[1]
        visited_partitions.append(R(G,T,e))
        
        #R(emoval) is the function that returns the partition you get by removing e from T
    ###
    
    ##Statistics from output from tree tools:
    visited_partitions_node_format = subgraph_to_node(visited_partitions)
    histogram = make_histogram(A, visited_partitions_node_format)
    total_variation = 0
    for k in histogram.keys():
        total_variation += np.abs( histogram[k] - 1 / len(A))
    print("total variation", total_variation)
    return [histogram, A, visited_partitions]
     
def TV(p,q):
#test comment
    total_variation = 0
    for k in p.keys():
        total_variation += np.abs(p[k] - q[k])
    return total_variation
#h1, A, partitions = test([2,3], 3)
    
def tree_walk(grid_size, k_part, steps = 100, MH = True, how_many = 'one', demand_equi = False):
    G = nx.grid_graph(grid_size)
    ##Tree walks:
    T = random_spanning_tree(G)
    e = list(T.edges())[0:k_part - 1]
    visited_partitions = []
    for i in range(steps):
        if demand_equi == True:
            #new = Equi_Step(G,T,e, False, MH)
            print("You haven't set this up yet!")
        if demand_equi == False:
            new = MH_step(G, T, e, False, MH)
        #This is the step that takes in the graph G, the spanning tree T, 
        #and the list of edges e that we are currently planning to delete.
        T = new[0]
        e = new[1]
        if how_many == 1:
            visited_partitions.append(R(G,T))
        if how_many == 'all':
            T = random_spanning_tree(G)
            visited_partitions += R_all(G,T, k_part)
        if (how_many != 1) and (type(how_many) == int):
            T = random_spanning_tree(G)
            visited_partitions += R_sample(G,T, k_part, how_many)
        
    return visited_partitions

def random_equi_partition_trees(graph, k_part, number = 100):
    equi_partitions = [] 
    counter = 0
    while len(equi_partitions) < number:
        counter += 1
        T = random_spanning_tree(graph)
        e = equi_partition(T, k_part)
        if e != False:
            print(len(equi_partitions), "waiting time:", counter, file=sys.stderr)
            equi_partitions.append( R(graph, T, e) )
            counter = 0
    return equi_partitions

def equi_partition(T, num_blocks):
    #Returns the partition if T can give an equi partition in num_blocks,
    #Else return false

    #Currently this is hard coded for 4 partitions -- but there shold be a good
    #and scalable algorithm
    if num_blocks == 4:
        e = equi_split(T)
        if e == False:
            return False
        if e != False:
            T.remove_edges_from([e])
            components = list(nx.weakly_connected_component_subgraphs(T))
            T.add_edges_from([e])
            e1 = equi_split(components[0])
            if e1 == False:
                return False
            e2 = equi_split(components[1])
            if e2 == False:
                return False
    else:
        print("you didn't set up functionality for more than 4 partitions yet!")
    return [e, e1, e2]


def equi_split(T):
    # Returns the partition if T can be split evenly in two
    # Else returns False
    label_weights(T)
    edge, weight = choose_best_weight(T)

    if weight == len(T) // 2:
        return edge

    return False


def label_weights(graph):
    """Label nodes of of a directed, rooted tree by their weights.

    :graph: NetworkX DiGraph.
    :returns: Nothing.

    The "weight" of a node is the size of the subtree rooted at itself.

    """
    def _label_weights(node):
        in_edges = graph.in_edges(node)

        if not in_edges:
            graph.nodes[node]["weight"] = 1
            return 1

        child_weights = [_label_weights(child) for child, _ in in_edges]

        weight = sum(child_weights) + 1
        graph.nodes[node]["weight"] = weight

        return weight

    root = [node for node in graph if graph.out_degree(node) == 0][0]
    _label_weights(root)


def choose_best_weight(graph):
    """Choose out edge from node with weight closest to n_nodes / 2.

    :graph: NetworkX Graph labeled by :func:`~label_weights`.
    :returns: Tuple (edge, weight).

    """

    best_node = None
    best_difference = float("inf")

    for node in graph:
        diff = abs(len(graph) / 2 - graph.nodes[node]["weight"])
        if diff < best_difference:
            best_node = node
            best_difference = diff

    edge = list(graph.edges(best_node))[0]
    weight = graph.nodes[best_node]["weight"]

    return (edge, weight)
