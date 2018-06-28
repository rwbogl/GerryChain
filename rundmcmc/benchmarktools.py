# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 10:41:37 2018

@author: MGGG
"""

from itertools import product
import numpy as np
import numpy as np

from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.decomposition import PCA

from benchmark_tests import chain_test
from scipy.stats import gaussian_kde

from treetools import test
import treetools
from matplotlib import pyplot as plt

import benchmark_tests
import metagraph_distances as md
import chain_to_benchmark as cb

def common_refinement(d1,d2):
    if set(d1.keys()) != set(d2.keys()):
        return "Keys do not match!"
    
    keys = d1.keys()

    d = dict()
    i = 0
    coord_to_values_dict = dict()

    for value1, value2 in product(set(d1.values()), set(d2.values()) ):
        coord_to_values_dict[(value1,value2)] = i
        i += 1 


    for node in keys:
        d[node] = (d1[node], d2[node])

    translated_dict = {node: coord_to_values_dict[d[node]] for node in keys}
    
    return translated_dict

def partition_entropy(d):
    """Returns the entropy of a partition"""    
    total = len(d)
    
    prob = {}
    
    for value in list(d.values()):
        prob[value] = sum(1 for x in d.keys() if d[x] == value)/total
        
    H = sum(- prob[key] * np.log(prob[key]) for key in prob.keys())
    return H


def mutual_information(dX,dY):
    d = common_refinement(dX,dY)
    HX = partition_entropy(dX)
    HY = partition_entropy(dY)
    HXY = partition_entropy(d)
    
    I = HX + HY - HXY
    return I
    
def mi_metric(d1,d2, normalised = False):
    d = common_refinement(d1,d2)
    H = partition_entropy(d)
    
    
    I = mutual_information(d1,d2)
    
    if normalised is True:
        return (H-I)/H
    else:
        return H-I
    
def lp_distance(p1, p2, p):
    if set(p1.keys()) != set(p2.keys()):
        return "Keys do not match!"
    
    keys = p1.keys()
    difference_vector = []
    for k in keys:
            difference_vector.append( p1[k] - p2[k])
    return np.linalg.norm(difference_vector, p)
        
    
def build_distance_matrix(partitions, metric = "information", p = 1):
    #Partition should be stored as dictionaries 
    n = len(partitions)
    M = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if metric == "information":
                M[i][j] = mi_metric(partitions[i], partitions[j])
            if metric == "Lp":
                #p = 1 is hamming
                M[i][j] = lp_distance(partitions[i], partitions[j],p)
    return M
 
    
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


############################################################

##Multi-Dimensional scaling

def make_mds_histogram(A, M_A, h1, dim = 2):
    #A is the ground set of partitions
    #M_A is the distance matrix
    #h1 is the histogram (labeled by A_node_lists -- below) that tells us the heat map
    color_list = []
    A_node_lists = [ set([frozenset(g.nodes()) for g in x]) for x in A]
    for x in A_node_lists:
        color_list.append(h1[str(x)])
    z = gaussian_kde(color_list)(color_list)
    
    similariaties = np.exp(M_A)
    seed = np.random.RandomState(seed=3)
    
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                       dissimilarity="precomputed", n_jobs=1)
    
    mds.fit(similariaties)
    pos = mds.embedding_
    # Rotate the data
    clf = PCA(n_components=2)    
    pos = clf.fit_transform(pos)
    
    s = 100
    
    from matplotlib import pyplot as plt
    plt.scatter(pos[:, 0], pos[:, 1], c=z, s=s, lw=0)    
    plt.show()
    

def make_mds_path(A, M_A, path, dim = 2):
    #A is the ground set of partitions
    #M_A is the distance matrix
    #path is path of indices, organized ot correspond to the indices of hte columns of M_A
    similariaties = np.exp(M_A)
    seed = np.random.RandomState(seed=3)
    
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                       dissimilarity="precomputed", n_jobs=1)
    
    mds.fit(similariaties)
    pos = mds.embedding_
    # Rotate the data
    clf = PCA(n_components=2)    
    pos = clf.fit_transform(pos)
    
    s = 100
    
    plt.scatter(pos[:, 0], pos[:, 1], s=s, lw=0)
    
    for t in range(len(path) -1):
        #print( pos[path[t]], pos[path[t+1]])
        plt.plot( [pos[path[t]][0], pos[path[t+1]][0]],
                 [pos[path[t]][1], pos[path[t+1]][1]], 'r', lw = M_A[ path[t]][ path[t+1]]  )
    plt.show()
    
def stress_test(A, M_A, h1, low_dim = 2, high_dim = 10):
    stress_list = []
    similariaties = M_A
    seed = np.random.RandomState(seed=3)
    for k in range(low_dim, high_dim):
        mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                           dissimilarity="precomputed", n_jobs=1)        
        mds.fit(similariaties)
        print(mds.stress_)
        stress_list.append(mds.stress_)
    from matplotlib import pyplot as plt
    plt.plot(stress_list)
    plt.show()
    
def testmds_histogram():
    #To see the stress function -- chcek out stress test. It doesn't seem to improve by increaseing dimension.
    #Probably this means that this space cannot be easily embedded. Maybe try hamming?
    from treetools import test
    h1, A, partitions = test([3,3], 4,300)
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    make_mds_histogram(A, M_A, h1)
   # stress_test(A, M_A, h1, 2, 25)
   
def testmds_boundary_histogram():
    h1, A, partitions = chain_test((3,3), 4,300)
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    #spectral_plot(A, M_A, h1, 3)
    make_mds_histogram(A, M_A, h1)
    
def testmds_walk():
    from treetools import test
    h1, A, partitions = test([3,3], 4,300)
    path_indices = make_path_indices(treetools.subgraph_to_node(A), treetools.subgraph_to_node(partitions))
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    make_mds_path(A, M_A, path_indices)
    
def testmds_boundary_walk():
    from treetools import test
    h1, A, partitions = chain_test([3,3], 4,300)
    path_indices = make_path_indices(treetools.subgraph_to_node(A), benchmark_tests.dictionary_list_to_node_set(partitions))
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    make_mds_path(A, M_A, path_indices)
    
def make_path_indices(ground_set, path):
    #Takes a walk over an indexed ground set (list), and returns the list of indices...
    indices = []
    for x in path:
        indices.append(ground_set.index(x))
    return indices
    
##SHould the conclusion from this be that the mds embedding is so bad that it's pointless ot use it? 
#testmds()
    
##########################################################
#Spectral Embeddings

def spectral_plot(A, M_A, h1, n = 2):
    
    color_list = []
    A_node_lists = [ set([frozenset(g.nodes()) for g in x]) for x in A]
    for x in A_node_lists:
        color_list.append(h1[str(x)])
    z = gaussian_kde(color_list)(color_list)
    
    se = manifold.SpectralEmbedding(n_components = 3, affinity = 'precomputed')
    Y = se.fit_transform(M_A)
    if n == 2:
        import matplotlib.pyplot as plt
        plt.scatter(Y[:, 0], Y[:, 1], c=z)
        plt.show()
    if n == 3:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(Y[:,0], Y[:,1], Y[:,2],  c = z)
        
def spectral_plot_walk(A, M_A, path, n = 2):
    
    color_list = []
    A_node_lists = [ set([frozenset(g.nodes()) for g in x]) for x in A]
    for x in A_node_lists:
        color_list.append(h1[str(x)])
    z = gaussian_kde(color_list)(color_list)
    
    se = manifold.SpectralEmbedding(n_components = 3, affinity = 'precomputed')
    Y = se.fit_transform(M_A)
    if n == 2:
        import matplotlib.pyplot as plt
        plt.scatter(Y[:, 0], Y[:, 1], c=z)
        plt.show()
        pos = Y
        for t in range(len(path) -1):
            #print( pos[path[t]], pos[path[t+1]])
            plt.plot( [pos[path[t]][0], pos[path[t+1]][0]],
                     [pos[path[t]][1], pos[path[t+1]][1]], 'r', lw = M_A[ path[t]][ path[t+1]]  )
    if n == 3:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(Y[:,0], Y[:,1], Y[:,2],  c = z)
        pos = Y
        for t in range(len(path) -1):
            #print( pos[path[t]], pos[path[t+1]])
            ax.plot( [pos[path[t]][0], pos[path[t+1]][0]],
                     [pos[path[t]][1], pos[path[t+1]][1]], [pos[path[t]][2], pos[path[t+1]][2]], 'r', lw = M_A[ path[t]][ path[t+1]]  )


def test_spectral():
    #it would be great to have a tool that lets us see which partitions are which by hovering over them
    h1, A, partitions = test([3,3], 4,1000)
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    #spectral_plot(A, M_A, h1, 3)
    spectral_plot(A, M_A, h1, 3)
    
def test_spectral_walk():
    from treetools import test
    h1, A, partitions = test([3,3], 4,300)
    path_indices = make_path_indices(treetools.subgraph_to_node(A), treetools.subgraph_to_node(partitions))
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    spectral_plot_walk(A, M_A, path_indices)
    
def test_spectral_walk_boundary():
    from treetools import test
    h1, A, partitions = chain_test([3,3], 4,300)
    path_indices = make_path_indices(treetools.subgraph_to_node(A), benchmark_tests.dictionary_list_to_node_set(partitions))
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    spectral_plot_walk(A, M_A, path_indices,3)
    #non-empty??
    #How is the support changing as the number of steps increases..
    #Flagify the complex
    
def test_spectral_boundary():
    h1, A, parts = chain_test((3,3),4, 1000)
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    #spectral_plot(A, M_A, h1, 3)
    spectral_plot(A, M_A, h1, 3)
    
def compare():
    #Need to start hard coding the set of partitions!
    h1, A, partitions_1 = chain_test((3,3),4, 200)
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    spectral_plot(A, M_A, h1, 3) 
    
    h2, A, partitions_2 = test([3,3], 4,200)
    dlist_A = partition_list_to_dictionary_list(A)
    spectral_plot(A, M_A, h1, 3)
    S1 = support(h1)
    S2 = support(h2)
    print(len(A), len(S1), len(S2), len(symdif(S1, S2)))

##Size of support tools:
    
def support(histogram):
    support = []
    for x in histogram.keys():
        if histogram[x] != 0:
            support.append(x)
    return support

def symdif(set1, set2):
    ourset = set([])
    for x in set1:
        if x not in set2:
            ourset.add(x)
    for x in set2:
        if x not in set1:
            ourset.add(x)
    return ourset
#test_spectral_boundary()
#test_spectral() 
    
#########CLUSTERING:
    
##Spectral partitioning
    


#K medioids
    
###Distances
    

###Comparing step sizes...
def compare_jump_histograms():
    h1, A, partitions = test([3,3], 4,300)
    path_indices_tree = make_path_indices(treetools.subgraph_to_node(A), treetools.subgraph_to_node(partitions))
    dlist_A = partition_list_to_dictionary_list(A)
    M_A = build_distance_matrix(dlist_A)
    h1, B, partitions = chain_test([3,3], 4,300)
    path_indices_boundary = make_path_indices(treetools.subgraph_to_node(A), benchmark_tests.dictionary_list_to_node_set(partitions))
    boundary_step_sizes = step_length_tally(M_A, path_indices_boundary)
    tree_step_sizes = step_length_tally(M_A, path_indices_tree)
    tree_cleaned = [x for x in tree_step_sizes if x != 0]
    boundary_cleaned = [x for x in boundary_step_sizes if x != 0]
    plt.hist(tree_cleaned)
    plt.hist(boundary_cleaned)
    
def step_length_tally(matrix, path):
    step_sizes = []
    for t in range(len(path) - 1):
        step_sizes.append( matrix[path[t]][path[t+1]])
    return step_sizes

#But maybe tree takes a bunch of big steps around the same general spot as boundary
###Big cluster vs. many bulbs with thin paths connecting them...

#So try computing matrix of pairwise distances of some subsample of the founds ... but this needs to be done
    #on a MUCH larger space. in 3x3 all the partitions are enumerated and its small of large numbers
    #it doesn't require enumeration tough, s o
    
#
def blobs():
    tree_partitions = treetools.subgraph_to_node(treetools.tree_walk([20,20], 4, 3000, False))
    tree_partitions_cleaned = list(set([frozenset(x) for x in tree_partitions]))
    print(len(tree_partitions_cleaned))
    boundary_partitions = benchmark_tests.dictionary_list_to_node_set(benchmark_tests.chain_walk((20,20), 4, 3000))
    boundary_partitions_cleaned = list(set([frozenset(x) for x in boundary_partitions]))
    print("building")
    D = cb.partitions_to_distance(tree_partitions_cleaned, md.shared_information_distance)
    print("next")
    E = cb.partitions_to_distance(boundary_partitions_cleaned, md.shared_information_distance)
    print("q")
    plt.hist(D.flatten(),color = 'r')
    plt.hist(E.flatten(),color = 'b')
