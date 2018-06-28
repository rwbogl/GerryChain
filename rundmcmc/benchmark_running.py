
from rundmcmc.grid import Grid

from rundmcmc.chain import MarkovChain
from rundmcmc.proposals import propose_random_flip
from rundmcmc.validity import Validator, contiguous, no_vanishing_districts
from rundmcmc.accept import always_accept
import numpy as np

def read_chain(graph, iterations):
  is_valid = Validator([contiguous,no_vanishing_districts])
  chain = MarkovChain(propose_random_flip, is_valid, always_accept, graph, total_steps = iterations)
  partitions = []
  
  
  for step in chain: 
    #print('parent = ')
    #print(step.parent.assignment)
    #print('current assignment')
    #print(step.assignment)
    #if step.flips:
      #if not (list(step.flips.keys())[0][0] == list(step.flips.keys())[0][1]):
  #  print(step)
    partitions.append(step.assignment)
    #print('Keys')
    #print(step.flips)
    #print(list(step.flips.keys())[0])
    #print(list(step.flips.keys())[0][0] == list(step.flips.keys())[0][1])

  #print(partitions)  
  newlist = [dict(s) for s in set(frozenset(d.items()) for d in partitions)]
  return partitions
  print(newlist)
  #distance_matrix = bmt.build_distance_matrix(newlist)
  return distance_matrix

def main():
  grid = Grid((3,3))
  print(grid)
  distance_matrix = read_chain(grid, 10)
  print(distance_matrix)
