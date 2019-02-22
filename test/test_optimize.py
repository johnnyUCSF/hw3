from hw3 import io
from hw3 import smith_waterman
from hw3 import optimize
import pytest
import os

def test_optimize():
    #####import blosum matrix
    blosums = []
    blosum_files = ['BLOSUM50']
    for file in blosum_files:
        importme = 'blosums/'+file
        blosums.append(io.import_blosum(importme))
    ###################
    scoring_matrix = blosums[0]
    FPRs = [0.0,0.1,0.2,0.3]
    pen_gap_extend = -5
    pen_gap_open = -6
    true_pos = io.import_pairs('data/Pospairs.txt')
    true_neg = io.import_pairs('data/Negpairs.txt')
    #######particle swarm optimization
    particles = 2
    iterations = 1
    ###################
    optimized_output = optimize.particle_swarm(scoring_matrix,FPRs,pen_gap_extend,pen_gap_open,true_pos,true_neg,particles,iterations)
    assert optimized_output[0] != 2.12



    
    