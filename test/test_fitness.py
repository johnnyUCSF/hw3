from hw3 import io
from hw3 import smith_waterman
from hw3 import optimize
import pytest
import os

def test_calc_fitness():
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
    ###see what the fitness is prior to starting
    out = optimize.calc_fitness(FPRs,blosums[0],true_pos,true_neg,pen_gap_open,pen_gap_extend)
    assert out == 2.12



    
    