from hw3 import io
from hw3 import smith_waterman
from hw3 import optimize
import pytest
import os

def test_smith_waterman():
    #####import unique blosum matrices
    blosums = []
    blosum_files = ['BLOSUM50']
    for file in blosum_files:
        importme = 'blosums/'+file
        scoring_matrix = (io.import_blosum(importme))
    ####import true positive and negative data
    true_pos = io.import_pairs('data/Pospairs.txt')
    ####define parameters
    pen_gap_extend = -5
    pen_gap_open = -6
    output = smith_waterman.smith_waterman(true_pos[0][0],true_pos[0][1],scoring_matrix,pen_gap_open,pen_gap_extend)
    out_score = output[2]

    assert out_score == 65.0




    
    