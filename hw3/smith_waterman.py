import numpy as np
from Bio import SeqIO
########
def smith_waterman(seq1,seq2,scoring_matrix,pen_gap_open,pen_gap_extend):
    #### first call a function to score the alignment matrix
    final_scores, final_dirs = calc_matrix(seq1,seq2,scoring_matrix,pen_gap_open,pen_gap_extend)
    #### then dynamically find the best path backwards, using the max num as the start point
    return(find_shortest(seq1,seq2,final_scores,final_dirs))

def calc_matrix(seq1,seq2,scoring_matrix,pen_gap_open,pen_gap_extend):
    #####final_scores contains the alignment scores for each cell
    final_scores = []
    #####final_dirs = diagonal if match (1), up/down (2), left/right(3), no direction (0)
    final_dirs = []
    ####use this to identify the location in the scoring matrix
    identities = list(scoring_matrix[0])
    ####initialize
    for i in range(len(seq1)+1):
        row = []
        row_dirs = []
        for j in range(len(seq2)+1):
            row_dirs.append(0)
            if j == 0 or i == 0:
                row.append(0)
            else:
                row.append(-100)
        final_scores.append(row)
        final_dirs.append(row_dirs)
    ####loop through again and calculate scores; indices are for scoring matrix
    for i in range(len(seq1)+1):
        ####
        row = []
        for j in range(len(seq2)+1):
            if j == 0 or i == 0:
                continue
            else:
                seq1_char = seq1[i-1]
                seq2_char = seq2[j-1]
                ###first find index using char in sequence; this index used to lookup in blosum
                ######check to see if letter exists, if not assign to *
                if seq1_char not in identities:
                    seq1_char = '*'
                ###
                if seq2_char not in identities:
                    seq2_char = '*'
                ###
                blosum_seq1 = identities.index(seq1_char)+1
                blosum_seq2 = identities.index(seq2_char)
                blosum_val = scoring_matrix[blosum_seq1][blosum_seq2]
                ###get old scores
                diag = final_scores[i-1][j-1]
                topright = final_scores[i-1][j]
                botleft = final_scores[i][j-1]
                ###define new scores
                new_diag = diag + blosum_val
                if final_dirs[i-1][j] == 2:
                    ###if going in the same direction, then use the extension penalty
                    new_topright = topright + pen_gap_extend
                    ###otherwise opening penalty is used
                else:
                    new_topright = topright + pen_gap_open
                if final_dirs[i][j-1] == 3:
                    new_botleft = botleft + pen_gap_extend
                else:
                    new_botleft = botleft + pen_gap_open
                ####evaluate which direction presents largest score
                if new_diag > new_botleft and new_diag > new_topright:
                    newscore = new_diag
                    newdir = 1
                else:
                    if new_botleft > new_topright:
                        newscore = new_botleft
                        newdir = 3
                    else:
                        newscore = new_topright
                        newdir = 2
                #####allows for local alignment; feature of smith waterman vs needleman wunsch
                if newscore < 0:
                    newscore = 0
                #####write in score and direction
                final_scores[i][j] = newscore
                final_dirs[i][j] = newdir
    ######
    return(final_scores,final_dirs)

def find_shortest(seq1,seq2,final_scores,final_dirs):
    #####
    seq1_align = ''
    seq2_align = ''
    #####
    max_score, i, j = find_max_indices(final_scores)
    #####dynamic programming
    while final_dirs[i][j] > 0:
        if final_dirs[i][j] == 1:
            seq1_align += seq1[i-1]
            seq2_align += seq2[j-1]
            i -= 1
            j -= 1
        else: 
            if final_dirs[i][j] == 2:
                seq1_align += seq1[i-1]
                seq2_align += '-'
                i -= 1
            if final_dirs[i][j] == 3:
                seq1_align += '-'
                seq2_align += seq2[j-1]
                j -= 1
    #####
    seq1_align = seq1_align[::-1]
    seq2_align = seq2_align[::-1]
    #####
    return seq1_align,seq2_align,max_score

def find_max_indices(final_scores):
    ####
    max_val = -100000
    max_i = 0
    max_j = 0
    ####
    for i in range(len(final_scores)):
        for j in range(len(final_scores[i])):
            if final_scores[i][j] > max_val:
                max_val = final_scores[i][j]
                max_i = i
                max_j = j
    ####
    return max_val, max_i, max_j