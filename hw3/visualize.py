from hw3 import io
from hw3 import smith_waterman
from hw3 import optimize
import copy
import matplotlib.pyplot as plt
#################

def byFirst(elem):
    return elem[0]

###########
def calc_ROC(true_pos, true_neg, total_scores):
    ####init 
    TPR = []
    FPR = []
    AUC = []
    stepsize = float(1)/float(len(true_neg))
    AUC_counter = 0.0
    #####
    for FP_limit in range(len(true_neg)+1):
        TP = 0.0
        FP = 0.0
        P = float(len(true_pos))
        N = float(len(true_neg))
        ####sort biggest to smallest
        total_scores.sort(key=byFirst,reverse=True)
        ####scroll through and count TPs and FPs
        if FP_limit == 0:
            for alignment in total_scores:
                if alignment[1] in true_pos:
                    TP += 1
                else:
                    break
        else:
            for alignment in total_scores:
                if FP >= FP_limit:
                    break
                if alignment[1] in true_pos:
                    TP += 1
                else:
                    FP += 1
        ####save the FPR, TPR, AUC
        TPR.append(TP/P)
        FPR.append(FP/N)
        AUC_counter += float(stepsize) * float(TP/P)
        AUC.append(AUC_counter)
        ####
    return(FPR,TPR,AUC)
##################
def plot_ROCs(true_pos,true_neg, scores_array,legend_array):
    ####
    i = 0
    plt.figure(figsize=(10,10))
    plt.plot([0, 1], [0, 1],'r--')
    ####
    for scores in scores_array:
        FPR, TPR, AUC = calc_ROC(true_pos,true_neg,scores)
        #####
        printme = legend_array[i] + ' AUC: ' + str(AUC[-1])
        #####
        plt.plot(FPR,TPR,label=printme)
        ##
        i+=1
    ####
    plt.legend(loc=0)
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.show()
########
def plot_ROC(true_pos,true_neg, scores,legend):
    FPR, TPR, AUC = calc_ROC(true_pos,true_neg,scores)
    #####
    printme = legend + ' AUC: ' + str(AUC[-1])
    #####
    plt.figure(figsize=(10,10))
    plt.plot(FPR,TPR,label=printme)
    plt.plot([0, 1], [0, 1],'r--')
    plt.legend(loc=0)
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.show()
