#!/usr/bin/env python3
import random as rnd
import numpy as np
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
import graph_tool as gt
from graph_tool.all import *
import pandas as pd
import time
import scipy.special
import scipy.stats

#skripta sadrži fuknciju koja radi histogram od pvalues svih napravljenih grafova
#sadrži funkciju za histograme svih podataka (varijabli reshuffleanih i početnih grafova)
#sadrži funkciju za određivanje pvalue jednog grafa

def ks_test(array, unshuffled):
    return scipy.stats.ks_2samp(array,unshuffled)

def pval_histo(path, array, name, zeta):
    plt.hist(array, bins = 100)
    plt.title(r' $\zeta$ = %.3f' %(zeta))
    plt.xlabel(name)
    plt.ylabel("number of occurences")
    plt.savefig("{0}/{1}".format(path,name))
    #plt.show()
    plt.gcf().clear()
    plt.close()

def zstat_histo(path, avg_array, array, unshuffled, name, zeta):
    avg_array = [np.round(np.asarray(array).mean(axis = 2),2),np.round(np.asarray(array).std(axis = 2),2)]
    np.seterr(divide='ignore', invalid='ignore')
    zstat = np.divide(np.asarray(unshuffled).flatten() - avg_array[0].flatten(),avg_array[1].flatten())
    zstat = [x for x in zstat if ~np.isnan(x)]
    zstat = [x for x in zstat if ~np.isinf(x)]
    plt.hist(zstat, bins = 100)
    plt.title(r' $\zeta$ = %.2f, $\mu$ = %.2f, $\sigma$ = %.2f, kurt = %.2f' %(zeta, np.mean(zstat), np.std(zstat), scipy.stats.kurtosis(zstat)))
    plt.xlabel(name)
    plt.ylabel("number of occurences")
    plt.savefig("{0}/{1}".format(path,name))
#    plt.show()
    plt.gcf().clear()
    plt.close()

def all_data_hist(path, array, name, unshuffled, exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, zeta):
    #array = [x for x in array if x != 0 and x != 1]
    #unshuffled = [x for x in unshuffled if x != 0 and x != 1]
    array_max = max(array)
    array_min = min(array)
    unshuffled_max = max(unshuffled)
    unshuffled_min = min(unshuffled)
    bin_num = array_max-array_min
    unsh_bin_num = unshuffled_max-unshuffled_min
#    print(array_max,array_min,unshuffled_max,unshuffled_min)
    fig, ax1 = plt.subplots()
    color = 'tab:blue'
    ax1.set_xlabel(name)
    ax1.set_ylabel('number of x occurences', color=color)
    ax1.hist(array, bins = max(bin_num,30), alpha = 0.7, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('number of x_p occurences', color=color) 
    ax2.hist(unshuffled, bins = max(unsh_bin_num,30), alpha = 0.7, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.title(r'N_vertices = %d N_graphs = %d N_processes = %d N_rnd = %d' %(N_vertices, N_random_graphs, N_processes,N_random_shuffles))
    fig.tight_layout()  
    plt.savefig("{0}/{1}.{2}".format(path, name,"png"))
    #plt.show()
    plt.gcf().clear()
    plt.close()
    out_to_tex(array,unshuffled,"{0}/{1}.{2}".format(path, name,"txt"))

def out_to_png(graph_on, array, name, unshuffled, exp_deg, N_vertices, N_random_shuffles, zeta):
    pval = 0
    histogram = plt.hist(array, bins = 100)
    plt.title(r'<k> = %d N_vertices = %d iter = %d $\zeta$ = %.3f' %(exp_deg, N_vertices, N_random_shuffles, zeta))
    plt.xlabel(name)
    plt.ylabel("number of occurences")
    surface = np.sum(histogram[0])
    plt.axvline(x = unshuffled,color = 'r')
    limit = 0
    for index,element in enumerate(histogram[1]):
        if(unshuffled <= element):
            limit = index - 11
            break
        limit = index
    part = np.sum(histogram[0][limit:])
    pval =part/surface
    if(graph_on == True):
        plt.savefig(name)
        #plt.show()
    plt.gcf().clear()
    plt.close()
    return pval

def out_to_tex_stat(path, array, name, stat, time_stops):
    text_file = open("{0}/{1},{2}".format(path, stat,name), "w")
    l = 0
    for part in array:
        text_file.write("network percentage = {0} %\n".format(round(100*(l+1)/time_stops,2)))
        text_file.write("[")
        for p in part:  
            #text_file.write('p_mean = {0}, sigma_p = {1}\n'.format(np.round(p[0],5),np.round(p[1],5))
            text_file.write("[")
            for n in p:
                text_file.write("{0} ".format(n))
            text_file.write("]")
        text_file.write("]\n")
        l += 1
    text_file.close()

def out_to_tex(array, unshuffled, name):
    text_file = open(name, "w")
    text_file.write("shuffled networks\n")
    np.asarray(array).flatten()
    text_file.write("{0}\n".format(array))
    text_file.write("original process\n")
    np.asarray(unshuffled).flatten()
    text_file.write("{0}\n".format(unshuffled))
    text_file.close()
    
def statistics(graph_on, array, variable, word, exp_deg, N_vertices, N_random_shuffles, zeta):
    pval = out_to_png(graph_on, array, word, variable, exp_deg, N_vertices, N_random_shuffles, zeta)
#    print(word, pval, "pval")
    return(pval)
