#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 12:03:42 2018

@author: tester
"""
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import graph_tool as gt
from graph_tool.all import *
import pandas as pd
import time
import scipy.special
import erdos_renyi_default_poisson #Poisson proces ili obični 
import set_and_reshuffle_poisson
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


#glavni dio programa u kojem mijenjamo parametre erdos_renyi grafa 
#sadrži petlju koja vrti po svim parametrima alpha (vjerojatnost slučajnog propadanja)
#unutar je druga petlja koja vrti po broju različitih kreiranih grafova
#u svakoj iteraciji za neki graf pozivamo funkciju process koja na njemu provede alpha beta i alphabeta 
#procese te napravi statistiku (largest component, onepath, twopath, threepath)
#nakon iteracija po svim zadanim grafovima, za dani alpha parametar napravimo još graf pvalues i graf
#na kojem su svi skupljeni podaci

def parse_args():
    """Parse input arguments
    """

    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    #TODO

    parser.add_argument('--exp_deg',
                        help=('Finetuning scenario'),
                        type = float,
                        default=2.0)

    parser.add_argument('--zeta',
                        help=('Finetuning scenario'),
                        type=float,
                        default=1)

    args = parser.parse_args()
    return args

def run_process(args):
	N_random_graphs = 100 #broj različitih grafova
	N_processes = 10
	N_random_shuffles = 100 #broj shuflleova na svakom grafu
	N_vertices = 1000 #broj čvorova
	time_stops = 20 #na koliko dijelova dijelimo ukupno vrijeme propadanja
	exp_deg = args.exp_deg
	p = exp_deg/(2*N_vertices)

	startTime = time.time()
	zeta = args.zeta # alpha = 1 - čisti alfa proces, alpha = 0 - čisti beta proces
	
	ks_lc = [[] for i in range(0,time_stops)]
	ks_op = [[] for i in range(0,time_stops)]
	ks_tp = [[] for i in range(0,time_stops)]
	ks_thp = [[] for i in range(0,time_stops)]
	mean_std_stat = [[] for i in range(0,time_stops)]
	p_lc = [[[] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	p_op = [[[] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	p_tp = [[[] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	p_thp = [[[] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	lc = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	op = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	tp = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	thp = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	largest = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	onepath = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	twopath = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
	threepath = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]

	path = "/share/storage/IRENA_STORAGE/data/k={0},zeta={1}".format(exp_deg,zeta)

	print(round(zeta,3), "z")
	k = 0
	#    bridovi = []
	#    čvorovi = []
	while(k < N_random_graphs):

	    print(k)
	    er = erdos_renyi_default_poisson.Erdos_Renyi(N_vertices,exp_deg,p)

	#        bridovi.append(er.g.num_edges())
	#        čvorovi.append(er.g.num_vertices())

	    graph_on = False
    
	    i = 0    

	    while(i < N_processes):
    
	        set_and_reshuffle_poisson.process(graph_on, i, k, er, zeta, largest, onepath, twopath, threepath, lc, op, tp, thp, p_lc, p_op, p_tp, p_thp, exp_deg, N_vertices, N_random_shuffles, time_stops)
		#print ('The script took {0} second !'.format(time.time() - startTime),"jedan k")

	        i += 1
    
	    k += 1
    
	#    print(bridovi)
	#    print(np.mean(bridovi))
	#    print(čvorovi)
	#    print(np.mean(čvorovi))
	#    deg = []
	#    g = 0
	#    while(g < N_random_graphs):
	#        deg.append(round(2*bridovi[g]/čvorovi[g],2))
	#        g += 1
	#    print(deg)
	#    print(np.mean(deg))
	print ('The script took {0} second !'.format(time.time() - startTime),"konačno")    
    


	set_and_reshuffle_poisson.signif_tests(largest, onepath, twopath, threepath, lc, op, tp, thp, ks_lc, ks_op, ks_tp, ks_thp, time_stops)

	set_and_reshuffle_poisson.stats_to_tex(path, ks_lc, ks_op, ks_tp, ks_thp, p_lc, p_op, p_tp, p_thp, time_stops, zeta, exp_deg)
      
	pval = True
	all_data = True
	zstat = True
	set_and_reshuffle_poisson.process_hist(mean_std_stat, zstat, path, pval, all_data, zeta, largest, onepath, twopath, threepath, lc, op, tp, thp, p_lc, p_op, p_tp, p_thp, exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, time_stops)        
        

	print ('The script took {0} second !'.format(time.time() - startTime),"konačno")

def main():
    args = parse_args()
    print(args)
    print(args.exp_deg, args.zeta)
    run_process(args)


if __name__ == '__main__':
    main()
    
