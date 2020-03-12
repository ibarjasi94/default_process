#!/usr/bin/env python3.5
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
import analytics
import stats_and_histo


#glavni dio programa u kojem mijenjamo parametre erdos_renyi grafa 
#sadrži petlju koja vrti po svim parametrima alpha (vjerojatnost slučajnog propadanja)
#unutar je druga petlja koja vrti po broju različitih kreiranih grafova
#u svakoj iteraciji za neki graf pozivamo funkciju process koja na njemu provede alpha beta i alphabeta 
#procese te napravi statistiku (largest component, onepath, twopath, threepath)
#nakon iteracija po svim zadanim grafovima, za dani alpha parametar napravimo još graf pvalues i graf
#na kojem su svi skupljeni podaci


N_random_graphs = 2 #broj različitih grafova
N_processes = 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
N_random_shuffles = 2 #broj shuflleova na svakom grafu
N_vertices = 1000 #broj čvorova
time_stops = 20 #na koliko dijelova dijelimo ukupno vrijeme propadanja
exp_degs = [2]
n_model = "Poisson" #"Poisson", "k-regular""scale free"
dynamics = "VM" # "SI", "VM"

for exp_deg in exp_degs:
    p = exp_deg/(2*N_vertices)
    kmin = 0
    kmax = 2*exp_deg+4
    exp_deg_lim = [kmin,kmax,kmin,kmax]
    
    startTime = time.time()
    zetas = [1]# alpha = 1 - čisti alfa proces, alpha = 0 - čisti beta proces
    for zeta in zetas:
        ks_lc = [[] for i in range(0,time_stops)]
        ks_op = [[] for i in range(0,time_stops)]
        ks_tp = [[] for i in range(0,time_stops)]
        ks_thp = [[] for i in range(0,time_stops)]
        average_degree = [[] for i in range(0,time_stops)]
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
        degree = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(0,time_stops)]
        potentialb = []
        alpha_t = []
        beta_t = []
        tot_t = []
        L_a_p_t = []
        L_b_p_t = []
        L_a_n_t = []
        L_b_n_t = []
        N_samci = []
        perc_time = []
        
        #timegrid = np.linspace(0,N_vertices/zeta,int(N_vertices/zeta*1000+1))
        
        print(round(zeta,3), "z")
        k = 0
        #bridovi = []
        #čvorovi = []
        while(k < N_random_graphs):
        
            print(k)
            er = erdos_renyi_default_poisson.Erdos_Renyi(N_vertices,exp_deg,p,zeta)
            N_samci.append(er.n_samci)
        
            #bridovi.append(er.g.num_edges())
            #čvorovi.append(er.g.num_vertices())
        
            graph_on = False
            
            i = 0    

            while(i < N_processes):
                
                potential, alpha_time, beta_time, tot_time, L = set_and_reshuffle_poisson.process(graph_on, i, k, er, zeta, degree, largest, onepath, twopath, threepath, lc, op, tp, thp, p_lc, p_op, p_tp, p_thp, exp_deg, N_vertices, N_random_shuffles, time_stops, perc_time, dynamics)
                potentialb.append(potential)
                alpha_t.append(alpha_time)
                beta_t.append(beta_time)
                tot_t.append(tot_time)
                L_a_p_t.append(L[0])
                L_b_p_t.append(L[1])
                L_a_n_t.append(L[2])
                L_b_n_t.append(L[3])
                #print ('The script took {0} second !'.format(time.time() - startTime),"jedan k")
                i += 1
            print ('The script took {0} second !'.format(time.time() - startTime),"konačno")    
            k += 1
            
        L_t = [L_a_p_t,L_b_p_t,L_a_n_t,L_b_n_t]
        #print(alpha_t)
        print(np.mean(N_samci))
        set_and_reshuffle_poisson.numeric_analytic(exp_deg,exp_deg_lim,zeta,N_vertices, \
                                                   potentialb, alpha_t, beta_t, tot_t, L_t,  \
                                                   onepath, op, twopath, tp, perc_time, n_model, dynamics)

#        print(bridovi)
#        print(np.mean(bridovi))
#        print(čvorovi)
#        print(np.mean(čvorovi))
#        deg = []
#        g = 0
#        while(g < N_random_graphs):
#            deg.append(round(2*bridovi[g]/čvorovi[g],2))
#            g += 1
#        print(deg)
#        print(np.mean(deg))
#        
#ovdje odkomentiraj
        print ('The script took {0} second !'.format(time.time() - startTime),"konačno")    
        
        
        set_and_reshuffle_poisson.signif_tests(largest, onepath, twopath, threepath, lc, op, tp, thp, ks_lc, ks_op, ks_tp, ks_thp, time_stops)
        
        set_and_reshuffle_poisson.stats_to_tex(ks_lc, ks_op, ks_tp, ks_thp, p_lc, p_op, p_tp, p_thp, time_stops, zeta, exp_deg)
        set_and_reshuffle_poisson.degree(average_degree, degree, time_stops, zeta, exp_deg)
        
             
        pval = True
        all_data = True
        zstat = True
        set_and_reshuffle_poisson.process_hist(mean_std_stat, zstat, pval, all_data, zeta, largest, onepath, twopath, threepath, lc, op, tp, thp, p_lc, p_op, p_tp, p_thp, exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, time_stops)        
        
        print ('The script took {0} second !'.format(time.time() - startTime),"konačno")
            
