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
import erdos_renyi_default_poisson_MH as erdos_renyi_default_poisson #Poisson proces ili obični 
import set_and_reshuffle_poisson_MH as set_and_reshuffle_poisson
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


#The main script that we use to set the parameters and iterate over them
#For each iteration a function is called to simulate the process on the network and collect all the necessary statistics from the desired time points od the default process
#All the arrays with the counts of statistics from the simulations are then saved as .npz files



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

    parser.add_argument('--dynamics',
			help=('Finetuning scenario'),
			type=str,
			default="SI")

    parser.add_argument('--N_vertices',
                        help=('Finetuning scenario'),
                        type=int,
                        default=1000)
    
    parser.add_argument('--N_samples',
			help=('Finetuning scenario'),
			type=str,
                        default="100,10,10")

    args = parser.parse_args()
    return args

def run_process(args):
    N_samples = [int(item) for item in args.N_samples.split(',')]
    N_random_graphs = N_samples[0] #number of different networks
    N_processes = N_samples[1]
    N_random_shuffles = N_samples[2] #number of shuffles on each network
    N_vertices = args.N_vertices #number of vertices
    time_frac_start = 5
    time_frac_stop = 10
    time_frac_tot = 20 #number of total time subdivisions
    exp_deg = args.exp_deg
    p = exp_deg/(2*N_vertices)
    n_model = "Poisson"
    dynamics = args.dynamics # "SI", "VM"

    kmin = 0
    kmax = 2*exp_deg+4
    exp_deg_lim = [kmin,kmax,kmin,kmax]

    zeta = args.zeta # alpha = 1 - pure exogenous, alpha = 0 - pure endogenous process (would not start without external contagion)
	
    lc = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    op = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    tpI = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    tpV = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    tpΛ = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]    
    tp = [tpI, tpV, tpΛ]
    thp = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    largest = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    onepath = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    twopathI = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    twopathV = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    twopathΛ = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]
    twopath = [twopathI,twopathV ,twopathΛ]
    threepath = [[[[] for k in range(0,N_processes)] for j in range(0,N_random_graphs)] for i in range(time_frac_start,time_frac_stop)]

    perc_time = []
    
    path = "/share/storage/IRENA_STORAGE/data/mahalanobis/k={0},zeta={1},dyn={2},N={3}".format(exp_deg,zeta,dynamics,N_vertices)
    
    print(round(zeta,3), "z")
    k = 0

    while(k < N_random_graphs):

        print(k)
        er = erdos_renyi_default_poisson.Erdos_Renyi(N_vertices,exp_deg,p)
        tot_subgraphs = [er.largest_component(),er.one_path(),er.two_path(),er.three_path()]



        graph_on = False
        
        i = 0    

        while(i < N_processes):
            
            set_and_reshuffle_poisson.process(graph_on, i, k, er, zeta, largest, onepath, twopath, threepath,\
                        lc, op, tp, thp, exp_deg, N_vertices, N_random_shuffles, time_frac_start, time_frac_stop,\
                        time_frac_tot, perc_time, dynamics,tot_subgraphs)

            i += 1
            
            
        k += 1
        
    add = True        
    set_and_reshuffle_poisson.process_save(add, path, zeta, largest, onepath, twopath, threepath, lc, op,\
                                            tp, thp, exp_deg, N_vertices, time_frac_start, time_frac_stop, \
                                                time_frac_tot, dynamics, perc_time)   

         
    
def main():
    args = parse_args()
    print(args)
    print(args.exp_deg, args.zeta,args.dynamics)
    run_process(args)


if __name__ == '__main__':
    main()
    
