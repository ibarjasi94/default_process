#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 15:38:02 2018

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
import stats_and_histo_MH as stats_and_histo

#The script contains a function that calls the functions that count all the test statistics (causal motifs and the largest component) at different percentages of network default. It then randomizes the default times, and collects the same data from the RRM.
#It contains a function that coordinates the saving of the data
def process(graph_on, i, k, er, zeta, largest, onepath, twopath, threepath, lc, op, tp, thp,exp_deg, N_vertices,\
            N_random_shuffles, time_frac_start, time_frac_stop, time_frac_tot,  perc_time, dynamics,tot_subgraphs):

    er.initial_properties(zeta,dynamics)
    frac_def = time_frac_stop/time_frac_tot
    potential = er.process(exp_deg,zeta,N_vertices,frac_def)
   
    
    one_tot = tot_subgraphs[1]
    two_totI,two_totV,two_totΛ = tot_subgraphs[2]
    three_tot = tot_subgraphs[3]
    largest_tot = tot_subgraphs[0]
   
    original_time = list(er.remember_time())
    l = time_frac_start
    while(l < time_frac_stop):
        p = l - time_frac_start
        gt.map_property_values(er.g_c.vertex_index, er.g_c.vp.time, lambda x: original_time[x])   
        perc_time.append(er.network_part(time_frac_tot,l))
        er.only_causal()
        er.g_c.set_edge_filter(er.g_c.ep.causal)
        er.set_directedness(False)
        largest[p][k][i].append(er.largest_component()/largest_tot)
        onepath[p][k][i].append(er.one_path()/one_tot)
        twoI,twoV,twoΛ = er.two_path() 
        twopath[0][p][k][i].append(twoI/two_totI)
        twopath[1][p][k][i].append(twoV/two_totV)
        twopath[2][p][k][i].append(twoΛ/two_totΛ)
        threepath[p][k][i].append(er.three_path()/three_tot)
        er.g_c.clear_filters()
        er.g_c.set_vertex_filter(er.g_c.vp.t_part)
        er.set_directedness(True)
        er.reset_edges()
        j = 0
        while(j < N_random_shuffles):
            er.shuffle_time()    
            er.only_causal()
            er.g_c.set_edge_filter(er.g_c.ep.causal)
            er.set_directedness(False)
            lc[p][k][i].append(er.largest_component()/largest_tot)
            op[p][k][i].append(er.one_path()/one_tot)
            tI,tV,tΛ = er.two_path() 
            tp[0][p][k][i].append(tI/two_totI)
            tp[1][p][k][i].append(tV/two_totV)
            tp[2][p][k][i].append(tΛ/two_totΛ) 
            thp[p][k][i].append(er.three_path()/three_tot)
            er.g_c.clear_filters()
            er.g_c.set_vertex_filter(er.g_c.vp.t_part)
            er.set_directedness(True)
            er.reset_edges()
            j += 1
        er.g_c.clear_filters()
         
        l += 1
    
    er.clear_graph()


def process_save(add, path, zeta, largest, onepath, twopath, threepath, lc, op, tp, thp,exp_deg, N_vertices,time_frac_start,\
                time_frac_stop,time_frac_tot,dynamics,perc_time):       

    
    np.save("{0}/perc_time,zeta={1},k={2},dyn={3}.npy".format(path,round(zeta,3), exp_deg, dynamics),perc_time)
    if(add == True):
        twopath_all = np.sum(twopath,axis = 0)
        tp_all = np.sum(tp,axis = 0) 

    l = time_frac_start
    while(l < time_frac_stop):
        p = l - time_frac_start
        stats_and_histo.out_to_tex(path,np.asarray(lc[p]).flatten(), np.asarray(largest[p]).flatten(), r'largest,zeta={0},network_default={1}%,k={2},dyn={3}'.format(round(zeta,3), round(100*(l+1)/time_frac_tot,2), exp_deg, dynamics))
        stats_and_histo.out_to_tex(path,np.asarray(op[p]).flatten(), np.asarray(onepath[p]).flatten(), r'onepath,zeta={0},network_default={1}%,k={2},dyn={3}'.format(round(zeta,3), round(100*(l+1)/time_frac_tot,2), exp_deg, dynamics))
        stats_and_histo.out_to_tex(path,np.asarray(tp[0][p]).flatten(), np.asarray(twopath[0][p]).flatten(), r'twopathI,zeta={0},network_default={1}%,k={2},dyn={3}'.format(round(zeta,3), round(100*(l+1)/time_frac_tot,2), exp_deg, dynamics))
        stats_and_histo.out_to_tex(path,np.asarray(tp[1][p]).flatten(), np.asarray(twopath[1][p]).flatten(), r'twopathV,zeta={0},network_default={1}%,k={2},dyn={3}'.format(round(zeta,3), round(100*(l+1)/time_frac_tot,2), exp_deg, dynamics))
        stats_and_histo.out_to_tex(path,np.asarray(tp[2][p]).flatten(), np.asarray(twopath[2][p]).flatten(), r'twopathΛ,zeta={0},network_default={1}%,k={2},dyn={3}'.format(round(zeta,3), round(100*(l+1)/time_frac_tot,2), exp_deg, dynamics))
        stats_and_histo.out_to_tex(path,np.asarray(tp_all[p]).flatten(), np.asarray(twopath_all[p]).flatten(), r'twopath,zeta={0},network_default={1}%,k={2},dyn={3}'.format(round(zeta,3), round(100*(l+1)/time_frac_tot,2), exp_deg, dynamics))
        stats_and_histo.out_to_tex(path,np.asarray(thp[p]).flatten(), np.asarray(threepath[p]).flatten(), r'threepath,zeta={0},network_default={1}%,k={2},dyn={3}'.format(round(zeta,3), round(100*(l+1)/time_frac_tot,2), exp_deg, dynamics))
        l += 1
