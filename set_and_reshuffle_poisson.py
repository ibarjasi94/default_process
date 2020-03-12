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
import stats_and_histo

#skripta sadrži funkciju koja provodi alpha ili beta proces propadanja na danom erdos renyi grafu,mjeri varijable
#te zatim radi reshuffleove na istom grafu pa opet mjeri
#također sadrži funkciju koja nakon svi generiranih grafova radi histograme od svih podataka
def process(graph_on, i, k, er, zeta, largest, onepath, twopath, threepath, lc, op, tp, thp, p_lc, p_op, p_tp, p_thp, exp_deg, N_vertices, N_random_shuffles, time_stops):
#    start = i * N_random_shuffles 
#    stop = (i+1) * N_random_shuffles 

    er.initial_properties(zeta)

    er.process()
    
#    graph_draw(er.g_c, vertex_text = er.g_c.vertex_index, vertex_halo = er.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
#    print("originalni graf")
#    er.set_directedness(False)
    
    one_tot = er.one_path()
    two_tot = er.two_path()
    three_tot = er.three_path()
    largest_tot = er.largest_component()
   
    original_time = list(er.remember_time())
    l = 0
    while(l < time_stops):
        gt.map_property_values(er.g_c.vertex_index, er.g_c.vp.time, lambda x: original_time[x])   
        er.network_part(time_stops,l)
#        graph_draw(er.g_c, vertex_text = er.g_c.vertex_index, vertex_halo = er.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
#        print(er.g_c.num_edges())
        er.only_causal()
        er.g_c.set_edge_filter(er.g_c.ep.causal)
#        graph_draw(er.g_c, vertex_text = er.g_c.vertex_index, vertex_halo = er.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
#        print(er.g_c.num_edges())
        er.set_directedness(False)
        largest[l][k][i].append(er.largest_component()/largest_tot)
        onepath[l][k][i].append(er.one_path()/one_tot)
        twopath[l][k][i].append(er.two_path()/two_tot)
        threepath[l][k][i].append(er.three_path()/three_tot)
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
            lc[l][k][i].append(er.largest_component()/largest_tot)
            op[l][k][i].append(er.one_path()/one_tot)
            tp[l][k][i].append(er.two_path()/two_tot)        
            thp[l][k][i].append(er.three_path()/three_tot)
            er.g_c.clear_filters()
            er.g_c.set_vertex_filter(er.g_c.vp.t_part)
            er.set_directedness(True)
            er.reset_edges()
            j += 1
        er.g_c.clear_filters()
        l += 1
    
    er.clear_graph()
    
#    print(onepath)
#    print(op)

    l = 0
    while(l < time_stops):
        p_lc[l][k].append(stats_and_histo.statistics(graph_on, lc[l][k][i], largest[l][k][i], r'largest zeta = {0},  network default = {1}%.txt'.format(round(zeta,3), round(100*(l+1)/time_stops,2)), exp_deg, N_vertices, N_random_shuffles, zeta))
        p_op[l][k].append(stats_and_histo.statistics(graph_on, op[l][k][i], onepath[l][k][i], r'onepath zeta = {0} network default = {1}%.txt'.format(round(zeta,3), round(100*(l+1)/time_stops,2)), exp_deg, N_vertices, N_random_shuffles, zeta))
        p_tp[l][k].append(stats_and_histo.statistics(graph_on, tp[l][k][i], twopath[l][k][i], r'twopath zeta = {0} network default = {1}%.txt'.format(round(zeta,3), round(100*(l+1)/time_stops,2)), exp_deg, N_vertices, N_random_shuffles, zeta))   
        p_thp[l][k].append(stats_and_histo.statistics(graph_on, thp[l][k][i], threepath[l][k][i], r'threepath zeta = {0} network default = {1}%.txt'.format(round(zeta,3), round(100*(l+1)/time_stops,2)), exp_deg, N_vertices, N_random_shuffles, zeta))
        l += 1

def process_hist(mean_std_stat, zstat, path, pval, all_data, zeta, largest, onepath, twopath, threepath, lc, op, tp, thp, p_lc, p_op, p_tp, p_thp, exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, time_stops):       
    
    if(pval == True):
        l = 0
        while(l < time_stops):
            stats_and_histo.pval_histo(path, np.asarray(p_lc[l]).flatten(), r'p-value_largest,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            stats_and_histo.pval_histo(path, np.asarray(p_op[l]).flatten(), r'p-value_onepath,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            stats_and_histo.pval_histo(path, np.asarray(p_tp[l]).flatten(), r'p-value_twopath,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            stats_and_histo.pval_histo(path, np.asarray(p_thp[l]).flatten(), r'p-value_threepath,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            l += 1
        
    if(all_data == True):
        l = 0
        while(l < time_stops):        
            stats_and_histo.all_data_hist(path, np.asarray(lc[l]).flatten(), r'largest,zeta={0},network_default={1}%,k={2},all_graphs'.format(round(zeta,3), round(100*(l+1)/time_stops,2), exp_deg), np.asarray(largest[l]).flatten(), exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, zeta)
            stats_and_histo.all_data_hist(path, np.asarray(op[l]).flatten(), r'onepath,zeta={0},network_default={1}%,k={2},all_graphs'.format(round(zeta,3), round(100*(l+1)/time_stops,2), exp_deg), np.asarray(onepath[l]).flatten(), exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, zeta)
            stats_and_histo.all_data_hist(path, np.asarray(tp[l]).flatten(), r'twopath,zeta={0},network_default={1}%,k={2},all_graphs'.format(round(zeta,3), round(100*(l+1)/time_stops,2), exp_deg), np.asarray(twopath[l]).flatten(), exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, zeta)
            stats_and_histo.all_data_hist(path, np.asarray(thp[l]).flatten(), r'threepath,zeta={0},network_default={1}%,k={2},all_graphs'.format(round(zeta,3), round(100*(l+1)/time_stops,2), exp_deg), np.asarray(threepath[l]).flatten(), exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, zeta)
            l += 1

    if(zstat == True):
        l = 0
        while(l < time_stops):
            stats_and_histo.zstat_histo(path, mean_std_stat, lc[l], largest[l], r'z-stat_largest,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            stats_and_histo.zstat_histo(path, mean_std_stat, op[l], onepath[l], r'z-stat_onepath,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            stats_and_histo.zstat_histo(path, mean_std_stat, tp[l], twopath[l], r'z-stat_twopath,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            stats_and_histo.zstat_histo(path, mean_std_stat, thp[l], threepath[l], r'z-stat_threepath,k={0},zeta={1},network_default={2}%.png'.format(exp_deg,round(zeta,3), round(100*(l+1)/time_stops,2)), zeta)
            l += 1
            
def signif_tests(largest, onepath, twopath, threepath, lc, op, tp, thp, ks_lc, ks_op, ks_tp, ks_thp, time_stops):
    
    l = 0
    while(l < time_stops):
        ks_lc[l].append(stats_and_histo.ks_test(np.asarray(lc[l]).flatten(), np.asarray(largest[l]).flatten()))
        ks_op[l].append(stats_and_histo.ks_test(np.asarray(op[l]).flatten(), np.asarray(onepath[l]).flatten()))
        ks_tp[l].append(stats_and_histo.ks_test(np.asarray(tp[l]).flatten(), np.asarray(twopath[l]).flatten()))
        ks_thp[l].append(stats_and_histo.ks_test(np.asarray(thp[l]).flatten(), np.asarray(threepath[l]).flatten()))
        l += 1
            
def stats_to_tex(path, ks_lc, ks_op, ks_tp, ks_thp, p_lc, p_op, p_tp, p_thp, time_stops, zeta, exp_deg):
    stats_and_histo.out_to_tex_stat(path, p_lc, r'largest_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "pval", time_stops)
    stats_and_histo.out_to_tex_stat(path, p_op, r'onepath_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "pval", time_stops)
    stats_and_histo.out_to_tex_stat(path, p_tp, r'twopath_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "pval", time_stops)
    stats_and_histo.out_to_tex_stat(path, p_thp, r'threepath_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "pval", time_stops)
    
    stats_and_histo.out_to_tex_stat(path, ks_lc, r'largest_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "ks", time_stops)
    stats_and_histo.out_to_tex_stat(path, ks_op, r'onepath_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "ks", time_stops)
    stats_and_histo.out_to_tex_stat(path, ks_tp, r'twopath_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "ks", time_stops)
    stats_and_histo.out_to_tex_stat(path, ks_thp, r'threepath_zeta={0},k={1}.txt'.format(round(zeta,3), exp_deg), "ks", time_stops)
    
    
    
