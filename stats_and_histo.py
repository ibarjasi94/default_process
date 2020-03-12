#!/usr/bin/env python3
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import graph_tool as gt
from graph_tool.all import *
import pandas as pd
import time  
import AME  
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy

  

#skripta sadrži fuknciju koja radi histogram od pvalues svih napravljenih grafova
#sadrži funkciju za histograme svih podataka (varijabli reshuffleanih i početnih grafova)
#sadrži funkciju za određivanje pvalue jednog grafa

def numeric_analytic_comparison(exp_deg, exp_deg_lim, zeta, N_vertices, potentialb,  alpha_t, beta_t, tot_t, L_t, onepath, op, twopath, tp, perc_time, name, n_model, dynamics):
    beta_potential_mean = np.mean(potentialb, axis = 0)
    #print(tot_t)
    #print(onepath)
    op_t = array_reshape(onepath,perc_time)
    op_s_t = array_reshape(op,perc_time)
    tp_t = array_reshape(twopath,perc_time)
    tp_s_t = array_reshape(tp,perc_time)
    t = np.linspace(0,max(beta_potential_mean[:,0]),200)
    #t = np.linspace(0,5,100)
    #sol1,sol2,L_analytic = analytics.potential_node_analytical(exp_deg,zeta,N_vertices,t)
    kin_min = exp_deg_lim[0]
    kin_max = exp_deg_lim[1]
    kout_min = exp_deg_lim[2]
    kout_max = exp_deg_lim[3]
    #print("pasmater")
    SIsol,rhosol,rhosol_a,rhosol_b,onesol,onesol_var,onesol_s,onesol_var_s,twosol,twosol_var,twosol_s,twosol_var_s = \
    AME.AME_solve(kin_min,kout_min,kin_max,kout_max,exp_deg,N_vertices,zeta,t,n_model,dynamics)
    #print ('The script took {0} second !'.format(time.time() - startTime),"ame")
    #print(SIsol)
    xlim = max(beta_potential_mean[:,0])
    #xlim = 20
    ylim = None
    plot_comparison(t,SIsol,potentialb,xlim,ylim,"t","N potencijalni")
    plot_comparison(t,rhosol,tot_t,xlim,ylim,"t","N propali")
    print("diff")
    print(onesol[-1]-onesol_s[-1])
    print(twosol[-1]-twosol_s[-1])
    plot_var_comp(t,[[onesol,onesol_var],[onesol_s,onesol_var_s]],[op_t,op_s_t],xlim,ylim,"t","onepath",name)
    plot_var_comp(t,[[twosol,twosol_var],[twosol_s,twosol_var_s]],[tp_t,tp_s_t],xlim,ylim,"t","twopath",name)
    plot_to_tex(t,[[onesol,onesol_var],[onesol_s,onesol_var_s]],[op_t,op_s_t],"onepath",name)
    plot_to_tex(t,[[twosol,twosol_var],[twosol_s,twosol_var_s]],[tp_t,tp_s_t],"twopath",name)
    print(alpha_t[1])
#    plt.errorbar(t,twosol,yerr = twosol_var)
#    plt.errorbar(t,twosol_s,yerr = twosol_var_s)
    #plot_comparison(t,rhosol_a,[alpha_t[0]],xlim,ylim,"t","N alpha propali")
    #plot_comparison(t,rhosol_b,[beta_t[0]],xlim,ylim,"t","N beta propali")
#    plot_comparison(t,sol1[:,0],potentialb,xlim,ylim,"t","N potencijalni")
#    plot_comparison(t,sol1[:,1],tot_t,xlim,ylim,"t","N propali")
    #print(onesol)
#    
#    plot_comparison(t,L_analytic[0],L_t[0],xlim,ylim,"t","N alpha+ potencijalni linkovi")
#    plot_comparison(t,L_analytic[1],L_t[1],xlim,ylim,"t","N beta+ potencijalni linkovi")
#    plot_comparison(t,L_analytic[2],L_t[2],xlim,ylim,"t","N alpha- potencijalni linkovi")
#    plot_comparison(t,L_analytic[3],L_t[3],xlim,ylim,"t","N beta- potencijalni linkovi")
#    
#    plot(L_a_t[0][:,0],L_a_t[1][:,1], xlim, ylim,"t","N alpha potencijalni linkovi")
#    plot(L_b_t[0][:,0],L_b_t[1][:,1], xlim, ylim,"t","N beta potencijalni linkovi")
    #plot(tot_t[0][:,0],tot_t[0][:,1], xlim, ylim,"t","N propali")
    
def array_reshape(statistic,perc_time):
    inst = np.asarray(statistic[0]).flatten()
    statistic_t = [[] for i in inst]
    for n,l in enumerate(statistic):
        #onepath_st = [round(np.mean(np.asarray(l).flatten()),2),round(np.std(np.asarray(l).flatten()),2)]
        #avg_onepath.append(onepath_st[0])
        for m,i in enumerate(np.asarray(l).flatten()):
            #print(m)
            statistic_t[m].append([perc_time[n],i])
    return statistic_t

def plot_to_tex(timegrid,solution,sim_array,ylab,name):
    text_file = open("{0},{1}.{2}".format(ylab,name,"txt"), "w")
    sol = np.asarray(solution[0][0]).flatten()
    sol_var = np.asarray(solution[0][1]).flatten()
    sol_s = np.asarray(solution[1][0]).flatten()
    sol_var_s = np.asarray(solution[1][1]).flatten()
    sim_array_mean = np.nanmean(sim_array[0], axis = 0)
    sim_array_std = np.nanstd(sim_array[0], axis = 0)
    sim_array_mean_s = np.nanmean(sim_array[1], axis = 0)
    sim_array_std_s = np.nanstd(sim_array[1], axis = 0)
    text_file.write("ANALYTIC\n")
    text_file.write("analytic time\n")
    text_file.write("{0}\n".format(timegrid))
    text_file.write("original process analytic mean\n")
    text_file.write("{0}\n".format(sol))
    text_file.write("original process analytic std\n")
    text_file.write("{0}\n".format(sol_var))
    text_file.write("shuffled process analytic mean\n")
    text_file.write("{0}\n".format(sol_s))
    text_file.write("shuffled process analytic std\n")
    text_file.write("{0}\n".format(sol_var_s))
    text_file.write("SIMULATION\n")
    text_file.write("simulation time\n")
    text_file.write("{0}\n".format(sim_array_mean[:,0]))
    text_file.write("original process simulated mean\n")
    text_file.write("{0}\n".format(sim_array_mean[:,1]))
    text_file.write("original process simulated std\n")
    text_file.write("{0}\n".format(sim_array_std[:,1]))
    text_file.write("shuffled process simulated mean\n")
    text_file.write("{0}\n".format(sim_array_mean_s[:,1]))
    text_file.write("shuffled process simulated std\n")
    text_file.write("{0}\n".format(sim_array_std_s[:,1]))
    text_file.close()
    
def plot_var_comp(timegrid,sol,sim_array,xlim,ylim,xlab,ylab,name):
    sim_array_mean = np.nanmean(sim_array[0], axis = 0)
    sim_array_std = np.nanstd(sim_array[0], axis = 0)
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    #plt.xscale("log")
    axes = plt.gca() 
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.grid()
    markers, caps, bars = plt.errorbar(timegrid, sol[0][0], yerr=sol[0][1])
    [bar.set_alpha(0.8) for bar in bars]
    markers, caps, bars_s = plt.errorbar(sim_array_mean[:,0], sim_array_mean[:,1], yerr = sim_array_std[:,1],xerr = sim_array_std[:,0],ms=2,fmt = 'o',markerfacecolor=(0, 0, 1, 0.5))
    [bar.set_alpha(0.8) for bar in bars_s]
    sim_array_mean_s = np.nanmean(sim_array[1], axis = 0)
    sim_array_std_s = np.nanstd(sim_array[1], axis = 0)
    markers, caps, bars_r = plt.errorbar(timegrid, sol[1][0], yerr=sol[1][1])
    [bar.set_alpha(0.8) for bar in bars_r]
    markers_s, caps_s, bars_r_s = plt.errorbar(sim_array_mean_s[:,0], sim_array_mean_s[:,1], yerr = sim_array_std_s[:,1],xerr = sim_array_std_s[:,0],ms=2,fmt = 'o',markerfacecolor=(1, 1, 1, 0.5))
    [bar.set_alpha(0.8) for bar in bars_r_s]
    axes.set_xlim([0,xlim])
    axes.set_ylim([0,ylim])  
#    axins = inset_axes(axes, width='40%', height='30%', loc=5)
#    axins.errorbar(timegrid, sol[0][0], yerr=sol[0][1])
#    axins.errorbar(sim_array_mean[:,0], sim_array_mean[:,1], yerr = sim_array_std[:,1],xerr = sim_array_std[:,0],ms=2,fmt = 'o',markerfacecolor=(0, 0, 1, 0.5))
#    axins.errorbar(timegrid, sol[1][0], yerr=sol[1][1])
#    axins.errorbar(sim_array_mean_s[:,0], sim_array_mean_s[:,1], yerr = sim_array_std_s[:,1],xerr = sim_array_std_s[:,0],ms=2,fmt = 'o',markerfacecolor=(1, 1, 1, 0.5))
#    axins.set_xlim([0,xlim/2])
#    plt.grid()
    axes.legend((bars,bars_s,bars_r,bars_r_s),("analytic-original","simulation-original","analytic-randomized","simulation-randomized"))
    plt.savefig("{0},{1}.{2}".format(ylab,name,"png"))
    plt.show()
    plt.clf()
    
def plot_comparison(timegrid,sol,sim_array,xlim,ylim,xlab,ylab):
    sim_array_mean = np.nanmean(sim_array, axis = 0)
    sim_array_std = np.nanstd(sim_array, axis = 0)
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(timegrid, sol, 'r')
    markers, caps, bars = plt.errorbar(sim_array_mean[:,0], sim_array_mean[:,1], yerr = sim_array_std[:,1],xerr = sim_array_std[:,0],ms=1,fmt = 'o',markerfacecolor=(0, 0, 1, 0.5))
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    axes = plt.gca()
    axes.set_xlim([None,xlim])
    axes.set_ylim([-2,ylim])
    plt.show()
    plt.clf()
    
def plot(array_time, array_value, xlim, ylim, xlab, ylab):
    plt.plot(array_time,array_value,'o',ms=1)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    axes = plt.gca()
    axes.set_xlim([0,xlim])
    axes.set_ylim([0,ylim])
    plt.show()
    plt.clf()
    
def ks_test(array, unshuffled):
    return scipy.stats.ks_2samp(array,unshuffled)

def mean_and_std(avg_array, array, name, stat, time_stops):
    for n,l in enumerate(array):
        st = [round(np.mean(np.asarray(l).flatten()),2),round(np.std(np.asarray(l).flatten()),2)]
        avg_array[n] = [st]
    print(avg_array)
    out_to_tex_stat(avg_array, name, stat, time_stops)
        
def zstat_histo(avg_array, array, unshuffled, name, zeta):
    avg_array = [np.round(np.asarray(array).mean(axis = 2),2),np.round(np.asarray(array).std(axis = 2),2)]
    print(avg_array)
    np.seterr(divide='ignore', invalid='ignore')
    zstat = np.divide(np.asarray(unshuffled).flatten() - avg_array[0].flatten(),avg_array[1].flatten())
    #zstat = (np.asarray(unshuffled).flatten() - avg_array[0].flatten())/avg_array[1].flatten()
    print(zstat)
    zstat = [x for x in zstat if ~np.isnan(x)]
    zstat = [x for x in zstat if ~np.isinf(x)]
    print(zstat)
    plt.hist(zstat, bins = 100)
    plt.title(r' $\zeta$ = %.2f, $\mu$ = %.2f, $\sigma$ = %.2f, kurt = %.2f' %(zeta, np.mean(zstat), np.std(zstat), scipy.stats.kurtosis(zstat)))
    plt.xlabel(name)
    plt.ylabel("number of occurences")
    plt.savefig(name)
    plt.show()
    plt.gcf().clear()
    
def pval_histo(array, name, zeta):
    plt.hist(array, bins = 100)
    plt.title(r' $\zeta$ = %.2f' %(zeta))
    plt.xlabel(name)
    plt.ylabel("number of occurences")
    plt.savefig(name)
    plt.show()
    plt.gcf().clear()

def all_data_hist(array, name, unshuffled, exp_deg, N_vertices, N_random_graphs, N_processes, N_random_shuffles, zeta):
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
    plt.savefig("{0}.{1}".format(name,"png"))
    plt.show()
    plt.gcf().clear()
    out_to_tex(array,unshuffled,"{0}.{1}".format(name,"txt"))

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
        plt.show()
    plt.gcf().clear()
    return pval

def out_to_tex_stat(array, name, stat, time_stops):
    text_file = open("{0},{1}".format(stat,name), "w")
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