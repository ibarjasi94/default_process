#!/usr/bin/env python3
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import graph_tool as gt
from graph_tool.all import *
import pandas as pd
import time
import scipy.special
import collections

#The script contains an Erdos-Renyi graph class, that has functions to simulate the process, do the time reshuffling and count different motifs


def subgraph_num(subgraph,graph):
    i = 0
    gen = gt.topology.subgraph_isomorphism(subgraph, graph, generator=True)
    for iter in gen:
        i += 1
    return i

class Erdos_Renyi():
    def __init__(self,N_vertices,exp_deg,p):
        self.g = gt.Graph()
        for i in range(0,N_vertices):
            self.g.add_vertex()
            for j in range(0,i):
                if(rnd.uniform(0,1) < p):
                    self.g.add_edge(i,j)
                if(rnd.uniform(0,1) < p):
                    self.g.add_edge(j,i)
        self.g.vp.stecaj = self.g.new_vertex_property("bool")
        self.g.vp.known = self.g.new_vertex_property("bool")
        self.g.vp.time = self.g.new_vertex_property("double")
        self.g.ep.time = self.g.new_edge_property("double")
        self.g.vp.t_part = self.g.new_vertex_property("bool")
        self.g.ep.causal = self.g.new_edge_property("bool")
        samci = []
        self.g_c = self.g.copy()

    def set_directedness(self,directedness):
        self.g_c.set_directed(directedness)

    def initial_properties(self, zeta, dynamics):
        self.g_c = self.g.copy()
        N = self.g_c.num_vertices()
        vert_times = self.poisson_interevent(N,zeta)
        rnd.shuffle(vert_times)
        self.index_time = dict()
        for v in self.g_c.vertices():
            self.index_time[vert_times[-1]] = [self.g_c.vertex_index[v],"alpha"]
            self.g_c.vp.known[v] = True
            vert_times = vert_times[:-1]
        edge_interevent = self.poisson_interevent(self.g_c.num_edges(),1)
        for e in self.g_c.edges():
            if(dynamics == "VM"):
                s,t = e
                kin = t.in_degree()
                self.g_c.ep.time[e] = edge_interevent[-1]*kin
            elif(dynamics == "SI"):
                self.g_c.ep.time[e] = edge_interevent[-1]
            edge_interevent = edge_interevent[:-1]
        for e in self.g_c.edges():
            self.g_c.ep.causal[e] = False

    def process(self,exp_deg,zeta,N_vertices,frac_def):
        index_time_sorted = collections.OrderedDict(sorted(self.index_time.items()))
        self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = True)
        nondefault = self.g_c.get_vertices()
        self.g_c.clear_filters()
        pot = []
        N_pot = 0
        pot.append([0,0])
        while((N_vertices - len(nondefault)) < frac_def*N_vertices):
            vert_time_list = list(index_time_sorted.keys())
            vert_time = vert_time_list[0]
            vert = index_time_sorted[vert_time][0]
            prcs = index_time_sorted[vert_time][1]
            while(self.g_c.vp.stecaj[vert] == True):
                del index_time_sorted[vert_time]
                index_time_list = list(index_time_sorted.items())
                vert_time, vert_process = index_time_list[0]
                vert,prcs = vert_process
            del index_time_sorted[vert_time]   
            self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = False)
            in_edges = self.g_c.get_in_edges(vert)
            self.g_c.clear_filters()
            self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = True)
            out_edges = self.g_c.get_out_edges(vert)
            N_pot = N_pot + len(out_edges) - len(in_edges)
            for edge in out_edges:
                s,t,i = edge
                edge_time = self.g_c.ep.time[edge]
                cascade_time = vert_time + edge_time
                index_time_sorted[cascade_time] = [t,"beta"]
            pot.append([vert_time,N_pot])
            self.g_c.vp.time[vert] = vert_time
            self.g_c.vp.stecaj[vert] = True
            nondefault = self.g_c.get_vertices()
            self.g_c.clear_filters()
            index_time_sorted = collections.OrderedDict(sorted(index_time_sorted.items()))
        time_list = self.g_c.vp.time.get_array()
        pot = np.asarray(pot)
        return pot

    def network_part(self,time_frac_tot,l):
        time_list_p = list(self.g_c.vp.time.get_array())
        time_list_p.sort()
        last = time_list_p[int((l+1)/time_frac_tot*self.g_c.num_vertices())-1]
        for v in self.g_c.vertices():
            if(self.g_c.vp.time[v] <= last):
                self.g_c.vp.t_part[v] = True
        self.g_c.set_vertex_filter(self.g_c.vp.t_part)
        return last        

    def clear_graph(self):
        for i in self.g_c.vertices():
            self.g_c.vp.stecaj[i] = False

    def largest_component(self):
        l = gt.topology.label_largest_component(self.g_c, directed = False)
        u = gt.GraphView(self.g_c, vfilt=l)   
        return u.num_vertices()
        
    def one_path(self):
        self.g_c.set_directed(False)
        onepath = 0
        onepath = self.g_c.num_edges()
        self.g_c.set_directed(True)
        return onepath
    
    def two_path(self):
        sub_twoI = gt.Graph(directed =True)
        sub_twoV = gt.Graph(directed =True)
        sub_twoΛ = gt.Graph(directed =True)
        for n in range(0,3):
            sub_twoI.add_vertex()
            sub_twoV.add_vertex()
            sub_twoΛ.add_vertex()
        edgesI = [[0,1],[1,2]]
        edgesV= [[1,0],[1,2]]
        edgesΛ = [[0,1],[2,1]]
        sub_twoI.add_edge_list(edgesI)
        sub_twoV.add_edge_list(edgesV)
        sub_twoΛ.add_edge_list(edgesΛ)
        sub_isoI= subgraph_num(sub_twoI,self.g_c)
        sub_isoV= subgraph_num(sub_twoV,self.g_c)
        sub_isoΛ= subgraph_num(sub_twoΛ,self.g_c)
        numI = int(sub_isoI)
        numV = int(sub_isoV)
        numΛ = int(sub_isoΛ)
        return numI, numV, numΛ
       
    def three_path(self):
        self.g_c.set_directed(False)
        sub_three1 = gt.Graph(directed = False)
        sub_three2 = gt.Graph(directed = False)
        sub_three3 = gt.Graph(directed = False)
        for n in range(0,4):
            sub_three1.add_vertex()
            sub_three2.add_vertex()
        for n in range(0,3):
            sub_three3.add_vertex()
        edges1 = [[0,1],[1,2],[2,3]]
        edges2 = [[0,1],[1,2],[1,3]]
        edges3 = [[0,1],[1,2],[2,0]]
        sub_three1.add_edge_list(edges1)
        sub_three2.add_edge_list(edges2)
        sub_three3.add_edge_list(edges3)
        sub_iso1 = subgraph_num(sub_three1,self.g_c)
        sub_iso2 = subgraph_num(sub_three2,self.g_c)
        sub_iso3 = subgraph_num(sub_three3 ,self.g_c)
        num = int(sub_iso1)/2 + int(sub_iso2)/6 + int(sub_iso3)/6
        self.g_c.set_directed(True)
        return num
 
    def reset_edges(self):
        for i in self.g_c.edges():
            self.g_c.ep.causal[i] = False
        
    def only_causal(self):
        for i in self.g_c.edges():
            first,second = i
            first_time = self.g_c.vp.time[self.g_c.vertex(first)]
            second_time = self.g_c.vp.time[self.g_c.vertex(second)]
            if(first_time < second_time):
                self.g_c.ep.causal[i] = True
        self.g_c.set_edge_filter(self.g_c.ep.causal)


    def poisson_interevent(self,N,lam):
        uniform = np.random.rand(N)
        interevent = [-np.log(x)/lam for x in uniform]
        return interevent
        
    def shuffle_time(self):
        time_list = self.g_c.vp.time.get_array()
        rnd.shuffle(time_list)
        gt.map_property_values(self.g_c.vertex_index, self.g_c.vp.time, lambda x: time_list[x])   

        
    def remember_time(self):
        original_time = self.g_c.vp.time.get_array()
        return original_time
