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

#modul stvara Erdos-Renyi graf te sadrži sve potrebno za propadanje čvorova u alpha i beta procesima
#reshuffle vremena tih propadanja te mjerenje varijabli na grafu (lc,op,tp,thp)


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
#        for v in self.g.get_vertices():
#            if(len(self.g.get_in_edges(v)) == 0):
#                samci.append(v)
#        self.g.remove_vertex(samci)

    def set_directedness(self,directedness):
        self.g_c.set_directed(directedness)

    def initial_properties(self, zeta):
        self.g_c = self.g.copy()
        N = self.g_c.num_vertices()
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_halo = self.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
        vert_times = self.poisson_interevent(N,zeta)
#        print(vert_interevent,"vert")
#        vert_times = [np.sum(vert_interevent[0:i]) for i in range(0,len(vert_interevent))]
        rnd.shuffle(vert_times)
#        print(len(vert_times))alpha,b
        self.index_time = dict()
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_halo = self.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
        for v in self.g_c.vertices():
#            if(len(self.g_c.get_in_edges(v)) == 0):
#                print(v,"samac")
            self.index_time[vert_times[-1]] = [self.g_c.vertex_index[v],"alpha"]
            self.g_c.vp.known[v] = True
            vert_times = vert_times[:-1]
        edge_interevent = self.poisson_interevent(self.g_c.num_edges(),1)
#        print(edge_interevent, "beta times")
        for e in self.g_c.edges():
            self.g_c.ep.time[e] = edge_interevent[-1]
            edge_interevent = edge_interevent[:-1]
#            print(self.g_c.ep.time[i])
        for e in self.g_c.edges():
            self.g_c.ep.causal[e] = False
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_halo = self.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
#        print(self.g_c.num_edges(),"edges",self.g_c.num_vertices(),"vertices")

    def process(self):
        index_time_sorted = collections.OrderedDict(sorted(self.index_time.items()))
#        print(index_time_sorted,"početna vremena")
        self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = True)
        nondefault = self.g_c.get_vertices()
        self.g_c.clear_filters()
#        alpha = 0
#        beta = 0
#        n = 0 
#        D = 450
        while(len(nondefault) != 0):
#            print(index_time_sorted,"update vremena")
            vert_time_list = list(index_time_sorted.keys())
            vert_time = vert_time_list[0]
            vert = index_time_sorted[vert_time][0]
            prcs = index_time_sorted[vert_time][1]
            while(self.g_c.vp.stecaj[vert] == True):
                del index_time_sorted[vert_time]
                index_time_list = list(index_time_sorted.items())
                vert_time, vert_process = index_time_list[0]
                vert,prcs = vert_process
#            print(vert,vert_time,"čvor i vrijeme")
            del index_time_sorted[vert_time]   
#            if(n < D):
#                if(prcs == "alpha"):
#                    alpha += 1
#                else:
#                    beta += 1
#            print(index_time_sorted,"obrisani propali")False
            self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = True)
            out_edges = self.g_c.get_out_edges(vert)
            for edge in out_edges:
                s,t,i = edge
                edge_time = self.g_c.ep.time[edge]
                cascade_time = vert_time + edge_time
                index_time_sorted[cascade_time] = [t,"beta"]
#            print("propao", vert, "u", vert_time)
            self.g_c.vp.time[vert] = vert_time
            self.g_c.vp.stecaj[vert] = True
            nondefault = self.g_c.get_vertices()
            self.g_c.clear_filters()
            index_time_sorted = collections.OrderedDict(sorted(index_time_sorted.items()))
#            n += 1
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_font_size=12, output_size=(600, 600))
        time_list = self.g_c.vp.time.get_array()
#        print(alpha, "alpha")
#        print(beta, "beta")
#        print(time_list,"original time list")

    def network_part(self,time_stops,l):
        time_list_p = list(self.g_c.vp.time.get_array())
        time_list_p.sort()
        last = time_list_p[int((l+1)/time_stops*self.g_c.num_vertices())-1]
        for v in self.g_c.vertices():
            if(self.g_c.vp.time[v] < last):
                self.g_c.vp.t_part[v] = True
        self.g_c.set_vertex_filter(self.g_c.vp.t_part)
        
    def clear_graph(self):
        for i in self.g_c.vertices():
            self.g_c.vp.stecaj[i] = False

    def largest_component(self):
#        graph_draw(self.g_c,  vertex_font_size=12, output_size=(600, 600))  
        l = gt.topology.label_largest_component(self.g_c, directed = False)
#        graph_draw(self.g_c,  vertex_halo = l,vertex_font_size=12, output_size=(600, 600))
        u = gt.GraphView(self.g_c, vfilt=l)   
        return u.num_vertices()
        
    def one_path(self):
        self.g_c.set_directed(False)
        onepath = 0
        onepath = self.g_c.num_edges()
        self.g_c.set_directed(True)
        return onepath
    
    def two_path(self):
        self.g_c.set_directed(False)
        sub_two = gt.Graph(directed = False)
        for n in range(0,3):
            sub_two.add_vertex()
        edges = [[0,1],[1,2]]
        sub_two.add_edge_list(edges)
        sub_iso = subgraph_num(sub_two,self.g_c)
        num = int(sub_iso/2)
        self.g_c.set_directed(True)
        return num
    
    def three_path(self):
        self.g_c.set_directed(False)
        sub_three1 = gt.Graph(directed = False)
        sub_three2 = gt.Graph(directed = False)
        for n in range(0,4):
            sub_three1.add_vertex()
            sub_three2.add_vertex()
        edges1 = [[0,1],[1,2],[2,3]]
        edges2 = [[0,1],[1,2],[1,3]]
        sub_three1.add_edge_list(edges1)
        sub_three2.add_edge_list(edges2)
        sub_iso1 = subgraph_num(sub_three1,self.g_c)
        sub_iso2 = subgraph_num(sub_three2,self.g_c)
        num = int(sub_iso1/2) + int(sub_iso2/6)
        self.g_c.set_directed(True)
        return num
        
    def reset_edges(self):
        for i in self.g_c.edges():
            self.g_c.ep.causal[i] = False
        
    def only_causal(self):
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_font_size=12, output_size=(600, 600))
#        print("all default links")
#        for v in self.g_c.vertices():
#            print(v,self.g_c.vp.time[v], "vertex time")
        for i in self.g_c.edges():
#            print("edge is", self.g_c.ep.causal[i])
            first,second = i
            first_time = self.g_c.vp.time[self.g_c.vertex(first)]
            second_time = self.g_c.vp.time[self.g_c.vertex(second)]
            if(first_time < second_time):
                self.g_c.ep.causal[i] = True
#                print("edge", i, "is", self.g_c.ep.causal[i])
        self.g_c.set_edge_filter(self.g_c.ep.causal)
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_font_size=12, output_size=(600, 600))
#        print("only causal links")

    def poisson_interevent(self,N,lam):
        uniform = np.random.rand(N)
        interevent = [-np.log(x)/lam for x in uniform]
        return interevent
        
    def shuffle_time(self):
        time_list = self.g_c.vp.time.get_array()
        rnd.shuffle(time_list)
        gt.map_property_values(self.g_c.vertex_index, self.g_c.vp.time, lambda x: time_list[x])   
#        print(time_list)
#        for i in self.g_c.vertices():
#            print(i, self.g_c.vp.time[i])
        
    def remember_time(self):
        original_time = self.g_c.vp.time.get_array()
#        print(original_time)
        return original_time
