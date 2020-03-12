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
from scipy.integrate import odeint


#modul stvara Erdos-Renyi graf te sadrži sve potrebno za propadanje čvorova u alpha i beta procesima
#reshuffle vremena tih propadanja te mjerenje varijabli na grafu (lc,op,tp,thp)


def subgraph_num(subgraph,graph):
    i = 0
    gen = gt.topology.subgraph_isomorphism(subgraph, graph, generator=True)
    for iter in gen:
        i += 1
    return i

class Erdos_Renyi():
    def __init__(self,N_vertices,exp_deg,p,zeta):
        #ERDOS RENYI
        #rnd.seed(12345)
        self.g = gt.Graph()
        for i in range(0,N_vertices):
            self.g.add_vertex()
            for j in range(0,i):
                if(rnd.uniform(0,1) < p):
                    self.g.add_edge(i,j)
                if(rnd.uniform(0,1) < p):
                    self.g.add_edge(j,i)
        print(self.g.num_edges())
        #LATTICE
#        row = 40
#        col = 25
#        self.g = gt.Graph()
#        self.g.add_vertex(N_vertices-1)
#        for r in range(0,row-1):
#            for c in range(0,col-1):
#                if(rnd.uniform(0,1) < 0.5):
#                    self.g.add_edge(c+r*col,c+r*col+1)
#                else:
#                    self.g.add_edge(c+r*col+1,c+r*col)
#                if(rnd.uniform(0,1) < 0.5):
#                    self.g.add_edge(c+r*col,c+(r+1)*col)
#                else:
#                    self.g.add_edge(c+(r+1)*col,c+r*col)
#        for c in range(0,col-1):
#            if(rnd.uniform(0,1) < 0.5):
#                self.g.add_edge(c+(row-1)*col,c+(row-1)*col+1)
#            else:
#                self.g.add_edge(c+(row-1)*col+1,c+(row-1)*col)
#        for r in range(0,row-1):
#            if(rnd.uniform(0,1) < 0.5):
#                self.g.add_edge(col-1+r*col,col-1+(r+1)*col)
#            else:
#                self.g.add_edge(col-1+(r+1)*col,col-1+r*col)
#        gt.draw.graph_draw(self.g, vertex_text = self.g.vertex_index,vertex_font_size=12, output_size=(600, 600))
        self.g.vp.stecaj = self.g.new_vertex_property("bool")
        self.g.vp.known = self.g.new_vertex_property("bool")
        self.g.vp.time = self.g.new_vertex_property("double")
        self.g.ep.time = self.g.new_edge_property("double")
        self.g.vp.t_part = self.g.new_vertex_property("bool")
        self.g.ep.causal = self.g.new_edge_property("bool")
        #self.timegrid = np.linspace(0,N_vertices/zeta,N_vertices/zeta*100+1)
        #gt.draw.graph_draw(self.g, vertex_text = self.g.vertex_index, vertex_halo = self.g.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
        samci = []
        for v in self.g.get_vertices():
            if(len(self.g.get_in_edges(v)) == 0):
                samci.append(v)
        self.n_samci = len(samci)
#        self.N_vert = self.g.num_vertices()
#        self.degree = (2*self.g.num_edges()/self.g.num_vertices())
#        self.g.remove_vertex(samci)

    def set_directedness(self,directedness):
        self.g_c.set_directed(directedness)

    def initial_properties(self, zeta, dynamics):
        self.g_c = self.g.copy()
        N = self.g_c.num_vertices()
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_halo = self.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
        vert_times = self.poisson_interevent(N,zeta)
        #print(vert_times)
#        print(vert_interevent,"vert")
        #vert_times = [np.sum(vert_interevent[0:i]) for i in range(0,len(vert_interevent))]
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
        #print(edge_interevent)
#        print(edge_interevent, "beta times")
        for e in self.g_c.edges():
            if(dynamics == "VM"):
                s,t = e
                kin = t.in_degree()
                self.g_c.ep.time[e] = edge_interevent[-1]*kin
            elif(dynamics == "SI"):
                self.g_c.ep.time[e] = edge_interevent[-1]
            edge_interevent = edge_interevent[:-1]
#            print(self.g_c.ep.time[i])
        for e in self.g_c.edges():
            self.g_c.ep.causal[e] = False
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_halo = self.g_c.vp.stecaj,vertex_font_size=12, output_size=(600, 600))
#        print(self.g_c.num_edges(),"edges",self.g_c.num_vertices(),"vertices")

    def process(self,exp_deg,zeta,N_vertices):
        index_time_sorted = collections.OrderedDict(sorted(self.index_time.items()))
#        print(index_time_sorted,"početna vremena")
        self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = True)
        nondefault = self.g_c.get_vertices()
        self.g_c.clear_filters()
        alpha = 0
        alpha_time = []
        beta = 0
        beta_time = []
        tot = 0
        tot_time = []
        pot = []
        N_pot = 0
        pot.append([0,0])
        L_a_p = []
        L_a_p.append([0,0])
        L_b_p = []
        L_b_p.append([0,0])
        L_a_n = []
        L_a_n.append([0,0])
        L_b_n = []
        L_b_n.append([0,0])
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
#            print(index_time_sorted,"obrisani propali")
            #bridovi koji potencijalno vode zarazu do čvora
            self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = False)
            in_edges = self.g_c.get_in_edges(vert)
            self.g_c.clear_filters()
            #nastavak koda
            self.g_c.set_vertex_filter(self.g_c.vp.stecaj, inverted = True)
            out_edges = self.g_c.get_out_edges(vert)
            if(prcs == "alpha"):
                alpha += 1
                tot += 1
                alpha_time.append([vert_time,alpha])
                L_a_p.append([vert_time,len(out_edges)])
                L_b_p.append([vert_time,np.nan])
                L_a_n.append([vert_time,len(in_edges)])
                L_b_n.append([vert_time,np.nan])
            else:
                beta += 1
                tot += 1
                beta_time.append([vert_time,beta])
                L_a_p.append([vert_time,np.nan])
                L_b_p.append([vert_time,len(out_edges)])
                L_a_n.append([vert_time,np.nan])
                L_b_n.append([vert_time,len(in_edges)])
            tot_time.append([vert_time,tot])
            N_pot = N_pot + len(out_edges) - len(in_edges)
            for edge in out_edges:
                s,t,i = edge
                edge_time = self.g_c.ep.time[edge]
                cascade_time = vert_time + edge_time
                index_time_sorted[cascade_time] = [t,"beta"]
            pot.append([vert_time,N_pot])
#            print("propao", vert, "u", vert_time)
            self.g_c.vp.time[vert] = vert_time
            self.g_c.vp.stecaj[vert] = True
            nondefault = self.g_c.get_vertices()
            self.g_c.clear_filters()
            index_time_sorted = collections.OrderedDict(sorted(index_time_sorted.items()))
#            n += 1
#        graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_font_size=12, output_size=(600, 600))
        time_list = self.g_c.vp.time.get_array()
        pot = np.asarray(pot)
        L_a_p = np.asarray(L_a_p)
        L_b_p = np.asarray(L_b_p)
        L_a_n = np.asarray(L_a_n)
        L_b_n = np.asarray(L_b_n)
        L = [L_a_p,L_b_p,L_a_n,L_b_n]
        alpha_time = np.asarray(alpha_time)
        beta_time = np.asarray(beta_time)
        tot_time = np.asarray(tot_time)
        return pot, alpha_time, beta_time, tot_time, L

    def network_part(self,time_stops,l):
        time_list_p = list(self.g_c.vp.time.get_array())
        time_list_p.sort()
        last = time_list_p[int((l+1)/time_stops*self.g_c.num_vertices())-1]
        for v in self.g_c.vertices():
            # < or <=??
            if(self.g_c.vp.time[v] <= last):
                self.g_c.vp.t_part[v] = True
        self.g_c.set_vertex_filter(self.g_c.vp.t_part)
        return last
        
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
        #gt.draw.graph_draw(self.g_c, vertex_text = self.g_c.vertex_index, vertex_font_size=12, output_size=(300, 300))
        #print("all default links")
#        for v in self.g_c.vertices():
#            print(v,self.g_    
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