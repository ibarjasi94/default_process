#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:46:26 2019

@author: tester
"""
import numpy as np
from scipy import integrate
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.special import factorial
from scipy.special import binom
import time


def dyn_dist(k,z,n_model):
    if(n_model == "Poisson"):
        k = np.asarray(k)
        return np.power(z,k)*np.exp(-z)/factorial(k)
    elif(n_model == "k-regular"):
        k = np.asarray(k)
        dist = []
        length = len(k)
        for n,k in enumerate(k):
            if(n != int((length-1)/2)):
                dist.append(0)
            else:
                dist.append(1)
        return dist 
    elif(n_model == "scale free"):
        k = np.asarray(k,dtype = float)
        k = k[1:]
        dist = np.power(k,-3)
        dist = np.insert(dist,0,0)
        dist = dist/np.sum(dist)
        return dist

def kmvec(kmin,kmax):
    kvec = []
    for k in range(kmin,kmax+1):
        for n in range(1,k+2):
            kvec.append(k)
    mvec = []
    for k in range(kmin,kmax+1):
        for n in range(0,k+1):
            mvec.append(n)
    return kvec,mvec

#def lower_ten(ten,kin_min,kout_min):
#    dims = ten.shape
#    for m in range(0,dims[0]):
#        for n in range(0,dims[1]):
#            for i in range(0,dims[2]):
#                for j in range(0,dims[3]):
#                    if(m+kin_min < n or i+kout_min < j):
#                        ten[m,n,i,j] = 0
#    return ten

def lower_ten(ten,kin_min,kout_min):
    dims = ten.shape
    for m in range(0,dims[0]):
        if(m == dims[0]-1):
            mdim = 0
        else:
            mdim = m + 1
        ten[m+kin_min,mdim,:,:] = 0
    for i in range(0,dims[2]):
        if(i == dims[2]-1):
            idim = 0
        else:
            idim = i + 1
        ten[:,:,i+kout_min,idim] = 0
    return ten

#def lower_ten5(ten,kin_min,kout_min):
#    dims = ten.shape
#    for m in range(0,dims[0]):
#        for n in range(0,dims[1]):
#            for i in range(0,dims[2]):
#                for j in range(0,dims[3]):
#                    for k in range(0,dims[4]):
#                        if(m+kin_min < n or i+kout_min < j or i+kout_min < k):
#                            ten[m,n,i,j,k] = 0
#    return ten

#def roll5(ten,kin_min,kout_min,axis):
#    dims = ten.shape
#    if(axis == 3):
#        ten_r = np.roll(ten,1,axis = axis)
##        for m in range(0,dims[0]):
##            if(m == dims[0]-1):
##                mdim = 0
##            else:
##                mdim = m + 1
##            ten_r[m+kin_min,mdim,:,:,:] = 0
#        for i in range(0,dims[2]):
#            if(i == dims[2]-1):
#                idim = 0
#            else:
#                idim = i + 1
#            ten_r[:,:,i+kout_min,idim,:] = 0
#    if(axis == 4):
#        ten_r = np.roll(ten,1,axis = axis)
##        for m in range(0,dims[0]):
##            if(m == dims[0]-1):
##                mdim = 0
##            else:
##                mdim = m + 1
##            ten_r[m+kin_min,mdim,:,:,:] = 0
#        for i in range(0,dims[2]):
#            if(i == dims[2]-1):
#                idim = 0
#            else:
#                idim = i + 1
#            ten_r[:,:,i+kout_min,:,idim] = 0        
#    return ten_r


def roll5(ten,kin_min,kout_min,axis):
    dims = ten.shape
    ten_r = ten.copy(order = 'F')
    if(axis == 3):
        for i in range(0,dims[2]):
            ten_r[:,:,i+kout_min,i,:] = 0
    if(axis == 4):
        for i in range(0,dims[2]):
            ten_r[:,:,i+kout_min,:,i] = 0    
    ten_r = np.roll(ten_r,1,axis = axis)
    return ten_r

def init(ten):
    dims = ten.shape
    for n in range(1,dims[1]):
        ten[:,n,:,:] = 0
    for j in range(1,dims[3]):
        ten[:,:,:,j] = 0
    return ten

def kten(kin_min,kout_min,kin_max,kout_max):
    kin = kin_max - kin_min + 1
    kout = kout_max - kout_min + 1
    kinten = np.zeros((kin,kin_max+1,kout,kout_max+1))
    koutten = np.zeros((kin,kin_max+1,kout,kout_max+1))
    minten = np.zeros((kin,kin_max+1,kout,kout_max+1))
    moutten = np.zeros((kin,kin_max+1,kout,kout_max+1))
    dims = kinten.shape
    for m in range(0,dims[0]):
        for n in range(0,dims[1]):
            for i in range(0,dims[2]):
                for j in range(0,dims[3]):
                    if(m+kin_min >= n and i+kout_min >= j):
                        kinten[m,n,i,j] = m + kin_min
                        koutten[m,n,i,j] = i + kout_min
                        minten[m,n,i,j] = n 
                        moutten[m,n,i,j] = j 
    return kinten,koutten,minten,moutten   

def kten5(kin_min,kout_min,kin_max,kout_max):
    kin = kin_max - kin_min + 1
    kout = kout_max - kout_min + 1
    kinten5 = np.zeros((kin,kin_max+1,kout,kout_max+1,kout_max+1))
    koutten5 = np.zeros((kin,kin_max+1,kout,kout_max+1,kout_max+1))
    minten5 = np.zeros((kin,kin_max+1,kout,kout_max+1,kout_max+1))
    moutten5 = np.zeros((kin,kin_max+1,kout,kout_max+1,kout_max+1))
    zoutten5 = np.zeros((kin,kin_max+1,kout,kout_max+1,kout_max+1))
    dims = kinten5.shape
    for m in range(0,dims[0]):
        for n in range(0,dims[1]):
            for i in range(0,dims[2]):
                for j in range(0,dims[3]):
                    for k in range(0,dims[4]):
                        if(m+kin_min >= n and i+kout_min >= j and i+kout_min >= k):
                            kinten5[m,n,i,j,k] = m + kin_min
                            koutten5[m,n,i,j,k] = i + kout_min
                            minten5[m,n,i,j,k] = n 
                            moutten5[m,n,i,j,k] = j 
                            zoutten5[m,n,i,j,k] = k
    return kinten5,koutten5,minten5,moutten5,zoutten5

def mout_reset(ten,kin,kout,kin_max,kout_max):
    dims = ten.shape
    ten_new = np.zeros(dims)
    for m in range(0,dims[0]):
        for n in range(0,dims[1]):
            for i in range(0,dims[2]):
                for j in range(0,dims[3]):
                    el = ten_new[m,n,i-j,0]
                    ten_new[m,n,i-j,0] = el + ten[m,n,i,j]
    return ten_new
 
def z_reset(ten):
    dims = ten.shape
    ten_new = np.zeros((dims[0],dims[1],dims[2],dims[3],dims[3]))
    for m in range(0,dims[0]):
        for n in range(0,dims[1]):
            for i in range(0,dims[2]):
                for j in range(0,dims[3]):
                    ten_new[m,n,i,j,0] = ten[m,n,i,j]
    return ten_new    

def roll_array(array,axis,kin_min,kout_min):
    array_roll = np.roll(array,1,axis = axis)
    array_roll = lower_ten(array_roll,kin_min,kout_min)
    return array_roll
    

def func(y,t,z,kin_min,kout_min,kin_max,kout_max,alpha,beta,kinten,koutten,minten,moutten,kinten5,koutten5,minten5,moutten5,zoutten5,kinvec,koutvec,n_model,dynamics):
    #startTime = time.time()
    kin = kin_max - kin_min + 1
    kout = kout_max - kout_min + 1
    arr2len = kin*(kin_max+1)*kout*(kout_max+1)*(kout_max+1)
    arrlen  = len(y) - arr2len
    n = 2
    s = y[0:np.int(arrlen/n)]
    i = y[np.int(arrlen/n):np.int(2*arrlen/n)]
#    i_a = y[np.int(2*arrlen/n):np.int(3*arrlen/n)]
#    i_b = y[np.int(3*arrlen/n):np.int(4*arrlen/n)]
    ic = y[-kin*(kin_max+1)*kout*(kout_max+1)*(kout_max+1):]
    s = np.reshape(s,(kin,kin_max+1,kout,kout_max+1))
    i = np.reshape(i,(kin,kin_max+1,kout,kout_max+1))
#    i_a = np.reshape(i_a,(kin,kin_max+1,kout,kout_max+1))
#    i_b = np.reshape(i_b,(kin,kin_max+1,kout,kout_max+1))
    ic = np.reshape(ic,(kin,kin_max+1,kout,kout_max+1,kout_max+1))
    kmminten = kinten - minten
    kmmoutten = koutten - moutten
    kmmoutten5 = koutten5 - moutten5
    if(dynamics == "SI"):
        Ften_in = alpha + beta*np.asarray(minten)
    elif(dynamics == "VM"):
        Ften_in = alpha + beta*np.asarray(np.divide(minten,kinten,out=np.zeros_like(minten), where=kinten!=0))
#    Ften_in_a = alpha
#    Ften_in_b = beta*np.asarray(minten)
    kimiy = np.multiply(kmminten,s)
    komoy = np.multiply(kmmoutten,s)
    moy = np.multiply(moutten,s)
    miy = np.multiply(minten,s)
    B_s_num = np.dot(np.dot(np.sum(np.multiply(komoy,Ften_in),axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    B_s_denom = np.dot(np.dot(np.sum(komoy,axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    C_s_num = np.dot(np.dot(np.sum(np.multiply(kimiy,Ften_in),axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    C_s_denom = np.dot(np.dot(np.sum(kimiy,axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    B_i_num = np.dot(np.dot(np.sum(np.multiply(moy,Ften_in),axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    B_i_denom = np.dot(np.dot(np.sum(moy,axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    C_i_num = np.dot(np.dot(np.sum(np.multiply(miy,Ften_in),axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    C_i_denom = np.dot(np.dot(np.sum(miy,axis = (1,3)),dyn_dist(koutvec,z,n_model)),dyn_dist(kinvec,z,n_model))
    if(B_s_denom == 0):
        B_s = 0
    else:
        B_s = B_s_num/B_s_denom
    if(C_s_denom == 0):
        C_s = 0
    else:
        C_s = C_s_num/C_s_denom
    if(B_i_denom == 0):
        B_i = 0
    else:
        B_i = B_i_num/B_i_denom
    if(C_i_denom == 0):
        C_i = 0
    else:
        C_i = C_i_num/C_i_denom
    sminp = roll_array(s,1,kin_min,kout_min)
    smoutp = roll_array(s,3,kin_min,kout_min)
    iminp = roll_array(i,1,kin_min,kout_min)
    imoutp = roll_array(i,3,kin_min,kout_min)
#    iminp_a = roll_array(i_a,1,kin_min,kout_min)
#    imoutp_a = roll_array(i_a,3,kin_min,kout_min)
#    iminp_b = roll_array(i_b,1,kin_min,kout_min)
#    imoutp_b = roll_array(i_b,3,kin_min,kout_min)
#    print("tu?")
#    icmoutp = np.roll(ic,1,axis = 3)
#    icmoutp[:,:,:,0,:] = 0
#    icmoutp = lower_ten5(icmoutp,kin_min,kout_min)
#    icmoutp = np.roll(icmoutp,1,axis = 4)
#    icmoutp[:,:,:,:,0] = 0
#    icmoutp = lower_ten5(icmoutp,kin_min,kout_min)
    icmoutp1 = roll5(ic,kin_min,kout_min,3)
    icmoutp = roll5(icmoutp1,kin_min,kout_min,4)
    dsdt = B_s*(np.multiply((kmminten+1),sminp)-np.multiply(kmminten,s)) + C_s*(np.multiply((kmmoutten+1),smoutp)-np.multiply(kmmoutten,s)) - Ften_in*s
    didt = B_i*(np.multiply((kmminten+1),iminp)-np.multiply(kmminten,i)) + C_i*(np.multiply((kmmoutten+1),imoutp)-np.multiply(kmmoutten,i)) + Ften_in*s
#    didt_a = B_i*(np.multiply((kmminten+1),iminp_a)-np.multiply(kmminten,i_a)) + C_i*(np.multiply((kmmoutten+1),imoutp_a)-np.multiply(kmmoutten,i_a)) + Ften_in_a*s
#    didt_b = B_i*(np.multiply((kmminten+1),iminp_b)-np.multiply(kmminten,i_b)) + C_i*(np.multiply((kmmoutten+1),imoutp_b)-np.multiply(kmmoutten,i_b)) + Ften_in_b*s
    dicdt = z_reset(np.multiply(Ften_in,s)) + C_i*(np.multiply((kmmoutten5+1),icmoutp)-np.multiply(kmmoutten5,ic))
    dsdt = np.reshape(dsdt,(kin*(kin_max+1)*kout*(kout_max+1)))
    didt = np.reshape(didt,(kin*(kin_max+1)*kout*(kout_max+1)))
#    didt_a = np.reshape(didt_a,(kin*(kin_max+1)*kout*(kout_max+1)))
#    didt_b = np.reshape(didt_b,(kin*(kin_max+1)*kout*(kout_max+1)))
    dicdt = np.reshape(dicdt,(kin*(kin_max+1)*kout*(kout_max+1)*(kout_max+1)))
    dydt = np.concatenate((dsdt,didt,dicdt),axis = 0)
    #print ('The script took {0} second !'.format(time.time() - startTime),"ame")
    return dydt

def odeint_vol2(dfunc, V0, times, args=()):
    r = integrate.ode(lambda t, X: dfunc(X, t, *args))
    r.set_integrator('vode', method='adams')
    r.set_initial_value(V0,times[0])
    V=[V0]
    for time in times[1:]:
        V.append(r.integrate(time))
    V = np.array(V)
    return V

def AME_solve(kin_min,kout_min,kin_max,kout_max,exp_deg,N_vertices,zeta,t,n_model,dynamics):
    startTime = time.time()
    kin = kin_max - kin_min + 1
    kout = kout_max - kout_min + 1
    kinvec = [n for n in range(kin_min,kin_max+1)]
    koutvec = [n for n in range(kout_min,kout_max+1)]
    alpha = zeta
    beta = 1
    exp_deg = 0.5*exp_deg
    kinten,koutten,minten,moutten = kten(kin_min,kout_min,kin_max,kout_max)
    kinten5,koutten5,minten5,moutten5,zoutten5 = kten5(kin_min,kout_min,kin_max,kout_max)
    s0 = np.full((kin,kin_max+1,kout,kout_max+1),1)
    s0 = init(s0)
    s0 = np.reshape(s0,((kin*(kin_max+1)*kout*(kout_max+1))))
    i0 = np.full((kin*(kin_max+1)*kout*(kout_max+1)),0)
#    i0_a = np.full((kin*(kin_max+1)*kout*(kout_max+1)),0)
#    i0_b = np.full((kin*(kin_max+1)*kout*(kout_max+1)),0)
    ic0 = np.full((kin*(kin_max+1)*kout*(kout_max+1)*(kout_max+1)),0)
    n = 2
    #y0 = np.concatenate((s0,i0,i0_a,i0_b,o0,ic0),axis = 0)
    y0 = np.concatenate((s0,i0,ic0),axis = 0)
    sol = odeint_vol2(func,y0,t,args=(exp_deg,kin_min,kout_min,kin_max,kout_max,\
                                    alpha,beta,kinten,koutten,minten,moutten,\
                                    kinten5,koutten5,minten5,moutten5,zoutten5,\
                                    kinvec,koutvec,n_model,dynamics))
    print ('The script took {0} second !'.format(time.time() - startTime),"ame")
    arr2len = kin*(kin_max+1)*kout*(kout_max+1)*(kout_max+1)
    arrlen = sol.shape[1] - arr2len
    ssol = sol[:,0:np.int(arrlen/n)]
    isol = sol[:,np.int(arrlen/n):np.int(2*arrlen/n)]
    #isol_a = sol[:,np.int(2*arrlen/n):np.int(3*arrlen/n)]
    #isol_b = sol[:,np.int(3*arrlen/n):np.int(4*arrlen/n)]
    icsol = sol[:,-arr2len:]
    rhosol = []
    rhosol_a = []
    rhosol_b = []  
    SIsol = []
    onepath = []
    onepath_var = []
    onepath_s = []
    onepath_var_s = []
    twopath = []
    twopath_var = []
    twopath_s = []
    twopath_var_s = []
    bin_in = binom(minten,2)
    bin_in[bin_in <=  0] = 0
    bin_out = binom(moutten,2)
    bin_out[bin_out <= 0] = 0
    onesum = 0
    twosum = 0
    for n in range(0,sol.shape[0]):
    #for n in range(1,2):
        #print(isol[n,:].shape)
        sreshape = np.reshape(ssol[n,:],(kin,kin_max+1,kout,kout_max+1))
        ireshape = np.reshape(isol[n,:],(kin,kin_max+1,kout,kout_max+1))
#        ireshape_a = np.reshape(isol_a[n,:],(kin,kin_max+1,kout,kout_max+1))
#        ireshape_b = np.reshape(isol_b[n,:],(kin,kin_max+1,kout,kout_max+1))
        icreshape = np.reshape(icsol[n,:],(kin,kin_max+1,kout,kout_max+1,kout_max+1))
        SIsum = np.dot(np.dot(np.sum(np.multiply(minten,sreshape),axis = (1,3)),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        rhosum = np.dot(np.dot(np.sum(icreshape,axis = (1,3,4)),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
#        rhosum_a = N_vertices*np.dot(np.dot(np.sum(ireshape_a,axis = (1,3)),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
#        rhosum_b = N_vertices*np.dot(np.dot(np.sum(ireshape_b,axis = (1,3)),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        one_rv = np.sum(np.multiply(zoutten5,icreshape),axis = (1,3,4))
#        one_in = np.sum(one_rv, axis = 1)
#        print(one_in)
#        print(np.mean(one_in))
#        one_in_mean = np.dot(one_in,dyn_dist(kinvec,exp_deg,n_model))
#        print(one_in_mean)
#        one_in_var = np.dot(np.power(one_in-one_in_mean,2),dyn_dist(kinvec,exp_deg,n_model))
#        print(one_in_var)
#        print(np.var(one_in))
        one_s_rv = 0.5*np.sum(np.multiply(moutten,np.sum(icreshape,axis = 4)),axis = (1,3))
        two_rv = np.sum(np.multiply(zoutten5,np.multiply(minten5,icreshape))+np.multiply(binom(minten5,2),icreshape)+np.multiply(binom(zoutten5,2),icreshape),axis = (1,3,4))
        two_s_rv = 0.25*np.sum(np.multiply(moutten,np.multiply(minten,ireshape))+np.multiply(bin_in,ireshape)+np.multiply(bin_out,ireshape),axis = (1,3))
        onesum = np.dot(np.dot(one_rv,dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        onevar = np.dot(np.dot(np.power((one_rv-onesum),2),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        onesum_s = np.dot(np.dot(one_s_rv,dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        onevar_s = np.dot(np.dot(np.power((one_s_rv-onesum_s),2),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        twosum = np.dot(np.dot(two_rv,dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        twovar = np.dot(np.dot(np.power((two_rv-twosum),2),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        twosum_s  = np.dot(np.dot(two_s_rv,dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        twovar_s = np.dot(np.dot(np.power((two_s_rv-twosum_s),2),dyn_dist(koutvec,exp_deg,n_model)),dyn_dist(kinvec,exp_deg,n_model))
        SIsol.append(N_vertices*SIsum)
        onepath.append(N_vertices*onesum)
        onepath_s.append(N_vertices*onesum_s)
        onepath_var.append(np.sqrt(N_vertices*onevar))
        onepath_var_s.append(np.sqrt(N_vertices*onevar_s))
        twopath.append(N_vertices*twosum)
        twopath_s.append(N_vertices*twosum_s)
        twopath_var.append(np.sqrt(N_vertices*twovar))
        twopath_var_s.append(np.sqrt(N_vertices*twovar_s))
        rhosol.append(N_vertices*rhosum)
#        rhosol_a.append(rhosum_a)
#        rhosol_b.append(rhosum_b)
#        print(twosum)
#        print(twovar)
#        print(onesum)
#        print(onevar)
    print(onesum_s*2)
    print(twosum_s*2)
    onepath = onepath/(N_vertices*onesum_s*2)
    onepath_var = onepath_var/(N_vertices*onesum_s*2)
    onepath_s = onepath_s/(N_vertices*onesum_s*2)
    onepath_var_s = onepath_var_s/(N_vertices*onesum_s*2)
    twopath = twopath/(N_vertices*twosum_s*4)
    twopath_var = twopath_var/(N_vertices*twosum_s*4)
    twopath_s = twopath_s/(N_vertices*twosum_s*4)
    twopath_var_s = twopath_var_s/(N_vertices*twosum_s*4)
    print ('The script took {0} second !'.format(time.time() - startTime),"end")
    return SIsol,rhosol,rhosol_a,rhosol_b,onepath,onepath_var,onepath_s,onepath_var_s,twopath,twopath_var,twopath_s,twopath_var_s
