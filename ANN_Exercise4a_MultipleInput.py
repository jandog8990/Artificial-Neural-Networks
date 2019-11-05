#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 00:44:22 2019

Multiple input firing neurons with post-synaptic firing

@author: alejandrogonzales
"""
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
from math import exp

# --------------------------------------------------------
# Functions for running integration and synapse calcs
# --------------------------------------------------------

# Calculat output synapse for given weights and times
def calculateSynapse(t, V, ws, Es, ts, tau_s):
    Ps = exp(-1*(t - ts)/tau_s)
    synOut = ws*(V - Es)*Ps
    return synOut

Ie = 0.0
El = -0.065
V_reset = El
V_th = -0.05
tau_m = 0.03
Rm = 9e7

# Synaptic weight array
ws = 0.9
ws1 = ws
ws2 = ws
ws_arr = [ws1, ws2]

# Excitatory variable array (turns on or off the synapsis)
Es = 0
Es1 = Es
Es2 = Es
Es_arr = [Es1, Es2]

# Time constants for synapses
tau_s = 0.01
tau_s1 = tau_s
tau_s2 = tau_s
tau_s_arr = [tau_s1, tau_s2]

# Input spike arrival times
nidx = 2
ts = 0
ts1 = ts
ts2 = ts
ts_arr = [ts1, ts2]

# Numerical example for time period
T = 20
dt = 10e-4
jmax = int(T/dt);

# Input risis
risis = 5
risis1 = risis
risis2 = risis
tisis = 1/risis
#dt_shift = tisis/20
dt_shift = dt*.4
k = int(tisis/dt)

# Initialize output firing rate voltages
#V = [0]*jmax
V = El
tspike = list()
tisi = list()
risi = list()
ts_arr = list()

# Initialize the average and std dev lists
risiAver = list()
risiStd = list()

# Outer loop trials for risis values
isisCount = 0
for i in np.arange(0, tisis, dt_shift):
    # Re-init the voltages
#    if (len(V) > 0):
#        V.clear
#        V = [0]*jmax
        
    # Re-initialize the integration arrays
    if (len(tspike) > 0):
        tspike.clear
        tisi.clear
        risi.clear
        ts_arr.clear
        
    # Reset voltages and spikes
#    V[0] = El
    V = El
    tspike.append(0)
    tisi.append(0)
    risi.append(0)
    ts_arr.append(0)
    
    count = 0
    told = 0
    ts = 0
    ts1 = 0
    ts2 = 0
#    ts_arr = [ts1, ts2]
    
    # Delta j shift
    dj_shift = i/dt
    
    #print(str(i) + " : " + str(dj_shift))
    
    # Integration time steps
    for j in range(1, int(jmax)):
        t = j*dt
        
        mt = dt/tau_m
        dj_diff = int(j - dj_shift)
        
        # Input spike 1 time check
        if (j % k == 0):
            ts1 = t
#            ts_arr[0] = ts1
            
        # Input spike 2 time check
        if (dj_diff % k == 0):
            ts2 = t
#            ts_arr[1] = ts2
            
#        print("ts_arr:")
#        print(str(ts_arr))
#        print("\n")
            
        # Will calculate output synapse firing using single values
#        synout1 = calculateSynapse(t, V[j], ws1, Es1, ts1, tau_s1)
#        synout2 = calculateSynapse(t, V[j], ws2, Es2, ts2, tau_s2)
        synout1 = calculateSynapse(t, V, ws1, Es1, ts1, tau_s1)
        synout2 = calculateSynapse(t, V, ws2, Es2, ts2, tau_s2)

        print(str(ts1) + " : " + str(ts2))
        
#        synSum = 0
#        for sn in range(nidx):
#            ws = ws_arr[sn]
#            Es = Es_arr[sn]
#            ts = ts_arr[sn]
#            tau_s = tau_s_arr[sn]
#            Ps = exp(-1*(t - ts)/tau_s)
#            synSum = synSum + ws*(V[j] - Es)*Ps
                        
        # ODE equations for output spikes and synapses
#        V[j+1] = V[j] + mt*(-V[j] + El + Rm*Ie - (synout1 + synout2))
        V = V + mt*(-V + El + Rm*Ie - (synout1 + synout2))
        if (V >= V_th):
            V = V_reset
#        if (V[j] >= V_th):
#            V[j+1] = V_reset
            isi = t - told
            told = t
            tspike.append(t)
            tisi.append(isi)
            risi.append(1/isi)
            count = count + 1
            
    # Store the current and average post synaptic spikes
    if (len(risi) > 1):
        risiAver.append(mean(risi))
        risiStd.append(np.std(risi))
        isisCount = isisCount + 1
        
# create the index array
timeShift = np.arange(0, 1, dt_shift/tisis)
plt.figure()
plt.plot(timeShift, risiAver)
plt.show()