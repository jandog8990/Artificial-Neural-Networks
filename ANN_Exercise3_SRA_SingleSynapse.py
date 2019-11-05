#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 22:56:18 2019

Spike Rate Adaptation and Single Synapse

@author: alejandrogonzales
"""
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
from math import exp

ws = 0.9
V_reset = -0.065
El = V_reset
Ek = -0.08; #-0.07
V_th = -0.05

# Excitatory/inhibitory
Es = 0.0
#Es = -0.08

# Time constants
tau_m = 0.03
tau_sra = 0.1
tau_s = 0.01
dgsra = 1.0

# Resistance current and voltage
Rm = 9e7
Ie = 0
IRme = Rm*Ie

# Time initialization and interspike interval synapsis
Nisis = 50 # number of trials for risis

# Integration time constants
T = 40
#dt = 1e-4
dt = 0.0001
jmax = int(T/dt)

# Risis input firing rate
risis = 20
tisis = 1/risis
k = tisis/dt


# Initialize output firing rate voltages
V = [0]*jmax
gsra = [0]*jmax
tspike = list()
tisi = list()
risi = list()

# Initialize the average and std dev lists
risiAver = list()
risiStd = list()

# Outer loop trials for risis values
isisCount = 0
for i in range(1, Nisis+1):
    # Re-init the voltages
    if (len(V) > 0):
        V.clear
        V = [0]*jmax
        
    # Re-initialize the integration arrays
    if (len(tspike) > 0):
        tspike.clear
        tisi.clear
        risi.clear
        
    # Input spike times and rates
    risis = i
    tisis = 1/risis
    k = int(tisis/dt)
    
    # Reset voltages and spikes
    V[0] = El
    tspike.append(0)
    tisi.append(0)
    risi.append(0)
    
    count = 0
    t = 0
    ts = 0
    told = 0
    gsra = 0
    
    # Loop for integration time steps
    for j in range(int(jmax)-1):
        t = j*dt
        
        # time constant calcs
        mt = dt/tau_m
        st = dt/tau_sra
        
        if (j % k == 0):
            ts = t
        #print("t = " + str(t) + ", ts = " + str(ts) + ", k = " + str(k) + ", risis = " + str(risis))
        
        # ODE equations for post synaptic spikes
        Ps = exp(-1*(t - ts)/tau_s)
        V[j+1] = V[j] + mt*(-V[j] + El + Rm*Ie - gsra*(V[j] - Ek) - ws*(V[j] - Es)*Ps)
        if (V[j] >= V_th):
            V[j+1] = V_reset
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
     
idx = range(len(risiAver))
plt.figure()
plt.plot(idx, risiAver, idx, risiStd)
plt.show()