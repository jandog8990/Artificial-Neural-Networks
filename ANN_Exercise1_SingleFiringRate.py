#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:16:48 2019

Basic Integrate and Fire Neuron

@author: alejandrogonzales
"""
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statistics import mean

dt = 0.001 # step size integration loop
N = 1000 
tau = 0.03

El = -0.065
V_reset = El
V_th = -0.05
Rm = 9e7
Ie = 0
dIe = 0.0002 # current increment for each trial
numTrials = 475

# Initialize output firing rate array
#ieArr = [0]*numTrials;
#risiAver = [0]*numTrials
ieArr = np.zeros(numTrials-1)
risiAver = np.zeros(numTrials-1)

# Init interspike arrays for integration
V = [0]*N
tspike = list()
tisi = list()
risi = list()

# Loop through outer trials for current
for i in range(numTrials-1):
    # Re-init the voltages
	if (len(V) > 0):
		V.clear
		V = [0]*N
    
    # Re-initialize the integration arrays
	if (len(tspike) > 0):
		tspike.clear
		tisi.clear
		risi.clear
   
    # Current current
	Ie = i*dIe;
    
    # Reset voltages and spikes
	V[0] = El
	tspike.append(0)
	tisi.append(0)
	risi.append(0)
    
	count = 0
	t = 0
	told = 0
    
    # Integration loop to calculte synapitc voltages
	for j in range(N-1):
		t = j*dt
		ct = dt/tau
        
		V[j+1] = V[j] + ct*(-V[j] + El + Ie)
		if (V[j] >= V_th):
			V[j+1] = V_reset
			isi = t - told
			told = t
			tspike.append(t)
			tisi.append(isi)
			risi.append(1/isi)
			count = count + 1
            
    # Store the current and average post synaptic spikes
	if (len(tisi) > 1):
		ieArr[i] = Ie
		risiAver[i] = mean(risi)
		print("Ie = " + str(Ie))
		print("Risi mean = " + str(risiAver[i]))
		print("\n")

# Plots for a single output synapse
#plt.figure()
plt.plot(ieArr, risiAver)
plt.show()

## Post synaptic voltages
#plt.figure()
#idx = range(N)
#plt.plot(idx, V)
#
## Spike times for integration
#plt.figure()
#plt.stem(tspike, np.ones(len(tspike)))
#
## Post synapitc firing rate
#plt.figure()
#plt.plot(tspike, risi)
