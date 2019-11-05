#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 19:52:02 2019

Rate Coding, Winner Takes All (Multi-node Network)

@author: alejandrogonzales
"""
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------
# Functions for running integration and synapse calcs
# --------------------------------------------------------

# Calculate the weight sum in the euler method voltages
# Input:
#   Y - Output voltages for nodes 1 to n
#   w - Weights for the current neuron k
#   n - total number of neurons
def weightSum(Y, w, n):
    #print("Weight Sum:")
    #print("Y = " + str(Y))
    #print("w = " + str(w))
    #print("\n")
    
    wsum = 0
    for j in range(n):
        Yj = Y[j]
        wj = w[j]
        wsum = wsum + wj*Yj
        
    return wsum
    
# --------------------------------------------------------

# Problem 1: 4 neuron network
n = 4
tau_s = 0.22
dt = 0.0001
T = 0.5
jmax = int(T/dt)

# Initialize the node voltages 
V = np.zeros([4, jmax])
Y = np.zeros([4, jmax])
t = np.zeros(jmax)
Y[0,:] = 1.0                # input Y1 in paper
Y[1,:] = 0.95               # input Y2 in paper
print("Voltage matrix size = " + str(V.shape))
print("Output matrix size = " + str(Y.shape))
print("Output Y0 size = " + str(len(Y[0,:])))
print(Y)
print("\n")

# Initialize the weights
w = np.zeros([4,4])
w[2,0] = 1
w[2,2] = 1
w[2,3] = -1.5
w[3,1] = 1
w[3,2] = -1.5
w[3,3] = 1


print("Weights:")
print(w)

# Non-linear initialization
a = 0.1
phi = list()
for v in np.arange(-2.0, 2.0, 0.1):
    if v >= 0:
        phiv = (v**2)/(a+v**2)
        phi.append(phiv)
v = 0

# Outer loop for kth neuron calculation
for k in range(2,n):
    for i in range(jmax-1):
        t[i+1] = i*dt
        
        # time constant calcs
        st = dt/tau_s
        
        # Where the Y[:,j] vector is column for all neurons at time i
        # and w[k,:] vector of weights for the current neuron k
        V[k,i+1] = V[k,i] + st*(-V[k,i] + weightSum(Y[:,i], w[k,:], n))
        
        # Check the voltage of the next step
        if (V[k,i+1] >= 0):
            v = V[k,i+1]
            Y[k,i+1] = v**2/(a+v**2)
        else:
            Y[k,i+1] = 0

print("Output Voltage Matrix:")
print(V)
print("\n")

print("Output Y Matrix:")
print(Y)
print("\n")

# Plot the output firing rates for two output neurons y2 (blue) and y4 (red)
plt.plot(t, Y[2,:], 'b', t, 1.3*Y[3,:], 'r')
plt.show()