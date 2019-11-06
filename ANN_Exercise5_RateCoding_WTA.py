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
y1 = 1.0
y2 = 0.5
Y[0,:] = y1                # input Y1 in paper
Y[1,:] = y2                # input Y2 in paper
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
print("\n")

# Non-linear initialization
a = 0.1
phi = list()
varr = list()
for v in np.arange(-2.0, 2.0, 0.1):
    if v >= 0:
        phiv = (v**2)/(a+v**2)
        phi.append(phiv)
        varr.append(v)

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

print("Varr len = " + str(len(varr)))
print("Phi V len = " + str(len(phi)))
plt.figure()
plt.plot(varr, phi)
plt.title("Phi(v) Output Function")
plt.ylabel("Phi(v)")
plt.xlabel("v (internal activation)")
plt.show()

# Plot the output firing rates for two output neurons y2 (blue) and y3 (red)
y3 = 1*Y[2,:]
y4 = 1.5*Y[3,:]
title1 = "WTA Output Firing Rates vs Time (y1 = {}, y2 = {})".format(y1, y2)
plt.figure()
plt.plot(t, y3, 'b', t, y4, 'r')
plt.title(title1)
plt.ylabel("Output risi y3(blue) & y4(red)")
plt.xlabel("Time (t)")
plt.show()

# Scatter plot for output firing rates (y3 & y4)
plt.figure()
plt.scatter(y3, y4, alpha=0.5, s=1)
plt.title("WTA Scatter Output Firing Rates (y3 & y4)")
plt.ylabel("Output risi y4")
plt.xlabel("Output risi y3")
plt.show()