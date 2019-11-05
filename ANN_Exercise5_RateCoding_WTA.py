#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 19:52:02 2019

Rate Coding, Winner Takes All (Multi-node Network)

@author: alejandrogonzales
"""
import numpy as np

# Problem 1: 4 neuron network
n = 4
tau_s = 0.1
dt = 0.001
T = 1
jmax = int(T/dt)
Ntest = 10

# Initialize the node voltages 
#V = np.zeros([4, jmax])
#Y = np.zeros([4, jmax])
V = np.zeros([4, Ntest])
Y = np.zeros([4, Ntest])
Y[0,:] = 1.0                # input Y1 in paper
Y[1,:] = 0.95               # input Y2 in paper
print("Voltage matrix size = " + str(V.shape))
print("Output matrix size = " + str(Y.shape))
print("Output Y0 size = " + str(len(Y[0,:])))
print(Y)
print("\n")

# Initialize inputs y1 and y2 neurons
y1 = 1.0
y2 = 0.95

# Initialize outputs y3 and y4 neurons
y3 = 0
y4 = 0

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
        
# Outer loop for kth neuron calculation
for k in range(0, n):
    for j in range(Ntest-1):
        t = j*dt
        
        # time constant calcs
        st = dt/tau_s
        
        V[k,j+1] = V[k,j] + 4
        Y[k,j+1] = Y[k,j] + 1
        
print("Output Voltage Matrix:")
print(V)
print("\n")

print("Output Y Matrix:")
print(Y)
print("\n")
        