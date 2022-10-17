#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from scipy import *
import time

N = 1000
L = 50.0
K = 20.0

a = (L*2)/N
b = 1/(a**2)

x = np.empty([N])
v = np.empty([N])
M = np.zeros((N,N))
En = np.empty([N])
Error = np.empty([N])

################ Constructing Matrix ###############
# X
for i in range(0, N):
    if (i == 0):
        x[i] = -L
    else:
        x[i] = x[i-1] + a
    
# potential
for i in range(0, N):
    if (x[i] < -K ) or (x[i] > K):
        v[i] = 1
    else:
        v[i] = 0

# Matrix
for i in range(0, N):
    
    M[i,i] = b + v[i]
    
    if (i != 0):
        M[i, i-1] = -0.5 * b
    
    if (i != (N-1)):
        M[i, i+1] = -0.5 * b
        
# Theoretical Energies
for i in range(0, N):
    En[i] = (((i+1)**2) *  ((np.pi)**2)) / (8 * ((2*L)**2))

###Eigenvalues and Eigenvectors of Matrix
eval, evec = np.linalg.eigh(M)

######### Plotting particle in a box #########
plt.plot(x, v, label = "Potential")
plt.title('Potential')
plt.xlabel('x')
plt.ylabel('Potential')
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.show()

############## WaveFunctions ##################

###Writing file for eigenvectors
f1 = open("Matrix2.txt", "w")

for i in range (0,10):
    for j in range (0,N): 
        f1.write(str(x[j]) + ' ' + str(evec[j,i]) + '\n')     
    f1.write('\n\n')
f1.close()


### PLotting Wavefunction 
data1 = np.loadtxt('Matrix2.txt')
for i in range(0, 5):
    plt.plot(x, v, label = "Potential")
    plt.plot(data1[(i*N):(N*(i+1)), 0], (data1[(i*N):(N*(i+1)), 1]), 'b', label = "Wavefunction")
    plt.ylim(-0.1, 0.1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Wavefunction')
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    plt.show()
    
########### Plotting Eigenvalues  ###############

### Plotting/ Table N vs. Eigenvalues
f2 = open("Eigenvalues2.txt", "w")
for i in range (0,N):
    f2.write(str(i+1) + ' ' + str((eval[i])/eval[0]) + '\n')   
f2.close()

f3 = open("Theoretical Energy.txt", "w")
for i in range (0,N):
    f3.write(str(i+1) + ' ' + str((En[i])/En[0]) + '\n') 
f3.close()

### Plotting Theoretical and Calculated Energies

data2 = np.loadtxt('Eigenvalues2.txt')
data3 = np.loadtxt('Theoretical Energy.txt')
num = 300

plt.plot(data2[:num, 0], data2[:num, 1], label = "Calculated")
plt.plot(data3[:num, 0], data3[:num, 1], label = "Theoretical")
plt.title('Theoretical and Calculated Energy')
plt.xlabel('x')
plt.ylabel('Energy')
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.show()

############## Error ###########################################

f4 = open("EnergyError.txt", "w")

for i in range (0,N):
    Error[i] = abs(100 * ((eval[i]/eval[0]) - (En[i]/En[0])) / (En[i]/En[0]))
    
    f4.write(str(i+1) + ' ' + str(Error[i]) + '\n')

f4.close()


data4 = np.loadtxt('EnergyError.txt')
    
for i in range(0, 300):
    plt.plot(data4[:i, 0], data4[:i, 1])
    plt.title('Error between Theoretical and Calculated Energy')
    plt.xlabel('N')
    plt.ylabel('Error')
plt.show()