#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz

# N = number of points
N = 400
# h = step size
h = 100.0 / N
# L = x goes from -L to L
L = 50
# K = potential barriers
K = 20
# C = constant for energies
C = (np.pi**2) / (8.0*(20.0**2))
# Constants
q = 1.0 / (h**2)
u = -0.5*q
sqrt = np.sqrt(1.0/K)
frac = np.pi/(2.0*K)

# Arrays
V = np.empty([N+1])
M = np.empty([N+1,N+1])
O = np.empty([N+1])
error1 = np.empty([N+1])

# Filling up values for x
x = np.linspace(-L,L,N+1)

# Filling up values for Potential
for i in range(0,N+1):
    if x[i] < -K or x[i] > K:
        V[i] = 1000000000
    else:
        V[i] = 0

# Filling up values for matrix M
for i in range (0,N+1):
    for j in range (0,N+1):
        if j == i:
            M[i][j] = V[i] + q
        elif j == i-1 or j == i+1:
            M[i][j] = u
        else:
            M[i][j] = 0

# Diagonalizing matrix M
spectrum,v=np.linalg.eigh(M)

# Defining function for calculated wavefunction
def calcpsi(j):
    for i in range(0,N+1):
        O[i] = v[i][j-1]
    return O

# Defining function for normalized calculated wavefunction
def psi(j):
    A = 1 / np.sqrt(trapz(calcpsi(j)**2,x))
    for i in range (0,N+1):
        O = A * calcpsi(j)
    return O

# Defining function for theoretical wavefunction
def wave(j):
    if j%2 == 0:
        for i in range(0,N+1):
            if x[i] < -K or x[i] > K:
                O[i] = 0
            else: 
                O[i] = sqrt*np.sin(j*frac * x[i])
    else:
        for i in range(0,N+1):
            if x[i] < -K or x[i] > K:
                O[i] = 0
            else: 
                O[i] = sqrt*np.cos(j*frac * x[i])
    return O


def error(j):
    for i in range (0,N+1):
        if wave(j)[i] == 0:
            O[i] = 0
        else:
            O[i] = 100*((wave(j)[i] - psi(j)[i]) / wave(j)[i])
    for k in range (0,N+1):
        if O[k] < 0:
            O[k] = -O[k]
    return O

#################### Plots ####################

# Plot for potential
plt.plot(x,V)
plt.xlabel('x')
plt.ylabel('V(x)')
plt.title('Plot of potential')
plt.show()

# Plot of psi-1
plt.plot(x,psi(1),color='b',label="calculated")
plt.plot(x,wave(1),color='r',label='theoretical',linestyle='--')
plt.xlabel('x')
plt.ylabel('\N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}(x)')
plt.title('Plot of \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}(x), N = '+str(N))
plt.legend()
plt.show()

# Plot of error for psi-1
plt.plot(x,error(1),color='b')
plt.xlabel('x')
plt.ylabel('% error')
plt.xlim(-K,K)
plt.ylim(0,1)
plt.title('Plot of error for \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}(x), N = '+str(N))
plt.show()

# Plot of psi-2
plt.plot(x,psi(2),color='b',label="calculated")
plt.plot(x,wave(2),color='r',label='theoretical',linestyle='--')
plt.xlabel('x')
plt.ylabel('\N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT TWO}(x)')
plt.title('Plot of \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT TWO}(x), N = '+str(N))
plt.legend()
plt.show()

# Plot of error for psi-2
plt.plot(x,error(2),color='b')
plt.xlabel('x')
plt.ylabel('% error')
plt.xlim(-K,K)
plt.ylim(0,1)
plt.title('Plot of error for \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT TWO}(x), N = '+str(N))
plt.show()

# Plot of psi-10
plt.plot(x,psi(10),color='b',label="calculated")
plt.plot(x,wave(10),color='r',label='theoretical',linestyle='--')
plt.xlabel('x')
plt.ylabel('\N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}\N{SUBSCRIPT ZERO}(x)')
plt.title('Plot of \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}\N{SUBSCRIPT ZERO}(x), N = '+str(N))
plt.legend()
plt.show()

# Plot of error for psi-10
plt.plot(x,error(10),color='b')
plt.xlabel('x')
plt.ylabel('% error')
plt.xlim(-K,K)
plt.ylim(0,1)
plt.title('Plot of error for \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}\N{SUBSCRIPT ZERO}(x), N = '+str(N))
plt.show()

# Plot of psi-100
plt.plot(x,psi(100),color='b',label="calculated")
plt.plot(x,wave(100),color='r',label='theoretical',linestyle='--')
plt.xlabel('x')
plt.ylabel('\N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}\N{SUBSCRIPT ZERO}\N{SUBSCRIPT ONE}(x)')
plt.title('Plot of \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}\N{SUBSCRIPT ZERO}\N{SUBSCRIPT ZERO}(x), N = '+str(N))
plt.legend()
plt.show()

# Plot of error for psi-100
plt.plot(x,error(100),color='b')
plt.xlabel('x')
plt.ylabel('% error')
plt.xlim(-K,K)
plt.ylim(0,1)
plt.title('Plot of error for \N{GREEK CAPITAL LETTER PSI}\N{SUBSCRIPT ONE}\N{SUBSCRIPT ZERO}\N{SUBSCRIPT ZERO}(x), N = '+str(N))
plt.show()