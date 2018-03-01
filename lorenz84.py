#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 17:39:00 2018

@author: ms5717
"""

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

# these are our constants

def Lorenz84(x,t):
  a=0.25
  b=4
  G=0.24
  F=8

  # compute state derivatives
  d = np.zeros(3)
  # first the 3 edge cases: i=1,2,N
  d[0] = -x[1]**2 - x[2]**2 - a*x[0] + a*F
  d[1] = x[0]*x[1] - b*x[0]*x[2]-x[1]+G
  d[2]=b*x[0]*x[1] + x[0]*x[2] - x[2]

  # return the state derivatives
  return d

x0 = np.ones(3) # initial state (equilibrium)

t = np.arange(0.0, 30.0, 0.01)

x = odeint(Lorenz84, x0, t)

# plot first three variables
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x[:,0],x[:,1],x[:,2])
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')
plt.show()