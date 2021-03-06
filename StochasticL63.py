

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 11:13:14 2018

@author: ms5717
"""

import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

#dW is the function that produces the Brownian motion
def dW(dt):
    return np.random.normal(loc = 0.0, scale = np.sqrt(dt))

s=10
b=8/3
r=28
sigma=0.3




t_i = 0
t_f = 30
N = 10000
dt = float(t_f - t_i)/N

ts = np.arange(t_i, t_f, dt)
#x,y,z  will save the solution
x   = np.zeros(N)
y   = np.zeros(N)
z   = np.zeros(N)

#Initial condition. Any other suggestion?
x[1]=0
y[1]=1
z[1]=30

for i in range(1,N-1):
    dww=dW(dt) #At each time we compute a shift of the Brownian motion
    x[i+1]=x[i] + s*(y[i]-x[i])*dt + sigma*x[i]*dww
    y[i+1]=y[i] + (r*x[i]-y[i]-x[i]*z[i])*dt + sigma*y[i]*dww
    z[i+1]=z[i] + (-b*z[i] + x[i]*y[i])*dt + sigma*z[i]*dww
    




#plt.hist2d(y, z, bins=50, norm=LogNorm())
#Histogram_xz=np.histogram2d(x,z,bins=50,normed=True)
    
    ###########
    

    
eig=np.linalg.eigvals(a)
plt.plot(eig)