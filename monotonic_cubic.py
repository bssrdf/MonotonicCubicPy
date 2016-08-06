# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 23:40:17 2016

@author: bzhao
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def cubic_mono(xi, yi, xnew, mono=False):

   ynew=np.zeros(xnew.shape)
   for n in range(0, xnew.shape[0]):
        i = np.searchsorted(xi, xnew[n])-1
        t = (xnew[n]-xi[i])/(xi[i+1]-xi[i])
        a0=yi[i]
        #print n, xnew[n], i, i+1, i+2
        if i-1 < 0:
                s1=(yi[i+1]-yi[i])/(xi[i+1]-xi[i])
        else:
                s1=0.5*((yi[i+1]-yi[i])/(xi[i+1]-xi[i])+ (yi[i]-yi[i-1])/(xi[i]-xi[i-1]))
        if i == xi.shape[0]-2:
                s2 = (yi[i+1]-yi[i])/(xi[i+1]-xi[i])
        else:
                s2=0.5*((yi[i+2]-yi[i+1])/(xi[i+2]-xi[i+1])+(yi[i+1]-yi[i])/(xi[i+1]-xi[i]))
        if mono:
            dk = (yi[i+1]-yi[i])/(xi[i+1]-xi[i])
            if dk*s1 <= 0.0:
                s1 = 0.0
            if dk*s2 <= 0.0:
                s2 = 0.0
            alpha = s1 / dk
            beta  = s2 / dk
            if alpha > 3.0:
                s1 = 3.0 * dk
            if beta > 3.0:
                s2 = 3.0 * dk                
        t2 = t*t
        t3 = t2*t
        h00 = 2.0*t3-3.0*t2+1.0
        h10 = t3-2.0*t2+t
        h01 = -2.0*t3+3.0*t2
        h11 = t3 - t2
        ynew[n]= h00*yi[i] + h10*(xi[i+1]-xi[i])*s1 + h01*yi[i+1] + h11*(xi[i+1]-xi[i])*s2
        #print n, t, h00, h01, h10, h11, ynew[n]
        print n, xnew[n], ynew[n], s1, s2
   return ynew


xi=np.array([7.99, 8.09, 8.19, 8.7, 9.2, 10.0, 11.0, 12, 15, 20])
yi=np.array([0.0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.8497, 0.998636, 
             0.999919, 0.999994])

f1 = interpolate.interp1d(xi, yi, kind='cubic')
f2 = interpolate.interp1d(xi, yi)
f3 = interpolate.PchipInterpolator(xi, yi)
xnew=np.arange(8, 20, 0.1)
y1new=f1(xnew)
y2new=f2(xnew)
y3new=f3(xnew)

y4new= cubic_mono(xi, yi, xnew)
y5new= cubic_mono(xi, yi, xnew, mono=True)

plt.figure()

plt.plot(xi, yi,'ko')
#plt.plot(xnew, y1new)
plt.plot(xnew, y2new,'b')
plt.plot(xnew, y3new,'g')
#plt.plot(xnew, y4new,'r')
plt.plot(xnew, y5new,'m')
plt.ylim(0, 1.1)
plt.xlim(7.9, 20)
plt.show()
