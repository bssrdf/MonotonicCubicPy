# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 01:13:51 2016

@author: bzhao
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 23:40:17 2016

@author: bzhao
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def bicubic_mono(xi, yi, zi, xnew, ynew, mono=False):

   assert len(xnew.shape) == 2
   assert len(ynew.shape) == 2
      
   znew=np.zeros(xnew.shape)
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
   return znew


xi=np.arange(0, 5.2, 0.2)
yi=np.arange(0, 5.2, 0.2)
xin=np.arange(0.1, 5.0, 0.05)
yin=np.arange(0.1, 5.0, 0.1)



x,y=np.meshgrid(xi, yi)
xn,yn=np.meshgrid(xin, yin)

zi = np.sin(x)*np.cos(y)

f1=interpolate.RectBivariateSpline(xi, yi, zi)
f2=interpolate.interp2d(xi, yi, zi, kind='cubic')
z1 = f1(xin, yin)
#z1 = z1.T
z2 = f2(xin, yin)


plt.figure()
plt.subplot(2,2,1)
plt.pcolormesh(x, y, zi)
plt.subplot(2,2,2)
plt.pcolormesh(xn, yn, z1)
plt.subplot(2,2,3)
plt.pcolormesh(xn, yn, z2)
plt.show()
