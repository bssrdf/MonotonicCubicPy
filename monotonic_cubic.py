# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 23:40:17 2016

@author: bzhao
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

xi=np.array([7.99, 8.09, 8.19, 8.7, 9.2, 10.0, 12, 15, 20])
yi=np.array([0.0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 
             0.999919, 0.999994])

f1 = interpolate.interp1d(xi, yi, kind='cubic')
f2 = interpolate.interp1d(xi, yi)
f3 = interpolate.PchipInterpolator(xi, yi)
xnew=np.arange(8, 20, 0.1)
y1new=f1(xnew)
y2new=f2(xnew)
y3new=f3(xnew)

plt.figure(1)
plt.plot(xi, yi,'ko')
#plt.plot(xnew, y1new)
plt.plot(xnew, y2new)
plt.plot(xnew, y3new)
plt.ylim(0, 1.5)
plt.xlim(7.9, 20)
plt.show()