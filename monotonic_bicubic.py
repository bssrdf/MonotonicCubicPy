# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 01:13:51 2016

@author: bzhao
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt


def basis_func(t):
    t2 = t*t
    t3 = t2*t
    h00 = 2.0*t3-3.0*t2+1.0
    h10 = t3-2.0*t2+t
    h01 = -2.0*t3+3.0*t2
    h11 = t3 - t2
    return h00,h10,h01,h11

def limit_tangent(dk1, dk2, s1, s2):
   if dk1*s1 <= 0.0:
       t1 = 0.0
   else:
       alpha = s1 / dk1           
       if alpha > 3.0:
           t1 = 3.0 * dk1
       else:
           t1 = s1
   if dk2*s2 <= 0.0:
       t2 = 0.0
   else:
       beta  = s2 / dk2    
       if beta > 3.0:
           t2 = 3.0 * dk2               
       else:
           t2 = s2
   return t1, t2
   
def limit_tangent_hyman(dk1, dk2, s):
   if dk1*dk2 <= 0.0:
       t = 0.0
   else:
       if s >= 0.0:
           if np.abs(dk1) < np.abs(dk2):
               t = 3.0*np.abs(dk1)
           else:
               t = 3.0*np.abs(dk2)
           if s < t:
               t = s
       else:
           if np.abs(dk1) < np.abs(dk2):
               t = -3.0*np.abs(dk1)
           else:
               t = -3.0*np.abs(dk2)
           if s > t:
               t = s
   return t   
   

def interp_1d(t, i, xi, yi, mono=False):
    if i-1 < 0:
        s1=(yi[i+1]-yi[i])/(xi[i+1]-xi[i])
    else:
        s1=0.5*((yi[i+1]-yi[i])/(xi[i+1]-xi[i])+ (yi[i]-yi[i-1])/(xi[i]-xi[i-1]))
    if i == xi.shape[0]-2:
        s2 = (yi[i+1]-yi[i])/(xi[i+1]-xi[i])
    else:
        s2=0.5*((yi[i+2]-yi[i+1])/(xi[i+2]-xi[i+1])+(yi[i+1]-yi[i])/(xi[i+1]-xi[i]))
    if mono:
       if i-1 < 0:
            dk1=(yi[i+1]-yi[i])/(xi[i+1]-xi[i])
       else:            
            dk1 = (yi[i]-yi[i-1])/(xi[i]-xi[i-1])
       dk2 = (yi[i+1]-yi[i])/(xi[i+1]-xi[i])
       t1=limit_tangent_hyman(dk1, dk2, s1)    
       dk1= (yi[i+1]-yi[i])/(xi[i+1]-xi[i])
       if i == xi.shape[0]-2:
           dk2 = (yi[i+1]-yi[i])/(xi[i+1]-xi[i])
       else:
           dk2 = (yi[i+2]-yi[i+1])/(xi[i+2]-xi[i+1])
       #t1,t2=limit_tangent(dk1, dk1, s1, s2)
       t2=limit_tangent_hyman(dk1, dk2, s2)    
    else:
       t1 = s1
       t2 = s2
    h00,h10,h01,h11=basis_func(t)
    ynew = h00*yi[i] + h10*(xi[i+1]-xi[i])*t1 + h01*yi[i+1] + h11*(xi[i+1]-xi[i])*t2
        #print n, t, h00, h01, h10, h11, ynew[n]
        #print n, xnew[n], ynew[n], s1, s2
    return ynew    

def bicubic_mono(xi, yi, zi, xnew, ynew, mono=False):
   assert len(xnew.shape) == 2
   assert len(ynew.shape) == 2
   assert xnew.shape==ynew.shape   
   znew=np.zeros(xnew.shape)
   for n in range(0, xnew.shape[1]):       
       for m in range(0, ynew.shape[0]):
           i = np.searchsorted(xi, xnew[m,n])-1
           j = np.searchsorted(yi, ynew[m,n])-1
           t = (xnew[m,n]-xi[i])/(xi[i+1]-xi[i])
           tmpj = interp_1d(t, i, xi, zi[j,:], mono)
           tmpjp1 = interp_1d(t, i, xi, zi[j+1,:], mono)
           if j-1 >= 0:
               tmpjm1 = interp_1d(t, i, xi, zi[j-1,:], mono)
           if j+2 < yi.shape[0]:
               tmpjp2 = interp_1d(t, i, xi, zi[j+2,:], mono)
           r = (ynew[m,n]-yi[j])/(yi[j+1]-yi[j])
           if j-1 < 0:           
               s1=(tmpjp1-tmpj)/(yi[j+1]-yi[j])
           else:           
               s1=0.5*((tmpjp1-tmpj)/(yi[j+1]-yi[j])+ (tmpj-tmpjm1)/(yi[j]-yi[j-1]))
           if j == yi.shape[0]-2:           
               s2 = (tmpjp1-tmpj)/(yi[j+1]-yi[j])
           else:           
               s2=0.5*((tmpjp2-tmpjp1)/(yi[j+2]-yi[j+1])+ (tmpjp1-tmpj)/(yi[j+1]-yi[j]))
           if mono:               
               if j-1 < 0:           
                   dk1 = (tmpjp1-tmpj)/(yi[j+1]-yi[j])
               else:
                   dk1 = (tmpj-tmpjm1)/(yi[j]-yi[j-1])
               dk2 = (tmpjp1-tmpj)/(yi[j+1]-yi[j])    
               t1=limit_tangent_hyman(dk1, dk2, s1)                   
               dk1 = (tmpjp1-tmpj)/(yi[j+1]-yi[j])
               if j == yi.shape[0]-2:           
                   dk2 = (tmpjp1-tmpj)/(yi[j+1]-yi[j])
               else:
                   dk2 = (tmpjp2-tmpjp1)/(yi[j+2]-yi[j+1])
               t2=limit_tangent_hyman(dk1, dk2, s2)                   
               #t1,t2=limit_tangent(dk1, dk1, s1, s2)               
           else:
               t1 = s1
               t2 = s2
           h00,h10,h01,h11=basis_func(r)           
           znew[m,n]= h00*tmpj + h10*(yi[j+1]-yi[j])*t1 + h01*tmpjp1 + h11*(yi[j+1]-yi[j])*t2                          
           #if abs(znew[m,n]-0.994465924832) < 1.e-6:
            #   print m, n, t, r, znew[m,n], h00, h10, h01, h11, tmpj, tmpjp1
             #  print m, n, t, r, xnew[m,n], xi[i], ynew[m,n], yi[j]
              # print m, n, j, i, zi[j+1,i], zi[j+1,i+1]
   return znew


xi=np.arange(0, 5.2, 0.2)
yi=np.arange(0.6, 5.2, 0.2)
xin=np.arange(0.1, 5.0, 0.05)
yin=np.arange(1.1, 5.0, 0.1)



x,y=np.meshgrid(xi, yi)
xn,yn=np.meshgrid(xin, yin)

zi = np.sin(x)*np.cos(y)

#f1=interpolate.RectBivariateSpline(xi, yi, zi)
f2=interpolate.interp2d(xi, yi, zi, kind='cubic')
#z1 = f1(xin, yin)
#z1 = z1.T
z2 = f2(xin, yin)

z3 = bicubic_mono(xi, yi, zi, xn, yn)
z4 = bicubic_mono(xi, yi, zi, xn, yn, mono=True)

print 'original min/max: ', zi.min(), '/', zi.max()
print 'scipy min/max: ', z2.min(), '/', z2.max()
print 'bicubic min/max: ', z3.min(), '/', z3.max()
print 'mono bicubic min/max: ', z4.min(), '/', z4.max()

plt.figure()
plt.subplot(2,2,1)
plt.pcolormesh(x, y, zi)
plt.title('original data', fontsize=15)
plt.subplot(2,2,2)
plt.pcolormesh(xn, yn, z2)
plt.title('scipy interp2d', fontsize=15)
plt.subplot(2,2,3)
plt.pcolormesh(xn, yn, z3)
plt.title('bicubic', fontsize=15)
plt.subplot(2,2,4)
plt.pcolormesh(xn, yn, z4)
plt.title('monotonic bicubic', fontsize=15)
plt.show()

plt.figure()
plt.subplot(2,2,1)
plt.pcolormesh(xn, yn, z2-z4)
plt.clim(-1.e-5, 1.e-5)
plt.subplot(2,2,2)
plt.pcolormesh(xn, yn, z3-z4)
plt.clim(-1.e-5, 1.e-5)
plt.show()

plt.figure()
plt.plot(z2[:,30],'k-')
plt.plot(z3[:,30],'g-')
plt.plot(z4[:,30],'r-')
plt.show()

plt.figure()
plt.plot(z2[:,30]-z4[:,30],'k-')
plt.plot(z3[:,30]-z4[:,30],'g-')
plt.show()



