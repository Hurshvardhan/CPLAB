# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 17:05:40 2016

@author: Sajeev
"""

import scipy,scipy.optimize
import matplotlib.pyplot as plt
import numpy as np
from scipy import array

R=8.314
T=array([250,300,350,400])
k1=array([4.1*10**-2,8.03*10**-2,1.36*10**-1,1.9*10**-1])
k1err=array([5.63*10**-2,1.55*10**-1,4.42*10**-1,6.51*10**-1])
k2=array([7.17*10**-4,1.85*10**-03,4.29*10**-03,1.13*10**-02])
k2err=array([1.83*10**-04,4.20*10**-04,1.43*10**-03,6.07*10**-03])

def curve(x,A,B):
    return A*np.exp(-B/(R*T))
    
def get_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2

pguess=[0.04,400]
res=scipy.optimize.curve_fit(curve,T,k1,[0.04,400],k1err)

''' for k1'''
P=res[0]
print P
pcov=res[1]
print pcov
print pcov.diagonal()**.5
[A,B]=P
resid = (k1-curve(T,A,B))/k1err
chisq=scipy.sum(resid**2)
chisqred = chisq/(len(resid)-len(P))
print chisqred
[A,B]=P
ycalc=curve(T,P[0],P[1])
r2=get_r2(T,k1,ycalc)
fig=plt.figure();
ax=fig.add_subplot(111)
ax.scatter(T,k1,color="red", label = "k1 experimental")
ax.plot(T,ycalc, color = "red", label = "k1 calculated")
ax.title.set_text('k1 vs T graph (R2=%f)'%(r2))
ax.xaxis.label.set_text("Temperature in degree celcius")
ax.yaxis.label.set_text("kinetic constant k1")
ax.legend()
fig.canvas.draw()
plt.show()

''' for k2'''
pguess=[0.04,400]
res=scipy.optimize.curve_fit(curve,T,k2,[0.04,400],k2err)

P=res[0]
print P
pcov=res[1]
print pcov
print pcov.diagonal()**.5
[A,B]=P
resid = (k2-curve(T,A,B))/k2err
chisq=scipy.sum(resid**2)
chisqred = chisq/(len(resid)-len(P))
print chisqred
[A,B]=P
ycalc=curve(T,P[0],P[1])
r2=get_r2(T,k2,ycalc)
fig=plt.figure();
ax=fig.add_subplot(111)
ax.scatter(T,k2,color="blue", label = "k2 experimental")
ax.plot(T,ycalc, color = "blue", label = "k2 calculated")  
ax.title.set_text('k2 vs T graph (R2=%f)'%(r2))
ax.xaxis.label.set_text("Temperature in degree celcius")
ax.yaxis.label.set_text("kinetic constant k2")
ax.legend()
fig.canvas.draw()
plt.show()