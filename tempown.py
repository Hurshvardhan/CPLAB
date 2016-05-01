# -*- coding: utf-8 -*-
"""
Created on Sun May 01 18:40:09 2016

@author: Hurshvardhai
"""
import numpy as np
from scipy.integrate import quad
from scipy.integrate import odeint
import matplotlib.pyplot as plt
k=(300*3.142*0.1)/(1*42)
Th=30
Tc=40
L=1
n=10
def temp(T,x):
            
        dTh =k*(T[1]-T[0])
        dTc =k*(T[1]-T[0])
        return [dTh,dTc]
T0=[Th,Tc]
x=np.linspace(0,L,n)
tem=odeint(temp,T0,x)
print tem[(n-1),1]
while(tem[(n-1),1]<80):
    Tc=Tc+0.01*(tem[(n-1),1]-Tc)
    T0=[Th,Tc]
    x=np.linspace(0,L,n)
    tem=odeint(temp,T0,x)
    print tem[(n-1),1]

i=0
x=0
while(i<n):
    
    plt.plot(x,tem[i,0],'ro')
    plt.plot(x,tem[i,1],'bo')
    x=x+0.1
    i=i+1
plt.show()
