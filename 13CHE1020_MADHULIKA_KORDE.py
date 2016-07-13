# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 20:25:28 2016

@author: Madhu
"""

from scipy import array
import scipy,scipy.optimize
import matplotlib.pyplot as plt
import numpy as np
class ExperimentalRun:
    def __init__(self, df, T):
        err1=array([5.63*10**-2,0.55*10**-1,0.72*10**-1,0.81*10**-1])
        err2=array([1.83*10**-04,4.20*10**-04,1.43*10**-03,6.07*10**-03])
        get_errk1 = err1+str(T)
        get_errk2 = err2+str(T)
def error_exp(kinetic_constants, exprun):
    exprun.k = kinetic_constants
    exprun.solve()
    get_errk1=array([5.63*10**-2,0.55*10**-1,0.72*10**-1,0.81*10**-1])
    get_errk2=array([1.83*10**-04,4.20*10**-04,1.43*10**-03,6.07*10**-03])
    err = scipy.concatenate((get_errk1, get_errk2))
    return err
        
    class ExperimentalRunFit(ExperimentalRun):
        def __init__(self, df, T):
            ExperimentalRun.__init__(self, df, T) 
        def fit(self):         
            (kopt, kcov, infodict, errmsg, ier) = scipy.optimize.leastsq(error_exp, self.k, args = (self,), full_output = 1)
            self.k = kopt
            self.kcov = kcov
            self.kerr = scipy.sqrt(scipy.diag(kcov))

        def get_confidence_intervals(self):
            self.solve()
            listdk = []
            for i in xrange(len(self.k)):
                k = self.k[i]
                dk = abs(k)/1e6 + 1e-20
                self.k[i] = self.k[i] + dk
                self.solve()
                listdk.append(k)
                self.k[i] = self.k[i] - dk
            err1 = get_errk1(listdk, self.kcov)
            err2 = get_errk2(listdk, self.kcov)
            self.solve()
            return 1.96*err1, 1.96*err2
R=8.314
T=array([250,300,350,400])
#considering the values of k1, k2,err1 and err2 from the given ipython program
#with 95% confidence intervals 
k1=array([4.1*10**-2,8.03*10**-2,1.36*10**-1,1.9*10**-1])
err1=array([5.63*10**-2,0.55*10**-1,0.72*10**-1,0.81*10**-1])
k2=array([7.17*10**-4,1.85*10**-03,4.29*10**-03,1.13*10**-02])
err2=array([1.83*10**-04,4.20*10**-04,1.43*10**-03,6.07*10**-03])
print k1
def curve(T,A,Ea):
   return A*np.exp(-Ea/(R*T))
    
def calc_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2
pguess=[0.05,400]
res=scipy.optimize.curve_fit(curve,T,k1,[0.05,400],err1)
p=res[0]
[a,b]=p
ycalc=curve(T,p[0],p[1])
r2=calc_r2(T,k1,ycalc)
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(T,k1,'*')
ax.plot(T,ycalc,'r')
ax.fill_between(T, k1-(err1*k1), k1+(err1*k1), color="#ff0000",alpha = 0.2)    
ax.title.set_text('R2=%f'%(r2))
fig.canvas.draw()
plt.show()