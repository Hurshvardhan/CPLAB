

import scipy,scipy.optimize
import matplotlib.pyplot as plt
import numpy as np
from scipy import array
R=8.314
x=array([250,300,350,400])
k1=array([4.1*10**-2,8.03*10**-2,1.36*10**-1,1.9*10**-1])
k1err=array([5.63*10**-2,1.55*10**-1,4.42*10**-1,6.51*10**-1])
k2=array([7.17*10**-4,1.85*10**-03,4.29*10**-03,1.13*10**-02])
k2err=array([1.83*10**-04,4.20*10**-04,1.43*10**-03,6.07*10**-03])
print k1
def curve(x,A,B):
    #[A,B]=p
    
    return A*np.exp(-B/(R*x))
    
def get_r2(x,y,ycalc):
    ymean=scipy.average(y)
    dymean2=(y-ymean)**2
    dycalc2=(y-ycalc)**2
    r2=1-sum(dycalc2)/sum(dymean2)
    return r2
pguess=[0.04,400]
res=scipy.optimize.curve_fit(curve,x,k1,[0.04,400],k1err)
print res
P=res[0]
print P
pcov=res[1]
print pcov
print pcov.diagonal()**.5
[A,B]=P
resid = (k1-curve(x,A,B))/k1err
chisq=scipy.sum(resid**2)
chisqred = chisq/(len(resid)-len(P))
print chisqred
[A,B]=P
ycalc=curve(x,P[0],P[1])
r2=get_r2(x,k1,ycalc)
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(x,k1,'*')
ax.plot(x,ycalc,'b')    
ax.title.set_text('R2=%f'%(r2))
fig.canvas.draw()
plt.show()