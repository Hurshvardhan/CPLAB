from scipy import array
import scipy,scipy.optimize
import matplotlib.pyplot as plt
import numpy as np
R=8.314
T=array([250,300,350,400])
k1=array([4.1*10**-2,8.03*10**-2,1.36*10**-1,1.9*10**-1])
err1=array([5.63*10**-2,1.55*10**-1,4.42*10**-1,6.51*10**-1])
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
pguess=[0.04,400]
res=scipy.optimize.curve_fit(curve,T,k1,[0.04,400],err1)
print res
p=res[0]
print p
pcov=res[1]
print pcov
print pcov.diagonal()**.5
[a,b]=p
residual = (k1-curve(T,a,b))/err1
chisq=scipy.sum(residual**2)
chisqred = chisq/(len(residual)-len(p))
print chisqred
[a,b]=p
ycalc=curve(T,p[0],p[1])
r2=calc_r2(T,k1,ycalc)
fig=plt.figure();
ax=fig.add_subplot(111)
ax.plot(T,k1,'*')
ax.plot(T,ycalc,'b')    
ax.title.set_text('R2=%f'%(r2))
fig.canvas.draw()
plt.show()