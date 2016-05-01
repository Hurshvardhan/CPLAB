# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 22:19:43 2016

@author: Aayushi
"""

from scipy.interpolate import spline
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import leastsq
import win32com.client
from scipy.optimize import fsolve
import statsmodels.stats.stattools as stools
"""This part of the code obtains data from excel"""
xl= win32com.client.gencache.EnsureDispatch('Excel.Application')
wb=xl.Workbooks('rates.xlsx')
sheet=wb.Sheets('Sheet1')
def assign(x, p):
    A,B = p
    return (A*x+B)
def residuals(p, y, x):
    A,B = p
    err = (((y-(A*x+B))**2)**0.5)
    return err   
def getdata(sheet, Range):
    data= sheet.Range(Range).Value
    data=scipy.array(data)
    data=data.reshape((1,len(data)))[0]
    return data
Cdata=getdata(sheet,"C4:C11")
tdata=getdata(sheet,"B4:B11")

"""This part fits a curve for the data, finds its derivative and 
makes a log log plot with the derivative to find k and n in rate
equation."""
x = tdata 
y = Cdata
# calculate polynomial 
z = np.polyfit(x, y, 10) 
f = np.poly1d(z) 
# calculate new x's and y's 
x_new = np.linspace(x[0], x[-1], 50) 
y_new = f(x_new) 
plt.plot(x,y,'o', x_new, y_new) 

plt.show()
t=0.0
x=getdata(sheet,"C4:C10")
x=scipy.log(x)
i=0
g=[]
while t<0.7:
    g.append(scipy.log(-(f(t+0.0001)-f(t))/0.0001))#taking derivative
    t=t+0.1
    i=i+1
print x
print g


guessAB = [1.0, 10.0]
fitting = leastsq(residuals, guessAB, args=(g, x))
popt=fitting[0]
A=popt[0]
B=popt[1]
print("n,k=")

print popt
print("k has volume units in l and time units in hrs as is in excel sheet")
greal = (A*x+B)
Q=1#l/hr
n1=input("Enter the number of reactors you want to insert:")
i=0
"""here after asking the user how many reactors he wants to insert, we can
create a cooresponding number of objects of our class, each giveing an outlet
concn."""
#Ca0=Cdata[0]
Ca0=1
c=[]
c.append(Ca0)
n=popt[0]
k=popt[1]
print("Initial Conc=",Ca0)
plt.title('Plot of C vs t from excel')
plt.plot(i,Ca0,'ro')
xnew = np.linspace(0,n1,300)
seti=[]
seti.append(0)
while(i<n1):
    
    seti.append((i+1))
    V1=input("Reactor Volume in l:")
    V=V1
    typ=input("ENTER +ve no FOR CSTR, -ve no FOR PFR:")
        
    
    y=1.0
    
    
    tow=V/Q
    
    """rate equation for cstr"""
    def cstr(conc):
        return((tow*k*(conc**n))+conc-Ca0)
        
        
    """rate equation for pfr"""
    def pfr(conc):
        if(n==1):
            b=(scipy.log(Ca0/conc))-(k*tow)
        else:
            b=-((Ca0**(1-n))/(1-n))+((conc**(1-n))/(1-n))-(k*tow)
        return b
        
    
    if typ>0:
        tow=V/Q
        l=fsolve(cstr,0.5)
        
            #solving for c with guess value 1
    else:
        tow=V/Q
        
        if(n==1):
            l=fsolve(pfr,0.5)
        else:
            l=((1-n)*(((Ca0**(1-n))/(1-n))-(k*tow)))**(1/(1-n))
        
    Ca0=l
    print("Conc=",l)
    c.append(l)
    
    i=i+1
xnew = np.linspace(0,n1,300)

power_smooth = spline(seti,c,xnew)
plt.title('Concn vs no of reactors')
plt.plot(xnew,power_smooth)
plt.show()        
        
        