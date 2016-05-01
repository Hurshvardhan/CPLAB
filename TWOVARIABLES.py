# -*- coding: utf-8 -*-
"""
Created on Sun May 01 11:01:34 2016

@author: Hurshvardhai
"""


import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import statsmodels.stats.stattools as stools
import win32com.client
xl= win32com.client.gencache.EnsureDispatch('Excel.Application')
wb=xl.Workbooks('tpdata.xlsx')
sheet=wb.Sheets('Sheet1')
   
def getdata(sheet, Range):
    data= sheet.Range(Range).Value
    data=scipy.array(data)
    data=data.reshape((1,len(data)))[0]
    return data
Ydata=getdata(sheet,"B3:B6")
xdata=getdata(sheet,"C3:C6")
ydata=getdata(sheet,"D3:D6")
def assign(x,Y, p):
    A,B,C = p
    return (A*x+B*Y+C)
def residuals(p, y,Y, x):
    A,B,C = p
    err = (((y-(A*x+B*Y+C))**2)**0.5)
    return err

guessAB = [1,1,1]
fitting = leastsq(residuals, guessAB, args=(ydata,Ydata, xdata))
popt=fitting[0]
pcov=fitting[1]
print("*******************************")
print("Optimized A, B AND C in AX+BY+C")
print popt[0],popt[1],popt[2]

