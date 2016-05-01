# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 21:06:07 2016

@author: admin
"""

import scipy, numpy
import scipy.optimize, scipy.stats
import numpy.random
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels
import statsmodels.stats
import statsmodels.stats.stattools as stools
import win32com.client

def fitdata(f,Xdata,Ydata,Errdata,pguess,ax=False,ax2=False):

    def error(p,Xdata,Ydata,Errdata):
        Y=f(Xdata,p)
        residuals=(Y-Ydata)/Errdata
        return residuals
    res=scipy.optimize.leastsq(error,pguess,args=(Xdata,Ydata,Errdata),full_output=1)    
    (popt,pcov,infodict,errmsg,ier)=res
    perr=scipy.sqrt(scipy.diag(pcov))
    M=len(Ydata)
    N=len(popt)
    Y=f(Xdata,popt)
    residuals=(Y-Ydata)/Errdata
    meanY=scipy.mean(Ydata)
    squares=(Y-meanY)/Errdata
    squaresT=(Ydata-meanY)/Errdata
    
    SSM=sum(squares**2)
    SSE=sum(residuals**2)
    SST=sum(squaresT**2)
    
    DFM=N-1
    DFE=M-N
    DFT=M-1
    
    MSM=SSM/DFM
    MSE=SSE/DFE
    MSM=SST/DFT
    
    '''R2'''
    R2=SSM/SST
    R2_adj=1-(1-R2)*(M-1)/(M-N-1)
    
    '''t-test'''
    t_stat=popt/perr
    t_stat=t_stat.real
    p_p=1.0-scipy.stats.t.cdf(t_stat,DFE)
    z=scipy.stats.t(M-N).ppf(0.95)
    p95=perr*z
    
    '''chi-square'''
    chisqred=sum(residuals**2)
    degfrdm=M-N
    chisqred_red=chisqred/degfrdm
    p_chi2=1.0-scipy.stats.chi2.cdf(chisqred,degfrdm)
    stderr_reg=scipy.sqrt(chisqred_red)
    chisqre=(p_chi2,chisqred,chisqred_red,degfrdm,R2,R2_adj)
    
    '''shapiro-wilk test'''
    w,p_shapiro=scipy.stats.shapiro(residuals)
    mean_res=scipy.mean(residuals)
    stddev_res=scipy.sqrt(scipy.var(residuals))
    t_res=mean_res/stddev_res
    p_res=1.0-scipy.stats.t.cdf(t_res,M-1)
    
    '''F-test'''
    F=MSM/MSE
    p_F=1-scipy.stats.f.cdf(F,DFM,DFE)
    
    '''durbin-watson'''
    dw=stools.durbin_watson(residuals)
    
    resanal=(p_shapiro,w,mean_res,p_res,F,p_F,dw)

    return popt,pcov,perr,p95,p_p,chisqre,resanal
    
  
def f(Xdata,p):
    Y=p[0]*Xdata+p[1]
    return Y
def error_fit(f,Xdata,popt,pcov):
    Y=f(Xdata,popt)
    dY=[]
    for i in xrange(len(popt)):
        p=popt[i]
        dp=abs(p)/1e6+1e-20
        popt[i]+=dp
        Yi=f(Xdata,popt)
        dy=(Yi-Y)/dp
        dY.append(dy)
        popt[i]-=dp
    dY=scipy.array(dY)  
    A=scipy.dot(dY.T,pcov)
    B=scipy.dot(A,dY)
    sigma2=B.diagonal()
    mean_sigma2=scipy.mean(sigma2)
    M=len(Xdata)
    N=len(popt)
    avg_stddev_data=scipy.sqrt(M*mean_sigma2/N)
    sigma=scipy.sqrt(sigma2)
    return sigma,avg_stddev_data
    
xl = win32com.client.gencache.EnsureDispatch("Excel.Application")
wb = xl.Workbooks('expt19.xlsx')
sheet = wb.Sheets('sheet1')
def getdata(sheet,Range):
    data = sheet.Range(Range).Value
    data = scipy.array(data)
    data = data.reshape((1, len(data)))[0]
    return data
Xdata =getdata(sheet,"F48:F58")
Ydata = getdata(sheet,"L48:L58")
Errdata = getdata(sheet,"M48:M58")
    
                        
N=2
pguess=N*[0.0]
    
popt,pcov,perr,p95,p_p,chisqre,resanal=fitdata(f,Xdata,Ydata,Errdata,pguess)
sigma,avg_stddev_data=error_fit(f,Xdata,popt,pcov)
X1=scipy.linspace(Xdata[0],Xdata[len(Ydata)-1],100)
Y1=popt[0]*X1+popt[1]
Y2=f(Xdata,popt)
plt.plot(Xdata,Ydata,'ro')
plt.plot(X1,Y1,'g')
print 'the parameters are ''\n',popt


