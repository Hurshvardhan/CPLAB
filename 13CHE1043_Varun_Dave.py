# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import scipy
import pandas as pd
headers = ['t_250', 'C_A_250', 'C_D_250', 't_300', 'C_A_300', 'C_D_300', 't_350', 'C_A_350', 'C_D_350', 't_400', 'C_A_400', 'C_D_400']
df = pd.read_csv("ExamProblemData.csv",  
                 header = 0,             
                 names = headers         
                 )  
import matplotlib.pyplot as plt
fig = plt.figure(); ax = fig.add_subplot(111)
ax.scatter(df.t_250, df.C_A_250, color="red")
ax.scatter(df.t_250, df.C_D_250, color = "blue")
ax.xaxis.label.set_text("time in seconds")
ax.yaxis.label.set_text("Concentration in kmol/m3")
fig.canvas.draw()

def reaction_model(N, t, kinetic_constants, R, rho):
    [NA, NB, NC, ND, V] = N  

    [k1, k2, k3] = kinetic_constants
    
    CA, CB, CC, CD = NA/V, NB/V, NC/V, ND/V
    
    r1 = k1*CA*CB
    r2 = k2*CA*CC
    r3 = k3*CB*CB
    
    dNAbydt = V*(-r1 - r2)
    dNBbydt = V*(-r1 - r3) + R
    dNCbydt = V*( r1 - r2)
    dNDbydt = V*( r2)
    dVbydt  = R/rho
    
    return [dNAbydt, dNBbydt, dNCbydt, dNDbydt, dVbydt]

import scipy.integrate

class ExperimentalRun:
    def __init__(self, df, T):
        time = 't_'+str(T)
        C_A = 'C_A_'+str(T)
        C_D = 'C_D_'+str(T)  
        self.t = df[time]
        self.CA = df[C_A]
        self.CD = df[C_D]

        ## Data ##
        self.NA0 = 100.0 #kmol
        self.R = 1.0 #kmol/s
        self.rho = 50.0 #kmol/s
        self.V0 = self.NA0/self.rho
        
        ## Guess ##
        self.k = [0.05, 0.005, 0.0005]   
        
    def solve(self, t = -1, N0 = -1):
        if t == -1:
            t = self.t
            N0 = [self.NA0, 0, 0, 0, self.V0]
        
            
        self.N = scipy.integrate.odeint(reaction_model, N0, t, args=(self.k, self.R, self.rho))
        self.CA_calc = self.N[:,0]/self.N[:,-1]
        self.CB_calc = self.N[:,1]/self.N[:,-1]
        self.CC_calc = self.N[:,2]/self.N[:,-1]
        self.CD_calc = self.N[:,3]/self.N[:,-1]
    
    def plot(self):
        fig = plt.figure(); ax = fig.add_subplot(111)
        ax.scatter(self.t, self.CA, color="red", label = "C_A experimental")
        ax.plot(self.t, self.CA_calc, color = "red", label = "C_A calculated")
        
        ax.scatter(self.t, self.CD, color="blue", label = "C_D experimental")
        ax.plot(self.t, self.CD_calc, color = "blue", label = "C_D calculated")
        
        ax.xaxis.label.set_text("time in seconds")
        ax.yaxis.label.set_text("concentration in kmol/m3")
        ax.legend()
        fig.canvas.draw()
        self.ax = ax
        
import scipy.optimize

def error_exp(kinetic_constants, exprun):
    exprun.k = kinetic_constants
    exprun.solve()
    errA = exprun.CA_calc - exprun.CA
    errD = exprun.CD_calc - exprun.CD
    err = scipy.concatenate((errA, errD))
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
        CA, CB, CC, CD = self.CA_calc, self.CB_calc, self.CC_calc, self.CD_calc
        listdCA, listdCB, listdCC, listdCD = [], [], [], []
        for i in xrange(len(self.k)):
            k = self.k[i]
            dk = abs(k)/1e6 + 1e-20
            self.k[i] = self.k[i] + dk
            self.solve()
            CAi, CBi, CCi, CDi = self.CA_calc, self.CB_calc, self.CC_calc, self.CD_calc
            dCA = (CAi - CA)/dk
            dCB = (CBi - CB)/dk
            dCC = (CCi - CC)/dk
            dCD = (CDi - CD)/dk
            listdCA.append(dCA)
            listdCB.append(dCB)
            listdCC.append(dCC)
            listdCD.append(dCD)
            self.k[i] = self.k[i] - dk
        errA = get_errY(listdCA, self.kcov)
        errB = get_errY(listdCB, self.kcov)
        errC = get_errY(listdCC, self.kcov)
        errD = get_errY(listdCD, self.kcov)
        self.solve()
        return 1.96*errA, 1.96*errB, 1.96*errC, 1.96*errD
    
    def plot_error(self):
        self.plot()
        ax = self.ax
        errA, errB, errC, errD = self.get_confidence_intervals()
        ax.fill_between(self.t, self.CA_calc-errA, self.CA_calc+errA, color="#ff0000",alpha = 0.2)
        ax.fill_between(self.t, self.CD_calc-errD, self.CD_calc+errD, color="#0000ff",alpha = 0.2)
        ax.figure.canvas.draw()
        print "k1 = %.2e (%.2e), k2 = %.2e (%.2e), k3 = %.2e (%.2e)"%(self.k[0],self.kerr[0]*1.96, self.k[1],self.kerr[1]*1.96, self.k[2],self.kerr[2]*1.96)
        
                
def get_errY(listdY, pcov):      
    listdY = scipy.array(listdY)
    left = scipy.dot(listdY.T, pcov)
    right = scipy.dot(left, listdY)
    sigma2Y = right.diagonal()
    sigmaY = scipy.sqrt(sigma2Y)
    errY = 1.96*sigmaY
    return errY     

def error_exp(kinetic_constants, exprun):
    exprun.k = [kinetic_constants[0], kinetic_constants[1], 0.0]  #k3 is always zero.
    exprun.solve()
    errA = exprun.CA_calc - exprun.CA
    errD = exprun.CD_calc - exprun.CD
    err = scipy.concatenate((errA, errD))
    return err

class ExperimentalRunNewModel(ExperimentalRunFit):
    def __init__(self, df, T):
        ExperimentalRunFit.__init__(self,df, T)
        self.Temp = T
    def fit(self):      
        k = [self.k[0], self.k[1]]
        (kopt, kcov, infodict, errmsg, ier) = scipy.optimize.leastsq(error_exp, k, args = (self,), full_output = 1)
        self.k = [kopt[0], kopt[1], 0.0]
        self.kcov = kcov
        self.kerr = scipy.sqrt(scipy.diag(kcov))   
        
    def get_confidence_intervals(self):
        self.solve()
        CA, CB, CC, CD = self.CA_calc, self.CB_calc, self.CC_calc, self.CD_calc
        listdCA, listdCB, listdCC, listdCD = [], [], [], []
        for i in [0,1]:
            k = self.k[i]
            dk = abs(k)/1e6 + 1e-20
            self.k[i] = self.k[i] + dk
            self.solve()
            CAi, CBi, CCi, CDi = self.CA_calc, self.CB_calc, self.CC_calc, self.CD_calc
            dCA = (CAi - CA)/dk
            dCB = (CBi - CB)/dk
            dCC = (CCi - CC)/dk
            dCD = (CDi - CD)/dk
            listdCA.append(dCA)
            listdCB.append(dCB)
            listdCC.append(dCC)
            listdCD.append(dCD)
            self.k[i] = self.k[i] - dk
        errA = get_errY(listdCA, self.kcov)
        errB = get_errY(listdCB, self.kcov)
        errC = get_errY(listdCC, self.kcov)
        errD = get_errY(listdCD, self.kcov)
        self.solve()
        return 1.96*errA, 1.96*errB, 1.96*errC, 1.96*errD
    def plot_error(self):
        self.plot()
        ax = self.ax
        errA, errB, errC, errD = self.get_confidence_intervals()
        ax.fill_between(self.t, self.CA_calc-errA, self.CA_calc+errA, color="#ff0000",alpha = 0.2)
        ax.fill_between(self.t, self.CD_calc-errD, self.CD_calc+errD, color="#0000ff",alpha = 0.2)
        ax.title.set_text("Plots for T = %s K"%(self.Temp))
        ax.figure.canvas.draw()
        print "k1 = %.2e (%.2e), k2 = %.2e (%.2e)"%(self.k[0],self.kerr[0]*1.96, self.k[1],self.kerr[1]*1.96)
    
exp400 = ExperimentalRunNewModel(df, '400')
exp400.k = [0.05, 0.0006, 0.000]
exp400.fit()
exp400.plot_error()