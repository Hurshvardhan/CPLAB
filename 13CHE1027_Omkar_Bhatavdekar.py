# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:19:01 2016

@author: Omkar Bhatavdekar
"""

import scipy
import pandas as pd

headers = ['t_250', 'C_A_250', 'C_D_250', 't_300', 'C_A_300', 'C_D_300', 't_350', 'C_A_350', 'C_D_350', 't_400', 'C_A_400', 'C_D_400']
df = pd.read_csv("ExamProblemData.csv",  #The name of the csv file
                 header = 0,             #We are saying that the first line of the data is the headers
                 names = headers         #We are supplying our own headers (not using the names in the file)
                 )    

import matplotlib.pyplot as plt
#matplotlib inline  
fig = plt.figure(); ax = fig.add_subplot(111)
ax.scatter(df.t_250, df.C_A_250, color="red")
ax.scatter(df.t_250, df.C_D_250, color = "blue")
ax.xaxis.label.set_text("time in seconds")
ax.yaxis.label.set_text("Concentration in kmol/m3")
fig.canvas.draw()

def reaction_model(N, t, kinetic_constants, R, rho):
    [NA, NB, NC, ND, V] = N  #note the presence of V in the list

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
        C_D = 'C_D_'+str(T)  #See the trick here?
        self.t = df[time]
        self.CA = df[C_A]
        self.CD = df[C_D]

        ## Data ##
        self.NA0 = 100.0 #kmol
        self.R = 1.0 #kmol/s
        self.rho = 50.0 #kmol/s
        self.V0 = self.NA0/self.rho
        
        ## Guess ##
        self.k = [0.05, 0.005, 0.0005]  #Just some random guesses (Not really!  After some trial and error :) ).  
        
    def solve(self, t = 0, N0 = -1):
        if t == 0:
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

exp250 = ExperimentalRun(df, '250')
exp250.solve()
exp250.plot()   # Just marvel at the convenience!

import scipy.optimize  #Our favourite fitting library

def error_exp(kinetic_constants, exprun):
    exprun.k = kinetic_constants
    exprun.solve()
    errA = exprun.CA_calc - exprun.CA
    errD = exprun.CD_calc - exprun.CD
    err = scipy.concatenate((errA, errD))
    return err
    
error_exp(exp250.k, exp250)

class ExperimentalRunFit(ExperimentalRun):  #Heads up!  The class declaration just took an argument!
    def __init__(self, df, T):
        ExperimentalRun.__init__(self, df, T) #We now have everything that the class ExperimenalRun would have. See below.
    def fit(self): #We are adding a method attribute that does not exist in ExperimentalRun class.         
        (kopt, kcov, infodict, errmsg, ier) = scipy.optimize.leastsq(error_exp, self.k, args = (self,), full_output = 1)
        self.k = kopt
        self.kcov = kcov
        self.kerr = scipy.sqrt(scipy.diag(kcov))
        print errmsg

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

exp250 = ExperimentalRunFit(df, '250')
exp250.k = [0.05, 0.0006, 0.000]
exp250.fit()
exp250.plot_error()

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
    def fit(self): #We are modifying a method attribute of ExperimentalRunFit class    
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

exp250 = ExperimentalRunNewModel(df, '250')
exp250.k = [0.05, 0.0006, 0.000]
exp250.fit()
exp250.plot_error()

exp300 = ExperimentalRunNewModel(df, '300')
exp300.k = [0.05, 0.0006, 0.000]
exp300.fit()
exp300.plot_error()

exp350 = ExperimentalRunNewModel(df, '350')
exp350.k = [0.05, 0.0006, 0.000]
exp350.fit()
exp350.plot_error()

exp400 = ExperimentalRunNewModel(df, '400')
exp400.k = [0.05, 0.0006, 0.000]
exp400.fit()
exp400.plot_error()



# Homework


    
class ExperimentalRun_k:
    def __init__(self):

        ## Data ##
        self.R = 8.314 #Pa-m3/K
        
        ## Guess ##
        self.c = [0.1, 7000.0]  #Just some random guesses (Not really!  After some trial and error :) ).  
        self.A0=0.05
        self.Ea0= 8000.0
        self.k1_T=scipy.array([exp250.k[0], exp300.k[0], exp350.k[0], exp400.k[0]])
        self.T= [250,300,350,400]       

    def Arrhenius_model(self):

        [self.A, self.Ea] = self.c
        
        self.k1_250_calc= self.A*scipy.exp(-self.Ea/self.R/250.0)
        self.k1_300_calc= self.A*scipy.exp(-self.Ea/self.R/300.0)
        self.k1_350_calc= self.A*scipy.exp(-self.Ea/self.R/350.0)
        self.k1_400_calc= self.A*scipy.exp(-self.Ea/self.R/400.0)
        
        self.k1_T_calc=[self.k1_250_calc, self.k1_300_calc, self.k1_350_calc, self.k1_400_calc]
    
        
    

exp_k1=ExperimentalRun_k()
exp_k1.Arrhenius_model()


def error_exp_k1(constants, exprun):
    exprun.c = constants
    exprun.Arrhenius_model()
    errk1 = exprun.k1_T_calc - exprun.k1_T
    return errk1

error_exp_k1(exp_k1.c, exp_k1)

class ExperimentalRunFit_k(ExperimentalRun_k):  #Heads up!  The class declaration just took an argument!
    def __init__(self):
        ExperimentalRun_k.__init__(self) #We now have everything that the class ExperimenalRun would have. See below.
    def fit(self): #We are adding a method attribute that does not exist in ExperimentalRun class.         
        (copt, ccov, infodict, errmsg, ier) = scipy.optimize.leastsq(error_exp_k1, self.c, args = (self,), full_output = 1)
        self.c = copt
        self.ccov = ccov
        self.cerr = scipy.sqrt(scipy.diag(ccov))

    def get_confidence_intervals(self):
        self.Arrhenius_model()
        k1_250, k1_300, k1_350, k1_400 = self.k1_250_calc, self.k1_300_calc, self.k1_350_calc, self.k1_400_calc
        listdk1_T = []
        for i in xrange(len(self.c)):
            c = self.c[i]
            dc = abs(c)/1e12 + 1e-30
            self.c[i] = self.c[i] + dc
            self.Arrhenius_model()
            k1_250i, k1_300i, k1_350i, k1_400i = self.k1_250_calc, self.k1_300_calc, self.k1_350_calc, self.k1_400_calc
            dk1_250 = (k1_250i - k1_250)/dc
            dk1_300 = (k1_300i - k1_300)/dc
            dk1_350 = (k1_350i - k1_350)/dc
            dk1_400 = (k1_400i - k1_400)/dc
            dk1_T=scipy.array([dk1_250,dk1_300,dk1_350,dk1_400])
            listdk1_T.append(dk1_T)
            self.c[i] = self.c[i] - dc
            
                        
        err_k1_T = get_errY(listdk1_T, self.ccov)
        
        self.Arrhenius_model()
        return 1.96*err_k1_T
        
    def plot(self):
        fig = plt.figure(); ax = fig.add_subplot(111)
        
        self.T1=scipy.arange(250.0,450.0,10.0)
        self.k1_T_calc1= self.c[0]*scipy.exp(-self.c[1]/self.R/self.T1)
                
        
        ax.scatter(self.T, self.k1_T, color="red", label = "k1_A optimised")
        ax.plot(self.T1, self.k1_T_calc1, color = "red", label = "k1_A calculated")
                
        ax.xaxis.label.set_text("T in degree celsius")
        ax.yaxis.label.set_text("k1 in Pa.m3/degree celsius")
        
        ax.legend()
        fig.canvas.draw()
        self.ax = ax
        
        

exp_k1 = ExperimentalRunFit_k()
exp_k1.c = [ 0.1, 7000.0]
exp_k1.fit()
exp_k1.plot()

print 'A=', exp_k1.c[0] ,'Pa-m3/K'
print 'Ea=', exp_k1.c[1] ,'J/mol'

print 'error_k1(confidence interval)', exp_k1.get_confidence_intervals()