# -*- coding: utf-8 -*-

import scipy
import pandas as pd

headers = ['t_250', 'C_A_250', 'C_D_250', 't_300', 'C_A_300', 'C_D_300', 't_350', 'C_A_350', 'C_D_350', 't_400', 'C_A_400', 'C_D_400']
df = pd.read_csv("C:\Users\Nandini Bhosale\Desktop\Python\ExamProblemData.csv",  #The name of the csv file
                 header = 0,             #We are saying that the first line of the data is the headers
                 names = headers         #We are supplying our own headers (not using the names in the file)
                 )

import matplotlib.pyplot as plt 



def reaction_model(N, t, kinetic_constants, R, rho):
    [NA, NB, NC, ND, V] = N  #note the presence of V in the list

    [k1, k2, k3] = kinetic_constants
    
    CA, CB, CC, CD = NA/V, NB/V, NC/V, ND/V
    
    r1 = k1*CA*CB
    r2 = k2*CA*CC
    r3 = k3*CB*CB
    
    dNAdt = V*(-r1 - r2)
    dNBdt = V*(-r1 - r3) + R
    dNCdt = V*( r1 - r2)
    dNDdt = V*( r2)
    dVdt  = R/rho
    
    return [dNAdt, dNBdt, dNCdt, dNDdt, dVdt]
    
import scipy.integrate

class ExperimentalRun:
    def __init__(self, df, T):
        time = 't_'+str(T)
        C_A = 'C_A_'+str(T)
        C_D = 'C_D_'+str(T)  #See the trick here? : Yes
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
    
class ExperimentalRunFit(ExperimentalRun):  #Heads up!  The class declaration just took an argument!
    def __init__(self, df, T):
        ExperimentalRun.__init__(self, df, T) #We now have everything that the class ExperimenalRun would have. See below.
    def fit(self): #We are adding a method attribute that does not exist in ExperimentalRun class.         
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
    def fit(self): #We are modifying a method attribute of ExperimentalRunFit class    
        k = [self.k[0], self.k[1]]
        (kopt, kcov, infodict, errmsg, ier) = scipy.optimize.leastsq(error_exp, k, args = (self,), full_output = 1)
        self.k = [kopt[0], kopt[1], 0.0]
        self.kcov = kcov
        self.kerr = scipy.sqrt(scipy.diag(kcov))
        return self.k[0],self.k[1]
        
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
        #print "k1 = %.2e (%.2e), k2 = %.2e (%.2e)"%(self.k[0],self.kerr[0]*1.96, self.k[1],self.kerr[1]*1.96)
        
k1s=[]        
k2s=[]
exp250 = ExperimentalRunNewModel(df, '250')
exp250.k = [0.05, 0.0006, 0.000]
exp250.fit()
exp250.plot_error()
#plt.show()

exp300 = ExperimentalRunNewModel(df, '300')
exp300.k = [0.05, 0.0006, 0.000]
exp300.fit()
exp300.plot_error()
#plt.show()

exp350 = ExperimentalRunNewModel(df, '350')
exp350.k = [0.05, 0.0006, 0.000]
exp350.fit()
exp350.plot_error()
#plt.show()

exp400 = ExperimentalRunNewModel(df, '400')
exp400.k = [0.05, 0.0006, 0.000]
exp400.fit()
exp400.plot_error()
#plt.show()

k1s.append((exp250.fit())[0])
k1s.append((exp300.fit())[0])
k1s.append((exp350.fit())[0])
k1s.append((exp400.fit())[0])

k2s.append((exp250.fit())[1])
k2s.append((exp300.fit())[1])
k2s.append((exp350.fit())[1])
k2s.append((exp400.fit())[1])

Ts=[250,300,350,400]

k1exp=k1s
k2exp=k2s
Texp=Ts

Tinv=[1.0/250.0,1.0/300.0,1.0/350.0,1.0/400.0]
from scipy import exp
from scipy import log
def eqn(parameters,kexp,Texp):
    [A,Ea]=parameters
    err=log(kexp)-log(A*exp(-Ea/8.314/Texp))
    return err

parameters=[0.1,1000]
(popt, pcov1, infodict, errmsg, ier) = scipy.optimize.leastsq(eqn, parameters, args = (k1exp,Texp), full_output = 1)
k1err=1.96*scipy.sqrt(scipy.diag(pcov1))

print "Value of A for 1st reaction =",popt[0], " with error", k1err[0]
print "Value of Ea for 1st reaction =",popt[1], " with error", k1err[1]
[A1,Ea1]=popt
(popt, pcov2, infodict, errmsg, ier) = scipy.optimize.leastsq(eqn, parameters, args = (k2exp,Texp), full_output = 1)
k2err=1.96*scipy.sqrt(scipy.diag(pcov2))
print "Value of A for 2nd reaction =",popt[0], " with error", k2err[0]
print "Value of Ea for 2nd reaction =",popt[1], " with error", k2err[1]
[A2,Ea2]=popt
k1cal=A1*exp(-Ea1/8.314/Texp)
k2cal=A2*exp(-Ea2/8.314/Texp)



"""
This did not work
def k(A,Ea):
    k=(A*exp(-Ea/8.314/Texp))
    return k
    
da1=A1/10**14
dEa1=Ea1/10**14
da2=A2/10**14
dEa2=Ea2/10**14


conf1=[]
conf1.append((k(A1+da1,Ea1)-k(A1,Ea1))/da1)
conf1.append((k(A1,Ea1+dEa1)-k(A1,Ea1))/dEa1)

conf2=[]
conf2.append((k(A2+da2,Ea2)-k(A2,Ea2))/da2)
conf2.append((k(A2,Ea2+dEa2)-k(A2,Ea2))/dEa2)

errk1 = 1.96*get_errY(conf1, pcov1)
errk2 = 1.96*get_errY(conf2, pcov2)
print errk1,errk2
"""
fig = plt.figure(); bx = fig.add_subplot(111)
bx.scatter(Tinv,log(k1exp), color="red", label ="Ln K1exp")
bx.plot(Tinv, log(A1*exp(-Ea1/8.314/Texp)), color = "red", label = "Ln K1cal")

bx.scatter(Tinv, log(k2exp), color="blue", label = "Ln K2exp")
bx.plot(Tinv, log(A2*exp(-Ea2/8.314/Texp)), color = "blue", label = "Ln K2cal")
        
bx.xaxis.label.set_text("1/Temperature")
bx.yaxis.label.set_text("Ln K")
bx.legend()
fig.canvas.draw()
plt.show()

