# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\admin\.spyder2\.temp.py
"""
import numpy as np
from scipy.integrate import quad
from scipy.integrate import odeint
import matplotlib.pyplot as plt
class absorption:
    CO2=1
    CH4=2
    H2O=3
    
    y1=.5
    y2=.5
    y3=0
    G0=100 #mol/hr
    G10=y1*G0 #mol/hr
    G20=y2*G0#mol/hr
    G30=y3*G0 #mol/hr
    
    x3=1
    x1=0
    x2=0
    LH=150 #mol/hr
    L1H=x1*LH #mol/hr
    L2H=x2*LH #mol/hr
    L3H=x3*LH #mol/hr
    
    H1= 16000*18*101325/1000 #atm
    H2= 37037*18*100000/1000 #bar
    Psat= .0031699*1000000 #MPa at 25degree celcius
    #Kl1=.163/3600 #m/hr
    #Kl2=.2/3600 #m/hr
    #kg3=1.22/3600 #m/hr
    Kl1=.01
    Kl2=.005
    kg3=2.452*(10**-6)
    H=3 #m
    A=.5 #m2
    a=1 #m-1
    
    P=101325 #atm
    T=25 #degreecelcius
    
    def mol(self):
        G0=self.G0   
        G10=self.G10
        G20=self.G20
        G30=self.G30
        y1=self.y1
        y2=self.y2
        y3=self.y3
        
        LH=self.LH
        L1H=self.L1H
        L2H=self.L2H
        L3H=self.L3H
        
        x1=self.x1
        x2=self.x2
        x3=self.x3
        H1=self.H1
        H2=self.H2
        Psat=self.Psat
        Kl1=self.Kl1
        Kl2=self.Kl2
        kg3=self.kg3
        H=self.H
        A=self.A
        a=self.a
        L1=50
        L2=50
        L3=130
        P=self.P
        def diff(G,z):
            
            

            
        
            
            
            
            
            VL=G[5]*.018/1000
            dG1=Kl1*a*A*((((G[0]*P)/(G[0]+G[1]+G[2]))/H1)-(G[3]/VL))
            dG2=Kl2*a*A*((((G[1]*P)/(G[0]+G[1]+G[2]))/H2)-(G[4]/VL))
            dG3=kg3*a*A*(Psat-((G[2]/(G[2]+G[1]+G[0]))*P))
            dL1=dG1
            dL2=dG2
            dL3=dG3
            
            
            
            
            
        
            return [dG1,dG2,dG3,dL1,dL2,dL3]
            
        
        Gi=[G10,G20,G30,L1,L2,L3]
        
        z=np.linspace(0,H,10)
        y=odeint(diff,Gi,z)
        print z
        print y
        dx1=0.1
        #dx2=.00001
        #dx3=.000001
        f=0
        err1=1
        err=1
        err2=1
    
        while (abs(err)>.01) or (abs(err1)>.01):
            Gi=[G10,G20,G30,L1,L2,L3] 
            a1=L1
            a2=L2
            a3=L3
            y=odeint(diff,Gi,z)
            err=y[9,5]-L3H
            err1=y[9,4]-L2H
            err2=y[9,3]-L1H
            print err,err1,err2    
            r=max(abs(err),abs(err1),abs(err2))
            
            if r==abs(err):
                L3=a3+dx1*r
            elif r==abs(err1):
                L2=a2-dx1*r
            elif r==abs(err2):
                L1=a1-dx1*r
                
        t1=y[0,3]+y[9,0]
        t2=y[0,4]+y[9,1]
        t3=y[0,5]+y[9,2]
        print y,t1,t2,t3
        
            
            
            
                
            
       
        return (err) 
           
            
            
            
        

