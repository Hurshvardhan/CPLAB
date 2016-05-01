# -*- coding: utf-8 -*-
"""
Created on Sun May 01 09:39:02 2016

@author: Hurshvardhai
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jan 09 19:11:50 2016

@author: Hurshvardhai
"" References
1.For Henrys constant values:
http://www.atmos-chem-phys.net/15/4399/2015/acp-15-4399-2015.pdf
2.kla For CO2 absorption:
Reporting of well stirred scrubber results: 
scrubbing of SO2 and CO2 by caustic 
solutions at atmospheric pressure 
Final Report 
Yinghui Liu, Dunyu Liu, and Terry Wall 
Chemical Engineering, the University of Newcastle, NSW, Australia
3.For kga of water: 
http://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html
4.For kla of Methane:
http://r.search.yahoo.com/_ylt=A2oKmM95TZVWwiQA_XHnHgx.;_ylu=X3oDMTByMjN2cHBuBHNlYwNzcgRwb3MDNwRjb2xvA3NnMwR2dGlkAw--/RV=2/RE=1452654073/RO=10/RU=http%3a%2f%2fwww.researchgate.net%2fpublication%2f225314921_Unsteady-State_Methane_Absorption_by_Water_before_Hydrate_Formation/RK=0/RS=QH1IkoCsorqa.nLpXaizNSXt5Go-

5.We shall take standard NTP conditions for the opration
Thus, T=298K P=1atm=101325Pa"""
import scipy
import numpy
import matplotlib.pyplot as plt
class scrubber:
    
    "5"
    P=101325.0#Pascals
    T=298.0#K
    xgc=0.5#mol fraction of CO2 in gas
    xgm=0.5#mol fraction of methane in gas
    xgw=0.0#mol fraction of water in gas
    n=20.0#step size but this can be changed
    "1"
    Hmlit=0.00014
    "1"
    Hclit=0.00034#in mol/(m3*Pa)
    Hm=(1/Hmlit)#Henrys const for methane in (Pa*m3/mol)
    Hc=(1/Hclit)#Henrys const for CO2 in (Pa*m3/mol)
    "2"
    klaclit=0.5249502#hr-1
    klac=klaclit/3600#kla for CO2 in s-1
    "4"
    klamlit=0.02#hr-1
    klam=(klamlit/3600)#kla for CH4 in s-1
    "3"
    kgawlit=158.89*(10**-6)#kg/m3*hr*Pa
    kgaw=kgawlit/(0.018*3600)#kga for evaporation of water in (mol/(m3*s*Pa))
    xLc=0.0#mol fraction of CO2 in liquid
    xLm=0.0
    xLw=1.0
    H=7.0#Height of scrubber
    A=5.0#Area for mass transfer assumed as 5m^2(a 5 m tall scrubber with 1m dia)
    """We define a class which will help the user use the concept of abstraction
    by interacting with the program via the object"""
    def sol(self):
        P=self.P;T=self.T
        xgc=self.xgc;xgm=self.xgm;xgw=self.xgw
        H=self.H;A=self.A
        n=self.n
        Hm=self.Hm;Hc=self.Hc
        klac=self.klac;klam=self.klam;kgaw=self.kgaw
        xLc=self.xLc;xlm=self.xLm;xLw=self.xLw
        G=0.002778*10#gas flow rate should be as low as possible
        L=0.002778#10m3/hr liquid flow rate and this is nearly constant throughout the operation
        Ln=(L*(10**6)/18)
        
        #0-water
        #1-CO2
        #2-methane
        i=0
        e=100
        h=H/n
        ng =numpy.array((n,3)) 
        nL =numpy.array((n,3))
        Gn=(P*G)/(8.314*T)#Total gas molar flow rate
        ngw=xgw*Gn
        ngc=xgc*Gn
        ngm=xgm*Gn
        nLw=150.0
        nLc=0.5
        nLm=0.1
        
        Gv1=nLw
        Gv2=nLc
        Gv3=nLm
        
        """The initial values of gas phase molar flow rate are fixed as per the
        conditions selected and the liquid phase molar flow rates are guessed
        assuming nearly 100% absorption of both CO2 and methane given that is the
        maximum possible concn at the bottom.
        Since the values are going to be changed in further iterations, it is important
        to store these numbers as Gv1,Gv2 and Gv3"""
        fig1 = plt.figure(1)
    
    
        while e>0.01:
            x=0
            e1=abs((nLw-Ln)*0.1) #water condition can be lax as compared to other two
            
            """The error of water is calculated based on the fact that we should have
            a flow rate of Ln at the top(calculated above). 0.1 is multiplied as
            this error can be more relaxed in magnitude as compared to the other
            two errors"""
            
            e2=abs(nLc)
            
            e3=abs(nLm)
            """Molar flow rate of CO2 AND CH4 at the top should be 0"""
            maxi=e1
            if(e2>maxi):
                maxi=e2
            if(e3>maxi):
                maxi=e3
            e=maxi
            #This was to treat the max of the errors first by conjugate gradient method
            """Conjugate gradient method. It makes no sense changing the guess values
            for two components if they have converged and only the one for the
            third component has to be manipulated. Thus we select the max error
            and attempt to work on it first :)"""
            j=0
            if(e<0):
                e=-e
             #This was to take absolute value   
            if(i>0):
                if (e==e1):
                    nLw=Gv1+((Ln-nLw)*0.1)
                    Gv1=nLw
                    nlc=Gv2
                    nlm=Gv3
                    """Shooting method! There need not be modifications in the other two guess
                    values while manipulating that of one component."""
                if(e==e2):
                    nLc=Gv2-0.1*nLc
                    Gv2=nLc
                    nlw=Gv1
                    nlm=Gv3
                if(e==e3):
                    nLm=Gv3-0.1*nLm
                    Gv3=nLm
                    nlc=Gv2
                    nlw=Gv1
            Gn=(P*G)/(8.314*T)
            Ln=(L*(10**6)/18)
            ngw=0.0*Gn;www=ngw#these are saved for later use
            ngc=0.5*Gn;ccc=ngc
            ngm=0.5*Gn;mmm=ngm
            """We should store these values of molar gas
            flow rate at the bottom since they will be modified"""
               
            while x<H:
                Gn=ngw+ngc+ngm
                
                k1=kgaw*A*(3115-(P*(ngw/Gn)))#psat=3115 at 1atm
                k2=kgaw*A*(3115-(P*((ngw+0.5*k1*h)/Gn)))
                k3=kgaw*A*(3115-(P*((ngw+0.5*k2*h)/Gn)))
                k4=kgaw*A*(3115-(P*((ngw+k3*h)/Gn)))
                """Range Kutta-4 method for ode soln"""
                f1=klac*A*(((P*(ngc/Gn))/Hc)-(nLc/(nLw*0.000018)))
                f2=klac*A*(((P*((ngc+0.5*f1*h)/Gn))/Hc)-((nLc)/(nLw*0.000018)))
                f3=klac*A*(((P*((ngc+0.5*f2*h)/Gn))/Hc)-((nLc)/(nLw*0.000018)))
                f4=klac*A*(((P*((ngc+f3*h)/Gn))/Hc)-((nLc)/(nLw*0.000018)))
                
                m1=klam*A*(((P*(ngm/Gn))/Hm)-(nLm/(nLw*0.000018)))
                m2=klam*A*(((P*((ngm+0.5*m1*h)/Gn))/Hm)-((nLm)/(nLw*0.000018)))
                m3=klam*A*(((P*((ngm+0.5*m2*h)/Gn))/Hm)-((nLm)/(nLw*0.000018)))
                m4=klam*A*(((P*((ngm+m3*h)/Gn))/Hm)-((nLm)/(nLw*0.000018)))
                
                ngw=ngw+((k1+2*k2+2*k3+k4)/6)*h
                ngc=ngc-((f1+2*f2+2*f3+f4)/6)*h
                ngm=ngm-((m1+2*m2+2*m3+m4)/6)*h
                
                                
                
                nLw=nLw+((k1+2*k2+2*k3+k4)/6)*h
                nLc=nLc-((f1+2*f2+2*f3+f4)/6)*h
                nLm=nLm-((m1+2*m2+2*m3+m4)/6)*h
                #this is rk4
                x=x+h
                j=j+1
            i=i+1
 
        
        f43=(0.5*Gn)-ngc
        return f43

i=0.5    
while i<=10.0:      
    sc1=scrubber()
    sc1.H=i
    f=sc1.sol()
    plt.plot(i,f,'bo')
    
    sc1.sol()
    i=i+0.5
    plt.xlabel('Height of column')
    plt.title('ABSORBED CO2 FOR DIFFERENT HEIGHTS')
        
    plt.ylabel('CO2 ABSORBED')
plt.show()
           
