# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 18:03:27 2016

@author: Hurshvardhai
"""
        
import numpy
#from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
import random
def sol():
    n=5
    
    c=numpy.zeros((n,n))
    n=n+2
    a=numpy.zeros((n,n))
    b=numpy.zeros((n,n))
    i=0
    j=0
    grass=0;her=0;carn=0;tp=0
    while(i<n):
        j=0
        while(j<n):
            if((i==0) or (j==0) or (i==(n-1)) or (j==(n-1))):
                a[i,j]=0
                b[i,j]=a[i,j]
            else:
                print("ENTER ELEMENT",(i),"th row ",j,"th column")
                a[i][j]=raw_input("num:")
                b[i][j]=a[i][j]
                c[(i-1)][(j-1)]=a[i][j]
            j=j+1
        i=i+1
    
    print c
    i=1
    no=input("No of times you want to run:")
    #colormap=colors.ListedColormap(["Brown","Green","Blue","Red"])
    plt.title("BARE LAND- BLACK, GRASS-RED, HERB-YELLOW, CARNIVORE-WHITE")
    plt.imshow(c,interpolation='nearest',cmap='hot',vmin=0,vmax=3)
    plt.show(block=False)
    plt.pause(3)
    m=0
    
    while(m<no):
        i=1
        cgr=0
        cher=0
        ccar=0
        while(i<(n-1)):
            j=1
            while(j<(n-1)):
                k=(i-1)
                her=0
                carn=0
                grass=0
                while(k<(i+2)):
                    l=(j-1)
                    while(l<(j+2)):
                        if((i==k) and (j==l)):
                            tp=0   
                        else:    
                            if(a[k][l]==1):
                                grass=grass+1
                            elif(a[k][l]==2):
                                her=her+1
                            elif(a[k][l]==3):
                                carn=carn+1
                        l=l+1
                    k=k+1
                if(a[i][j]==0):
                    if(grass>0):
                        b[i][j]=1#grass around
                    elif((her>0) and (grass>0)):
                        b[i][j]=2#grass+herbivore
                    elif((carn>0) and (her>0)):
                        b[i][j]=3#herbivore+carnivore
                elif(a[i][j]==1):
                    cgr=cgr+1
                
                elif(a[i][j]==2):
                    if(grass<1):
                        b[i][j]=0
                    #elif(carn>0):
                     #   b[i][j]=0
                    elif(her>2):
                        b[i][j]=0
                    cher=cher+1
                elif(a[i][j]==3):
                    if(her<1):
                        b[i][j]=0
                    elif(carn>2):
                        b[i][j]=0
                    ccar=ccar+1
                c[(i-1)][(j-1)]=b[i][j]                                            
                j=j+1
            i=i+1
        
        a=b
        """img=Image.new("1",(n,n))
        pixels=img.load()
        for i in range(img.size[0]):
            for j in range(img.size[1]):
                pixels[i,j]=a[i][j]
                img.show()"""
        
        print(c)
        print("No of grass",cgr)
        print("No of herbivores",cher)
        print("No of carnivores",ccar)
        plt.plot(m,ccar,'bo')
        plt.imshow(c,interpolation='nearest',cmap='hot',vmin=0,vmax=3)
        plt.show(block=False)
        
        plt.pause(0.3)
        plt.title("BARE LAND- BLACK, GRASS-RED, HERB-YELLOW, CARNIVORE-WHITE")
        m=m+1
        
        
    return her,carn
    plt.show()
    
    