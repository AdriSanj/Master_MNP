# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 16:51:28 2020

@author: adria
"""
### Programa de desintegraciones nucleares

from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative

def diff_prog(f,x0,h):      ### formula de orden 1 a introducir en richardson
    yp=(f(x0+h)-f(x0))/h
    return yp


def richardson_gen_N(N1,f,x0,h,j):
    if j==1:
        yp=2*N1(f,x0,h/2)-N1(f,x0,h)
    elif j>1:
        N=np.zeros(j+1)     ### si orden es j, necesitamos j+1 terminos de 
                            ### orden 1 para el NJ y j+2 para yp 
        for l in range(0,j+1):
            N[l]=N1(f,x0,h/2**l)
        ### A partir de aqui empieza el bucle
        for i in range(2,j+1):
            N2=np.zeros(j-i+2)    ### Vector de resultados asociado al paso i
            for k in range(0,len(N2)):
                N2[k]=(N[k+1]-N[k])/(2**(i-1)-1)+N[k+1]
            N=np.copy(N2)   ###necesario para el siguiente paso
        ### Aqui la derivada de richardson
        yp=((2**j)*N2[1]-N2[0])/(2**j-1)
    return yp


def Exponents(L,t):
    L=np.reshape(L,(len(L),1))
    E=np.exp(-L*t)
    return E

def Acoef(L):
    #### Funcion que calcula los coeficientes que multiplican a cada exponencial 
    l=len(L)
    A=np.zeros((l,l))
    for i in range(0,l):
        for j in range(0,l):
            if j>=i:
                num=1
                den=1
                for k in range(0,j+1):
                    if k<j:
                        num=L[k]*num
                    if i!=k:
                        den=(L[k]-L[i])*den
                A[j,i]=num/den
    return A

def desintegrations(A,E,t,L):
    return np.matmul(A,E(L,t))



Lambdas=np.array([1.,1.4,5.,3.]) 

MA=Acoef(Lambdas)
print(MA)
        ### vector de constantes de decaimiento (inverso de la vida media)

time=np.linspace(0,10,100)

#### en esta parte estamos poniendo las graficas del numero de desintegraciones

N=desintegrations(MA,Exponents,time,Lambdas)
N=10000*N    #### multiplicamos por el numero inicial de particulas
for k in range(len(Lambdas)):
    plt.plot(time,N[k])

plt.show()


