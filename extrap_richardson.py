# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:55:35 2020

@author: adria
"""

from math import *
from scipy.misc import derivative
import numpy as np
import matplotlib.pyplot as plt

### Programa de extrapolacion de Richardson

### Definicion de formulas de derivacion de orden 1.
### Todas son para la primera derivada.

def diff_prog(f,x0,h):
    yp=(f(x0+h)-f(x0))/h
    return yp

def diff_reg(f,x0,h):
    yp=(f(x0)-f(x0-h))/h
    return yp

### Definicion de formulas de derivacion de orden 2

def diff_cent(f,x0,h):
    ### Formula para derivada primera centrada en x0
    yp=(f(x0+h)-f(x0-h))/(2*h)
    return yp

#### Definicion de las formulas de Richardson, en orden ascendiente.
#### N1 es una formula de derivacion de orden 1, mientras que N2
#### es de orden dos (por ahora solo estoy tomando primeras derivadas)
def Richardson_1(N1,f,x0,h):
    yp=2*N1(f,x0,h/2)-N1(f,x0,h)
    return yp

def Richardson_2(N2,f,x0,h):
    yp=(4*N2(f,x0,h/2)-N2(f,x0,h))/3
    return yp

#### Forumlas generales de Richardson

### f_p(x0)=(2**i*N_i(f,x0,h/2)-N_i)/(2**i-1)
### N_i=(N_(i-1)(h/2)-N_(i-1)(h))/(2**(i-1)-1)+N_(i-1)


#### esta funcion nos genera las N_(j) bajo el paso h
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
    return N2,yp


def deriv_richardson(Njh,Njh2,j):
    yp=((2**j)*Njh2-Njh)/(2**j - 1)
    return yp


h=0.1
N3h,yp1=richardson_gen_N(diff_prog,sin,1,h,3)
N3h2,yp2=richardson_gen_N(diff_prog,sin,1,h/2,3)
print(N3h)
print(N3h2)

deriv=deriv_richardson(N3h[0],N3h2[0],3)
print(deriv)
print(yp1)