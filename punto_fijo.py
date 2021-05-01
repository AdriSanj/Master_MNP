# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 22:00:44 2021

@author: adria
"""
import numpy as np


def punto_fijo(a0,g,it):
    for i in range(it):
        a=g(a0)
        a0=a
    return a

def dicotomia(a,b,f):
    c=(a+b)/2.
    fa=f(a)
    fb=f(b)
    fc=f(c)
    for i in range(100):
        if fa*fc<0:
            fc=fb
            b=c
            c=(a+b)/2
            fc=f(c)
        elif fc*fb<0:
            fc=fa
            a=c
            c=(a+b)/2
            fc=f(c)
    return c
        
        
        
        
def func(a):
    acos=np.arccos(np.exp(a)/4)
    return np.log(-acos**2-a**2)-np.log(np.sin(acos))+np.log(0.2*acos)

def func2(a):
    acos=np.arccos(np.exp(a)/4)
    b=np.log((-acos**2-a**2)*np.sin(acos)/(0.2*acos))
    return a-b

def func3(a):
    acos=acos=np.arccos(np.exp(a)/4)
    b=np.sqrt(-0.2*np.exp(a)*acos/np.sin(acos)+acos**2)
    return a-b

print(func3(0.5))
print(func3(0.8))

c=dicotomia(0.5,0.8,func3)

print(c)

print(func3(c))

print(np.sqrt(np.arccos(np.exp(c)/4)**2)+c**2)
