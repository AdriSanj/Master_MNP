t# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 16:40:39 2020

@author: adria
"""

import numpy as np
from math import *

##### Metodo de Jacobi con criterio de parada 
### en primer lugar, con error absoluto |x(k+1)-x(k)|<eps
## en este caso tomando un bucle for de bastantes iteraciones
## pero lo idoneo es un while. Ya que no piloto whiles, hago for.

### Definimos una matriz A = D -E -F. D diagonal, E triang inf, F triang sup

A=np.array([[5.,-2.,3],[-3.,9.,1.],[2.,-1.,-7]])

b=np.array([-1.,2.,.3])

#print(np.linalg.solve(A,b))            #### solucion exacta

xk=np.array([0.,0.,0.])         ###vector de soluciones. Si es de ceros, todos tienen que tener 0.

#print(A)
#print(np.linalg.det(A))         ### comprobacion de que es no singular

D=np.zeros((len(b),len(b)))         ### matriz diagonal
M=np.zeros((len(b),len(b)))         ### triangular inf y triang sup

for i in range(len(b)):
    for j in range(len(b)): 
        if i==j:
            D[i,j]=A[i,j]
        elif i>j or i<j:
            M[i,j]=A[i,j]
            

#print(D)
#print(M)


for k in range(30):         
    c=np.copy(xk)
    for i in range(len(b)):
        sum1=0          ### sumamos los elementos de la fila i con el vector copia c
        for j in range(0,len(b)):
            sum1+=M[i,j]*c[j]       ### en este metodo tenemos M asi porque ya se anula solo el elemento ii
        xk[i]=(-sum1+b[i])/D[i,i]
        #print(xk[i])
    if abs(np.dot(xk,xk)-np.dot(c,c))<1e-12:
        print("Iteracion final: "+str(k))
        break

print(xk)


