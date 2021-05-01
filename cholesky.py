# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:12:31 2020

@author: adria
"""
from math import *
import numpy as np

def cholesky2(A):       ### esta es la mia
    n=len(A)
    L=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            suma=0
            for k in range(j):
                suma+=L[i][k]*L[j][k]
            if i == j:
                L[i][i]=sqrt(A[i][i]-suma)
            else:
                L[i][j]=(A[i][j]-suma)/L[j][j]
    return L

def cholesky(A):
    """Performs a Cholesky decomposition of A, which must 
    be a symmetric and positive definite matrix. The function
    returns the lower variant triangular matrix, L."""
    n = len(A)

    # Create zero matrix for L
    L = [[0.0] * n for i in range(n)]

    # Perform the Cholesky decomposition
    for i in range(n):
        for k in range(i+1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in range(k))
            
            if (i == k): # Diagonal elements
                # LaTeX: l_{kk} = \sqrt{ a_{kk} - \sum^{k-1}_{j=1} l^2_{kj}}
                L[i][k] = sqrt(A[i][i] - tmp_sum)
            else:
                # LaTeX: l_{ik} = \frac{1}{l_{kk}} \left( a_{ik} - \sum^{k-1}_{j=1} l_{ij} l_{kj} \right)
                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L
 
A = [[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]]
L = cholesky2(A)
L2=np.asarray(cholesky(A))

print(A)

L=cholesky(A)

print(np.matmul(L,np.transpose(L)))