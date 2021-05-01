# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:11:18 2020

@author: adria
"""

#### IMPORTANTE --------------------
#### Con esto, debería declarar funciones y así poder copiar y pegar codigo para el futuro



from math import *
import numpy as np

##Metodos de resolucion de sistemas de ecuaciones lineales. Todos son sistemas compatibles determinados

## para comprobar que es un sistema compatible determinado, la matriz ha de tener inversa (np.linalg.inv())

### Supongamos un sistema diagonal (mas sencillo que cualquier otro)

#Diag=np.array([[2.,0,0],[0,-2.,0],[0,0,-5.]])   ### matriz de coeficientes

#b1=np.array([2.,-4.,10.])          ###array del sistema A x = b

### invertimos la matriz Diag, multiplicamos por b1 y sale. Este es el metodo rapido.
##Dinv=np.linalg.inv(Diag)
##x1=np.dot(Dinv,b1)
##print(x1)

### Pero como todo en la vida, vamos a hacer un esfuerzo para aprender.
def diagres(Diag,b1):
    x1=[]     ##creamos la lista de resultados. La convertiremos mas adelante en array
    for i in range(len(b1)):        ## tomamos los elementos de la diagonal y se los dividimos a cada b[i]
        x1.append(b1[i]/Diag[i,i])
    return x1
    
      #conversion a array



###### triangular superior
b2=np.array([4.,2.,-3.])
Tsup=np.array([[2.,-3.,1.],[0,-2.,5.],[0,0,-5.]])

def triangsup(Tsup,b2):
    b=len(b2)-1 ### ya que lo vamos a usar, denotamos b como la longitud del vector de coeficientes y asi 
                ### no tenemos que escribir len(b2) todo el tiempo. Sin embargo, antes lo ponemos porque
                ### si no, nos crea un array de soluciones x2 de menos dimension
    x2=np.ones(len(b2))
    x2[b]=b2[b]/Tsup[b,b]

    for i in range(b-1,-1,-1):        ##bucle inverso
        coef=0                         ### elementos que vamos a restar en el siguiente paso para obtener x[i]
        for j in range(i+1,b+1):            ## desde el elemento i+1 hasta el elemento b, que se obtiene poniendo de paso final b+1 (cosas de python)
            coef+=Tsup[i,j]*x2[j]            ### en cada paso añadimos otro mas en función de cuantos tengammos de antes)
            x2[i]=(b2[i]-coef)/Tsup[i,i]        ##una vez obtenidos todos los coeficientes, se los restamos al elemento del vector b2
    return x2




##### Triangular inferior

### voy a trabajar con la misma por pura pereza, pero era bueno comprobar una 5 X 5

b3=np.array([-3.,2.,4.])

Tinf=np.array([[-5.,0,0],[5.,-2.,0],[1.,-3.,2.]])

def trianginf(Tinf,b3):
    x3=np.ones(len(b3))
    x3[0]=b3[0]/Tinf[0,0]
    for i in range(1,len(b3)):
        coef=0
        for j in range(0,i):            ### recordemos que el primer elemento es el 0-esimo
            coef+=Tinf[i,j]*x3[j]       ##analogo al anterior
            x3[i]=(b3[i]-coef)/Tinf[i,i]
    return x3


#### Metodo de reduccion de Gauss

## Y lo facil que sería hacerlo a mano
### como no dan muchas señales de ser ultra necesario, lo hare en el futuro si eso
## centrandome primero en el metodo de LU


#### Metodo de descomposicion LU


LU=np.array([[2.,3.,-2],[0,-2.,6.],[-4.,1.,0]])      ### matriz a factorizar

VecLU=np.array([2.,-1.,4.])


#print(LU)

### vamos a sustituir la matriz inicial por el producto de las matrices L y U
### l incluye en la diagonal unos, y un triangulo inferior. U  es triangula superior
### a la vez que obtenemos las L_ij y las U_ij, las vamos sustituyendo en LU

def facLU(LU,VecLU):
    for k in range(1,len(VecLU)):       #primer paso, para la columna 0
        LU[k,0]=LU[k,0]/LU[0,0]
    
    print("Tras la primera iteracion")
    print(LU)
    
    for i in range(1,len(VecLU)):
        for j in range(1,len(VecLU)):
            if i<=j:
                coefi=0
                for l in range(0,i):        ##recuerda siempre que el rango llega al i-1
                    coefi+=LU[i,l]*LU[l,j]      ##esto nos mete los coeficientes superiores y diagonales
                LU[i,j]=LU[i,j]-coefi       #### elementos correspondientes a U
            elif i>j:
                coefj=0
                for m in range(0,j):        ##aqui tambien, esto es j-1
                    coefj+=LU[i,m]*LU[m,j]      ##coeficientes inferiores
                LU[i,j]=(LU[i,j]-coefj)/LU[j,j]         ## elementos correspondientes a L

    print("Iteracion final")
    print(LU)

        ###la virgueria es separar la matriz LU en dos, L y U
            ### automaticamente invertir L con b, habiendo definido y=Ux, asi obtenemmos y
                ### luego, invertimos U y obtenemos x=U**(-1)y

    U=np.zeros((len(VecLU),len(VecLU)))     ##vamos a obtener la matriz U copiando los elementos de la matriz LU
    for i in range(len(VecLU)):
        for k in range(len(VecLU)):
            if i<k or i==k:             ### no se por que pero pensaba que la sintaxis era al reves
                U[i,k]=LU[i,k]
            
    invU=np.linalg.inv(U)

    #print(invU)

    L=np.identity(len(VecLU))

    for i in range(len(VecLU)):
        for k in range(len(VecLU)):
            if i>k:
                L[i,k]=LU[i,k]

    invL=np.linalg.inv(L)

    #print(invL)        ### es casi la misma, cambia el signo en la diagonal inferior

    Y=trianginf(invL,VecLU)

    Xfin=triangsup(invU,Y)
    
    return Xfin


###### Metodo de factorizacion QR ---- HACER EN FUNCION
    
## Q es una matriz ortogonal y R una triangular superior


def mat(v):
    m=np.zeros((len(v),len(v)))
    for i in range(len(v)):
        for j in range(len(v)):
            m[i,j]=v[i]*v[j]
    return m
    
#Aini=np.array([[1.,3.,-2],[-2.,4.,-1.],[-1.,-1.,3.]])
#bini=np.array([1.,0,-2])

Aini=np.array([[0,0,-1.,-2.],[0,0,0,1.],[0,1.,2.,3.],[-1.,-2.,-3.,-4.]])
bini=np.array([2.,-1.,3.,5.])
print(Aini[0::,0])
#print(Aini[1::,1])
        ###recuerda generalizarlo para el resto de componentes de la matriz, cambiando el 0 por la i

m=[]        ## en esta lista almacenaremos las matrices H ortogonales donde aplicamos las transformaciones

for i in range(len(bini)-1):
    x=np.zeros(len(bini))
    x[i::]=Aini[i::,i]       ### funciona porque funciona
    print("vector x en la iteracion "+str(i))
    print(x)
    
    normx=np.linalg.norm(x) 
    print(normx)     
    y=np.zeros(len(x))
    signo=np.sign(x[i])
    if signo==0:            #### asumimos que el signo es positivo en caso de que sea cero
        signo=1
    y[i]=normx*signo

    print("vector y en la iteracion " +str(i))
    print(y)
    norma=np.linalg.norm(x-y)
    if norma==0:
        continue
    v=(x-y)/norma
    #print("vector v en la iteracion " +str(i))
    #print(v)
    matriz=mat(v)
    
    H=np.identity(len(v))-2*matriz
    m.append(H)
    Aini=np.matmul(H,Aini)
    for j in range(len(bini)):          ## los valores muy pequeños casi cero los hacemos cero
        for k in range(len(bini)):
            if abs(Aini[j,k])<1e-12:
                Aini[j,k]=0
    ##Ht=np.transpose(H)
    ###print(np.matmul(H,H))            #### es la matriz identidad, H=Ht
    print("La matriz diagonal R es, en la iteracion " + str(i) +":")        #### si la matriz es tocha, comentar estas tres lineas, que no hacen falta
    print(Aini)
    print("+---------------------+")

Q=m[0]          ### hacemos el producto de todas las matrices H empezando por la primera

for i in range(1,len(m)):
    Q=np.matmul(Q,m[i])

for j in range(len(bini)):          ## los valores muy pequeños casi cero los hacemos cero
        for k in range(len(bini)):
            if abs(Q[j,k])<1e-12:
                Q[j,k]=0
print("La matriz ortogonal Q es:")      
print(Q)
Qt=Q.transpose()
#print(Qt)
###print(np.linalg.inv(Q))          #### con esto comprobamos que Q es ortogonal

Yfin=np.dot(Qt,bini)
print(Yfin)
Xfin=triangsup(Aini,Yfin)
###solucion por fin
print(Xfin)


