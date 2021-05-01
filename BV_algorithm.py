# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 16:41:38 2020

@author: adria
"""
from qiskit import *

#s=1010001       ### cadena 

#x1=1000000           #### input que asignamos a la cadena 
#x2=0100000          ### asi sucesivamente hasta siete veces, que es el numero de bits
### de esta forma, podemos designar los caracteres de s

### pero estamos en computacion cuantica, no?

### hacemos un estado inicial 0 (up) y el output sera el estado -

##aplicamos puerta hadamar al input

##luego oracle

##otra vez hadamar

### medir la salida. Empezamos con dos qubits

### comenzamos con la cadena secreta

#s=101011

qc=QuantumCircuit(6+1,6)        ## 6+1 qubits y 6 bits clasicos
qc.h([0,1,2,3,4,5])
qc.x(6)
qc.h(6)

qc.barrier()        ### separa operaciones entre elementos del circuito cuantico

qc.cx(5,6)      ### lo colocamos en las zonas donde haya unos con el seis
qc.cx(3,6)
qc.cx(1,6)
qc.cx(0,6)

qc.barrier()

qc.h([0,1,2,3,4,5])


qc.measure([0,1,2,3,4,5],[0,1,2,3,4,5])

print(qc.draw())        ### los 6 bits clasicos nos dan el resultado de la medida

### despues de la medida, obtenemos bits clasicos en nuestro registrador c:6

### vamos a ejecutar este circuito para ver si obtenemos s

simulation = Aer.get_backend('qasm_simulator')
result = execute(qc, backend = simulation, shots = 1).result()
counts = result.get_counts()

print(counts)