# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import math
import numpy as np
import sympy as sy
import cmath
from scipy.integrate import quad 
#import pandas as pd 
from mpl_toolkits.mplot3d import Axes3D 
n= 0.01 #(w0*(2*pi/L)=2*delta(q) )
m =1
j = cmath.sqrt(-1)
w_0 = 1
m0 =0.00001
T =0.00000001
N = 100
q = np.arange(-np.pi,np.pi,2*np.pi/N)


def K(q):
    return 4*(w_0**2)*(np.sin(q/2)**2)

def w_td(q):
    return cmath.sqrt(K(q)+m0)

## Using zero temperature limit

def iota_D_kxx(q,I,J):
    return (1/N)*(np.exp(j*q*(I-J))/(m*w_td(q)))

def iota_D_kpp(q,I,J):
    return (1/N)*(np.exp(j*q*(I-J))*(m*w_td(q)))




def iota_Ds_kxx(I,J):
    total_sum = 0
    for qi in q:
        total_sum = total_sum + iota_D_kxx(qi,I,J) 
    return total_sum

def iota_Ds_kpp(I,J):
    total_sum = 0
    for qi in q:
        total_sum = total_sum + iota_D_kpp(qi,I,J) 
    return total_sum 

def arrayx(La):
    two_d_arrayx = []    
    for I in range (La):
        row = []
        for J in range (La):
            value = iota_Ds_kxx(I,J)
            row.append(value)
        two_d_arrayx.append(row) 
    return np.array(two_d_arrayx)

def arrayp(La):
    two_d_arrayp = []    
    for I in range (La):
        row = []
        for J in range (La):
            value = iota_Ds_kpp(I,J)
            row.append(value)
        two_d_arrayp.append(row)
    return np.array(two_d_arrayp)
    
def S(La):
    print ("\n We are calculating for La =",La)
    prod_matrix = np.dot(arrayx(La),arrayp(La))
    
    print("\n The XP matrix is =",prod_matrix.real)
    eigenvalues = (np.linalg.eigvals(prod_matrix)).real
    print("\n The eigen values of XP are",eigenvalues)
    prod_eigenval= np.prod(np.array(eigenvalues))
    return 0.5*np.sum(np.log(eigenvalues))
    #return (0.5*np.log(prod_eigenval))


print(arrayx(2).real)
print(arrayp(2).real)
#print(prod_matrix(2).real)
#print("\n Eigvals of X",np.linalg.eigvals(arrayx(5)))
    
#("Value of S",S(50)) 
print("The entanglement entropy",S(10))
#La_values =[1,2,3,4,5,15,25,35,45,55,65,75,85,95,96,97,98,99,100]
La_values = np.arange(3,100,3)
#La_values = [5,10]
SA = []

for La in La_values:
    SA.append(S(La).real)
    
plt.scatter(La_values,SA)
plt.xlabel("La")
plt.ylabel("SA")
plt.show()    