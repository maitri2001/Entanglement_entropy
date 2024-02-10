# -*- coding: utf-8 -*-


"""
Analysis on the Propagators without measurements
"""

import matplotlib.pyplot as plt
import math
import numpy as np
import sympy as sy
import cmath
from scipy.integrate import quad 
#import pandas as pd 
from mpl_toolkits.mplot3d import Axes3D 
n= 0.01 #(w0*(2*pi/L)=2*delta(q) )
T =0.000000001
m =1
j = cmath.sqrt(-1)
w_0 = 1

## Define the functions

def K(q):
    return 4*(w_0**2)*(np.sin(q/2)**2)
def A(q,w):
    return m*((w+j*n)**2)-m*((w_0)**2)*K(q)

def B(q,w):
    return m*((-w+j*n)**2)-m*((w_0)**2)*K(q)

def C(q,w):
    return -4*j*n*w*(np.tanh(w/(2*T))**(-1))

#def C(q,w):
    #return (A(q,w)-B(q,w))*(np.tanh(w/(2*T))**(-1))


def G_k(q,w):
    return np.divide( -C(q,w),(A(q,w)*B(q,w)))

def G_R(q,w):
    return 1/B(q,w)

def G_A(q,w):
    return 1/A(q,w)
#def G_k2(q,w):
    #return (np.tanh(w/(2*T)))**(-1)*(G_R(q, w)-G_A(q, w))
#def G_k3(q,w):
    return -(np.tanh(w/(2*T)))**(-1)*(A(q,w)-B(q,w))/(A(q,w)*B(q,W))


q = np.arange(-np.pi,np.pi,0.01)
w = np.arange(-2,2,0.01)

#############################################
# absolute value of the Keldysh part
plt.plot(q,abs(G_k(q,5)), label = "Keldysh Green's function (w=5) ")
plt.plot(q,abs(G_k(q,2)), label = "Keldysh Green's function (w=2) ")
plt.plot(q,abs(G_k(q,2.1)), label = "Keldysh Green's function (w=2.1) ")
plt.plot(q,abs(G_k(q,0.5)), label = "Keldysh Green's function (w=0.5) ")
plt.xlabel('q')
plt.ylabel('G_k(q)')
plt.title("Keldysh Greens's Function")
plt.legend()
#plt.show()

plt.plot(w,abs(G_k((np.pi/3),w)),label = "Keldysh Green's function (q=pi/3)")
plt.plot(w,abs(G_k((0),w)),label = "Keldysh Green's function (q=0)")
plt.plot(w,abs(G_k((np.pi/2),w)),label = "Keldysh Green's function (q=pi/2)")
plt.plot(w,abs(G_k((np.pi),w)),label = "Keldysh Green's function (q=pi)") 
plt.xlabel('w')
plt.ylabel('G_k(w)')
plt.title("Keldysh Greens's Function")
plt.legend()
plt.show()

# absolute value of the Retarded part
plt.plot(q,abs(G_R(q,5)), label = "Retarded Green's function (w=5) ")
plt.plot(q,abs(G_R(q,2)), label = "Retarded Green's function (w=2) ")
plt.plot(q,abs(G_R(q,2.1)), label = "Retarded Green's function (w=2.1) ")
plt.plot(q,abs(G_R(q,0.5)), label = "Retarded Green's function (w=0.5) ")
plt.xlabel('q')
plt.ylabel('G_R(q)')
plt.title("Retarded Greens's Function")
plt.legend()
#plt.show()

plt.plot(w,abs(G_R((np.pi/3),w)),label = "Retarded Green's function (q=pi/3)")
plt.plot(w,abs(G_R((0),w)),label = "Retarded Green's function (q=0)")
plt.plot(w,abs(G_R((np.pi/2),w)),label = "Retarded Green's function (q=pi/2)")
plt.plot(w,abs(G_R((np.pi),w)),label = "Retarded Green's function (q=pi)") 
plt.xlabel('w')
plt.ylabel('G_R(w)')
plt.title("Retarded Greens's Function")
plt.legend()
#plt.show()

# absolute value of the Advanced part
plt.plot(q,abs(G_A(q,5)), label = "Advanced Green's function (w=5) ")
plt.plot(q,abs(G_A(q,2)), label = "Advanced Green's function (w=2) ")
plt.plot(q,abs(G_A(q,2.1)), label = "Advanced Green's function (w=2.1) ")
plt.plot(q,abs(G_A(q,0.5)), label = "Advanced Green's function (w=0.5) ")
plt.xlabel('q')
plt.ylabel('G_A(q)')
plt.title("Advanced Greens's Function")
plt.legend()
#plt.show()

plt.plot(w,abs(G_A((np.pi/3),w)),label = "Advanced Green's function (q=pi/3)")
plt.plot(w,abs(G_A((0),w)),label = "Advanced Green's function (q=0)")
plt.plot(w,abs(G_A((np.pi/2),w)),label = "Advanced Green's function (q=pi/2)")
plt.plot(w,abs(G_A((np.pi),w)),label = "Advanced Green's function (q=pi)") 
plt.xlabel('w')
plt.ylabel('G_A(w)')
plt.title("Advanced Greens's Function")
plt.legend()
#plt.show()

#############################################################
#Comparison of the real parts of Advanced and Retarded Green's Function
plt.plot(w,(G_R((np.pi/3),w).real),label = "Real part of Retarded Green's function (q=pi/3)")
plt.plot(w,(G_A((np.pi/3),w).real),label = "Real part of Advanced Green's function (q=pi/3)")
plt.xlabel('w')
plt.ylabel('real part')
plt.legend()
plt.show()

#Comparison of the imaginary parts of Advanced and Retarded Green's Function
plt.plot(w,(G_R((np.pi/3),w).imag),label = "Imaginary part of Retarded Green's function (q=pi/3)")
plt.plot(w,(G_A((np.pi/3),w).imag),label = " Imaginary part of Advanced Green's function (q=pi/3)")
plt.plot(w,(G_k((np.pi/3),w).imag),label = " Imaginary part of keldysh Green's function (q=pi/3)")
plt.xlabel('w')
plt.ylabel('Imaginary part')
plt.legend()
plt.show()

###############################################################################
# 2D intensity plot of Imaginary part of retarded Green's Function as a function of (q,w)
fig = plt.figure() 
ax = plt.axes(projection="3d")
q = np.arange(-np.pi,np.pi,0.01)
w = np.arange(0,2,0.01)
Q,W = np.meshgrid(q,w)
Z1 = (G_R(Q,W).imag)
surf = ax.plot_surface(Q, W, Z1)
#plt.colorbar
#ax.view_init(90,0) 
plt.xlabel('Q')
plt.ylabel('W')
plt.title("2D intensity plot of Imaginary part of G_R(q,w)")
plt.show()
##########################################################################
# 2D intensity plot of Imaginary part of Keldysh Green's Function as a function of (q,w)
fig = plt.figure() 
ax = plt.axes(projection="3d")
q = np.arange(-np.pi,np.pi,0.01)
w = np.arange(0,2,0.01)
Q,W = np.meshgrid(q,w)
Z2 = (G_k(Q,W).imag)
surf = ax.plot_surface(Q, W, Z2)
#plt.colorbar
#ax.view_init(45,0) 
plt.xlabel('Q')
plt.ylabel('W')
plt.title("2D intensity plot of Imaginary part of G_k(q,w)")
plt.show()
#############################################################################
# 2D intensity plot comparison between the imaginary part of the Keldysh and Retarded Green's Function
fig = plt.figure() 
ax = plt.axes(projection="3d")
q = np.arange(-np.pi,np.pi,0.01)
w = np.arange(0,2,0.01)
Q,W = np.meshgrid(q,w)
Z1 = 10000*(G_R(Q,W).imag)
Z2 = (G_k(Q,W).imag)
surf = ax.plot_surface(Q, W, Z1)
surf = ax.plot_surface(Q, W, Z2)
ax.view_init(50,45)
plt.xlabel('Q')
plt.ylabel('W')
plt.legend()
#plt.title("")
plt.show()
##################################################################
fig = plt.figure() 
ax = plt.axes(projection="3d")
q = np.arange(-np.pi,np.pi,0.01)
w = np.arange(0,10,0.01)
Q,W = np.meshgrid(q,w)
R1 = np.divide(G_k(Q,W),(2*j*(G_R(Q,W).imag)))
#R1 = np.divide(G_k3(Q,W),(G_R(Q, W))-(G_A(Q,W)))
surf = ax.plot_surface(Q, W, R1)
#plt.colorbar
ax.view_init(0,0) 
plt.xlabel('Q')
plt.ylabel('W')
plt.show()
###########################################################################

def D_k_p(I,J):
    total_sum = 0
    for q in np.arange(-np.pi,np.pi,0.01):
        resultReal = quad(lambda w: G_k(q,w).real,-5,4.99)[0]
        resultImag = quad(lambda w: G_k(q,w).imag,-5,4.99)[0]
        Result = resultReal + j*resultImag # after integration wrt w on D_k(q,w)
        total_sum = total_sum +  np.exp(j*q*(I-J))*Result*(1/(2*np.pi))
    return total_sum
    
#print (D_k_p(2,10))
    
    
    

     




########################################################
'''# Analysis on the 4 by 4 propagator
def inverse(arr):
    return 1/(arr[0,0]*arr[1,1]-arr[0,1]*arr[1,0])*(np.array([[arr[1,1],-arr[0,1]],[-arr[1,0],arr[0,0]]]))
def K2(q):
    return 4*(sy.sin(q/2)**2)

def D_A_in(q,w):
    return np.array( [[-1,(-w-j*n)],[(w-j*n),-w_0**2*K2(q)]])
#print (D_A_in)
def D_R_in(q,w):
    return np.array([[-1,(-w + j*n)],[(w + j*n),-w_0**2*K2(q)]])
  
#print(np.linalg.inv(D_A_in))
def D_A(q,w):
    return inverse(D_A_in(q,w))
    
def D_R(q,w):
    return inverse(D_R_in(q,w))
    
def D_k(q,w):
    return ((sy.tanh(w/(2*T))**(-1)))*np.subtract(D_R(q,w),D_A(q,w))

q, w = sy.symbols('q w')

fig = plt.figure() 
ax = plt.axes(projection="3d")
q = np.arange(-np.pi,np.pi,0.01)
w = np.arange(0,10,0.01)
Q,W = np.meshgrid(q,w)
R2 = np.divide(D_k(Q,W),(2*j*(D_R(Q,W).imag)))
#R1 = np.divide(G_k3(Q,W),(G_R(Q, W))-(G_A(Q,W)))
surf = ax.plot_surface(Q, W, R2)
#plt.colorbar
ax.view_init(0,0) 
plt.xlabel('Q')
plt.ylabel('W')
plt.show()    
    
#print(D_k(q,w))
def D_kpp(q,w):
    return D_k(q,w)[0,0]
#print(D_kpp(q,w))
def D_kxx(q,w):
    return D_k(q,w)[1,1]


def D_kpp_position(I,J):
    total_sum = 0
    for q in np.arange(-np.pi,np.pi,0.01):
        resultReal = quad(lambda w: D_kpp(q,w).real,-5,4.99)[0]
        resultImag = quad(lambda w: D_kpp(q,w).imag,-5,4.99)[0]
        Result = resultReal + j*resultImag # after integration wrt w on D_k(q,w)
        total_sum = total_sum +  np.exp(j*q*(I-J))*Result*(1/(2*np.pi))
    return total_sum


print (D_k(7,8)[0,0])'''

####################################################################
# Analysis on the 4 by 4 propagator
m0 =0.00001
T =0.00000001
N = 100
q = np.arange(-np.pi,np.pi,2*np.pi/N)

def w_td(q):
    return cmath.sqrt(K(q)+m0)

def iota_D_kxx(q,I,J):
    return (1)*(np.exp(j*q*(I-J))/(m*w_td(q)))

def iota_D_kpp(q,I,J):
    return (1)*(np.exp(j*q*(I-J))*(m*w_td(q)))




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


'''two_d_arrayx = []    
for I in range (10):
    row = []
    for J in range (10):
        value = iota_Ds_kxx(I,J)
        row.append(value)
    two_d_arrayx.append(row)    
    
two_d_arrayp = []    
for I in range (10):
    row = []
    for J in range (10):
        value = iota_Ds_kpp(I,J)
        row.append(value)
    two_d_arrayp.append(row)
arrayx = np.array(two_d_arrayx)
arrayp = np.array(two_d_arrayp)
prod_matrix = np.dot(arrayx,arrayp) 
eigenvalues = np.linalg.eigvals(prod_matrix)
prod_eigenval= np.prod(np.array(eigenvalues))
print(0.5*np.log(prod_eigenval))'''  

'''S = []
for La in np.arange(5,25,1):
    two_d_arrayx = []    
    for I in range (La):
        row = []
        for J in range (La):
            value = iota_Ds_kxx(I,J)
            row.append(value)
        two_d_arrayx.append(row)    
        
    two_d_arrayp = []    
    for I in range (La):
        row = []
        for J in range (La):
            value = iota_Ds_kpp(I,J)
            row.append(value)
        two_d_arrayp.append(row)
    arrayx = np.array(two_d_arrayx)
    arrayp = np.array(two_d_arrayp)
    prod_matrix = np.dot(arrayx,arrayp) 
    U, Sing_val, VT = np.linalg.svd(prod_matrix)
    eigenvalues = np.linalg.eigvals(prod_matrix)
    prod_eigenval= np.prod(np.array(eigenvalues))
    prod_Sing_val= np.prod(Sing_val)
    S.append (0.5*np.log( prod_Sing_val))

 
Sub_size = np.arange(5,25,1)
plt.scatter(Sub_size,S)
plt.xlabel("La")
plt.ylabel("SA")
plt.show()'''

'''def S(La):
    two_d_arrayx = []    
    for I in range (La):
        row = []
        for J in range (La):
            value = iota_Ds_kxx(I,J)
            row.append(value)
        two_d_arrayx.append(row)    
        
    two_d_arrayp = []    
    for I in range (La):
        row = []
        for J in range (La):
            value = iota_Ds_kpp(I,J)
            row.append(value)
        two_d_arrayp.append(row)
    def arrayx(La):
        return np.array(two_d_arrayx)
    def eigval_x(La):
        return np.linalg.eigvals(arrayx(La))
    def arrayp(La):
        return np.array(two_d_arrayp)
    prod_matrix = np.dot(arrayx(La),arrayp(La)) 
    eigenvalues = np.linalg.eigvals(prod_matrix)
    prod_eigenval= np.prod(np.array(eigenvalues))
    return (0.5*np.log(prod_eigenval))'''

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
    U, Sing_val, VT = np.linalg.svd(prod_matrix)
    #print("Prod matrix =",prod_matrix)
    Single_values = Sing_val[2:]
    eigenvalues = (1/N**2)*(np.linalg.eigvals(prod_matrix)).real
    print("\n The eigen values of XP are",eigenvalues)
    #print("\n The singular values are",Sing_val)
    prod_eigenval= np.prod(np.array(eigenvalues))
    prod_Sing_val= np.prod(Single_values)
    return np.sum(np.log(eigenvalues))
    #return (0.5*np.log(prod_eigenval))


#print("\n Eigvals of X",np.linalg.eigvals(arrayx(5)))
    
#("Value of S",S(50)) 

#La_values =[1,2,3,4,5,15,25,35,45,55,65,75,85,95,96,97,98,99,100]
La_values = [5,10]
SA = []

for La in La_values:
    SA.append(S(La).real)
    
plt.scatter(La_values,SA)
plt.xlabel("La")
plt.ylabel("SA")
plt.show()