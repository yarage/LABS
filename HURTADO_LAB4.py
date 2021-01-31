#! /usr/bin/python3

from numpy import *

import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import random

import math as math


##########
# P1
##########

def p1(T,theta):
    return (theta/T)**2*(3*8.314472)*exp(theta/T)/(exp(theta/T)-1)**2


# Las constantes serian estas, pero por el orden definimos:
# w1 = 6*1e12
# hbar = 1.054571*1e-34
# kb = 1.380649*1e-23


w1 = 6
hbar = 1.054571
kb = 1.380649

# Consideraremos las potencias de 10 por separado
#operando manualmente el exponente es 10^1


#definimos la temp de einst:

te1 = w1*hbar*10/kb

xgr1 = linspace(0,250,10000)
ygr1 = []
for i0 in xgr1:
    ygr1.append(p1(i0,te1))

plt.plot(xgr1,ygr1,color = 'red', label = 'frecuencia = 6 x 10^12 Hz')
plt.grid()
plt.legend()
plt.xlabel('Temperatura')
plt.ylabel('c_v')
plt.title('Modelo einstein - cv de solidos')
plt.show()


##########
# P2
##########

def Br(j,x):
    return 1/j*((j+1/2)*(1/math.tanh((j+1/2)*x))-1/2*(1/math.tanh(x/2)))

def im(n):
    return 2*9.274*1e-24*5/2*Br(5/2,n)

aux2 = linspace(-3500,3500,1000)
xgr2a = []
xgr2b = []
xgr2c = []
for i1 in aux2:
    T1 = 300
    xgr2a.append(2*9.274*1e-1*i1/(kb*T1))
    T2 = 500
    xgr2b.append(2*9.274*1e-1*i1/(kb*T2))
    T3 = 700
    xgr2c.append(2*9.274*1e-1*i1/(kb*T3))
ygr2a = []
ygr2b = []
ygr2c = []
for i in xgr2a:
    ygr2a.append(im(i))
for i in xgr2b:
    ygr2b.append(im(i))
for i in xgr2c:
    ygr2c.append(im(i))


plt.plot(aux2,ygr2a,color = 'red', label = 'T = 300K')
plt.plot(aux2,ygr2b,color = 'blue', label = 'T = 500K')
plt.plot(aux2,ygr2c,color = 'orange', label = 'T = 700K')
plt.grid()
plt.legend()
plt.xlabel('Campo H')
plt.ylabel('Imanación')
plt.title('Imanación vs H')
plt.show()

plt.plot(xgr2a,ygr2a,color = 'red', label = 'T = 300K')
plt.plot(xgr2b,ygr2b,color = 'blue', label = 'T = 500K')
plt.plot(xgr2c,ygr2c,color = 'orange', label = 'T = 700K')
plt.grid()
plt.legend()
plt.xlabel('eta')
plt.ylabel('Imanación')
plt.title('Imanación vs eta')
plt.show()

##########
# P4
##########

# Despejando la ecuacion de estado tenemos:

def pres(v,T):
    return 8*T/(3*v-1) - 3/v**2

def gr3(xgr,T):
    ygr3 = []
    for i2 in xgr:
        ygr3.append(pres(i2,T))
    return ygr3

xgr3 = linspace(0.5,4,1000)


ygr3a = gr3(xgr3,0.8)
ygr3b = gr3(xgr3,0.9)
ygr3c = gr3(xgr3,1)
ygr3d = gr3(xgr3,1.1)
ygr3e = gr3(xgr3,1.2)
ygr3f = gr3(xgr3,1.3)
    
plt.plot(xgr3,ygr3a,color = 'red', label = 'T = 0.8')
plt.plot(xgr3,ygr3b,color = 'blue', label = 'T = 0.9')
plt.plot(xgr3,ygr3c,color = 'orange', label = 'T = 1')
plt.plot(xgr3,ygr3d,color = 'green', label = 'T = 1.1')
plt.plot(xgr3,ygr3e,color = 'black', label = 'T = 1.2')
plt.plot(xgr3,ygr3f,color = 'purple', label = 'T = 1.3')
plt.grid()
plt.legend()
plt.xlabel('Volumen')
plt.xlim(0.5,4)
plt.ylim(0,4)
plt.ylabel('Presion')
plt.title('Presion vs Volumen (Ec. Van der Waals)')
plt.show()





xgr33 = linspace(0.5,8,10000)
ygr33 = gr3(xgr33,0.9)

#hallamos los puntos críticos
xc1 = 0.718597188
yc1 = pres(0.718597188,0.9)


xc2 = 1.528504964
yc2 = pres(1.528504964,0.9)


xf = 7
yf = pres(7,0.9)


plt.plot(xgr33,ygr33,color = 'blue', label = 'T = 0.9')
plt.plot(xc1,yc1, 'o', color = 'green',label = 'P min')
plt.plot(xc2,yc2, 'o', color = 'black', label = 'P max')
plt.plot(xf,yf,'o',color = 'red',label = 'Punto fijo escogido \nV = 7')
plt.plot(0.55,pres(0.55,0.9),'*',color = 'black',label = 'Pcalc')
plt.grid()
plt.legend()
plt.xlabel('Volumen')
plt.xlim(0.5,7.5)
plt.ylim(0,1.7)
plt.ylabel('Presion')
plt.title('Presion vs Volumen (Ec. Van der Waals)')
plt.show()




#----------------------------
# Para integrar con la cuadratura de Gauss-Legendre:
#----------------------------

GLr1 = zeros([1,2],float64)
GLr1[0,0:2] = [0, 2]

GLr2 = zeros([2,2],float64)
GLr2[0,0:2] = [sqrt(1/3),1]
GLr2[1,0:2]= [-sqrt(1/3),1]

GLr3 = zeros([3,2],float64)
GLr3[0,0:2] = [0,8/9]
GLr3[1,0:2] = [sqrt(3/5),5/9]
GLr3[2,0:2]= [-sqrt(3/5),5/9]

GLr4 = zeros([4,2],float64)
GLr4[0,0:2] = [sqrt(3/7 - 2/7 * sqrt(6/5)), (18+sqrt(30))/36]
GLr4[1,0:2]= [-sqrt(3/7 - 2/7 * sqrt(6/5)), (18+sqrt(30))/36]
GLr4[2,0:2] = [sqrt(3/7 + 2/7 * sqrt(6/5)), (18-sqrt(30))/36]
GLr4[3,0:2]= [-sqrt(3/7 + 2/7 * sqrt(6/5)), (18-sqrt(30))/36]

GLr5 = zeros([5,2],float64)
GLr5[0,0:2] = [0,128/225]
GLr5[1,0:2] = [1/3*sqrt(5-2*sqrt(10/7)), (322+13*sqrt(70))/900]
GLr5[2,0:2]= [-1/3*sqrt(5-2*sqrt(10/7)), (322+13*sqrt(70))/900]
GLr5[3,0:2] = [1/3*sqrt(5+2*sqrt(10/7)), (322-13*sqrt(70))/900]
GLr5[4,0:2]= [-1/3*sqrt(5+2*sqrt(10/7)), (322-13*sqrt(70))/900]



def cGL(ff,N,c,lim1,lim2):
    #limite N = 5
    if N == 1:
        GLr = GLr1
    
    if N == 2:
        GLr = GLr2
    
    if N == 3:
        GLr = GLr3
    
    if N == 4:
        GLr = GLr4
    
    if N == 5:
        GLr = GLr5

    
    if c == 1:
        #print(GLr)
        sum1 = 0
        for i1 in range(N):
            a = ff(GLr[i1,0])
            b = GLr[i1,1]
            sum1 = sum1 + a*b
        res = sum1
    else:
        if c == 2:
            sum2 = 0
            for i2 in range(N):
                a1 = ff((lim2-lim1)*(GLr[i2,0])/2 + (lim1+lim2)/2)
                b1 = GLr[i2,1]
                sum2 = sum2 + a1*b1
            res = (lim2-lim1)*sum2/2
        else:
            print("Condicion para los límites de integracion mal definidos")
            print("c = 1 -> Límites de integracion -1 y 1")
            print("c = 2 -> Límites de integracion deben ser especificados")
            print("         como lim1 y lim2 (lim1<lim2)")

    return res



def prese(v):
    return 8*0.9/(3*v-1) - 3/v**2

def calcg(v):
    aux0 = prese(v)
    aux1 = cGL(prese,3,2,v,7)
    aux2 = aux1 - 0.2987755102040816*(7-v)
    aux3 = aux2 + (aux0-0.2987755102040816)*v
    return aux3

#volumen
vgrf = linspace(0.55,7,1000)
pgrf = []
for i in vgrf:
    pgrf.append(prese(i))

plt.plot(pgrf,calcg(vgrf),color = 'red', label = 'g vs p')
plt.plot(yc1,calcg(xc1),'o',color = 'blue',label = 'P min')
plt.plot(yc2,calcg(xc2),'o',color = 'green',label = 'P max')
plt.plot(prese(7),0,'x',color = 'black',label = 'P referencia')
plt.plot(prese(0.55),calcg(0.55),'*',color = 'black', label = 'P calc' )
plt.title('g vs p\n(Referencia V = 7)')
plt.xlabel('Presion')
plt.ylabel('g')
plt.legend()
plt.grid()
plt.show()





