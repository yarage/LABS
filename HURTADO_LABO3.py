#! /usr/bin/python3

from numpy import *

import matplotlib.pyplot as plt
from pylab import plot,show

import warnings
warnings.filterwarnings("ignore")



def hmac10(niv,Et,npar):
    cn = len(niv)#cantidad de niveles de energía
    # este programa solo funciona cuando los niveles de energía son 10

    #reordenemos de menor a mayor los niveles
    niv.sort()

    ct1 = []


    ct1 = []
    for i1 in range(cn):
        if niv[i1] != 0:
            ct1.append(Et//niv[i1]+1)
        else:
            ct1.append(0)
            aux1 = i1
    ct1[aux1] = npar+1
    

    c1 = []
    c2 = []
    c3 = []
    c4 = []
    c5 = []
    c6 = []
    c7 = []
    c8 = []
    c9 = []
    c10 = []


    for i2 in range(ct1[0]):
        for i3 in range(ct1[1]):
            for i4 in range(ct1[2]):
                for i5 in range(ct1[3]):
                    for i6 in range(ct1[4]):
                        for i7 in range(ct1[5]):
                            for i8 in range(ct1[6]):
                                for i9 in range(ct1[7]):
                                    for i10 in range(ct1[8]):
                                        for i11 in range(ct1[9]):
                                            aux2 = i2*niv[0] + i3*niv[1] + i4*niv[2] + i5*niv[3] + i6*niv[4] + i7*niv[5] + i8*niv[6] + i9*niv[7] + i10*niv[8] + i11*niv[9]
                                            aux3 = i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10 + i11
                                            if aux2 == Et and aux3 == npar:
                                                c1.append(i2)
                                                c2.append(i3)
                                                c3.append(i4)
                                                c4.append(i5)
                                                c5.append(i6)
                                                c6.append(i7)
                                                c7.append(i8)
                                                c8.append(i9)
                                                c9.append(i10)
                                                c10.append(i11)
    #c1 es de E = 0e
    #c10 es de E = 9e
    aux4 = len(c1)

    Mstes = zeros([10,aux4],int)

    Mstes[0,0:aux4] = c1
    Mstes[1,0:aux4] = c2
    Mstes[2,0:aux4] = c3
    Mstes[3,0:aux4] = c4
    Mstes[4,0:aux4] = c5
    Mstes[5,0:aux4] = c6
    Mstes[6,0:aux4] = c7
    Mstes[7,0:aux4] = c8
    Mstes[8,0:aux4] = c9
    Mstes[9,0:aux4] = c10


    Mstes = transpose(Mstes)
    
    return Mstes


def ft(n):
    #factorial
    fact = 1

    if n == 0:
       fact = 1
    else:
       for i in range(1,n + 1):
           fact = fact*i
    return fact


def nummic(nlvl,MAC,npar):
    c = nlvl
    cc = len(MAC)
    nm = zeros([cc,1],int)
    for ix in range(cc):#fila
        ax1 = 1
        for iy in range(c):#columna
            ax1 = ax1*ft(MAC[ix,iy])
        nm[ix] = ft(npar)/ax1
    
    return nm



#Para ordenar de menor a mayor la cantidad de microestados
def or1(A):
    ctt = len(A)
    C = zeros([ctt,1],int)
    for i in range(ctt):
        C[i,0] = A[i,0]
    B = transpose(C)
    B.sort()
    B = transpose(B)
    return B

#Para ordenar y distribuir alrededor de los valores más altos de una lista
#inicialmente ordenada de menor a mayor
def or2(M):
    tm = len(M)
    ax1 = []
    ax2 = []
    p1 = 0
    while p1 < tm:
        if p1/2 == p1//2:
            ax1.append(p1)#pares
        else:
            ax2.append(p1)#impares
        p1 = p1 + 1
    #invertimos la lista de los impares
    ax2.sort(reverse=True)

    for i in ax2:
        ax1.append(i)

    N = []
    for i1 in ax1:
        N.append(M[i1,0])
    return N




def beta10(om,E,b,N):
    return log(om)-9*E/b-N*log(1+exp(-b*E)+exp(-2*b*E)+exp(-3*b*E)+exp(-4*b*E)+exp(-5*b*E)+exp(-6*b*E)+exp(-7*b*E)+exp(-8*b*E)+exp(-9*b*E))


def beta3(om,E,b,N):
    return log(om)-6*E/b-N*log(1+exp(-b*E)+exp(-2*b*E))


#Metodo de la falsa posicion para encontrar la raiz de una ecuacion de forma numerica
def fposi(f,om,E,N,xl,xu,er):
    #print("----------------------------------------------------------")
    #print("=========      Método de la Falsa Posicion      ==========")
    #print("----------------------------------------------------------")

    aux1 = 0
    aux2 = 0
    it = 0
    imax = 5000
    A = zeros(imax-1,float)
    xr = xu - (f(om,E,xu,N)*(xl-xu)/(f(om,E,xl,N)-f(om,E,xu,N)))
    ea = er + 1 #para que se cumpla al inicio la condicion del while
    esc = 1
    #print("Raiz / iteracion / error")
    while ea >= er and it < imax and esc == 1:
        it = it + 1
        xrold = xr
        xr = xu - (f(om,E,xu,N)*(xl-xu)/(f(om,E,xl,N)-f(om,E,xu,N)))
        if xr !=  0 and xr != xrold:
            ea = abs((xr - xrold)/xr)*100
        test = f(om,E,xl,N)*f(om,E,xr,N)
        if test < 0:
            aux1 = xr
            xu = aux1
        else:
            if test > 0:
                aux2 = xr
                xl = aux2
            else:
                print("the root is : ",xr)
                ea = 0
                esc = 2
        A[it-1] = ea
        #print("{0:.6f}".format(xr),it,ea)
    B = zeros(it,float)
    for i in range(it):
        B[i] = A[i]
    C = range(1,it+1)
    #plt.yscale('log')

    #quitar de comentario si se desea ver la evolucion del error
    #plot(C, B, color = 'green', label = 'Falsa posicion')
        

    print("Beta : ","{0:.6f}".format(xr), "  /  Error: ", ea)  


    return xr


#ahora calculare la probabilidad de que una particula tenga cierta energia
def prob10(k,M,lvl,N):#k : macroestados ; M : nro micro ; lvl : niveles energia
    # N : cantidad de particulas
    nt = sum(M)
    t1 = len(M)
    pr = zeros([10])
    for i1 in range(10):
        p = 0
        for i2 in range(t1):
            p = p + k[i2,i1]/N * M[i2,0]/nt
        pr[i1] = p
    return pr






########################################################################
########################################################################


print("\n Pregunta 1")
print("\n")

#niveles de energia
nvels = [0,1,2,3,4,5,6,7,8,9]

print("Los niveles de energia en esta pregunta son:")
print(nvels)


epsilon = 0.1

########################################################################
########################################################################


print("\n")
print("-> Con 6 partículas")
k1 = hmac10(nvels,9,6)

#Para hallar microestados
M1 = nummic(10,k1,6)

print("\n Macroestados* ")
print("\n")
print(k1)
print("\n Microestados(omega) de cada macroestado*")
print("\n")
print(M1)

M11 = or1(M1)
of1 = or2(M11)
x1 = range(1,len(of1)+1)

#numero total de microestados:
nt1 = sum(M11)


plt.plot(x1,of1,color = 'orange',label = '6 particulas')
plt.plot(x1,of1,'x',color = 'red')
plt.grid()
plt.legend()
plt.xlabel('N')
plt.ylabel('Numero de microestados (Omega)')
plt.show()


print("\n Ahora para calcular la temperatura:")



#Primero estimamos en que intervalo cae el resultado, es decir, el punto donde corta el eje x
x=linspace(0,10,100)
plt.plot(x,beta10(sum(M1),epsilon,x,6),label = 'Ecuacion a resolver')
plt.grid()
plt.title('Grafica de la funcion auxiliar para encontrar β')
plt.xlabel("β")
plt.ylabel("f(Ɛ,β)")
plt.show()

bet1 = fposi(beta10,sum(M1),epsilon,6,0.2,10,1e-05)

print("\n Entonces, la temperatura es 1/(k_b*beta), entonces:")
print("T = ",1/bet1,"/ k_b")





########################################################################
########################################################################


print("\n")
print("-> Con 12 partículas")
k2 = hmac10(nvels,9,12)
M2 = nummic(10,k2,12)

print("\n Macroestados* ")
print("\n")
print(k2)
print("\n Microestados(omega) de cada macroestado*")
print("\n")
print(M2)

M22 = or1(M2)
of2 = or2(M22)
x2 = range(1,len(of2)+1)

plt.plot(x2,of2,color = 'red',label = '12 particulas')
plt.plot(x2,of2,'x',color = 'blue')
plt.grid()
plt.legend()
plt.xlabel('N')
plt.ylabel('Numero de microestados (Omega)')
plt.show()



print("\n Ahora para calcular la temperatura:")



#Primero estimamos en que intervalo cae el resultado, es decir, el punto donde corta el eje x
x=linspace(0,10,100)
plt.plot(x,beta10(sum(M2),epsilon,x,12),label = 'Ecuacion a resolver')
plt.grid()
plt.title('Grafica de la funcion auxiliar para encontrar β')
plt.xlabel("β")
plt.ylabel("f(Ɛ,β)")
plt.show()

bet2 = fposi(beta10,sum(M2),epsilon,12,0.2,2,1e-05)

print("\n Entonces, la temperatura es 1/(k_b*beta), entonces:")
print("T = ",1/bet2,"/ k_b")




########################################################################
########################################################################


print("\n")
print("-> Con 18 partículas")
k3 = hmac10(nvels,9,18)
M3 = nummic(10,k3,18)

print("\n Macroestados* ")
print("\n")
print(k3)
print("\n Microestados(omega) de cada macroestado*")
print("\n")
print(M3)

M33 = or1(M3)
of3 = or2(M33)
x3 = range(1,len(of3)+1)

plt.plot(x3,of3,color = 'blue',label = '18 particulas')
plt.plot(x3,of3,'x',color = 'black')
plt.grid()
plt.legend()
plt.xlabel('N')
plt.ylabel('Numero de microestados (Omega)')
plt.show()



print("\n Ahora para calcular la temperatura:")



#Primero estimamos en que intervalo cae el resultado, es decir, el punto donde corta el eje x
x=linspace(0,10,100)
plt.plot(x,beta10(sum(M3),epsilon,x,18),label = 'Ecuacion a resolver')
plt.grid()
plt.title('Grafica de la funcion auxiliar para encontrar β')
plt.xlabel("β")
plt.ylabel("f(Ɛ,β)")
plt.show()

bet3 = fposi(beta10,sum(M3),epsilon,18,0.2,2,1e-05)

print("\n Entonces, la temperatura es 1/(k_b*beta), entonces:")
print("T = ",1/bet3,"/ k_b")




#Probabilidad de que una particula tenga cierta energia

PR1 = prob10(k1,M1,nvels,6)
PR2 = prob10(k2,M2,nvels,12)
PR3 = prob10(k3,M3,nvels,18)


plt.plot(nvels,PR1,color ='orange',label='6 particulas')
plt.plot(nvels,PR2,color = 'red',label='12 particulas')
plt.plot(nvels,PR3,color = 'black',label='18 particulas')
plt.grid()
plt.xticks([0,1,2,3,4,5,6,7,8,9])
plt.legend()
plt.xlabel('Nivel de energía')
plt.ylabel('Probabilidad')
plt.show()

########################################################################
########################################################################



def hmac3(niv,Et,npar):
    cn = len(niv)#cantidad de niveles de energía
    # este programa solo funciona cuando los niveles de energía son 3

    #reordenemos de menor a mayor los niveles
    niv.sort()

    ct1 = [npar+1,npar+1,npar+1]


    c1 = []
    c2 = []
    c3 = []
    
    for i2 in range(ct1[0]):
        for i3 in range(ct1[1]):
            for i4 in range(ct1[2]):
                aux2 = i2*niv[0] + i3*niv[1] + i4*niv[2]
                #print(aux2)
                aux3 = i2 + i3 + i4 
                #print(aux3)
                if aux2 == Et and aux3 == npar:
                    c1.append(i2)
                    c2.append(i3)
                    c3.append(i4)

    #c1 es de E = 0e
    #c10 es de E = 2e
    aux4 = len(c1)

    Mstes = zeros([3,aux4],int)

    Mstes[0,0:aux4] = c1
    Mstes[1,0:aux4] = c2
    Mstes[2,0:aux4] = c3



    Mstes = transpose(Mstes)
    
    return Mstes


def temp(M):
    #Primero sumamos el numero de microestados
    dum = sum(M)
    


    
    return T


def prob3(k,M,lvl,N):#k : macroestados ; M : nro micro ; lvl : niveles energia
    # N : cantidad de particulas
    nt = sum(M)
    t1 = len(M)
    pr = zeros([3])
    for i1 in range(3):
        p = 0
        for i2 in range(t1):
            p = p + k[i2,i1]/N * M[i2,0]/nt
        pr[i1] = p
    return pr







print("\n Pregunta 2")
print("\n")

nvels2 = [0,1,2]

print("Los niveles de energia en esta pregunta son:")
print(nvels2)

epsilon1 = 0.1

#####################################################################
#####################################################################

print("\n")
print("-> Con 18 partículas")
k5 = hmac3(nvels2,6,18)
M5 = nummic(3,k5,18)


print("\n Macroestados* ")
print("\n")
print(k5)
print("\n Microestados(omega) de cada macroestado*")
print("\n")
print(M5)

M55 = or1(M5)
of5 = or2(M55)
x5 = range(1,len(of5)+1)


plt.plot(x5,of5,color = 'orange',label = '18 particulas')
plt.plot(x5,of5,'x',color = 'red')
plt.grid()
plt.legend()
plt.xlabel('N')
plt.ylabel('Numero de microestados (Omega)')
plt.show()





print("\n Ahora para calcular la temperatura:")

#Primero estimamos en que intervalo cae el resultado, es decir, el punto donde corta el eje x
x=linspace(0,30,1000)
plt.plot(x,beta3(sum(M5),epsilon1,x,18),label = 'Ecuacion a resolver')
plt.grid()
plt.title('Grafica de la funcion auxiliar para encontrar β')
plt.xlabel("β")
plt.ylabel("f(Ɛ,β)")
plt.show()

bet5 = fposi(beta3,sum(M5),epsilon1,18,0.2,1000,1e-05)

print("\n Entonces, la temperatura es 1/(k_b*beta), entonces:")
print("T = ",1/bet5,"/ k_b")





########################################################################
########################################################################







print("\n")
print("-> Con 36 partículas")
k6 = hmac3(nvels2,6,36)
M6 = nummic(3,k6,36)


print("\n Macroestados* ")
print("\n")
print(k6)
print("\n Microestados(omega) de cada macroestado*")
print("\n")
print(M6)

M66 = or1(M6)
of6 = or2(M66)
x6 = range(1,len(of6)+1)


plt.plot(x6,of6,color = 'orange',label = '36 particulas')
plt.plot(x6,of6,'x',color = 'red')
plt.grid()
plt.legend()
plt.xlabel('N')
plt.ylabel('Numero de microestados (Omega)')
plt.show()




print("\n Ahora para calcular la temperatura:")

#Primero estimamos en que intervalo cae el resultado, es decir, el punto donde corta el eje x
x=linspace(0,30,1000)
plt.plot(x,beta3(sum(M6),epsilon1,x,36),label = 'Ecuacion a resolver')
plt.grid()
plt.title('Grafica de la funcion auxiliar para encontrar β')
plt.xlabel("β")
plt.ylabel("f(Ɛ,β)")
plt.show()

bet6 = fposi(beta3,sum(M6),epsilon1,36,0.2,1500,1e-05)

print("\n Entonces, la temperatura es 1/(k_b*beta), entonces:")
print("T = ",1/bet6,"/ k_b")





########################################################################
########################################################################



print("\n")
print("-> Con 54 partículas")
k7 = hmac3(nvels2,6,54)
M7 = nummic(3,k7,54)

print("\n Macroestados* ")
print("\n")
print(k7)
print("\n Microestados(omega) de cada macroestado*")
print("\n")
print(M7)

M77 = or1(M7)
of7 = or2(M77)
x7 = range(1,len(of7)+1)


plt.plot(x7,of7,color = 'orange',label = '54 particulas')
plt.plot(x7,of7,'x',color = 'red')
plt.grid()
plt.legend()
plt.xlabel('N')
plt.ylabel('Numero de microestados (Omega)')
plt.show()





print("\n Ahora para calcular la temperatura:")

#Primero estimamos en que intervalo cae el resultado, es decir, el punto donde corta el eje x
x=linspace(0,30,1000)
plt.plot(x,beta3(sum(M7),epsilon1,x,54),label = 'Ecuacion a resolver')
plt.grid()
plt.title('Grafica de la funcion auxiliar para encontrar β')
plt.xlabel("β")
plt.ylabel("f(Ɛ,β)")
plt.show()

bet7 = fposi(beta3,sum(M7),epsilon1,54,0.2,1500,1e-05)

print("\n Entonces, la temperatura es 1/(k_b*beta), entonces:")
print("T = ",1/bet7,"/ k_b")



#Probabilidad de que una particula se encuentre en cierto nivel
# de energia

PR5 = prob3(k5,M5,nvels2,18)
PR6 = prob3(k6,M6,nvels2,36)
PR7 = prob3(k7,M7,nvels2,54)

plt.plot(nvels2,PR5,color = 'orange',label='18 particulas')
plt.plot(nvels2,PR7,color = 'red',label='36 particulas')
plt.plot(nvels2,PR6,color = 'black',label='54 particulas')
plt.grid()
plt.xticks([0,1,2])
plt.legend()
plt.xlabel('Nivel de energía')
plt.ylabel('Probabilidad')
plt.show()


