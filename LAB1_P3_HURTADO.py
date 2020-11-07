#! /usr/bin/python3

from numpy import *
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import random


###########################################################
#######   TERCERA PREGUNTA:      ##########################
###########################################################


print(" ---------------------")
print("   Tercera pregunta: ")
print(" ---------------------")

#generacion de data:
N = [100,1000,100000]
ctt1 = []
ctt2 = []
ctt3 = []

for ii in N:
    print("\n-> Usando N = ",ii)
    ct_s = ii
    d_c = []
    for i1 in range(0,ct_s):
        n = random.randint(1,6)#números aleatorios del 1-6
        d_c.append(n)
    #Si queremos ver los datos generados:    print(d_c)

    ######## Parte 1: A = {salir número impar}
    cant1 = 0
    for i2 in range(ct_s):
        if d_c[i2] == 1 or d_c[i2] == 3 or d_c[i2] == 5:
            cant1 = cant1 + 1
    print("La cantidad de impares es: ", cant1)
    ctt1.append(cant1/ii)

    ######## Parte 2: B = {salir un número inferior a 4}
    cant2 = 0
    for i3 in range(ct_s):
        if d_c[i3] == 1 or d_c[i3] == 2 or d_c[i3] == 3:
            cant2 = cant2 + 1
    print("La cantidad de menores a 4 es: ", cant2)
    ctt2.append(cant2/ii)

    ######### Parte 3: A ∩ B
    cant3 = 0
    for i4 in range(ct_s):
        if d_c[i4] == 5:
            cant3 = cant3 + 1
    print("Suceso A intersección (no B): ", cant3)
    ctt3.append(cant3/ii)

    
    print("\n")



print("Mostraremos los resultados en [100,1000,100000] lanzamientos:")
print("(a) La frecuencia relativa de A: ",ctt1)
print("(b) La frecuencia relativa de B: ",ctt2)
print("(c) La frecuencia relativa de A ∩ (no B): ",ctt3)

print("\n")

print(" Para visualizar mejor el comportamiento:")
print(" Simularemos [100 500 900 ... 99300 99700] lanzamientos")

N = range(100,100000,400)

ctt1 = []
ctt2 = []
ctt3 = []

for ii in N:
    ct_s = ii
    d_c = []
    for i1 in range(0,ct_s):
        n = random.randint(1,6)#números aleatorios del 1-6
        d_c.append(n)
    #Si queremos ver los datos generados:    print(d_c)

    ######## Parte 1: A = {salir número impar}
    cant1 = 0
    for i2 in range(ct_s):
        if d_c[i2] == 1 or d_c[i2] == 3 or d_c[i2] == 5:
            cant1 = cant1 + 1
    #print("La cantidad de impares es: ", cant1)
    ctt1.append(cant1/ii)

    ######## Parte 2: B = {salir un número inferior a 4}
    cant2 = 0
    for i3 in range(ct_s):
        if d_c[i3] == 1 or d_c[i3] == 2 or d_c[i3] == 3:
            cant2 = cant2 + 1
    #print("La cantidad de menores a 4 es: ", cant2)
    ctt2.append(cant2/ii)

    ######### Parte 3: A ∩ B
    cant3 = 0
    for i4 in range(ct_s):
        if d_c[i4] == 5:
            cant3 = cant3 + 1
    #print("Suceso A intersección (no B): ", cant3)
    ctt3.append(cant3/ii)

plt.plot(N,ctt2, color = 'red',label = 'Frecuencia relativa de B')
plt.plot(N,ctt1, color = 'orange',label = 'Frecuencia relativa de A')
plt.plot(N,ctt3, color = 'blue',label = 'Frecuencia relativa de A ∩ (no B)')

plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.axhline(y=0.5,linestyle = 'dashed',color = 'gray',label='Frec. relativa = 0.5')
plt.axhline(y=1/6,linestyle = 'dashed',color = 'black',label='Frec. relativa = 1/6 = 0.1666...')

plt.grid()
plt.legend()
plt.title('Convergencia de las frecuencias relativas')
plt.xlabel('N lanzamientos - Simulaciones')
plt.ylabel('Frecuencia relativa')
plt.ylim(0,1)
plt.show()

print("\n")
