#! /usr/bin/python3

#Librerías usadas:

from numpy import *
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import random


###########################################################
#######   PRIMERA PREGUNTA:      ##########################
###########################################################

print(" ---------------------")
print("   Primera pregunta: ")
print(" ---------------------")

ct_s = 1000 #Cantidad de simulaciones:

aps = [] # Lista donde guardo los datos de las simulaciones
         # 1 -> Dos o más personas cumplen años el mismo día 
         # 0 -> Ninguno cumple el mismo día que otro

for ii in range(0,ct_s):
    d_c = []
    for i1 in range(50):
        n = random.randint(1,365) #Generación de números aleatorios
        d_c.append(n)
    
    ct1 = 0
    for i2 in range(0,49):
        for i3 in range(i2+1,50):
            if d_c[i2] == d_c[i3]:
                ct1 = ct1 + 1

    if ct1 > 0:
        aps.append(1)

    else:
        aps.append(0)
    
print(aps)

si = 0 # si: cantidad de simulaciones en las que si repiten fecha de cumpleaños
for i4 in range(ct_s):
    if aps[i4] == 1:
        si = si + 1

print("\n ------------------------------------------")
print("  Se repiten    : ", si)
print("  No se repiten : ", ct_s - si)
print(" ------------------------------------------")    

print("\n Entonces la probabilidad de ganar la apuesta en")
print(" porcentaje es aproximadamente: ",(ct_s - si)/ct_s *100)


