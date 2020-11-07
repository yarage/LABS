#! /usr/bin/python3

#Librerías usadas:

from numpy import *
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import random


###########################################################
#######   SEGUNDA PREGUNTA:      ##########################
###########################################################

print(" ---------------------")
print("   Segunda pregunta: ")
print(" ---------------------")

#Funciones a utilizar posteriormente:

def reslists(A,B):
    #Función para restar dos listas
    C = []
    for i in range(len(A)): 
        C.append(A[i]-B[i])
    return C
def clist1(V,cant):
    A = []
    for i in range(cant):
        A.append(V)
    return A

    
def desv(DT,CD,prom):
    # DT   : Veces que salió el mimsmo número de caras en nuestras 10k simulaciones
    # CD   : Cantidades de caras obtenidas en cada simulacion
    # prom : Promedio de caras en las 10k simulaciones

    suma1 = 0
    for i in range(len(CD)):
        aux = DT[i]*(CD[i] - prom)**2

        suma1 = suma1 + aux
    dest = sqrt(suma1/sum(DT))
    return dest

### Desarrollo del problema: ###


N = [10,100,1000]

# Las siguientes 3 listas almacenan los datos generados.
# Guardan el número de caras generadas en cada simulacion.
ctt1 = []
ctt2 = []
ctt3 = []

aux1 = 0

csimul = 10000

for ii in N:
    ct_s = ii
    aux1 = aux1 + 1
    
    for iii in range(csimul):
        d_c = []
        for i1 in range(ct_s):
            n = random.randint(0,1)#consideraremos 0 -> cara; 1 -> sello
            d_c.append(n)
        #print(d_c)
        cant1 = 0
        for i2 in d_c:
            if i2 == 0:
                cant1 = cant1 + 1
        
        if aux1 == 1:
            ctt1.append(cant1)
        if aux1 == 2:
            ctt2.append(cant1)
        if aux1 == 3:
            ctt3.append(cant1)
           
#print(ctt1)#de 10
#print(ctt2)#de 100
#print(ctt3)#de 1000

hct1 = zeros([11],int)
hct11 = range(11)

hct2 = zeros([101],int)
hct22 = range(101)

hct3 = zeros([1001],int)
hct33 = range(1001)

# Los siguientes bucles, ordenan los datos, contando cuantas veces
# se repite cada nro. de caras en el total de simulaciones
# (Haciendo esto, tenemos listos los datos a graficar)

for i3 in hct11:
    for i4 in ctt1:
        if i4 == i3:
            hct1[i3] = hct1[i3] + 1
hct1 = hct1/csimul
#print(hct1)

for i3 in hct22:
    for i4 in ctt2:
        if i4 == i3:
            hct2[i3] = hct2[i3] + 1
hct2 = hct2/csimul
#print(hct2)

for i3 in hct33:
    for i4 in ctt3:
        if i4 == i3:
            hct3[i3] = hct3[i3] + 1
hct3 = hct3/csimul
#print(hct3)


#Graficando los histogramas uno a uno:

fig = plt.figure(u'Gráfica de barras') # Figure
ax = fig.add_subplot(111) # Axes

ax.bar(range(11), hct1, width=0.8, align='center')
ax.set_xticks(range(11))
ax.set_xticklabels(range(11))
ax.set_title('Histograma \n(10Mil simulaciones con 10 monedas)')
ax.set_ylabel('Frecuencia relativa')
ax.set_xlabel('Número de caras')
ax.grid()
plt.show()



fig = plt.figure(u'Gráfica de barras') # Figure
ax = fig.add_subplot(111) # Axes

ax.bar(range(101), hct2, width=0.8, align='center')
ax.set_xticks(range(0,101,5))
ax.set_xticklabels(range(0,101,5))
ax.set_title('Histograma \n(10Mil simulaciones con 100 monedas)')
ax.set_ylabel('Frecuencia relativa')
ax.set_xlabel('Número de caras')
ax.grid()
plt.show()



fig = plt.figure(u'Gráfica de barras') # Figure
ax = fig.add_subplot(111) # Axes

ax.bar(range(1001), hct3, width=0.8, align='center')
ax.set_xticks(range(0,1001,50))
ax.set_xticklabels(range(0,1001,50))
ax.set_title('Histograma \n(10Mil simulaciones con 1000 monedas)')
ax.set_ylabel('Frecuencia relativa')
ax.set_xlabel('Número de caras')
ax.grid()
plt.show()



#Ahora, lo que sigue, será para calcular el promedio de caras:

ccar1 = hct1*csimul
ccar2 = hct2*csimul
ccar3 = hct3*csimul

#Entonces ahora, los promedios:

# Con 10 monedas, simulando 'csimul' veces:
prom1 = sum(ccar1*hct11)/csimul

# Con 100 monedas, simulando 'csimul' veces:
prom2 = sum(ccar2*hct22)/csimul

# Con 1000 monedas, simulando 'csimul' veces:
prom3 = sum(ccar3*hct33)/csimul

print("\n---------------------------\n")

print(" Promedio del nro. de caras usando [10, 100, 1000] monedas \n simulando diez mil veces por caso :\n")
print(prom1,prom2,prom3)

print("\n---------------------------\n")

#Ahora para calcular la desviación estándar:

dves1 = desv(ccar1,hct11,prom1)
print('Desviacion Estandar con 10 monedas : ',dves1)
print("Fluctuaciones relativas            : ", dves1/prom1)

print("\n---------------------------\n")

dves2 = desv(ccar2,hct22,prom2)
print('Desviacion Estandar con 100 monedas : ',dves2)
print("Fluctuaciones relativas             : ", dves2/prom2)

print("\n---------------------------\n")

dves3 = desv(ccar3,hct33,prom3)
print('Desviacion Estandar con 1000 monedas : ',dves3)
print("Fluctuaciones relativas              : ", dves3/prom3)

print("\n---------------------------\n")
