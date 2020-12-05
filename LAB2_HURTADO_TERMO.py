#! /usr/bin/python3

from numpy import *
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import random


def fact(n):

    fact = 1

    if n == 0:
       fact = 1
    else:
       for i in range(1,n + 1):
           fact = fact*i
    return fact



######
#pregunta 1 a - completed



B=[]


N = 10**23
T = 300
k = 1.380649 * 10**(-23)
p_E = 3/2 * N * k * T


print("Por la forma de la distribución, podemos deducir que será una gaussiana ")
print("centrada en p_E, y este valor es:\n")
print(p_E)


sigma = sqrt(2/(3*N))*p_E

def f1(E,N,p_E,sigma):
    d_E = E - p_E
    R = 1/(sqrt(2*pi)*sigma)*exp(-(d_E)**(2)/(2*sigma**2))
    return R


A = linspace(621.292,621.2921,10000)



for i1 in A:
    aux1 = f1(i1,N,p_E,sigma)
    B.append(aux1)
#print(B)

plt.plot(A,B, color = 'orange')
plt.title('Densidad de Probabilidad')
plt.grid()
plt.show()


#################################################
#pregunta 2 finished:
#################################################



k = 1.380649 * 10**(-23)
#k=1
tm1 = 60

SAA = [1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]

SAev = zeros([tm1,24],float64)

s1t = []
s2t = []
sst = []

tm_gr = range(1,tm1 + 1)


print("\n El cuadro con la evolución paso a paso se muestra a continuación:\n")

for i1 in range(tm1):
    SAev[i1,0:24] = SAA[0:24]
    print(SAA)
    n1 = random.randint(0,23) 
    n2 = n1
    while n2 == n1 or SAA[n2] == SAA[n1]:
        n2 = random.randint(0,23)

    
    SAA[n1] = -SAA[n1]
    SAA[n2] = -SAA[n2]

    

    #La entropía en cada 'tiempo':
    # S = k*ln(omega)

    #sistema A:
    up1 = 0
    dw1 = 0
    for i2 in range(8):
        if SAA[i2] == 1:
            up1 = up1 + 1
        else:
            dw1 = dw1 + 1

    s1 = k*log(fact(8)/(fact(up1)*fact(dw1)))
    s1t.append(s1)

    
    #sistema A':
    up2 = 0
    dw2 = 0
    for i3 in range(8,24):
        if SAA[i3] == 1:
            up2 = up2 + 1
        else:
            dw2 = dw2 + 1

    s2 = k*log(fact(16)/(fact(up2)*fact(dw2)))
    s2t.append(s2)

    #S total:
    ss = s1+s2
    sst.append(ss)

#print(SAA)

SA = SAA[0:8]
SAp = SAA[8:24]
print("\n Sistema A:")
print(SA)
print("\n Sistema A':")
print(SAp)

print("\n Entropia ( A / A' / total )")
print(s1,s2,ss)


plt.plot(tm_gr,sst,linestyle='--',marker='.',color = 'blue',label = "Entropia de A'+A")
plt.title("Stotal vs 'tiempo' ")
plt.xlabel(' "tiempo" ')
plt.grid()
plt.legend()
plt.show()

plt.plot(tm_gr,s1t,linestyle='--',marker='.',color = 'orange',label = "Entropia de A")
plt.title("S1 vs 'tiempo' ")
plt.xlabel(' "tiempo" ')
plt.grid()
plt.legend()
plt.show()

plt.plot(tm_gr,s2t,linestyle='--',marker='.',color = 'green',label = "Entropia de A' ")
plt.title("S2 vs 'tiempo' ")
plt.xlabel(' "tiempo" ')
plt.grid()
plt.legend()
plt.show()




###############################################
# PREGUNTA 3: completed
###############################################


qA = [0,1,2,3]
omA = [1,3,6,10]
omB = [10,6,3,1]
omT = [10,18,18,10]


plt.plot(qA,omA,linestyle='--',marker='o',color = 'blue',label = "Ω_A")
plt.plot(qA,omB,linestyle='--',marker='o',color = 'orange',label = "Ω_B")
plt.plot(qA,omT,linestyle='--',marker='o',color = 'red',label = "Ω_Total")
plt.grid()
plt.xticks([0,1,2,3])
plt.yticks([1,3,6,10,18])
plt.legend()
plt.title('Ω vs q_A')
plt.xlabel('q_A')
plt.ylabel('Ω')
plt.show()








