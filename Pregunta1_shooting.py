"""
METODO SHOOTING LINEAL
PREGUNTA 1
"""
# RESOLVER EDO de la forma: y''=p(x)*y'+q(x)*y+r(x)
import numpy as np
import math
import matplotlib.pyplot as plt
#Definimos los valores iniciales para resolver la EDO con valores de frontera:
a=0
b=1
N=1000
alpha=1#y(a)=alpha
beta=-0.5#y(b)=beta
h=(b-a)/N #Este es el 'paso' con el que trabajaremos 
#La ecuacion diferencial a resolver es la siguiente
#     y''=f(x,y,y')=-2*x*y
def p(x):
    return 0
def q(x):
    return -2*x
def r(x):
    return 0
#Definimos los valores que usaremos para aplicar Runge Kutta de 4to orden:
u_1=np.zeros(N+1)
u_2=np.zeros(N+1)
v_1=np.zeros(N+1)
v_2=np.zeros(N+1)
#CONDICIONES INICIALES DE LA EDO
u_1[0]=alpha
u_2[0]=0
v_1[0]=0
v_2[0]=1
for i in range(N):
    #Ahora resolvemos usando RK4
    x=a+i*h
    k1_1=h*u_2[i]
    k1_2=h*(p(x)*u_2[i]+q(x)*u_1[i]+r(x))
    k2_1=h*(u_2[i]+k1_2/2)
    k2_2=h*(p(x+h/2)*(u_2[i])+q(x+h/2)*(u_1[i]+k1_1/2)+r(x+h/2))
    k3_1=h*(u_2[i]+k2_2/2)
    k3_2=h*(p(x+h/2)*(u_2[i]+k2_2/2)+q(x+h/2)*(u_1[i]+k2_1/2)+r(x+h/2))
    k4_1=h*(u_2[i]+k3_2)
    k4_2=h*(p(x+h)*(u_2[i]+k3_2)+q(x+h)*(u_1[i]+k3_1)+r(x+h))
    u_1[i+1]=u_1[i]+(k1_1+2*k2_1+2*k3_1+k4_1)/6
    v_2[i+1]=v_2[i]+(k1_2+2*k2_2+2*k3_2+k4_2)/6
    #La segunda parte de la iteración es correspondiente a v1 y v2:
    k1_1_prima=h*v_2[i]
    k1_2_prima=h*(p(x)*v_2[i]+q(x)*v_1[i])
    k2_1_prima=h*(v_2[i]+k1_2_prima/2)
    k2_2_prima=h*(p(x+h/2)*(v_2[i]+k1_2_prima/2)+q(x+h/2)*(v_1[i]+k1_1_prima/2))
    k3_1_prima=h*(v_2[i]+k2_2_prima/2)
    k3_2_prima=h*(p(x+h/2)*(v_2[i]+k2_2_prima/2)+q(x+h/2)*(v_1[i]+k2_1_prima/2))
    k4_1_prima=h*(v_2[i]+k3_2_prima)
    k4_2_prima=h*(p(x+h)*(v_2[i]+k3_2_prima)+q(x+h)*(v_1[i]+k3_1_prima))
    v_1[i+1]=v_1[i]+(k1_1_prima+2*k2_1_prima+2*k3_1_prima+k4_1_prima)
    v_2[i+1]=v_2[i]+(k1_2_prima+2*k2_2_prima+2*k3_2_prima+k4_2_prima)    
#Definimos nuevas variables para trabajar:
w1=np.zeros(N+1)
w2=np.zeros(N+1)
w1[0]=alpha    
w2[0]=(beta-u_1[N])/v_1[N]    
#Empezamos la nueva iteración para hallar w1  
for i in range(N):    
    w1[i+1]=u_1[i]+w2[0]*v_1[i]
    w2[i+1]=u_2[i]+w2[0]*v_2[i]
    x=a+i*h
w1=np.append(w1,b)#Agregamos el ultimo valor del intervalo
print('Tenemos que para x=a+i*h con h={} la solucion w1 es:\n{}'.format(h,w1))

#Ahora pasamos a graficar la solucion con el paso indicado
X_dom=np.zeros(N+1)#En esta matriz contendrá el dominio de X
X_dom[0]=a
for i in range(N):
    X_dom[i+1]=X_dom[i]+h
plt.plot(X_dom,w2)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Grafica de la Solucion Numerica-Shooting')












