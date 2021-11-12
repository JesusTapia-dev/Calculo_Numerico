^{"""
Metodo de disparo no lineal
PREGUNTA 2
"""
# RESOLVER EDO de la forma: y''=f(x,y,y')=-2*x*y^2
import numpy as np
import matplotlib.pyplot as plt

#Definimos los valores iniciales para resolver la EDO con valores de frontera:
a=0
b=1
alpha=1
beta=1.5
N=200
M=100#Maximo numero de iteraciones
TOL=10**-5#Indica la tolerancia para deterne el
#Aproximamos w1_i a y(xi) and w2_i a y'(xi)
#Definimos la funcio f que usaremos
def f(x,y,y_prima):
    return -2*x*y**2
def fy(x,y,y_prima):#Representa la derivada parcial con respecto a y
    return -4*x*y
def fy_prima(x,y,y_prima):#Representa la derivada parcial con respecto a y'
    return 0
#Definimos:
h=(b-a)/N
k=1
TK=(beta-alpha)/(b-a)
w1=np.zeros(N+1)
w2=np.zeros(N+1)

while k<=M:
    w1[0]=alpha
    w2[0]=TK
    u1=0
    u2=1
    for i in range(N):#Metodo de RK4
        i=i+1
        x=a+(i-1)*h
        k1_1=h*w2[i-1]
        k1_2=h*f(x,w1[i-1],w2[i-1])
        k2_1=h*(w2[i-1]+k1_2/2)
        k2_2=h*f(x+h/2,w1[i-1]+k1_1/2,w2[i-1]+k1_2/2)
        k3_1=h*(w2[i-1]+k2_2/2)
        k3_2=h*f(x+h/2,w1[i-1]+k2_1/2,w2[i-1]+k2_2/2)
        k4_1=h*(w2[i-1]+k3_2)
        k4_2=h*f(x+h,w1[i-1]+k3_1,w2[i-1]+k3_2)
        w1[i]=w1[i-1]+(k1_1+2*k2_1+2*k3_1+k4_1)/6
        w2[i]=w2[i-1]+(k1_2+2*k2_2+2*k3_2+k4_2)/6
        k1_1_prima=h*u2
        k1_2_prima=h*(fy(x,w1[i-1],w2[i-1])*u1++fy_prima(x,w1[i-1],w2[i-1])*u2)
        k2_1_prima=h*(u2+k1_2_prima/2)
        k2_2_prima=h*(fy(x+h/2,w1[i-1],w2[i-1])*(u1+k1_1_prima/2)+fy_prima(x+h/2,w1[i-1],w2[i-1])*(u2++k1_2_prima/2))
        k3_1_prima=h*(u2+k2_2_prima/2)
        k3_2_prima=h*(fy(x+h/2,w1[i-1],w2[i-1])*(u1+k2_1_prima/2)+fy_prima(x+h/2,w1[i-1],w2[i-1])*(u2+k2_2_prima/2))
        k4_1_prima=h*(u2+k3_2_prima)
        k4_2_prima=h*(fy(x+h,w1[i-1],w2[i-1])*(u1+k3_1_prima)+fy_prima(x+h,w1[i-1],w2[i-1])*(u2+k3_2_prima))
        u1=u1+(k1_1_prima+2*k2_1_prima+2*k3_1_prima+k4_1_prima)/6
        u2=u2+(k1_2_prima+2*k2_2_prima+2*k3_2_prima+k4_2_prima)/6
    TK=TK-(w1[N]-beta)/u1
    k=k+1
print('Los valores de y en el intervalo pedido son:\n {}'.format(w1))
X_dom=np.zeros(N+1)#En esta matriz contendrá el dominio de X
X_dom[0]=a
for i in range(N):
    X_dom[i+1]=X_dom[i]+h
plt.plot(X_dom,w1)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Grafica de la Solucion Numerica-Shooting No Lineal')

