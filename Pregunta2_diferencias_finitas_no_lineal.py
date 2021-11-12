"""
Metodo diferencias finitas no lineal Metodo de Newton
"""
import numpy as np
import matplotlib.pyplot as plt
#La ED con valores de frontera a resolver es la siguiente:
## RESOLVER EDO de la forma: y''=f(x,y,y')=-2*x*y^2 y(0)=1 y(1)=1.5
#Definimos los parametros de entrada:
a=0
b=1
alpha=1
beta=1.5
N=199
M=100#Nos da la condicion de maximo numero de iteraciones antes de detenernos
h=(b-a)/(N+1)
def f(x,y,y_prima):
    return -2*x*y**2
def fy(x,y,y_prima):
    return -4*x*y
def fy_prima(x,y,y_prima):
    return 0
w=np.zeros(N+2)
w[0]=alpha
w[N+1]=beta
a_m=np.zeros(N+1)
b_m=np.zeros(N+1)
c=np.zeros(N+1)
d=np.zeros(N+1)
#Las matrices que nos serviran para hallar la solucion del sistema de ecuaciones tridiagonal:
l=np.zeros(N+1)
u=np.zeros(N+1)
z=np.zeros(N+1)
v=np.zeros(N+1)
for i in range(1,N+1):
    w[i]=alpha+(beta-alpha)*h/(b-a)
k=1#inicializamos el valor con el que comprobaremos que no se supere M
while k<=M:
    x=a+h
    t=(w[2]-2*alpha)/h
    a_m[1]=2+fy(x,w[1],t)*h**2
    b_m[1]=-1+fy_prima(x,w[1],t)*h/2
    d[1]=-(2*w[1]-w[2]-alpha+f(x,w[1],t)*h**2)
    for i in range(2,N):
        x=a+i*h
        t=(w[i+1]-w[i-1])/(2*h)
        a_m[i]=2+fy(x,w[i],t)*h**2
        b_m[i]=-1+fy_prima(x,w[i],t)*h/2
        c[i]=-1-fy_prima(x,w[i],t)*h/2
        d[i]=-(2*w[i]-w[i+1]-w[i-1]+f(x,w[i],t)*h**2)
    x=b-h
    t=(beta-w[N-1])/(2*h)
    a_m[N]=2+fy(x,w[N],t)*h**2
    c[N]=-1-fy_prima(x,w[N],t)*h/2
    d[N]=-(2*w[N]-beta-w[N-1]+f(x,w[N],t)*h**2)
    l[1]=a_m[1]
    u[1]=b_m[1]/a_m[1]
    z[1]=d[1]/l[1]
    for i in range(2,N):
        l[i]=a_m[i]-c[i]*u[i-1]
        u[i]=b_m[i]/l[i]
        z[i]=(d[i]-c[i]*z[i-1])/l[i]
    l[N]=a_m[N]-c[N]*u[N-1]
    z[N]=(d[N]-c[N]*z[N-1])/l[N]
    v[N]=z[N]
    w[N]=w[N]+v[N]
    for i in range(N-1,0,-1):
        v[i]=z[i]-u[i]*v[i+1]
        w[i]=w[i]+v[i]
    k=k+1
print('Los valores de "y" en el intervalo pedido con un ancho={} son:\n {}'.format(h,w))
X_dom=np.zeros(N+2)#En esta matriz contendrá el dominio de X
X_dom[0]=a
for i in range(N+1):
    X_dom[i+1]=X_dom[i]+h
plt.plot(X_dom,w)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Grafica de la Solucion Numerica-Shooting No Lineal')


