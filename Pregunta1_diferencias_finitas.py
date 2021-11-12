"""
Diferencias Finitas Lineal
Problema 1
"""
import numpy as np
import math
import matplotlib.pyplot as plt
#La ecuacion diferencial a resolver es:
# y''(x)=p(x)*y'+q(x)*y+r(x)
#para a<=x<=b y(a)=alpha y(b)=beta
#La ecuacion a resolver es:
#y''=-2*x*y
#y(0)=1 y(1)=-0.5
def p(x):
    return 0
def q(x):
    return -2*x
def r(x):
    return 0
a=0
b=1
N=999
alpha=1
beta=-0.5
h=(b-a)/(N+1)
x=a+h
a_m=np.zeros(N+1)
b_m=np.zeros(N+1)
c=np.zeros(N+1)
d=np.zeros(N+1)
a_m[1]=2+q(x)*h**2 
b_m[1]=-1+(h/2)*p(x)
d[1]=-r(x)*h**2+(1+h*p(x)/2)*alpha
for i in range(2,N):
    x=a+i*h
    a_m[i]=2+q(x)*h**2 
    b_m[i]=-1+(h/2)*p(x)
    c[i]=-1-h*p(x)/2
    d[i]=-r(x)*h**2
x=b-h
a_m[N]=2+q(x)*h**2
c[N]=-1-(h/2)*p(x)
d[N]=-r(x)*h**2+(1-h*p(x)/2)*beta
#Los siguientes pasos corresponden al desarrollo 
#de un sistema tridiagonal
l=np.zeros(N+1)
u=np.zeros(N+1)
z=np.zeros(N+1)
l[1]=a_m[1]
u[1]=b_m[1]/a_m[i]
z[1]=d[1]/l[1]
for i in range(2,N):
    l[i]=a_m[i]-c[i]*u[i-1]
    u[i]=b_m[i]/l[i]
    z[i]=(d[i]-c[i]*z[i-1])/l[i]
l[N]=a_m[N]-c[N]*u[N-1]
z[N]=(d[N]-c[N]*z[N-1])/l[N]    
#Definimos una nueva matriz para trabajar
w=np.zeros(N+2)
w[0]=alpha
w[N+1]=beta
w[N]=z[N]
for i in range(N-1,1,-1):
    w[i]=z[i]-u[i]*w[i+1]
#Finalmente imprimimos los datos:
print('Los valores de <y> con el paso h={} son: \n {}'.format(h,w))

#Graficamos la solución 
X_dom=np.zeros(N+2)
X_dom[0]=a
for i in range(N+1):
    X_dom[i+1]=X_dom[i]+h
plt.plot(X_dom,w)
plt.title('Grafica de la solucion numerica-Diferencias Finitas')
plt.xlabel('X')
plt.ylabel('Y')