import math

import sympy
from einsteinpy.symbolic import MetricTensor,ChristoffelSymbols,RicciTensor, RicciScalar,EinsteinTensor
from sympy import Function,exp

syms=sympy.symbols('t r theta phi') #definimos las variables que tendremos
G,M,c,m=sympy.symbols('G M c m ')  #para poder usar las ctes
A=Function('A')(syms[1])
B=Function('B')(syms[1])

#definimos el tensor métrico:
list2d=[[0 for i in range(4)] for i in range(4)]
list2d[0][0]=exp(2*A)
list2d[1][1]=-exp(2*B)
list2d[2][2]=-syms[1]**2
list2d[3][3]=-(syms[1]*sympy.sin(syms[2]))**2

#creamos el tensor métrico con su forma conveniente, y los símbolos de christoffel
metric=MetricTensor(list2d,syms)

metric_ch=ChristoffelSymbols.from_metric(metric)
simplified=metric_ch.simplify()

tensorricci=RicciTensor.from_metric(metric)
simplified2=tensorricci.simplify()

scalarricci=RicciScalar.from_metric(metric)
simplified3=scalarricci.simplify()

tensoreinstein=EinsteinTensor.from_metric(metric)
simplified4=tensoreinstein.simplify()

#mostremos la métrica
print('Tensor métrico:')
for a in metric.tensor():
    print(a)
print('\n') #salto de línea

#y ahora los símbolos de christoffel
print('Simbolos de Christoffel:')
for a in metric_ch.tensor():
    for b in a:
        print(b)
    print()
print()
#esta serie de prints nos mostrará, en 4d, 4 distintos arrays pq simbolos son 2 veces contra, una cov, conque irá mostrando,
#en cada array, primero la coordenada temporal, despues r, theta y phi, en ese orden.

#también nos la podemos apañar para mostrar el tensor de Ricci, el escalar de Ricci y el Tensor de Einstein:
print('Tensor Ricci:')
for a in tensorricci.tensor():
    print(a)
print('\n') #salto de línea

print('Tensor Einstein:')
for a in tensoreinstein.tensor():
    print(a)
print('\n') #salto de línea

print(f'Escalar de Ricci: \n {scalarricci.tensor()}') #salen con signos cambiados, es solo convenio
