import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

#definimos una función para borrar lo que haya en el archivo phi.txt
def borrar_contenido(nombre_archivo):
    with open(nombre_archivo, 'w') as archivo:
        pass  # No es necesario escribir nada

#obtenemos la ubicación del directorio actual
directorio =os.getcwd()

#crear un bucle if donde se pregunte si existe el archivo phi.txt, y en caso de que no exista, crearlo en el directorio actual
if os.path.exists(os.path.join(directorio, 'phi.txt')):
    print('El archivo "phi.txt" existe en el directorio.')
    #como existe, borramos el contenido
    borrar_contenido('phi.txt')
else:
    print('El archivo "phi.txt" no existe en el directorio.')
    #si no existe, creamos el archivo
    archivo = open('phi.txt', 'w')
    archivo.close()

# Declaremos ahora nuestras variables
nx = 1001  # número de puntos total
time = 0.0 #tiempo inicial
dx = 0.002 #espaciado espacial
dt = 0.0005 #espaciado temporal
inv2dx = 1.0 / (2.0 * dx) #variable para simplificar
A = 1.0 #constante
sigma = 0.1 #constante
x0 = 0.0 #constante

# Inicializamos las funciones: de esta forma las funciones van de 0 a nx-1, con nx elementos
xh = np.zeros(nx)
phi = np.zeros(nx)
ppi = np.zeros(nx)
psi = np.zeros(nx)
ppi1 = np.zeros(nx)
psi1 = np.zeros(nx)
ppi2=np.zeros(nx)
psi2=np.zeros(nx)
ppi3=np.zeros(nx)
psi3=np.zeros(nx)

#Definimos la malla espacial, que va de -1 a 1 con un salto de 0.002
xh = np.linspace(-1.0, 1.0, nx)
#Datos iniciales del campo escalar
phi=A * np.exp(-(xh - x0)**2 / sigma**2)
psi=-2.0 * (xh - x0) * phi / sigma**2

#Escribimos el campo escalar en phi.txt para el instante inicial time=0
with open('phi.txt', 'a') as f:
    for i in range(nx):
        f.write(f'{time} {xh[i]} {phi[i]}\n')
    f.write('\n')

#Pidámosle al usuario que figura del texto quiere representar.
#En función de ello, se tomarán unos valores u otros de alpha, beta, y lastiter
opcion=int(input('Indique que figura del PDF quiere representar, indicando el número de la misma: '))
if opcion==1:
    alpha = 1.0
    beta = 0
    lastiter = 3000
elif opcion==2:
    alpha = 1.0
    beta = 1.0
    lastiter = 3000
elif opcion==3:
    alpha = np.zeros(nx)
    beta = 0.0
    #con lastiter=3000 se ve, pero se aprecia mejor la diferencia entre ambos pulsos para mayor tiempo
    lastiter = 4000
    for i in range(0, nx):
        alpha[i] = 0.25 * np.tanh(10 * xh[i]) + 0.75
elif opcion==4:
    alpha = 1.0
    beta = np.linspace(-1.0, 1.0, nx)
    lastiter = 3000
elif opcion==5:
    alpha = 1.0
    beta = np.linspace(-1.0, 1.0, nx)
    lastiter = 10000
else:
    print('Ha introducido un número no válido')


if opcion==1 or opcion==2:
    #Empezamos el algoritmo para la evolución temporal
    for iter in range(lastiter):
        time += dt

        #Runge-kutta (son las ecuaciones que aparecen en el paper, sin añadir sumando extras ni nada)
        #Primer paso para ppi1 y psi1, son las estrellas
        for i in range(1, nx-1):
            ppi1[i] = ppi[i] + dt * inv2dx * (alpha * (psi[i + 1] - psi[i - 1]) + beta * (ppi[i + 1] - ppi[i - 1]))
            psi1[i] = psi[i] + dt * inv2dx * (alpha * (ppi[i + 1] - ppi[i - 1]) + beta * (psi[i + 1] - psi[i - 1]))
        #Condiciones de frontera ppi1 y psi1
        ppi1[0] = ppi1[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi1[1] - ppi1[2])
        psi1[0] = ppi1[0]
        ppi1[nx-1] = ppi1[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (ppi1[nx - 2] - ppi1[nx - 3])
        psi1[nx-1] = -ppi1[nx - 1]

        #Segundo paso para las dos estrellas
        for i in range(1,nx-1):
            ppi2[i]=0.75*ppi[i]+0.25*ppi1[i]+0.25*dt*inv2dx*(alpha*(psi1[i+1]-psi1[i-1])+beta*(ppi1[i+1]-ppi1[i-1]))
            psi2[i]=0.75*psi[i]+0.25*psi1[i]+0.25*dt*inv2dx*(alpha*(ppi1[i+1]-ppi1[i-1])+beta*(psi1[i+1]-psi1[i-1]))
        #Condiciones de frontera para las dos estrellas
        ppi2[0]=ppi2[2]+(xh[0]-xh[2])/(xh[1]-xh[2])*(ppi2[1]-ppi2[2])
        psi2[0]=ppi2[0]
        ppi2[nx-1]=ppi2[nx-3]+(xh[nx-1]-xh[nx-3])/(xh[nx-2]-xh[nx-3])*(ppi1[nx-2]-ppi1[nx-3])
        psi2[nx-1]=-ppi2[nx-1]

        #Tercer paso
        for i in range(1,nx-1):
            ppi3[i]=(1/3)*ppi[i]+(2/3)*ppi2[i]+(2/3)*dt*inv2dx*(alpha*(psi2[i+1]-psi2[i-1])+beta*(ppi2[i+1]-ppi2[i-1]))
            psi3[i]=(1/3)*psi[i]+(2/3)*psi2[i]+(2/3)*dt*inv2dx*(alpha*(ppi2[i+1]-ppi2[i-1])+beta*(psi2[i+1]-psi2[i-1]))
        ppi3[0] = ppi3[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi3[1] - ppi3[2])
        psi3[0] = ppi3[0]
        ppi3[nx - 1] = ppi3[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (ppi3[nx - 2] - ppi3[nx - 3])
        psi3[nx - 1] = -ppi3[nx - 1]

        #La solucion final no será más que:
        phi=phi+dt*(alpha*ppi3+beta*psi3)
        ppi=ppi3
        psi=psi3

        #Vamos guardando los valores en phi.txt cada 50 iteraciones, y mostramos el tiempo y la iteración que vamos
        if iter % 50 == 0:
            with open('phi.txt', 'a') as f:
                for i in range(nx):
                    f.write(f'{time} {xh[i]} {phi[i]}\n')
                f.write('\n')
            print(f'Output time,it: {time}, {iter}')
elif opcion==3:
    #Empezamos el algoritmo para la evolución temporal
    for iter in range(lastiter):
        time += dt

        #Runge-kutta (ahora hay que tener en cuenta que alpha es un array)
        #Primer paso para ppi1 y psi1, son las estrellas
        for i in range(1, nx - 1):
            ppi1[i] = ppi[i] + dt * inv2dx * (
                        alpha[i] * (psi[i + 1] - psi[i - 1]) + (alpha[i + 1] - alpha[i - 1]) * psi[i] + beta * (
                            ppi[i + 1] - ppi[i - 1]))
            psi1[i] = psi[i] + dt * inv2dx * (
                        alpha[i] * (ppi[i + 1] - ppi[i - 1]) + (alpha[i + 1] - alpha[i - 1]) * ppi[i] + beta * (
                            psi[i + 1] - psi[i - 1]))
        #Condiciones de frontera ppi1 y psi1
        ppi1[0] = ppi1[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi1[1] - ppi1[2])
        psi1[0] = ppi1[0]
        ppi1[nx - 1] = ppi1[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (
                    ppi1[nx - 2] - ppi1[nx - 3])
        psi1[nx - 1] = -ppi1[nx - 1]

        #Segundo paso para las dos estrellas
        for i in range(1, nx - 1):
            ppi2[i] = 0.75 * ppi[i] + 0.25 * ppi1[i] + 0.25 * dt * inv2dx * (
                        alpha[i] * (psi1[i + 1] - psi1[i - 1]) + (alpha[i + 1] - alpha[i - 1]) * psi1[i] + beta * (
                            ppi1[i + 1] - ppi1[i - 1]))
            psi2[i] = 0.75 * psi[i] + 0.25 * psi1[i] + 0.25 * dt * inv2dx * (
                        alpha[i] * (ppi1[i + 1] - ppi1[i - 1]) + (alpha[i + 1] - alpha[i - 1]) * ppi1[i] + beta * (
                            psi1[i + 1] - psi1[i - 1]))
        #Condiciones de frontera para las dos estrellas
        ppi2[0] = ppi2[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi2[1] - ppi2[2])
        psi2[0] = ppi2[0]
        ppi2[nx - 1] = ppi2[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (
                    ppi1[nx - 2] - ppi1[nx - 3])
        psi2[nx - 1] = -ppi2[nx - 1]

        #Tercer paso
        for i in range(1, nx - 1):
            ppi3[i] = (1 / 3) * ppi[i] + (2 / 3) * ppi2[i] + (2 / 3) * dt * inv2dx * (
                        alpha[i] * (psi2[i + 1] - psi2[i - 1]) + (alpha[i + 1] - alpha[i - 1]) * psi2[i] + beta * (
                            ppi2[i + 1] - ppi2[i - 1]))
            psi3[i] = (1 / 3) * psi[i] + (2 / 3) * psi2[i] + (2 / 3) * dt * inv2dx * (
                        alpha[i] * (ppi2[i + 1] - ppi2[i - 1]) + (alpha[i + 1] - alpha[i - 1]) * ppi2[i] + beta * (
                            psi2[i + 1] - psi2[i - 1]))
        ppi3[0] = ppi3[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi3[1] - ppi3[2])
        psi3[0] = ppi3[0]
        ppi3[nx - 1] = ppi3[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (
                    ppi3[nx - 2] - ppi3[nx - 3])
        psi3[nx - 1] = -ppi3[nx - 1]

        # La solucion final no será más que:
        phi = phi + dt * (alpha * ppi3 + beta * psi3)
        ppi = ppi3
        psi = psi3

        # Vamos guardando los valores en phi.txt cada 50 iteraciones, y mostramos el tiempo y la iteración que vamos
        if iter % 50 == 0:
            with open('phi.txt', 'a') as f:
                for i in range(nx):
                    f.write(f'{time} {xh[i]} {phi[i]}\n')
                f.write('\n')
            print(f'Output time,it: {time}, {iter}')
elif opcion==4 or opcion==5:
    #Empezamos el algoritmo para la evolución temporal
    for iter in range(lastiter):
        time += dt

        #Runge-kutta (ahora hay que tener en cuenta que beta es un array)
        #Primer paso para ppi1 y psi1, son las estrellas
        for i in range(1, nx - 1):
            ppi1[i] = ppi[i] + dt * inv2dx * (
                        alpha * (psi[i + 1] - psi[i - 1]) + beta[i] * (ppi[i + 1] - ppi[i - 1]) + (
                            beta[i + 1] - beta[i - 1]) * ppi[i])
            psi1[i] = psi[i] + dt * inv2dx * (
                        alpha * (ppi[i + 1] - ppi[i - 1]) + beta[i] * (psi[i + 1] - psi[i - 1]) + (
                            beta[i + 1] - beta[i - 1]) * psi[i])
        #Condiciones de frontera ppi1 y psi1
        ppi1[0] = ppi1[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi1[1] - ppi1[2])
        psi1[0] = ppi1[0]
        ppi1[nx - 1] = ppi1[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (
                    ppi1[nx - 2] - ppi1[nx - 3])
        psi1[nx - 1] = -ppi1[nx - 1]

        #Segundo paso para las dos estrellas
        for i in range(1, nx - 1):
            ppi2[i] = 0.75 * ppi[i] + 0.25 * ppi1[i] + 0.25 * dt * inv2dx * (
                        alpha * (psi1[i + 1] - psi1[i - 1]) + beta[i] * (ppi1[i + 1] - ppi1[i - 1]) + (
                            beta[i + 1] - beta[i - 1]) * ppi1[i])
            psi2[i] = 0.75 * psi[i] + 0.25 * psi1[i] + 0.25 * dt * inv2dx * (
                        alpha * (ppi1[i + 1] - ppi1[i - 1]) + beta[i] * (psi1[i + 1] - psi1[i - 1]) + (
                            beta[i + 1] - beta[i - 1]) * psi1[i])
        #Condiciones de frontera para las dos estrellas
        ppi2[0] = ppi2[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi2[1] - ppi2[2])
        psi2[0] = ppi2[0]
        ppi2[nx - 1] = ppi2[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (
                    ppi1[nx - 2] - ppi1[nx - 3])
        psi2[nx - 1] = -ppi2[nx - 1]

        #Tercer paso
        for i in range(1, nx - 1):
            ppi3[i] = (1 / 3) * ppi[i] + (2 / 3) * ppi2[i] + (2 / 3) * dt * inv2dx * (
                        alpha * (psi2[i + 1] - psi2[i - 1]) + beta[i] * (ppi2[i + 1] - ppi2[i - 1]) + (
                            beta[i + 1] - beta[i - 1]) * ppi2[i])
            psi3[i] = (1 / 3) * psi[i] + (2 / 3) * psi2[i] + (2 / 3) * dt * inv2dx * (
                        alpha * (ppi2[i + 1] - ppi2[i - 1]) + beta[i] * (psi2[i + 1] - psi2[i - 1]) + (
                            beta[i + 1] - beta[i - 1]) * psi2[i])
        ppi3[0] = ppi3[2] + (xh[0] - xh[2]) / (xh[1] - xh[2]) * (ppi3[1] - ppi3[2])
        psi3[0] = ppi3[0]
        ppi3[nx - 1] = ppi3[nx - 3] + (xh[nx - 1] - xh[nx - 3]) / (xh[nx - 2] - xh[nx - 3]) * (
                    ppi3[nx - 2] - ppi3[nx - 3])
        psi3[nx - 1] = -ppi3[nx - 1]

        # La solucion final no será más que:
        phi = phi + dt * (alpha * ppi3 + beta * psi3)
        ppi = ppi3
        psi = psi3

        # Vamos guardando los valores en phi.txt cada 50 iteraciones, y mostramos el tiempo y la iteración que vamos
        if iter % 50 == 0:
            with open('phi.txt', 'a') as f:
                for i in range(nx):
                    f.write(f'{time} {xh[i]} {phi[i]}\n')
                f.write('\n')
            print(f'Output time,it: {time}, {iter}')

# Cargar los datos del archivo
data = np.loadtxt('phi.txt')

# Separar los datos en diferentes variables
tiempo = data[:, 0]
posicion = data[:, 1]
funcion = data[:, 2]

# Crear una figura 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Hacer el gráfico
ax.scatter(tiempo, posicion, funcion, c='black', s=3,alpha=0.1)

# Etiquetas de los ejes
ax.set_xlabel('Tiempo')
ax.set_ylabel('Posición')
ax.set_zlabel('Función')

# Mostrar el gráfico
plt.show()