import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import os
import csv

carpeta_principal="black holes comprimido"
carpetas = os.listdir(carpeta_principal)
print(carpetas)
for carpeta in carpetas:
    ruta_carpeta = os.path.join(carpeta_principal, carpeta)
    carpetas2=os.listdir(ruta_carpeta)
    for carpeta2 in carpetas2:
        ruta_carpeta2=os.path.join(ruta_carpeta,carpeta2)
        #Iterador que nos ira pasando por todos los elementos de la carpeta
        iterador = os.scandir(ruta_carpeta2)

        # Hacemos una lista de la energia porque se iran mostrando sus valores para cada l,m
        energia = []
        # Indice con el que iremos rellenando la lista energia
        indice = 0
        #f(u) que sera la suma de todas las contribuciones de las distintas f(u) para los valores de l,m
        funcion_u_total=0
        #donde guardaremos el tiempo, que es el mismo para un mismo BH, y que usaremos para representar funcion_u_total
        tiempo=[]

        # Vamos viendo los archivos uno a uno, y vemos sus elementos
        for archivo in iterador:
            # En todas las carpetas hay un archivo llamado Metadata que no nos interesa
            if archivo.name == 'Metadata':
                #Abrir el archivo en modo lectura
                archivo2 = open(archivo, 'r')

                #Leer el contenido del archivo línea por línea
                lineas = archivo2.readlines()

                #Extraer las etiquetas y valores y asignarlos al diccionario
                diccionario = {}
                for linea in lineas:
                    if '=' in linea and linea.strip() != '':
                        etiqueta, valor = linea.strip().split('=',1)
                        etiqueta = etiqueta.strip()
                        valor = valor.strip()
                        diccionario[etiqueta] = valor
                # Paso 7: Cerrar el archivo
                archivo2.close()

                # Mostrar el diccionario
                print(diccionario)  #vemos que nos lo hace bien, de aqui obtener los parametros de interior

                continue
                #si no lo consigues, era poner 'continue'
            if archivo.is_file():
                data = np.loadtxt(archivo)
                tiempo = data[:, 0]
                real_psi = data[:, 1]
                im_psi = data[:, 2]

                # con longitud nos podremos ir moviendo entre cada elemento especifico del archivo
                longitud = len(tiempo)

                # Creamos lista vacia para las funciones u_prima que se van a integrar, cuyos valores se guardaran aqui
                in_im_psi = np.zeros(longitud)
                in_real_psi = np.zeros(longitud)

                # se empieza a integrar desde el primer elemento ya que el elemento 0 es igual a 0 (se integra una linea)
                for i in range(1, longitud):
                    # integramos entre los puntos inicial y el final i
                    in_im_psi[i] = integrate.simps(im_psi[:i], tiempo[:i])
                    in_real_psi[i] = integrate.simps(real_psi[:i], tiempo[:i])

                funcion_u = -real_psi * in_im_psi + in_real_psi * im_psi
                energia.append(integrate.simps(funcion_u, tiempo))
                print(f"Nombre archivo: {archivo}, Energia*72*pi**2={energia[indice]}")
                indice += 1
                funcion_u_total+=funcion_u


        # La energia total no sera mas que la suma
        energia_total=sum(energia)
        print(f'Lo que aparece en el paper es {energia_total}')
        chern_pontryagin = sum(energia)/(72*np.pi**2)
        print(chern_pontryagin)
        diccionario['chern_pontryagin']=chern_pontryagin

        #mostramos la funcion f(u) total sumado para todos los l,m
        #plt.plot(tiempo,funcion_u_total)
        #plt.show()

        #diferenciamos entre los casos de Precissing, Non precissing, Aligned,
        if diccionario['system-type']=='Precessing':

            claves_seleccionadas = ['catalog-tag','system-type','chern_pontryagin','eccentricity','initial-ADM-energy','initial-orbital-angular-momentum','initial-separation', 'initial-mass1','initial-mass2','initial-total-mass','initial-bh-chi1x','initial-bh-chi1y','initial-bh-chi1z','initial-bh-chi2x','initial-bh-chi2y','initial-bh-chi2z','relaxed-mass-ratio-1-over-2','final-mass','final-chi','peak-luminosity-ergs-per-sec','final-kick']
            #añado lo siguiente para que la primera fila del txt sea de texto que no se puede introducir y sea un archivo numerico
            clave=['#']
            clave_total = clave + claves_seleccionadas
            claves_seleccionadas_bonita = []
            for i in range(0, len(claves_seleccionadas)):
                claves_seleccionadas_bonita.append(f'{i}:{claves_seleccionadas[i]}')
            clave_total_bonita = clave + claves_seleccionadas_bonita


            directorio =os.getcwd()
            if os.path.exists(os.path.join(directorio, 'metadata_precessing_eccentricity.txt')):
                pass
            else:
                print('El archivo "metadata_precessing_eccentricity.txt" no existe en el directorio.')
                #creamos el archivo entonces
                archivo = open('metadata_precessing_eccentricity.txt', 'w')
                writer = csv.writer(archivo, delimiter=' ')
                # Escribe las claves del diccionario en la primera fila
                writer.writerow(clave_total_bonita)
                archivo.close()

            #Una vez tenemos el archivo creado, añadamos lo que nos interesa
            with open('metadata_precessing_eccentricity.txt', 'a') as f:
                writer = csv.writer(f, delimiter=' ')
                writer.writerow([diccionario[clave] for clave in claves_seleccionadas])

                #y con lo siguiente eliminamos las filas en blanco del archivo
            with open('metadata_precessing_eccentricity.txt','r') as f:
                lineas = f.readlines()
            with open('metadata_precessing_eccentricity.txt', "w") as archivo:
                for linea in lineas:
                    if linea.strip():
                        archivo.write(linea)


        elif diccionario['system-type']=='Nonspinning':
            claves_seleccionadas = ['catalog-tag', 'system-type', 'chern_pontryagin','eccentricity', 'initial-ADM-energy',
                                    'initial-orbital-angular-momentum', 'initial-separation', 'initial-mass1',
                                    'initial-mass2', 'initial-bh-chi1z', 'initial-bh-chi2z',
                                    'relaxed-mass-ratio-1-over-2', 'final-mass', 'final-chi',
                                    'peak-luminosity-ergs-per-sec', 'final-kick']
            # añado lo siguiente para que la primera fila del txt sea de texto que no se puede introducir y sea un archivo numerico
            clave = ['#']
            clave_total = clave + claves_seleccionadas
            claves_seleccionadas_bonita = []
            for i in range(0, len(claves_seleccionadas)):
                claves_seleccionadas_bonita.append(f'{i}:{claves_seleccionadas[i]}')
            clave_total_bonita = clave + claves_seleccionadas_bonita

            directorio = os.getcwd()
            if os.path.exists(os.path.join(directorio, 'metadata_non_spinning.txt')):
                pass
            else:
                print('El archivo "metadata_non_spinning.txt" no existe en el directorio.')
                # creamos el archivo entonces
                archivo = open('metadata_non_spinning.txt', 'w')
                writer = csv.writer(archivo, delimiter=' ')
                # Escribe las claves del diccionario en la primera fila
                writer.writerow(clave_total_bonita)
                archivo.close()

            # Una vez tenemos el archivo creado, añadamos lo que nos interesa
            with open('metadata_non_spinning.txt', 'a') as f:
                writer = csv.writer(f, delimiter=' ')
                writer.writerow([diccionario[clave] for clave in claves_seleccionadas])

                # y con lo siguiente eliminamos las filas en blanco del archivo
            with open('metadata_non_spinning.txt', 'r') as f:
                lineas = f.readlines()
            with open('metadata_non_spinning.txt', "w") as archivo:
                for linea in lineas:
                    if linea.strip():
                        archivo.write(linea)


        elif diccionario['system-type']=='Aligned':
            claves_seleccionadas = ['catalog-tag', 'system-type','chern_pontryagin', 'eccentricity', 'initial-ADM-energy',
                                    'initial-orbital-angular-momentum', 'initial-separation', 'initial-mass1',
                                    'initial-mass2', 'initial-bh-chi1z', 'initial-bh-chi2z',
                                    'relaxed-mass-ratio-1-over-2', 'final-mass', 'final-chi',
                                    'peak-luminosity-ergs-per-sec', 'final-kick']
            # añado lo siguiente para que la primera fila del txt sea de texto que no se puede introducir y sea un archivo numerico
            clave = ['#']
            clave_total = clave + claves_seleccionadas
            claves_seleccionadas_bonita = []
            for i in range(0, len(claves_seleccionadas)):
                claves_seleccionadas_bonita.append(f'{i}:{claves_seleccionadas[i]}')
            clave_total_bonita = clave + claves_seleccionadas_bonita

            directorio = os.getcwd()
            if os.path.exists(os.path.join(directorio, 'metadata_aligned.txt')):
                pass
            else:
                print('El archivo "metadata_aligned.txt" no existe en el directorio.')
                # creamos el archivo entonces
                archivo = open( 'metadata_aligned.txt', 'w')
                writer = csv.writer(archivo, delimiter=' ')
                # Escribe las claves del diccionario en la primera fila
                writer.writerow(clave_total_bonita)
                archivo.close()

            # Una vez tenemos el archivo creado, añadamos lo que nos interesa
            with open( 'metadata_aligned.txt', 'a') as f:
                writer = csv.writer(f, delimiter=' ')
                writer.writerow([diccionario[clave] for clave in claves_seleccionadas])

                # y con lo siguiente eliminamos las filas en blanco del archivo
            with open( 'metadata_aligned.txt', 'r') as f:
                lineas = f.readlines()
            with open( 'metadata_aligned.txt', "w") as archivo:
                for linea in lineas:
                    if linea.strip():
                        archivo.write(linea)


