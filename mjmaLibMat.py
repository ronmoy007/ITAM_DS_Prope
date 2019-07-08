import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import scipy.linalg 
import pprint

def funcPrueba():
    print('Hola mundo')
    return

######################## Funciones Auxiliares
# Función que recibe 2 números y devuleve el punto medio entre ellos
def ObtenerPuntoMedio(nbrX1,nbrX2):
    nbrPuntoMedio=((nbrX2-nbrX1)/2)+nbrX1
    return nbrPuntoMedio

def CalcErrorRelat(p1,p2):
    return np.abs(p1-p2)/math.fabs(p2)

def GenerarDictGaussLeg():

    dictGeneral={}

    dictPesosN0={0:2}
    dictPesosN1={0:1, 1:1}
    dictPesosN2={0:5/9, 1:8/9, 2:5/9}
    dictPesosN3={0:0.347855, 1:0.652145, 2:0.652145, 3:0.347855}
    dictPesosN4={0:0.236927, 1:0.478629, 2:0.568889, 3:0.478629, 4:0.236927}

    dictNodosN0={0:0}
    dictNodosN1={0:-math.sqrt(1/3), 1:math.sqrt(1/3)}
    dictNodosN2={0:-math.sqrt(3/5), 1:0, 2:math.sqrt(3/5)}
    dictNodosN3={0:-0.861136, 1:-0.339981, 2:0.339981, 3:0.861136}
    dictNodosN4={0:-0.90618, 1:-0.538469, 2:0, 3:0.538469, 4:0.90618}

    dictGeneral[0]={'pesos':dictPesosN0, 'nodos':dictNodosN0}
    dictGeneral[1]={'pesos':dictPesosN1, 'nodos':dictNodosN1}
    dictGeneral[2]={'pesos':dictPesosN2, 'nodos':dictNodosN2}
    dictGeneral[3]={'pesos':dictPesosN3, 'nodos':dictNodosN3}
    dictGeneral[4]={'pesos':dictPesosN4, 'nodos':dictNodosN4}
    
    return dictGeneral

def GenerarDictGaussHer():
    
    pi=math.pi
    
    dictGeneral={}
    
    dictNodosN0={0:0}
    dictNodosN1={0:-.707107,  1:.707107}
    dictNodosN2={0:-1.224745, 1:0,        2:1.224745}
    dictNodosN3={0:-1.650680, 1:-.524648, 2:0.524648, 3:1.650680}
    dictNodosN4={0:-2.020183, 1:-0.958572, 2:0,        3:0.958572, 4:2.020183}
    dictNodosN5={0:-2.350605, 1:-1.335849,2:-.436077, 3:0.436077,4:1.335849, 5:2.350605}
    
    dictPesosN0={0:math.sqrt(pi)}
    dictPesosN1={0:0.886227,  1:0.886227}
    dictPesosN2={0:0.295409,  1:1.181636,       2:0.295409}
    dictPesosN3={0:0.813128*10**-1, 1:0.804914, 2:0.804914, 3:0.813128*10**-1}
    dictPesosN4={0:0.199532*10**-1, 1:0.393619, 2:0.945308, 3:0.393619, 4:0.199532*10**-1}
    dictPesosN5={0:0.453001*10**-2, 1:0.157067, 2:0.724629, 3:0.724629, 4:0.157067,   5:0.453001*10**-2}

    dictGeneral[0]={'pesos':dictPesosN0, 'nodos':dictNodosN0}
    dictGeneral[1]={'pesos':dictPesosN1, 'nodos':dictNodosN1}
    dictGeneral[2]={'pesos':dictPesosN2, 'nodos':dictNodosN2}
    dictGeneral[3]={'pesos':dictPesosN3, 'nodos':dictNodosN3}
    dictGeneral[4]={'pesos':dictPesosN4, 'nodos':dictNodosN4}
    dictGeneral[5]={'pesos':dictPesosN5, 'nodos':dictNodosN5}
    
    return dictGeneral

def VerificarMatricesEnPantalla(mtrxA, mtrxB, mtrxX):
    print('Verificamos que al multiplicar AX, obtengamos B')
    print('mtrxA:')
    pprint.pprint(mtrxA)
    print('mtrxX:')
    pprint.pprint(mtrxX)
    print('mtrxA@mtrxX:')
    pprint.pprint(mtrxA@mtrxX)
    print('mtrxB:')
    pprint.pprint(mtrxB)

######################## Funciones para calcular aproximaciones por diferencias 
def AproximarDifAd1(f,x,h):
    nbrRslt=(f(x+h)-f(x))/h
    return nbrRslt

def AproximarDifAd2(f,x,h):
    nbrRslt=(f(x+2*h)-2*f(x+h)+f(x))/h**2
    return nbrRslt

'''
def funcAproxDifAd(f,x,h):
    nbrRslt1=funcAproxDifAd1(f,x,h)
    nbrRslt2=funcAproxDifAd2(f,x,h)
    return nbrRslt1,nbrRslt2
'''

######################## Funciones para calcular integrales por diferenes métodos 
# Función que obtiene la integral de una función mediante el método del rectángulo
# Recibe como parámetros: 
    # func: la función a integrar
    # a y b: los rangos de la integral
    # n: la cantidad de intervalos que se usarán para generar los rectángulos
def IntegrarXRect(func, a, b, n):
    
    #Variables de la función
    arrIntervalos=[]
    nbrResultadoFuncion=0
    nbrIndice=0
    nbrAcumulado=0
    nbrPuntoMedio=0
    nbrDx=0
    
    #Generamos el arreglo de intervalos con los parámetros de entrada de la función
    arrIntervalos=np.linspace(a,b,n+1)
    
    #Obtenemos el valor de cada rango (dx). Puesto que los intervalos son equidistantes, bastaa con calcular dx una vez 
    nbrDx=arrIntervalos[1]-arrIntervalos[0]
    
    #Recorremos el intervalo dado para obtener los puntos medios y acumularlos
    for i in arrIntervalos:
        
        if nbrIndice>0:
            
            #Obtenemos el punto medio de cada rango
            nbrPuntoMedio=ObtenerPuntoMedio(arrIntervalos[nbrIndice-1],arrIntervalos[nbrIndice])

            #El área de un rectángulo es:
            # (base) X (altura)

            #Obtenemos el área del rectángulo. (nbrDx es la base y el resultado de la función la altura)
            nbrArea=func(nbrPuntoMedio)*nbrDx

            #Acumulamos el valor para aproximarnos a la integral
            nbrAcumulado=nbrAcumulado+nbrArea
        
        nbrIndice=nbrIndice+1
            
    return nbrAcumulado

# Función que obtiene la integral de una función mediante el método del trapezoide
# Recibe como parámetros: 
    # func: la función a integrar
    # a y b: los rangos de la integral
    # n: la cantidad de intervalos que se usarán para generar los rectángulos
def IntegrarXTrap(func, a, b, n):
    
    #Variables de la función
    arrIntervalos=[]
    nbrArea=0
    nbrIndice=0
    nbrAcumulado=0
    nbrBase1=0
    nbrBase2=0
    nbrDx=0
    
    #Generamos el arreglo de intervalos con los parámetros de entrada de la función
    arrIntervalos=np.linspace(a,b,n+1)
    
    #Obtenemos el valor de cada rango (dx). Puesto que los intervalos son equidistantes, bastaa con calcular dx una vez 
    nbrDx=arrIntervalos[1]-arrIntervalos[0]
    
    #Recorremos el intervalo dado para obtener los puntos medios y acumularlos
    for i in arrIntervalos:
        
        if nbrIndice>0:
            
            #La forma de obtener el área de un trapezoide es:
            # (base mayor + base menor) X (altura)
            # ----------------------------------
            #                 2

            # Obtenemos las bases
            nbrBase1=func(arrIntervalos[nbrIndice-1])
            nbrBase2=func(arrIntervalos[nbrIndice])

            #Obtenemos el área (nbrDx es la altura)
            nbrArea=(((nbrBase1+nbrBase2)*(nbrDx))/2)

            #Acumulamos el valor para aproximarnos a la integral
            nbrAcumulado=nbrAcumulado+nbrArea
        
        nbrIndice=nbrIndice+1
        
    return nbrAcumulado

def IntegrarXGaussLeg(func,a,b,n):

    dictGausLeg={}
    dictGausLeg=GenerarDictGaussLeg()
    dictGausLeg

    Wi=0
    Xi=0
    nbrParametro=0
    nbrAcumulado=0
    nbrResultado=0

    nbrIndice=0
    
    #Recorremos el intervalo dado
    for i in range(0,n+1):

        Wi=dictGausLeg[n]['pesos'][i]
        Xi=dictGausLeg[n]['nodos'][i]

        nbrParametro=(.5)*(((b-a)*(Xi))+a+b)
        nbrAcumulado=nbrAcumulado+(Wi*func(nbrParametro))

    nbrResultado=((b-a)/2)*(nbrAcumulado)
    
    return nbrResultado

def IntegrarXGaussHer(n,func):
    
    #Se genera el diccionario Gauss-Hermite con los nodos y pesos
    dictGausHer={}
    dictGausHer=GenerarDictGaussHer()

    #Generamos un arreglo de acuerdo al rango que vamos a barrer
    arrResultadosFunc=np.zeros(n+1)
    
    #Barremos el rango solicitado
    for i in range(n+1):

        #Buscamos el nodo correspondiente en el diccionario
        nbrXi=dictGausHer[n]['nodos'][i]
        
        #Ponemos el resultado de la función del nodo correspondiente en el arreglo temporal de resultados
        arrResultadosFunc[i]=func(nbrXi)
      
    #Convertimos el diccionario de pesos en un numpy array
    arrayPesos = np.array(list(dictGausHer[n]['pesos'].values()))
    
    #Calculamos el producto dot con el arreglo de pesos y los resultado de la función evaluados en los distintos nodos
    nbrProductoDot=np.dot(arrayPesos,arrResultadosFunc)
        
    return nbrProductoDot

######################## Funciones para graficar integrales por diferenes métodos
# Función que grafica la integral de una función mediante rectángulos
# Recibe como parámetros: 
    # func: la función a integrar
    # a y b: los rangos de la integral
    # n: la cantidad de intervalos que se usarán para generar los rectángulos
def GraficarIntXRect(func,a,b,n):
    
    #Generamos el arreglo de intervalos con los parámetros de entrada de la función
    arrIntervalos=np.linspace(a,b,n+1)
    
    #Obtenemos el valor de cada rango (dx). Puesto que los intervalos son equidistantes, basta con calcular dx una vez 
    nbrDx=arrIntervalos[1]-arrIntervalos[0]
    
    nbrIndice=0
    
    #Recorremos el intervalo dado para obtener los puntos medios y acumularlos
    for i in arrIntervalos:
            
        #Para generar el rectángulo, necesitamos 4 vertices:
        #  (xn-1,f(pm))----f(pm)----(xn,f(pm))
        #             |               |
        #             |               |
        #             |               |
        #             |               |
        #             |               |
        #             |               |
        #             |               |
        #             |               |
        #             |<---- dx ----->|
        #             |               |
        #             |               |
        #       (xn-1,0)    pm     (xn,0)
        
        #
        
        if nbrIndice>0:
           
            XnMenos1=arrIntervalos[nbrIndice-1]
            Xn=arrIntervalos[nbrIndice]
            
            #Obtenemos el punto medio de cada rango para determinar la altura que tendrá cada rectángulo
            nbrPuntoMedio=ObtenerPuntoMedio(arrIntervalos[nbrIndice-1],arrIntervalos[nbrIndice])
            nbrFuncPm=func(nbrPuntoMedio)

            #Generamos la línea de la izquierda
            x=[XnMenos1,XnMenos1]
            y=[0,nbrFuncPm]
            plt.plot(x,y)

            #Generamos la línea de arriba
            x=[XnMenos1,Xn]
            y=[nbrFuncPm,nbrFuncPm]
            plt.plot(x,y)

            #Generamos la línea de la derecha
            x=[Xn,Xn]
            y=[nbrFuncPm,0]
            plt.plot(x,y)
            
        nbrIndice=nbrIndice+1
    
    plt.grid(True)
    plt.show()
        
    return 

def GraficarIntXTrap(func,a,b,n):
    
    #Generamos el arreglo de intervalos con los parámetros de entrada de la función
    arrIntervalos=np.linspace(a,b,n+1)
    
    #Obtenemos el valor de cada rango (dx). Puesto que los intervalos son equidistantes, basta con calcular dx una vez 
    nbrDx=arrIntervalos[1]-arrIntervalos[0]
    
    nbrIndice=0
    
    #Recorremos el intervalo dado para obtener los puntos medios y acumularlos
    for i in arrIntervalos:
            
        #Para generar el rectángulo, necesitamos 4 vertices:
        #                      (xn,f(xn))
        #                     - |
        #                    -  |
        #                   -   |
        #                  -    |
        #                 -     |
        #                -      |
        #               -       |
        #(xn-1,f(xn-1))-        |
        #             |         |
        #             |         |
        #             |         |
        #       (xn-1,0)        (xn,0)
        
        #
        
        if nbrIndice>0:
            
            #Obtenemos el punto medio de cada rango
            XnMenos1=arrIntervalos[nbrIndice-1]
            Xn=arrIntervalos[nbrIndice]
            nbrFuncXnMenos1=func(XnMenos1)
            nbrFuncXn=func(Xn)

            #Generamos la línea de la izquierda
            x=[XnMenos1,XnMenos1]
            y=[0,nbrFuncXnMenos1]
            plt.plot(x,y)

            #Generamos la línea de arriba
            x=[XnMenos1,Xn]
            y=[nbrFuncXnMenos1,nbrFuncXn]
            plt.plot(x,y)

            #Generamos la línea de la derecha
            x=[Xn,Xn]
            y=[nbrFuncXn,0]
            plt.plot(x,y)
        
        nbrIndice=nbrIndice+1

    plt.grid(True)
    plt.show()
        
    return 

#######Funciones para resolver sistemas de ecuaciones matriciales
def ResolverMatrizXPLU(mtrxA, mtrxB):

    #Generamos las matricex Q y R con el módulo linalg de scipy:
    mtrxQ, mtrxR = scipy.linalg.qr(mtrxA)

    #Generamos la matriz p, utilizando la traspuesta de Q, y aplicando el producto punto con B
    mtrxP = np.dot(mtrxQ.T, mtrxB)


    #Generamos la matriz que contiene los vectores que resuelven el sistema A,B
    #Para ello, obtenemos el producto punto entre la inversa de R y P
    mtrxX = np.dot(scipy.linalg.inv(mtrxR), mtrxP) # Solve Rx=y 

    #Verificamos que al multiplicar AX, obtengamos B
    return mtrxX

####
def CalcularECM(arrY, arrYY):
    arrTmp=(arrY-arrYY)**2
    nbrResultado=np.sum(arrTmp)
    return nbrResultado