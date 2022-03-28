import Matrices_Vectores as matrivec
import numeros_complejos as numc
import numpy as np
import math
from matplotlib import pyplot as plt

"""Funcion que calcula la probabilidad de encontrar una particula en una posicion"""

def probabilidadPosicion(posicion, vector1):
    return (numc.moduloCplx(vector1[posicion][0]) ** 2) / matrivec.normaVector(vector1) ** 2


"""Funcion que calcula la amplitud al transitar de un vector a otro"""


def amplitudTransicion(v1, v2):
    bra = matrivec.adjunta(v2)
    amplitud = matrivec.productMatrix(bra, v1)[0][0]
    norma1 = matrivec.normaVector(v1)
    norma2 = matrivec.normaVector(v2)
    a = numc.divisionCplx(amplitud, numc.productoCplx((norma1, 0), (norma2, 0)))
    return a


"""Funcion que calcula la probabilidad, se implementa para sacar la probabilidad
de la amplitud al transitar de un vector a otro pero se puede usar cuando sea prudente"""


def probabilidad(c):
    return numc.moduloCplx(c) ** 2


"""Funcion que calcula el valor esperado (media) y la varianza"""


def valorEsperadoYVarianza(matriz, v):
    if matrivec.esHermitiana(matriz):
        bra = matrivec.adjunta(matrivec.productoMatrices(matriz, v))[0]
        for i in range(len(bra)):
            bra[i] = [(bra[i])]
        mu = matrivec.productoMatrices(matrivec.transpuesta(bra), v)[0][0]
        x = matrivec.multiplicacionEscalar_M(mu, matrivec.matrizIdentidad(len(matriz)))
        y = matrivec.adjunta(matriz, matrivec.inversoAditivo_M(x))
        cuadrado = matrivec.productoMatrices(y, y)
        var = matrivec.productoMatrices(matrivec.productoMatrices(matrivec.adjunta(v), cuadrado), v)[0][0]
        a = mu, var
    else:
        a = "Matriz no hermitiana"
    return a


"""Funcion que calcula los vectores normalizados y valores propios de un observable.
El observable se escribe con de fila a fila y el numero imaginario debe ir con un multiplo y j"""


def valoresPropios(observable):
    matriz = np.array(observable)
    propios = np.linalg.eig(matriz)
    valores, vectores = [], []
    for i in range(len(propios[0])):
        valores += [(propios[0][i].real, propios[0][i].imag)]

    for j in range(len(propios[1])):
        vect = []
        for k in range(len(propios[1][j])):
            vect += [[(propios[1][j][k].real, propios[1][j][k].imag)]]
        vectores += [vect]
    vectNormalizado = []
    for h in range(len(vectores)):
        vectNormalizado += normalizarVector(vectores[h])
    return valores, vectNormalizado


"""Funcion que normaliza un vector"""


def normalizarVector(vector):
    vectorNormalizado, v = [], []
    norma = matrivec.normaVector(vector)
    for i in range(len(vector)):
        v += [numc.productoCplx((1 / norma, 0), vector[i][0])]
    vectorNormalizado += [v]
    return vectorNormalizado


"""Funcion que calcula la probabilidad de que el sistema transite a un vector propio"""


def probabilidadPropios(vectorEstado, vectorPropio):
    prob = []
    bra = matrivec.adjunta(vectorEstado)[0]
    for i in range(len(bra)):
        bra[i] = [(bra[i])]
    for j in range(len(vectorPropio)):
        prob += [((numc.moduloCplx(matrivec.productoMatrices(matrivec.transpuesta(bra), vectorPropio)[0][0])) ** 2, 0)]
    return prob


"""Funcion que calcula el valor medio usando la probabilidad y los valores propios"""


def meanValue(probabilidad, valoresPropios):
    m, a = [], (0, 0)
    for i in range(len(probabilidad)):
        m += [numc.productoCplx(probabilidad[i], valoresPropios[i])]
    for j in range(len(m)):
        a = numc.sumaCplx(m[i], a)
    return a


"""Funcion para graficar las probabilidades"""


def grafica(v):
    prob = [i for i in range(len(v))]
    vector = [v[i][0] for i in range(len(v))]
    plt.bar(prob, vector, color="#FF0080")
    plt.show()


"""Funcion que calcula el estado final de una matriz unitaria (dinamica)"""


def dinamica(steps, matriz, estadoInicial):
    if matrivec.esUnitaria(matriz):
        if steps == 1:
            a = matrivec.productoMatrices(matriz, estadoInicial)[0]
        else:
            for i in range(steps):
                matriz = matrivec.productoMatrices(matriz, matriz)
            a = matrivec.productoMatrices(matriz, estadoInicial)[0]
    return a