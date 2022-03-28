import numeros_complejos
import math

def sumaVectores(v1, v2):

    result = []
    longitud = len(v1)

    for c in range(longitud):
        c1 = v1[c]
        c2 = v2[c]
        result.append((numeros_complejos.sumaCplx(c1, c2)))

    return result


def inversoAditivo(v):

    result = []
    longitud = len(v)

    for i in range(longitud):
        c1 = -1 * v[i][0]
        c2 = -1 * v[i][1]

        result.append((c1, c2))

    return result


def multiplicacionEscalar(escalar, v):

    result = []
    longitud = len(v)

    for i in range(longitud):
        producto = numeros_complejos.productoCplx(escalar, v[i])
        result.append(producto)

    return result


def adicionMatrices(m1, m2):

    result = []
    m = len(m1)

    for i in range(m):
        suma = sumaVectores(m1[i], m2[i])
        result.append(suma)

    return result


def inversoAditivo_M(m):

    result = []
    filas = len(m)

    for i in range(filas):
        result.append(inversoAditivo(m[i]))

    return result


def multiplicacionEscalar_M(escalar, matriz):

    filas = len(matriz)
    result = []

    for i in range(filas):
        producto = multiplicacionEscalar(escalar, matriz[i])
        result.append(producto)

    return result


def transpuesta(matriz):

    filas = len(matriz)
    columnas = len(matriz[0])

    if filas == columnas:
        result = [[0 for j in range(columnas)] for i in range(filas)]

    else:
        result = [[0 for j in range(filas)] for i in range(columnas)]

    for i in range(filas):
        for j in range(columnas):
            result[j][i] = matriz[i][j]

    return result


def conjugada(matriz):

    filas = len(matriz)
    columnas = len(matriz[0])

    if type(matriz[0]) is tuple:
        matriz2 = []

        for i in range(filas):
            matriz2.append(numeros_complejos.conjugadoCplx(matriz[i]))

    else:
        matriz2 = [[0 for j in range(columnas)] for i in range(filas)]

        for i in range(filas):
            for j in range(columnas):
                matriz2[i][j] = numeros_complejos.conjugadoCplx(matriz[i][j])

    return matriz2


def adjunta(matriz):

    result = conjugada(matriz)
    result = transpuesta(result)

    return result


def productoMatricesComp(m1, m2):

    filasM1 = len(m1)
    filasM2 = len(m2)
    columnasM1 = len(m1[0])
    columnasM2 = len(m2[0])

    assert columnasM1 == filasM2, "Las matrices no son compatibles para hacer un producto entre ellas"

    resultado = [[0 for j in range(columnasM2)] for i in range(filasM1)]

    for i in range(filasM1):
        for j in range(columnasM2):
            suma = (0, 0)
            for k in range(filasM2):
                suma = numeros_complejos.sumaCplx(suma, numeros_complejos.productoCplx(m1[i][k], m2[k][j]))
            resultado[i][j] = suma

    return resultado


def accionMatrizVector(matriz, vector):

    filasM = len(matriz)
    result = [[0 for j in range(1)] for i in range(filasM)]

    for i in range(filasM):
        for j in range(1):
            suma = (0, 0)
            for k in range(len(vector)):
                suma = numeros_complejos.sumaCplx(suma, numeros_complejos.productoCplx(matriz[i][k], vector[k]))
            result[i][j] = suma

    return result


def productoInternoVectores(v1, v2):

    v3 = conjugada(v1)

    suma = (0, 0)
    for i in range(len(v2)):
        producto = numeros_complejos.productoCplx(v3[i], v2[i])
        suma = numeros_complejos.sumaCplx(suma, producto)

    return suma


def normaVector(vector):

    result = 0
    for i in vector:
        result += i[0] ** 2 + i[1] ** 2
    return math.sqrt(result)


def distanciaVectores(v1, v2):

    result = 0
    for i in range(len(v1)):
        result += (v1[i][0] - v2[i][0]) ** 2 + (v1[i][1] - v2[i][1]) ** 2
    return math.sqrt(result)


def matrizIdentidad(long):

    mUnitaria = [[(0, 0) for j in range(long)] for i in range(long)]
    for i in range(long):
        mUnitaria[i][i] = (1, 0)
    return mUnitaria


def esUnitaria(a):

    if productoMatrices(a, adjunta(a)) == matrizIdentidad(len(a)):
        ans = True
    else:
        ans = False
    return ans


def esHermitiana(matriz):

    return matriz == adjunta(matriz)


def productoTensor(a, b):

    result = [[[[]] for j in range(len(a[0]) * len(b[0]))] for i in range(len(a) * len(b))]

    for i in range(len(a) * len(b)):
        for j in range(len(a[0]) * len(b[0])):
            x, y = i // len(b), j // len(b[0])
            res = multiplicacionEscalar_M(a[x][y], b)
            x1, y1 = i % len(b), j % len(b[0])
            result[i][j] = res[x1][y1]
    return result


def main():

    c1 = (-1 / math.sqrt(12), 1 / math.sqrt(12))
    print(numeros_complejos.moduloCplx(c1) ** 2)


"""
Funcion para obtener la transpuesta de una matriz o un vector
"""


def transpuesta(matriz):

    if isinstance(matriz[0], int):
        resultado = [[0 for j in range(1)] for i in range(len(matriz))]

        for i in range(len(resultado)):
            for j in range(1):
                resultado[i][j] = matriz[i]

        return resultado

    elif len(matriz[0]) == 1:
        resultado = []
        for i in range(len(matriz)):
            for j in range(1):
                resultado.append(matriz[i][j])
        return resultado

    else:
        filas = len(matriz)
        columnas = len(matriz[0])
        if filas == columnas:
            result = [[0 for j in range(columnas)] for i in range(filas)]

        else:
            result = [[0 for j in range(filas)] for i in range(columnas)]

        for i in range(filas):
            for j in range(columnas):
                result[j][i] = matriz[i][j]

        return result

"""
Funcion para obtener el producto entre dos matrices de tamaños compatibles
"""


def productoMatrices(m1, m2):

    filasM1 = len(m1)
    filasM2 = len(m2)
    columnasM1 = len(m1[0])
    columnasM2 = len(m2[0])

    assert columnasM1 == filasM2, "Las matrices no son compatibles para hacer un producto entre ellas"

    resultado = [[0 for j in range(columnasM2)] for i in range(filasM1)]

    for i in range(filasM1):
        for j in range(columnasM2):
            suma = 0
            for k in range(filasM2):
                suma += (m1[i][k] * m2[k][j])
            resultado[i][j] = suma

    return resultado


"""
Funcion para obtener el resultado de multiplicar 
"""


def escalarXMatriz(escalar, matriz):

    resultado = [[0 for j in range(len(matriz[0]))] for i in range(len(matriz))]

    for i in range(len(resultado)):
        for j in range(len(resultado[0])):
            resultado[i][j] = escalar * matriz[i][j]

    return resultado


"""
Funciones para crear la matrices probabilistica y cuantica
"""


def crearMatrizProb(nRendijas, nBlancos):

    total = nRendijas + nBlancos + 1
    probRendija = 1 / nRendijas
    probBlanco = 1 / 3
    matrizResultado = [[0 for j in range(total)] for i in range(total)]
    cont = 0
    contR = 1
    col = 1
    inicio = nRendijas + 1

    for i in range(1, nRendijas + 1):
        matrizResultado[i][0] = round(probRendija, 3)

    for i in range(nRendijas + 1, nRendijas + nBlancos + 1):
        for j in range(1, nRendijas + 1):
            if j == col:
                matrizResultado[i][j] = round(probBlanco, 3)
                cont += 1

            if cont == 3 and contR < nRendijas:
                col += 1
                cont = 0
                contR += 1

    # va sumando 1  en la posicion [i][i] despues de el numero de rendijas
    for i in range(total):
        for j in range(total):
            if i == inicio and j == inicio:
                matrizResultado[i][j] = 1
                inicio += 1

    return matrizResultado


def crearMatrizCuant(nRendijas, nBlancos):

    total = nRendijas + nBlancos + 1
    probRendija = (1 / round(math.sqrt(nRendijas), 3), 0)
    probBlanco = [(-1 / round(math.sqrt(6), 3), 1 / round(math.sqrt(6), 3)),
                  (-1 / round(math.sqrt(6), 3), -1 / round(math.sqrt(6), 3)),
                  (1 / round(math.sqrt(6), 3), -1 / round(math.sqrt(6), 3))]
    matrizResultado = [[(0, 0) for j in range(total)] for i in range(total)]
    cont = 0
    contR = 1
    contB = 0
    col = 1
    inicio = nRendijas + 1

    for i in range(1, nRendijas + 1):
        matrizResultado[i][0] = probRendija

    for i in range(nRendijas + 1, nRendijas + nBlancos + 1):
        for j in range(1, nRendijas + 1):
            if j == col:
                matrizResultado[i][j] = probBlanco[contB]
                contB += 1
                cont += 1

            if cont == 3 and contR < nRendijas:
                col += 1
                cont = 0
                contB = 0
                contR += 1

    for i in range(total):
        for j in range(total):
            if i == inicio and j == inicio:
                matrizResultado[i][j] = (1, 0)
                inicio += 1

    return matrizResultado


def imprimirMatriz(matriz):

    for i in matriz:
        print(i)
    print()


def main():
    m = crearMatrizProb(3, 5)
    imprimirMatriz(m)


if __name__ == '__main__':
    main()



