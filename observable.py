import Matrices_Vectores as matrivec
import numeros_complejos as numc
import probabilidad as pr

def simulacion4_1():
    v = [[(-3, -1)], [(0, -2)], [(0, 1)], [(2, 0)]]
    v1 = [[(1, 0)], [(0, -1)]]
    v2 = [[(0, 1)], [(1, 0)]]
    prob = pr.probabilidadPosicion(2, v)
    ampl = pr.amplitudTransicion(v1, v2)
    print("Probabilidad de encontrar la particula en la posicion 2: ", prob)
    print("Amplitud transitar de 1 a 2: ", numc.redondear(ampl, 2))
    print("Probabilidad de transitar de 1 a 2: ", pr.probabilidad(ampl))


def ejercicio4_3_1():
    s = [[0, 1], [1, 0]]
    x = [[(0, 0), (1, 0)], [(1, 0), (0, 0)]]
    inicial = [(1, 0), (0, 0)]
    valores, vectores = pr.valoresPropios(s)
    final = matrivec.accionMatrizVector(x, inicial)
    for i in range(len(final)):
        final[i] = [final[i]]
    print("Valores Propios: ", valores)
    print("Vectores Propios Normalizados: ", vectores)
    print("Estado Final: ", final)


def ejercicio4_3_2():
    s = [[0, 1], [1, 0]]
    inicial = [[(1, 0)], [(0, 0)]]
    valores, vectores = pr.valoresPropios(s)
    probabilidad = pr.probabilidadPropios(inicial, vectores)
    m = pr.meanValue(probabilidad, valores)
    print("Probabilidad: ", probabilidad)
    print("Valor medio: ", m)
    pr.grafica(probabilidad)


def ejercicio4_4_1():
    v1 = [[(0, 0), (1, 0)], [(1, 0), (0, 0)]]
    v2 = [[((2 ** 0.5) / 2, 0), ((2 ** 0.5) / 2, 0)], [((2 ** 0.5) / 2, 0), (-(2 ** 0.5) / 2, 0)]]
    if matrivec.esUnitaria(v1) and matrivec.esUnitaria(v2):
        productov1v2 = matrivec.productoMatrices(v1, v2)
        productov2v1 = matrivec.productoMatrices(v2, v1)
        if matrivec.esUnitaria(productov1v2) and matrivec.esUnitaria(productov2v1):
            print("v1, v2 y el producto entre ambos son matrices unitarias")
    else:
        print("No son unitarias")


def ejercicio4_4_2():
    mapaUnitario = [[(0, 0), (1 / (2 ** 0.5), 0), (1 / (2 ** 0.5), 0), (0, 0)],
                    [(0, 1 / (2 ** 0.5)), (0, 0), (0, 0), (1 / (2 ** 0.5), 0)],
                    [(1 / (2 ** 0.5), 0), (0, 0), (0, 0), (0, 1 / (2 ** 0.5))],
                    [(0, 0), (1 / (2 ** 0.5), 0), (-1 / (2 ** 0.5), 0), (0, 0)]]
    inicial = [[(1, 0)], [(0, 0)], [(0, 0)], [(0, 0)]]
    if cM.matrixUnitaria(mapaUnitario):
        dinamica = pr.dinamica(3, mapaUnitario, inicial)
        for i in range(len(dinamica)):
            dinamica[i] = [dinamica[i]]
        print("Estado final despues de 3 time steps: ", dinamica)
        print("Probabilidad de encontar la bola cuantica en la posicion 3: ", pr.probabilidadPosicion(3, dinamica))
    else:
        print("Mapa no unitario")


def main():
    simulacion4_1()
    ejercicio4_3_1()
    ejercicio4_3_2()
    ejercicio4_4_1()
    ejercicio4_4_2()


main()