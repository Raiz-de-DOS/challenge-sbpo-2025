
import pandas as pd
from pyscipopt import Model
import numpy as np

TOLERANCE = 10e-6
inf = 10e7
epsilon = 10e-3

class Instancia():
    def __init__(self):
        self.cant_ordenes = None
        self.cant_items = None
        self.cant_pasillos = None
        self.lb = None
        self.ub = None
        self.cant_elem_por_orden = [[]]
        self.cant_elem_por_pasillo = [[]]

    def leer_datos(self, nombre_archivo):
    #LEO LOS DATOS DE ALGUNA MANERA DEL ARCHIVO DE TEXTO

def cargar_instancia():

    nombre_archivo = "nombre.lp"#ahora lo pongo

    # Crea la instancia vacía
    instancia = Instancia()

    # Llena la instancia con los datos del archivo de entrada
    instancia.leer_datos(nombre_archivo)

    return instancia

def planteo_pasillos_fijos(prob, instancia, cant_pasillos_disponibles):

    dict_W = {}
    dict_A = {}
    a_star = cant_pasillos_disponibles

    #Si la orden pertenece a la wave
    for o in range (instancia.cant_ordenes):
        dict_W[o] = prob.addVar(vtype = 'B', name = f"W_{o}", lb = 0, ub = 1)

    #Si se usa el pasillo a
    for a in range(instancia.cant_pasillos):
        dict_A[a] = prob.addVar(vtype = 'B', name = f"A_{a}", lb=0, ub=1)

    u_orden = instancia.cant_elem_por_orden
    u_pasillos = instancia.cant_elem_por_pasillo

    #La cantidad de elementos de la wave está en rango
    for o in range(instancia.cant_ordenes):
        suma_u_o_i = 0
        for i in range(instancia.cant_items):
            suma_u_o_i += u_orden[o][i] * dict_W[o]
        prob.addCons(suma_u_o_i >= instancia.lb)
        prob.addCons(suma_u_o_i <= instancia.ub)

    #Todos los elementos que tomo de un pasillo tienen stock
    for i in range(instancia.cant_items):
        izq, der = 0, 0
        for o in range(instancia.cant_ordenes):
            izq += u_orden * dict_W[o]
        for a in range(instancia.cant_pasillos):
            der += u_pasillos[i][a] * dict_A[a]
        prob.addCons(izq <= der)

    #La cantidad de pasillos usados es A* (pasado por parámetro)
    suma_a = 0
    for a in range(instancia.cant_pasillos):
        suma_a += dict_A[a]
    prob.addCons(suma_a == a_star)

    prob.setObjective(suma_u_o_i, "maximize")

    return dict_W, dict_A
def planteo_busqueda_binaria(prob, instancia, k):

#Sea k en [lb/A ; ub] una constante fija que es el resultado del cociente (fijo)
#por lo que ahora la función objetivo es minimizar una diferencia que correspondía a la restricción (ver modelo)

    dict_W = {}
    dict_A = {}

    # Si la orden pertenece a la wave
    for o in range(instancia.cant_ordenes):
        dict_W[o] = prob.addVar(vtype='B', name=f"W_{o}", lb=0, ub=1)

    # Si se usa el pasillo a
    for a in range(instancia.cant_pasillos):
        dict_A[a] = prob.addVar(vtype='B', name=f"A_{a}", lb=0, ub=1)

    u_orden = instancia.cant_elem_por_orden
    u_pasillos = instancia.cant_elem_por_pasillo

    # La cantidad de elementos de la wave está en rango
    for o in range(instancia.cant_ordenes):
        suma_u_o_i = 0
        for i in range(instancia.cant_items):
            suma_u_o_i += u_orden[o][i] * dict_W[o]
        prob.addCons(suma_u_o_i >= instancia.lb)
        prob.addCons(suma_u_o_i <= instancia.ub)

    #Restricción de la forma f = k * g
    pasillos_utilizados = sum(dict_A[a] for a in range(instancia.cant_pasillos))
    prob.addCons(k * pasillos_utilizados - epsilon <= suma_u_o_i)
    prob.addCons(k * pasillos_utilizados + epsilon >= suma_u_o_i)

    #Función objetivo
    for i in range(instancia.cant_items):
        izq, der = 0, 0
        for o in range(instancia.cant_ordenes):
            izq += u_orden * dict_W[o]
        for a in range(instancia.cant_pasillos):
            der += u_pasillos[i][a] * dict_A[a]

    prob.setObjective(izq - der, "minimize")

    return dict_W, dict_A

def armar_lp_pasillos_fijos(prob, instancia, a_actual):
    # Agregar las variables
    planteo_pasillos_fijos(prob, instancia, a_actual)

    # Escribir el lp a archivo
    prob.writeProblem()

def armar_lp_busqueda_binaria(prob, instancia, k_actual):
    # Agregar las variables
    planteo_busqueda_binaria(prob, instancia, k_actual)

    # Escribir el lp a archivo
    prob.writeProblem()

def resolver_lp(prob):
    # Definir los parametros del solver
    prob.setParam("limits/gap", TOLERANCE)
    # Resolver el lp
    prob.optimize()
    print(f"Estado de optimización: {prob.getStatus()}")  # Debugging


def mostrar_resultados(prob, instancia, dict_W, dict_A):

    n = sum(dict_W[o] for o in instancia.cant_ordenes)
    print(n)
    for j in range(n):
        if dict_W[j]==1:
            print(j)
            print("\n")

    m = sum(dict_A[a] for a in instancia.cant_pasillos)
    print(m)
    for s in range(m):
        if dict_A[s]==1:
            print(s)
            print("\n")

def main():

    instancia = cargar_instancia() #Lectura de datos desde el archivo de entrada
    prob: Model = Model() #Definición del problema

    # DEFINICIÓN DEL MODELO Y SOLVER
    # Decidimos qué modelo vamos a usar (el que tenga que resolver menos PL
    rango_k = np.log(instancia.ub - (instancia.lb / instancia.cant_pasillos)) + np.log(epsilon**-1)
    if (instancia.cant_pasillos <= rango_k):
        maximo = -inf
        for a_actual in range(1, instancia.cant_pasillos + 1):  # incluyo ese máximo
            dict_W, dict_A = planteo_pasillos_fijos(prob, instancia, a_actual)
            resolver_lp(prob)
            valor_objetivo_actual = prob.getObjVal() /a_actual
            if (maximo < valor_objetivo_actual):
                maximo = valor_objetivo_actual
                dict_W_sol = dict_W
                dict_A_sol = dict_A
                valor_objetivo = valor_objetivo_actual
            dict_W.clear()
            dict_A.clear()
            prob.resetParams()
        print(valor_objetivo)
    else:
        l = instancia.lb, u = instancia.ub
        dict_W, dict_A = {} , {}
        while(u - l > epsilon):
            dict_W.clear()
            dict_A.clear()
            prob.resetParams()
            k = (u-l)/2
            dict_W, dict_A = planteo_busqueda_binaria(prob, instancia, k)
            resolver_lp(prob)
            valor_objetivo_actual = prob.getObjVal()
            if (valor_objetivo_actual > 0): #Era infactible con ese valor de k alto
                u = k
            else:
                l = k
        dict_W_sol = dict_W
        dict_A_sol = dict_A
        print(prob.getObjVal())

    # Mostrar los resultados
    mostrar_resultados(prob, instancia, dict_W_sol, dict_A_sol)

if __name__ == '__main__':
    main()
