# -*- coding: latin-1 -*-
"""
Clases para facilitar la entrada de datos en forma de tablas y listas de 
valores.
"""

from numarray import *
from matematicas import *

class Tabla:
    def __init__(self, v):
        """
        Crea una tabla 
    
        v: lista de elementos. La primera fila  contiene las coordenadas x de 
           cada columna menos el primer elemento de la fila, cuyo valor se 
           ignora, y el primer elemento de cada fila la coordenada de la fila. 
           Por ejemplo:
           
           sum123 = Tabla([
            [ None,     1,  2,  3 ],

            [ 1,        2,  3,  4 ],
            [ 2,        3,  4,  5 ],
            [ 3,        4,  5,  6 ]])

        """

        self.interpolador =  InterpoladorLineal( 
                [ v[0][1:], [ x[0] for x in v[1:] ] ],
                transpose(array([ x[1:] for x in v[1:]])))

    def __call__(self, x, y):
        """
        Calcula el valor interpolado en el punto x, y.
        """

        return self.interpolador([x,y])

class Lista:
    def __init__(self, *v):
        """
        Crea una lista a partir de un numero variable de arguementos 
        con el formato, por ejemplo:
        senos = Lista(
                -pi,    0,
                -pi/2, -1,
                 0,     0,
                 pi/2,  1,
                 pi,    0
            )
        Devuelve una instancia que se puede llamar y que utiliza 
        interpolación mediante splines cúbicas para calcular los 
        valores.
        """

        x = []
        y = []
        for i in range(len(v)/2):
            x.append(v[2*i])
            y.append(v[2*i + 1])

        self.interpolador = InterpoladorSpline(x, y)
    
    def __call__(self, x):
        """
        Calcula el valor interpolado en el punto x.
        """

        return self.interpolador(x)
