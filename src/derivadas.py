#!/usr/bin/env python
# -*- coding: latin 1 -*-

# Las derivadas aqui calculadas siguen el convenio del Padfield:
# Las fuerzas estan divididas por la masa
#   X' = X/M
#   Y' = Y/M
#   Z' = Z/M
# Los momentos por las inercias:
#                 -1
#     { L' }   [  Ixx   0   -Ixz ]  { L }   
#     { M' } = [   0   Iyy    0  ]  { M } = 
#     { N' }   [ -Ixz   0    Izz ]  { N }   
#              [ Izz/(IxxIzz - Ixz^2)     0     Ixz/(Ixx*Izz - Ixz^2) ] { L }
#              [           0        1/Iyy                         ]     { M }
#              [ Ixz/(IxxIzz - Ixz^2)     0     Ixx/(Ixz*Izz - Ixz^2) ] { N }


# Cargamos el modelo del helicoptero
import sys
import imp
nombre_modelo = sys.argv[1]
try:
    (file_modelo, path_modelo, desc_modelo) = imp.find_module(nombre_modelo)
except ImportError:
    sys.exit("Error: no se encontro el modelo")
modulo = imp.load_module(nombre_modelo, file_modelo, path_modelo, desc_modelo)
modelo = modulo.modelo
modelo.init()

modelo.viento_vel = 0
modelo.viento_dir = 0
modelo.P = 99461.0575
modelo.T = 283
modelo.motor.nmotores = 2
modelo.x_i = modelo.y_i = 0
modelo.z_i = 1000

from matematicas import *
from numarray import *
import numarray.linear_algebra as la

def extract(dic, lis):
    """ 
    Devuelve una lista con el valor que se encuentra en dic
    para cada elemento de lis
    """
    r = []
    for l in lis:
        r.append(dic[l])
    return r

class Derivadas:
    def __init__(self):
        # Variables a derivar
        self.a = ['X', 'Y', 'Z', 'L', 'M', 'N', 'Dq0', 
                'Dq1', 'Dq2', 'Dq3', 'DOm', 'DT2Q1', 'DTQ1']
        # Variables que derivan
        self.b = ['u', 'v', 'w', 'p', 'q', 'r', 'q0', 'q1', 'q2', 'q3', 
                'Om', 'DTQ1', 'Q1', 'th0', 'th1c', 'th1s', 'th0T']
        # Velocidades
        self.V = []
        # Las derivadas: la derivada X_u para todas las velocidades es una
        # lista que se encuentra en d['X'][u]
        self.d = {}
        for x in self.a:
            self.d[x] = {}
            for y in self.b:
                self.d[x][y] = []
    
    def asarray(self):
        """
        Devuelve el valor de las derivadas de estabilidad en forma de array:
            El primer indica indica la velocidad.
            El segundo indice la variable derivada.
            El tercer indice la variable que deriva.
        Por ejemplo:
            r = d_Modelo.asarray()
            print r[10, 2, 3]

        Como:
            d_Modelo.a[2] = 'Z'
            d_Modelo.b[3] = 'p'
        Suponiendo que fuese por ejemplo:
            d_Modelo.V[10] = 15

        Entonces mostraria por pantalla el valor de Z_p para la velocidad
        de 15 m/s
        """

        v = len(self.V)
        r = array(shape = (v, len(self.a), len(self.b)), type='Float64')
        for k in range(v):
            for i in range(len(self.a)):
                for j in range(len(self.b)):
                    r[k, i, j] = self.d[self.a[i]][self.b[j]][k]
        return r

    def append_b(self, n, l):
        """
        l es una lista que contiene las derivadas de todas las variables
        respecto a n, para una nueva velocidad. Esta funcion anade esta
        lista.
        """
        c = 0
        for a in self.a:
            self.d[a][n].append(l[c])
            c += 1
    
    def eigen(self, v):
        """ 
        Calcula los autovalores para la velocidad v
        """
        r = array(shape = (len(self.b)-4, len(self.a)), type='Float64')
        for i in range(len(self.a)):
            for j in range(len(self.b)):
                if self.b[j] not in ['th0', 'th1c', 'th1s', 'th0T']:
                    r[j, i] = self.d[self.a[i]][self.b[j]][v]
        return la.eigenvalues(r)
    

d_Modelo = Derivadas()

# Funcion a derivar
def f(x):
    modelo.trim = True
    modelo.controles.th0 =  x[16]
    modelo.controles.th1c = x[17]
    modelo.controles.th1s = x[18]
    modelo.controles.th0T = x[19]
    return modelo.paso(x[:16], 0)

# Velocidad en nudos, 1kn = 0.51444 m/s
V = concatenate([arange(0,2,0.2), arange(2, 162, 2.)])      
for x in V:
    # Trimamos el modelo
    tr = modelo.trimado(x*0.51444, 0.0, 0.0, 0.0, 1000, 1.225)
    x0 = extract(tr,[
         "u", "v", "w", "p", "q", "r",
         "q0", "q1", "q2", "q3",
         "Om", "DTQ1", "Q1"
         ])
    x0.extend([0, 0.0, 1000.0])
    x0.extend([tr['th0'], tr['th1c'], tr['th1s'], tr['th0T']])

    # Evitamos que actuen los mandos
    modelo.trim = True
    
    # Empezamos a calcular derivadas
    print "V= ", x
    d_Modelo.V.append(x)
    for i in range(len(d_Modelo.b)):
        if d_Modelo.b[i] in ['th0', 'th1c', 'th1s', 'th0T']:
            k = i + 3
        else:
            k = i
        d_Modelo.append_b(d_Modelo.b[i], dpart_1(f, x0, k))

def pop_prox(x, l):
    """ 
    Extrae de la lista l el elemento que se encuentre mas proximo
    a x.
    """
    def prox(a, b):
        return abs(b-a)
    c = 0
    p = prox(x, l[0])
    for i in range(1, len(l)):
        if prox(x, l[i])<p:
            c = i
            p = prox(x, l[i])
    return l.pop(c)

eigen_x = [[] for i in range(13)]
eigen_y = [[] for i in range(13)]

e0 = d_Modelo.eigen(0).tolist()
for v in range(len(e0)):
    eigen_x[v].append(complex(e0[v]).real)
    eigen_y[v].append(complex(e0[v]).imag)

for x in range(1, len(V)):
    e = d_Modelo.eigen(x).tolist()
    for v in range(len(e0)):
        b = pop_prox(e0[v], e)
        eigen_x[v].append(complex(b).real)
        eigen_y[v].append(complex(b).imag)
    e0 = [ eigen_x[i][x] + eigen_y[i][x]*1j for i in range(len(eigen_x)) ]

from pylab import *
rc('figure', figsize=(4.72, 3.15))
rc('text', usetex=True)

c = ['b', 'g', 'r', 'c', 'm', 'y', 'k' ]
c *= 3
for i in (5, 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12):
   plot([eigen_x[i][0]], [eigen_y[i][0]], c[i] + 'o')
   plot([eigen_x[i][-1]], [eigen_y[i][-1]], c[i] + '^')
   plot(eigen_x[i], eigen_y[i], c[i] + '-')
show()

def latex(s):
    """ 
    Transforma la cadena s en una cadena latex
    """
    try:
        r = {
                'th0':  r'\theta_0',
                'th1c': r'\theta_{1c}',
                'th1s': r'\theta_{1s}',
                'th0T': r'\theta_{0_T}'
            }[s]
    except KeyError:
        r = s
    return r

for a in ['X', 'Y', 'Z', 'L', 'M', 'N']:
    for b in ['u', 'v', 'w', 'p', 'q', 'r', 'th0', 'th1s', 'th1c', 'th0T']:
        plot(d_Modelo.V, d_Modelo.d[a][b])
        ylabel(r'$%s_{%s}$' % (a, latex(b)))
        savefig('%s_%s_%s.eps' % (modelo.nombre, a,b), orientation='landscape')
        close()
        
        print a,b
            
