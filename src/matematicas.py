# -*- coding: latin-1 -*-
"""
Funciones matematicas de todo tipo
"""

from math import asin, atan
import cmath
from numarray import *
import numarray.linear_algebra as la
import copy
import logging
import pylab

################################################################################
#               UTILIDADES                                                     #
################################################################################

def escalon(a=0.0, b=1.0, A = 0.0, B=1.0):
    """
    Devuelve una funcion escalon:
     B |   ........     
       |   .      .     
     A |....      ........
    --------------------------------
           a      b     
    
    Si b = None --> b = inf
    Si a = None --> a = inf

    Valores por defecto:
        a = A = 0.0
        b = B = 1.0
    """

    if a == None:
        def r(t):
            if t<=b:
                return B
            else:
                return A
    elif b == None:
        def r(t):
            if t>=a:
                return B
            else:
                return A
    else:
        def r(t):
            if t>=a and t<=b:
                return B
            else:
                return A
    return r

def sign(x):
    """ 
    Devuelve el signo de x (0 se interpreta con signo positivo)  
        x>=0 --> +1
        x <0 --> -1
    """

    if x<0:
        return -1
    elif x>=0:
        return +1

def cross(a, b):
    """ 
    Calcula el producto vectorial de dos vectores. 
    Los vectores pueden ser cualquier elemento iterable, pero el resultado 
    es un array.
    """

    return array([
         a[1]*b[2] - a[2]*b[1], 
         a[2]*b[0] - a[0]*b[2], 
         a[0]*b[1] - a[1]*b[0] 
    ])

def mod(x):
    """ 
    Devuelve el modulo del vector x
    """
    s = 0
    for c in x:
        s += c**2
    return sqrt(s)

def rlen(x):
    """ 
    Alias para range(len(x)) 
    """

    return range(len(x))

def reverse(x):
    """ 
    Ordena al reves una lista. 
    """

    y = []
    for i in rlen(x):
        y.append(x[-1 - i])
    return y

def rad(x):
    """ 
    Convierte de grados a radianes 
    """

    return x*pi/180.

def deg(x):
    """ 
    Convierte de radianes a grados 
    """

    return x*180./pi


def localiza(x, L):
    """ 
    Dado un contenedor L ordenado de menor a mayor , y un valor x, devuelve el 
    indice del array tal que L[i]<=x<=L[i+1] mediante el metodo de la 
    biseccion. Si el elemento no esta en la lista devuelve alguno de los
    extremos, de forma que el interpolador lineal se convierte automaticamente 
    en un extrapolador. 
    """

    a=0
    b=len(L)-1
    while not ( (L[a+(b-a-1)/2]<=x) and (L[a+(b-a-1)/2+1]>=x) ):
        if x<L[a+(b-a-1)/2]:
            b=a+(b-a-1)/2
        else:
            a=a+(b-a-1)/2+1
        if a==len(L)-1:
            return len(L)-2
        elif b==0:
            return 0

    return a+(b-a-1)/2

def ang(x, y):
    """ 
    Calcula el angulo que forma el punto (x,y).
    Util porque distingue porque vale para cualquier cuadrante y porque trata
    sin problemas los casos +-pi/2.
    """

    if x == 0:
        if y>0:
            a = pi/2
        elif y<0:
            a = -pi/2
        else:
            a = 0.0
    else:
        a = atan(y/x)
   
    if  x<0:
        if y<0:
            a = a - pi
        else:
            a = a + pi
    if a>=pi:
        a -= 2*pi
    elif a<=-pi:
        a += 2*pi
    return a

################################################################################
#               SUAVIZADO                                                      #
################################################################################
def smooth(x, cparada=None, cfijo=None):
    """
    Suaviza los datos x, hasta que se cumpla una condicion de parada.
    cparada(n, x): devuelve verdadero cuando se cumpla la condicion de parada, 
    para n pasadas y datos x
    cfijo(i) : devuelve verdadero si el dato en el indice i no es modificable
    """
    
    if cfijo == None:
        def cfijo(i):
            if i<2 or i>len(x)-3:
                return True
            else:
                return False

    if cparada == None:
        def cparada(n, y):
            if n==1:
                return True
            else:
                return False

    # Arrays temporales
    y = x.copy()
    z = x.copy()
    # numero de pasadas
    n = 0
    
    # kernel gaussiano para 4 vecinos
    p = 2
    c = [0.1174, 0.1975, 0.2349, 0.1975, 0.1174]
    
    # Aplicamos sucesivas paradas mientras no se la condicion de parada
    while not cparada(n, y):
        for i in rlen(y):
            if not cfijo(i):
                z[i] = 0.0
                for k in range(-p, p +1):
                    z[i] += c[k]*y[i+k]
        y, z = z, y
        n += 1
    return y

################################################################################
#               CALCULO                                                        #
################################################################################

def solveNewton(h, ix, Nmax, eps):
    """ 
    Resuelve por el metodo de Newton. 
    En realidad aplica el metodo iterativo:
        Xn+1 = Xn + h(Xn)
    Mientras sea:
        |Xn+1 - Xn|>eps y n<Nmax
    Se supone que la funcion h(x) representa:
        h(x) = -f(x)/f'(x)
    para que el metodo calcule f(x)=0
    """

    n, error = 0, eps + 1
    x0=ix
    while error>eps:
        x1=x0 + h(x0)
        error=abs(x1-x0)
        x0=x1
        if n>Nmax:
            logging.error("Metodo de Newton no converge")
            break
        n += 1

    return x1

def solveNewtonMulti(f, ix, Nmax, eps):
    """
    Metodo de Newton para la ecuacion vectorial f(x) = 0.
    ix es el punto inicial
    Nmax maximo de iteraciones
    eps error numerico
    """
    n, error, x0 = 0, eps + 1, ix
    s = len(ix)
    J = array(shape=(s,s), type='Float64')
    while error>eps:
        for i in range(s):
                J[i,:] = dpart_1(f, x0, i, eps)
        x1 = x0 - dot(la.inverse(J), f(x0))
        error = mod(x1 - x0)
        x0 = x1
        if n>Nmax:
            logging.error("Metodo de Newton(multi) no converge")
            break
        n += 1
    return x1


def solveSecante(f, xa, xb, Nmax, eps):
    """
    Resuelve f=0 mediante el metodo de la secante.

    xa y xb son los puntos iniciales para calcular la secante, de forma ideal 
    deberia ser signo(f(xa))!=signo(f(xb)) y  xa<xb
    """
    error = eps + 1
    n = 0
    xa = float(xa)
    xb = float(xb)
    ya = f(xa)
    yb = f(xb)
    while error>eps:
        # Si tenemos la mala suerte de una horizontal 
        # probamos con otro punto
        if yb == ya:
            xb = (xb + xa)/2.
            yb = f(xb)
            continue
        
        # Interseccion con cero de la recta que une A y B
        x = (xa*yb - ya*xb)/(yb - ya)
        y = f(x)
        
        # A lo mejor ha habido suerte
        if x == xa or x == xb:
            return x
        
        # Decidimos que punto se marcha, si A o B
        # Supongamos que se marcha B
        if y!=ya: 
            x1 = (xa*y - ya*x)/(y-ya) 
        else:
            # x1 es infinito
            x1 == None
        # Supongamos que se marcha A
        if y!=yb:
            x2 = (x*yb - y*xb)/(yb-y) 
        else:
            # x2 es infinito
            x2 = None
        
        # Alguno de los 2 puntos no es posible
        if x1 == None:
            xa = x
            ya = y
            continue
        
        if x2 == None:
            xb = x
            yb = y
            continue

        # Hay que decidir entre dos puntos
        if sign(ya)!=sign(yb):
            # La solucion esta entre xa y xb, nos aseguramos
            # de que el nuevo punto caiga en este intervalo
            # si x1==None o x2==None fallan los test de intervalo
            # como debe ser
            if xa<x1 and x1<xb:
                # 1 cae dentro
                if xa<x2 and x2<xb:
                    # Ambos caen dentro, nos quedamos con el mejor
                    y1 = f(x1)
                    y2 = f(x2)
                    if abs(y1)<abs(y2):
                        xb = x
                        yb = y
                    else:
                        xa = x
                        ya = y
                else:
                    # 1 cae dentro, pero 2 fuera
                    xb = x
                    yb = y
            else:
                # 1 cae fuera, al menos uno tiene que caer dentro
                # por ser distinto signo ya, yb, luego 2 cae dentro
                xa = x
                ya = y
        else:
            # La solucion no tiene por que estar dentro del intervalo
            # Nos limitamos a coger el mejor
            y1 = f(x1)
            y2 = f(x2)
            if abs(y1)<abs(y2):
                xb = x
                yb = y
            else:
                xa = x
                ya = y
            if xa>xb:
                xa, xb = xb, xa
                ya, yb = yb, ya

        error = abs(y)
        n += 1
        if n>Nmax:
            logging.error("Secante no converge")
            return None
    return x

def deriv_1(f, x, h=0.001):
    """ 
    Calcula la primera derivada de f en el punto x 
    """
    
    return (f(x+h) - f(x-h))/(2*h)

def dpart_1(f, x, p, h=0.001):
    """ 
    Calcula la derivada parcial de f en el punto x, en la coordenada p 
    """

    x1 = copy.deepcopy(x)
    x2 = copy.deepcopy(x)
    x1[p] -= h
    x2[p] += h
    return (f(x2) - f(x1))/(2*h)
    


################################################################################
#           AJUSTE POR MINIMOS CUADRADOS                                       #
################################################################################

def lsq_fit(x, y, f):
    """ 
    Calcula un ajuste lineal del tipo: 
        a0*f[0](x) + a1*f[1](x) + ... + aM-1*f[M-1](x)
    para los N datos "y" definidos en las N posiciones "x". 
    Devuelve los coeficientes "a" y la desviacion estimadas en cada punto.
    """

    M = len(f)
    b = array(shape=(M), type='Float64')
    c = array(shape=(M, M), type='Float64')
    
    for j in range(M):
        for k in range(M):
            c[k, j] = sum([dot(f[j](x[i]),f[k](x[i])) for i in rlen(x)])
        b[j] = sum([dot(y[i], f[j](x[i])) for i in rlen(x)])
    d = la.inverse(c)
    return dot(d, b), map(sqrt, diagonal(d))

def f_pol(n):
    """ 
    Devuelve una funcion polinomica de grado n 
    """
    
    return (lambda x: x**n)

def f_cos(n):
    """ 
    Devuelve una funcion coseno de frecuencia n 
    """

    return (lambda x: cos(n*x))

def f_sin(n):
    """ 
    Devuelve una funcion seno de frecuencia n 
    """

    return (lambda x: sin(n*x))

def f_lin(c, f):
    """ 
    Devuelve la funcion c[i]*f[i] 
    """
    
    def r(x):
        return sum([c[i]*f[i](x) for i in rlen(f)])
    
    return r

################################################################################
#           INTERPOLADORES                                                     #
################################################################################

class InterpoladorLineal:
    """ 
    Interpolador lineal para tablas multidimensionales 
    """
    
    def __init__(self, x, y):
        # nos aseguramos de que los tipos son correctos
        if not hasattr(x[0], '__getitem__'):
            self.x = [x]
        else:
            self.x = x
        self.y = array(y)
        
        # la dimension de cada punto de la tabla
        self.dim_x = len(self.x)
    
    def __call__(self, ix):
        """ 
        Interpola el valor en el punto ix 
        """
        
        if not  hasattr(ix, '__getitem__'):
            x = [ix]
        else:
            x = ix

        y0 = self.y
        for d in range(self.dim_x):
            i = localiza(x[d], self.x[d])
            x0, x1, x2 = x[d], self.x[d][i], self.x[d][i+1]
            y1 = ((x2 - x0)*y0[i] + (x0 -x1)*y0[i+1])/(x2-x1)
            y0 = y1
        return y1

class InterpoladorSpline:
    """ 
    Interpola mediante splines cubicas datos unidimensionale. 
    """
    
    def __init__(self, x, y, tipo = 'Clamp', t0 = None, t1 = None):
        """ 
        x: vector con las coordenadas x
        y: vector con las coordenadas y
        t0: tangente en x[0],y[0]
        t1: tangente en x[n],y[n] (len(x) == len(y) == n + 1)
        """
        
        self.x = x
        self.y = y

        if t0 == None:
            t0 = (y[1] - y[0])/(x[1] - x[0])
        
        if t1 == None:
            t1 = (y[-1] - y[-2])/(x[-1] - x[-2])

        # Construimos el sistema de ecuaciones que nos determinan los 
        # coeficientes de las splines
        n = len(x) - 1
        C = zeros( shape = (n+1, n+1), type='Float64' )
        b = array( shape =(n+1,), type='Float64')
        
        self.h = h = array( [x[i+1] - x[i] for i in range(n)], type='Float64')
        
        # Condiciones en los extremos
        if tipo == 'Clamp':
            C[0,:2] = [2, 1]
            b[0] = 6/(h[0]**2)*(y[1] - y[0] - t0*h[0])
            C[n,-2:] = [1,2]
            b[n] = 6/(h[n-1]**2)*(y[n-1] - y[n] + t1*h[n-1])
        else:
            C[0,0] = 1.0
            b[0] = 0.0
            C[n,n] = 1.0
            b[n] = 0.0

        # Grueso de la matriz
        for i in range(1,n):
            C[i][i-1:i+2] = [h[i-1], 2*(h[i-1] + h[i]), h[i]]
            b[i] = 6*((y[i+1] - y[i])/h[i] + (y[i-1] - y[i])/h[i-1])

        # Salvamos los coeficientes
        self.a = la.solve_linear_equations(C, b)
        
        # Para poder consultarlo
        self.t0 = t0
        self.t1 = t1
    
    def __call__(self, ix):
        """
        Calcula el valor interpolado en el punto ix.
        """
        
        x, y = self.x, self.y
        i = localiza(ix, x)
        dx0, dx1 = ix - x[i], x[i+1] -ix
        a0, a1 = self.a[i], self.a[i+1]
        h = self.h[i]
        
        return ( 
                (a0*dx1**3 + a1*dx0**3)/(6*h) 
                + (y[i] - a0*h**2/6)*dx1/h + 
                (y[i+1] - a1*h**2/6)*dx0/h 
        )




################################################################################
#           QUATERNIOS                                                         #
################################################################################

def mul_quat(qA, qB):
    """ 
    Multiplica dos cuaternios 
    """
    
    qAv = array(qA[1:])
    qBv = array(qB[1:])
    q0 = qA[0]*qB[0] - dot(qAv, qBv)
    q1, q2, q3 = qA[0]*qBv + qB[0]*qAv + cross(qAv, qBv)
    return q0, q1, q2, q3

class Quaternion:
    def __init__(self, q0=1., q1=0., q2=0., q3=0.):
        """ 
        Inicializa el quaternio como (q0, (q1,q2,q3)) 
        q0 es la parte escalar
        q1,q2,q3 la parte vectorial
        """
        
        self.q0, self.q1, self.q2, self.q3 = q0, q1, q2, q3
    
    def escalar(self):
        """ 
        Devuelve la parte escalar del quaternio. 
        """
        
        return self.q0
    
    def vector(self):
        """ 
        Devuelve la parte vectorial del quaternio.
        """
        
        return array([self.q1, self.q2, self.q3], type='Float64')

    @classmethod
    def rot(cls, ang, x, y, z):
        """ 
        Devuelve un quaternio de la rotacion ang por el vector 
        (x,y,z) 
        """

        q0 = cos(ang/2.)
        q1, q2, q3 = ( 
                sin(ang/2.)*array([x, y, z], type='Float64')/
                sqrt(x**2 + y**2 + z**2) )
        return Quaternion(q0, q1, q2, q3)

    @classmethod
    def euler(cls, ch, th, fi):
        """ 
        Devuelve el quaternio de rotacion equivalente a los angulos 
        de euler dados 
        """
        
        return Quaternion.rotacion(
                [ (ch, 0, 0, 1), (th, 0, 1, 0), (fi, 1, 0, 0) ])
    
    def conj(self):
        """ 
        Devuelve el quaternio conjugado 
        """
        
        return Quaternion(self.q0, -self.q1, -self.q2, -self.q3)
    
    def toMatrix(self):
        """ 
        Devuelve la matriz equivalente de rotacion
        """
        
        q0, q1, q2, q3 = self.q0, self.q1, self.q2, self.q3
        M = array(shape = (3, 3), type='Float64')
        
        return array([
            [ 
              1. - 2.*(q2**2 + q3**2),  
              2.*(-q0*q3 + q1*q2), 
              2.*(q0*q2 + q1*q3) 
            ],
            [
              2.*(q0*q3 + q1*q2),
              1. - 2.*(q1**2 + q3**2), 
              2.*(-q0*q1 + q2*q3)
            ],
            [
              2.*(-q0*q2 + q1*q3),
              2.*(q0*q1 + q2*q3), 
              1. - 2.*(q1**2 + q2**2)
            ]], type='Float64')
    
    def toEuler(self):
        """ 
        Devuelve los angulos de euler en la tupla (ch, th, fi).
        """

        q0, q1, q2, q3 = self.q0, self.q1, self.q2, self.q3
        
        # normalizamos
        m = sqrt(q0**2 + q1**2 + q2**2 + q3**2)
        q0 /= m
        q1 /= m
        q2 /= m
        q3 /= m

        Sth = -2*(q1*q3 - q0*q2)
        # Como -pi/2<=th<=pi/2 podemos hacer esto tranquilos
        # Para th = +-pi/2 es posible que por error numerico se Sth>1
        if Sth>1.0:
            Sth = 1.0
        th = asin(Sth)
        
        # si no hay singularidad (th != +-pi/2)
        if abs(Sth) != 1.0:
            ch = ang(1 - 2*(q2**2 + q3**2), 2*(q1*q2 + q0*q3))
            fi = ang(1 - 2*(q1**2 + q2**2), 2*(q0*q1 + q2*q3))
        
        # si no damos un valor cualquiera a ch, ya que hay varias posibilidades
        # para ch y th que nos dan los mismos ejes
        else:
            ch = 0
            fi = ang(q0*q2 + q1*q3, q1*q2 - q0*q3)
        
        return ch, th, fi
    
    def __add__(self, other):
        return Quaternion(
                self.q0 + other.q0, 
                self.q1 + other.q1, 
                self.q2 + other.q2, 
                self.q3 + other.q3)

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            q0 = self.q0*other.q0 - dot(self.vector(), other.vector())
            q1, q2, q3 = ( self.q0*other.vector() + other.q0*self.vector() 
                    + cross(self.vector(), other.vector()) )
        else:
            q0, q1, q2, q3 = ( other*self.q0, other*self.q1, 
                    other*self.q2, other*self.q3 )
        
        return Quaternion(q0, q1, q2, q3)
    
    def __rmul__(self, other):
        return self*other
 
    def __sub__(self, other):
        return Quaternion(
                self.q0 - other.q0, 
                self.q1 - other.q1, 
                self.q2 - other.q2, 
                self.q3 - other.q3 )

    def __repr__(self):
        return "[%f, (%f, %f, %f)]" % (self.q0, self.q1, self.q2, self.q3)
        
    @classmethod    
    def rotacion(cls, l_quat):
        """ 
        Devuelve el quaternio resultante de aplicar sucesivamente las 
        rotaciones contenidas en l_quat. Cada rotacion es una tupla de formato:
            (angulo, x, y, z)
        """
        
        return reduce(lambda x,y: x*y, 
                [Quaternion.rot(r[0], r[1], r[2], r[3])  for r in l_quat])

    def __iter__(self):
        return [self.q0, self.q1, self.q2, self.q3].__iter__()

################################################################################
#           POLINOMIOS                                                         #
################################################################################
def polyeval(p, x):
    """
    Evalua el polinomio en x
    """
    n = len(p) - 1
    r = p[n]
    for i in range(n-1, -1, -1):
        r = p[i] + r*x
    return r

def polydiv(p, a):
    """
    Divide el polinomio p por el monomio x-a
    """
    n = len(p) - 1
    q = array(shape=(n,), type='Complex64')
    
    r = p[n]
    for i in range(n-1, -1, -1):
        s = p[i]
        q[i] = r
        r = s + a*r
    return q

def psolve(p, eps = 0.001, x0 = 0., x1 = 1., x2 = 2., Nmax=100):
    """
    Encuentra una raiz mediante el metodo de Muller
    """

    P0 = polyeval(p, x0)
    P1 = polyeval(p, x1)
    P2 = polyeval(p, x2)
    
    error = eps + 1.
    n = 0
    while error>eps:
        if n>Nmax:
            print "polysolve: n=%d>Nmax=%d" % (n, Nmax)
            return None

        q = (x0 - x1)/(x1 - x2)
        A = q*P0 - q*(1 + q)*P1 + q*q*P2
        B = (2*q + 1)*P0 - (1+q)*(1+q)*P1 + q*q*P2
        C = (1+q)*P0
        # Si tenemos suerte es C = 0, en cuyo caso
        # interrumpimos YA, porque entonces es E = 0
        # y dividiriamos por cero
        if C == 0:
            return x0

        D = cmath.sqrt(B*B  - 4*A*C)
        E1 = B - D
        E2 = B + D
        if abs(E1)>abs(E2):
            E = E1
        else:
            E = E2
        x = x0 - (x0 - x1)*2*C/E
        
        error = abs(x-x0)
        x0, x1, x2 = x, x0, x1
        P0, P1, P2 = polyeval(p, x0), P0, P1
    
    return x0

def polysolve(p, eps = 0.001, x0 = 0., x1 = 1., x2 = 2., Nmax=100):
    """
    Encuentra todas las raices del polinomio, aplicando el metodo
    de Muller y dividiendo por la raiz resultante
    """

    r = []
    q = p
    while len(q)>1:
        s = psolve(q, eps, x0, x1, x2, Nmax)
        r.append(s)
        q = polydiv(q, s)
    return r

def estabilidad(a, f, xa=-4, ya=-4, xb=0.5, yb=4, Nx=100, Ny=100):
    """
    Dibuja la region de estabilidad.
    El polinomio de estabilidad tiene la forma:

    q(r) = sum( (a[j] - w*f[j](w))*r^(p-j), j, 0, p)
    donde p es el numero de pasos del esquema

    xa, ya, xb, yb forma el rectangulo donde se realiza el calculo.
    Nx, Ny son la resolucion que se utiliza para la malla.
    """

    def pol(x, y):
        """
        Calcula el polinomio característico en un punto generico x, y.
        """
        
        # Si un coeficiente no esta presente, vale cero
        la = len(a)
        lf = len(f)
        # lt - 1 = numero de pasos-->lt coeficientes
        # generalmente lf>la
        lt = max(la, lf)
        def geta(i):
            if i>=la:
                return 0.0
            else:
                return a[i]
        
        def getf(i):
            if i>=lf:
                return lambda x: 0.0 + 0.0j
            else:
                return f[i]
        
        # Punto del plano complejo
        w = x + y*1j
        # Lista con los coeficientes
        p = zeros(shape=(lt,), type='Complex64')
       
        # Finalmente construimos el polinomio
        for j in range(lt): 
            p[lt - 1 - j] = geta(j) - w*getf(j)(w)

        return p

    Z = array(shape=(Nx, Ny), type='Float64')
    dx = float(xb - xa)/(Nx-1)
    dy = float(yb - ya)/(Ny-1)
    X = arange(xa, xb + dx, dx) 
    Y = arange(ya, yb + dy, dy)

    mini = Nx - 1
    maxi = 0
    minj = Ny - 1
    maxj = 0
    for i in range(0, Nx):
        for j in range(0, Ny):
            p = pol(X[i], Y[j])
            # Los ceros del polinomio
            z = polysolve(p)
            # Nos quedamos con la raiz de mayor modulo
            r = max(map(abs, z))
            if r<=1:
                if i<mini:
                    mini = i
                if i>maxi:
                    maxi = i
                if j<minj:
                    minj = j
                if j>maxj:
                    maxj = j

            Z[j, i] = r
    
    # Damos un margen
    mini -= 3
    maxi += 3
    minj -= 3
    maxj += 3
    # Evitamos salirnos de limites
    if mini<0:
        mini = 0
    if maxi>Nx - 1:
        maxi = Nx - 1
    if minj<0:
        minj = 0
    if maxj>Ny - 1:
        maxj = Ny - 1

    # Por algun motivo solo funciona con copias
    X0 = X[mini:maxi].copy()
    Y0 = Y[minj:maxj].copy()
    Z0 = Z[minj:maxj, mini:maxi].copy()
    
    pylab.contour(X0, Y0, Z0, arange(0, 1.05, 0.05), colors='k')
    pylab.show()

    return X, Y, Z



if __name__=="__main__":
    def f(x):
        A = array([1, 0, 0], type='Float64')
        B = array([
            [1, 1, 1],
            [1, 1, 0],
            [1, 0, 0]], type='Float64')
        return A + dot(B, x)
    def g(x):
        a, b = x
        return array([a + b - 2, (a - 1)*(b + 2)], type='Float64')
    
    print solveNewtonMulti(f, [0, 0, 0], 100, 1e-5)
    print solveNewtonMulti(g, [0, 0], 100, 1e-5)

