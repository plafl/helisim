from ctypes import *

################################################################################
#       TIPOS ADICIONALES DE LA LIBRERIA EN C                                  #
################################################################################
class CPunto(Structure):
    _fields_ = [
            ( "x", c_double * 3),
            ( "y", c_double * 3),
            ( "D", c_long      ),
            ( "n", c_double    ),
            ( "a", c_double    ),
            ( "b", c_double    ),
            ( "c", c_double    ),
            ( "d", c_double    ),
            ( "e", c_double    ),
            ( "f", c_double    ) ]


class CTren(Structure):
    _fields_ = [
            ( "M",          c_double                 ),
            ( "Ixx",        c_double                 ),
            ( "Iyy",        c_double                 ),
            ( "Izz",        c_double                 ),
            ( "Ixz",        c_double                 ),
            ( "FC",         c_double * 3             ),
            ( "MC",         c_double * 3             ),
            ( "FB",         c_double * 3             ),
            ( "MB",         c_double * 3             ),
            ( "x",          c_double * 13            ),
            ( "t",          c_double                 ),
            ( "T0",         c_double                 ),
            ( "P",          POINTER(POINTER(CPunto)) ),
            ( "nP",         c_long                   ),
            ( "colision",   c_long                   ),
            ( "x0",         POINTER(c_double)        ),
            ( "x1",         POINTER(c_double)        ),
            ( "x2",         POINTER(c_double)        ),
            ( "xa",         POINTER(c_double)        ),
            ( "xi",         POINTER(c_double)        ),
            ( "F0",         POINTER(c_double)        ),
            ( "F1",         POINTER(c_double)        ),
            ( "F2",         POINTER(c_double)        ),
            ( "Fi",         POINTER(c_double)        ),
            ( "Fa",         POINTER(c_double)        ),
            ( "ti",         c_double                 ),
            ( "t0",         c_double                 ),
            ( "t1",         c_double                 ),
            ( "t2",         c_double                 ),
            ( "dt1",        c_double                 ),
            ( "dt2",        c_double                 ) ]

################################################################################
#       INICIALIZACION DE LA LIBRERIA EN C                                     #
################################################################################
c_colision = cdll.LoadLibrary("c_colision.dll")

# struct Tren * malloc_tren(long)
c_colision.malloc_tren.arg_types = [c_long]
c_colision.malloc_tren.restype = POINTER(CTren)

# void free_tren(struct Tren *)
c_colision.free_tren.arg_tpes = [POINTER(CTren)]
c_colision.free_tren.restype = None

# void step_tren(struct Tren *, double)
c_colision.step_tren.arg_tpes = [POINTER(CTren), c_double]
c_colision.step_tren.restype = None

################################################################################
#       FUNCIONES DE CONVERSION DE TIPOS                                       #
################################################################################
def PyC_Punto(P):
    """
    A partir de un diccionario que representa un punto P devuelve una 
    estructura de C.
    """
    C = CPunto()

    C.x[0] = P['x'][0]
    C.x[1] = P['x'][1]
    C.x[2] = P['x'][2]
    
    C.y[0] = P['y'][0]
    C.y[1] = P['y'][1]
    C.y[2] = P['y'][2]

    C.D = P['D']
    
    C.n = P['n']
    C.a, C.b, C.c = P['a'], P['b'], P['c']
    C.d, C.e, C.f = P['d'], P['e'], P['f']

    return C

    
def CPy_Punto(C, P):
    """
    A partir de la estructura en C de un punto asigna los valores en el  
    diccionario de python que representa un punto.
    """
    
    P['x'][0] = C.x[0]
    P['x'][1] = C.x[1]
    P['x'][2] = C.x[2]

    P['y'][0] = C.y[0]
    P['y'][1] = C.y[1]
    P['y'][2] = C.y[2]

    P['D'] = C.D
    P['a'], P['b'], P['c'] = C.a, C.b, C.c
    P['d'], P['e'], P['f'] = C.d, C.e, C.f

################################################################################
#       CLASE TREN                                                             #
################################################################################
class Tren:
    def __init__(self, modelo):
        """
        modelo contiene los datos masicos:
            M, Ixx, Iyy, Izz, Ixz
        """
        self.modelo = modelo
        self.puntos = []
        
        def ceros(n):
            return [0.0 for x in range(n)]

        self.x  = ceros(13)
        self.F  = ceros(13)
        self.FC = ceros(3)
        self.MC = ceros(3)
        self.FB = ceros(3)
        self.MB = ceros(3)

    def init(self):
        """
        A partir de los puntos reserva la memoria necesaria
        """
        # pC es un puntero a la imagen en C
        nP = len(self.puntos)
        self.pC = c_colision.malloc_tren(nP)
        # Estas magnitudes se inicializan en la imagen y
        # ya no se tocan desde python
        for i in range(nP):
            self.pC.contents.P[i]= pointer(PyC_Punto(self.puntos[i]))
        self.colision = self.pC.contents.colision = 0

    def __del__(self):
        # Liberamos la memoria reservada dentro de C
        c_colision.free_tren(self.pC)

    def reset(self):
        self.colision = self.pC.contents.colision = 0

    def punto(self, x, n, a, b, c, d, e, f):
        self.puntos.append({
            'x':    x,
            'y':    [0., 0., 0.],
            'D':    0,
            'n':    n,
            'a':    a,
            'b':    b,
            'c':    c,
            'd':    d,
            'e':    e,
            'f':    f})
    
    def step(self, t):
        def copyPC(a, b, n):
            for i in range(n):
                b[i] = c_double(a[i])
        
        def copyCP(a, b, n):
            for i in range(n):
                b[i] = a[i]

        
            
        # De Pyton a C
        self.pC.contents.M   = c_double(self.modelo.M  )
        self.pC.contents.Ixx = c_double(self.modelo.Ixx)
        self.pC.contents.Iyy = c_double(self.modelo.Iyy)
        self.pC.contents.Izz = c_double(self.modelo.Izz)
        self.pC.contents.Ixz = c_double(self.modelo.Ixz)
         
        copyPC(self.FC, self.pC.contents.FC, 3)
        copyPC(self.MC, self.pC.contents.MC, 3)
        copyPC(self.FB, self.pC.contents.FB, 3)
        copyPC(self.MB, self.pC.contents.MB, 3)
         
        copyPC(self.x, self.pC.contents.x, 13)
        self.pC.contents.t = c_double(self.t)
    
        self.pC.contents.T0 = c_double(self.T0)
         
        # Llamada
        c_colision.step_tren(self.pC, c_double(t))
         
        # De C a Python
        copyCP(self.pC.contents.x,  self.x, 13)
        copyCP(self.pC.contents.F0, self.F, 13)
        self.t = self.pC.contents.t
        self.colision = self.pC.contents.colision
        
