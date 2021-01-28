# -*- coding: latin-1 -*-
"""
Contiene diversas utilidades para la gestion de valores que se obtienen en 
forma de series temporales:
"""

import time

class Historial:
    """ 
    La clase historial se compone de una lista cuyo primer elemento ha sido el 
    ultimo en añadirse, su segundo elemento el penultimo, ... asi hasta el 
    n-esimo elemento a partir del cual no se han guardado mas datos.
    """

    def __init__(self, longitud):
        """
        longitud: Tamaño del buffer para el historial, es decir, guarda tantos
        valores pasados como se especifica en longitud.
        """

        # Creamos la lista con la lonigtud especificada
        self.longitud = longitud
        self.valores  = [None for x in range(longitud)]
        
        # señala la posicion del ultimo elemento añadido en la lista
        self.zero = -1
    
    def __getitem__(self, k):
        """
        Elemento k-veces anterior, es decir 0 es el mas reciente, 1 el anterior,
        etc... asi hasta tantos elementos como se especificaron con el 
        parámetro longitud al crear la instancia.
        """

        return self.valores[( self.zero - k )%self.longitud]
    
    def __setitem__(self, k, v):
        """
        Asigna el valor v al elemento k-veces anterior.
        """

        self.valores[( zero - k )%self.longitud] = v
    
    def append(self, v):
        """ 
        Añade un nuevo elemento a la lista. 
        """
        self.zero = (self.zero + 1) % self.longitud
        self.valores[self.zero] = v

class SerieTemporal(Historial):
    """ 
    Extiende historial para llevar un control del instante en que se añadió un
    elemento al historial, asi como obtener mediante interpolacion el valor en 
    un instante arbitrario. Para que esto ultimo sea posible es obvio que el 
    valor debe de ser numérico.
    """
    
    def __init__(self, longitud, timeFunc=time.time):
        """
        timeFunc: funcion que cada vez que se la llama devuelve el tiempo en ese
        mismo momento. Puede ser el tiempo real o no. Por defecto si no se
        espefica nada consulta el reloj del sistema. Cada vez que se añade un
        elemento con append se llama a esta funcion para asignarle un valor de
        tiempo al valor numérico.
        """
        
        Historial.__init__(self, longitud)
        self.timeFunc = timeFunc
    
    def __setitem__(self, k, v):
        """ 
        Sin sentido, no hace nada. Si se quiere añadir un valor nuevo, ver
        append.
        """

        pass
    
    def __getitem__(self, k):
        """ 
        Obtiene el valor del k-esimo elemento. 
        """

        it = Historial.__getitem__(self, k)
        if it != None:
            return it[1]
        else:
            return it
    
    def append(self, v, t=None):
        if t==None:
            Historial.append(self, (self.timeFunc(), v))
        else:
            Historial.append(self, (t, v))

    def time(self, k):
        """ 
        Obtiene el tiempo del k-esimo elemento (por el final). 
        """
        
        it = Historial.__getitem__(self, k)
        if it != None:
            return it[0]
        else:
            return None
    
    def times(self):
        """ 
        Devuelve un iterador por la lista de tiempos. 
        """
        
        class IterTimes:
            def __init__(self, st):
                self.c = 0
                self.st = st
            
            def __iter__(self):
                return self
            
            def next(self):
                t = self.st.time(self.c)
                if t == None or self.c == self.st.longitud:
                    raise StopIteration
                else:
                    self.c += 1
                    return t
        
        return IterTimes(self)
    
    def separa(self):
        t = []
        x = []
        for i in self.valores:
            if i == None:
                break
            t.append(i[0])
            x.append(i[1])
        return t, x

    def __call__(self, t):
        """ 
        Interpola linealmente el valor en el instante t 
        """
        
        c = c1 = c2 = 0
        d = d1 = d2 = None
        while d1==None or d2==None or d<=d2:
            _t = self.time(c)
            if _t == None:
                break
            else:
                d = abs(t - _t)
            if d1 == None or d<d1:
                c1, d1 = c, d
            elif d2 == None or d<d2:
                c2, d2 = c, d
            c += 1
        t1, t2 = map(self.time, [c1, c2])
        v1, v2 = map(self.__getitem__, [c1, c2])
        # a lo mejor la lista solo contiene un elemento
        if t1 == t2:
            return v1
        else:
            # si contiene por lo menos dos, pues interpolacion lineal
            return (v2*(t - t1) - v1*(t - t2))/(t2 - t1)
        
class Reloj:
    """ 
    Permite controlar el tiempo mediante los metodos arrancar y parar. 
    """
    
    def __init__(self, timeFunc=time.time):
        """
        timeFunc: funcion que cada vez que se la llama devuelve el tiempo en ese
        mismo momento. Puede ser el tiempo real o no. Por defecto si no se
        espefica nada consulta el reloj del sistema.
        """

        self.timeFunc = timeFunc
        self.t0 = None
        self.retraso = 0
        self.funciona = False
    
    def arranca(self):
        """ 
        Arranca el reloj, si ya estaba arrancado no pasa nada. 
        """
        
        if self.t0 == None:
            self.t0 = self.timeFunc()
            self.funciona = True
        if not self.funciona:
            self.retraso += self.timeFunc() - self.parada
            self.funciona = True
    
    def para(self):
        """ 
        Para el reloj. 
        """

        if self.funciona:
            self.parada = self.timeFunc()
            self.funciona = False

    def __call__(self):
        """ 
        Consulta el tiempo actual. 
        """
        
        if self.funciona:
            t = self.timeFunc() - self.t0 - self.retraso
        elif self.t0!=None:
            t = self.parada - self.t0 - self.retraso
        else:
            t = None
        return t

################################################################################
#       TESTING                                                                #
################################################################################
if __name__=='__main__':
    st = SerieTemporal(10, time.time)
    st.append(1)
    time.sleep(0.1)
    st.append(2)
    for t in st.times():
        print "%f" % t
    t1, t2 = st.time(0), st.time(1)
    t = (t1 + t2)/2
    print st(t)
