from historial import *
from numarray import arange
from matematicas import *
        
class Euler:
    # Esta se llama en input.py
    def __init__(self):
        pass

    # Esta en FlightGear.reset()
    def init(self, F, ti, xi):
        """
        Coloca las condiciones iniciales
        """
        self.F, self.ti, self.xi = F, ti, xi
        
        self.x0 = xi
        self.t0 = ti
        self.F0 = self.F(xi, ti)
        
    # En FlightGear.loop()
    def __call__(self, t):
        """
        Calcula la respuesta hasta el instante t
        """
        
        if t<self.t0:
            return 
        
        self.x1 = self.x0
        self.t1 = self.t0
        self.F1 = self.F0

        self.x0 = self.x1 + self.F1*(t - self.t1)
        self.t0 = t
        self.F0 = self.F(self.x0, self.t0)

class RungeKutta4:
    def __init__(self):
        pass
    
    def init(self, F, ti, xi):
        """
        Coloca las condiciones iniciales
        """
        self.F, self.ti, self.xi = F, ti, xi
        
        self.x0 = xi
        self.t0 = ti
        self.F0 = self.F(xi, ti)


    def __call__(self, t):
        # La primera vez obtenemos el punto inicial
        # de las condiciones iniciales
       
        if t<self.t0:
            return

        self.t1 = self.t0
        self.x1 = self.x0
        self.F1 = self.F0

        # paso
        h = t - self.t1
        
        # Estimaciones de la derivada
        k1 = self.F1
        k2 = self.F(self.x1 + h*k1/2., self.t1 + h/2.)
        k3 = self.F(self.x1 + h*k2/2., self.t1 + h/2.)
        # Calculamos el nuevo punto
        k4 = self.F(self.x1 + h*k3, self.t1 + h)
        
        self.x0 = self.x1 + h/6.*(k1 + 2*k2 + 2*k3 + k4)
        self.t0 = t
        self.F0 = self.F(self.x0, self.t0)

class AdamsBashforth2:
    """
    Esquema de Adams-Bashforth de orden 2, de paso variable.
    Nota:
        No calcula el paso para integrar hasta el instante
        t, solo da un paso, por lo que es responsabilidad
        de otra funcion elegir el paso.
    """

    def __init__(self):
        pass
    
    def init(self, F, ti, xi):
        """
        Coloca las condiciones iniciales:
            
            x(ti) = xi
        
        y define F(x, t), tal que:
            
            dx
            -- = F(x, t)
            dt
        """
        self.F, self.ti, self.xi = F, ti, xi
        # x0, t0 guardan los ultimos resultados, de momento
        # esos son las condiciones iniciales
        self.x0 = xi
        self.t0 = ti
        self.F0 = self.F(xi, ti)

        self.ci = 0
    
    def __call__(self, t):
        """
        Integra hasta el instante de tiempo t
        El resultado e instante final se guardan en:
            x0
            t0 
        Si el esquemas es multipaso de orden p, los anteriores
        pasos se encuenran en:
            x1, x2, ..., xp
            t1, t2, ..., tp
        """
        dtmin = 1e-8
        # Si pedimos un tiempo equivocado, ignoramos
        if t <= self.ti or ( self.ci != 0 and t<=self.t0 ):
            return 
        
        # Hemos arrancado el esquema?
        if self.ci == 0:
            self.ci = 1
            
            self.x1 = self.x0
            self.t1 = self.t0
            self.t0 = self.t1 + dtmin

            self.dt1 = t - self.t0
            self.F1 = self.F0
            
            # Un Euler para la primera aproximacion
            self.x0 = self.x1 + self.dt1*self.F1
            self.F0 = self.F(self.x0, self.t0)
       
        self.t2 = self.t1
        self.t1 = self.t0
        self.t0 = t
            
        self.dt2 = self.dt1
        self.dt1 = t - self.t1
            
        self.x2 = self.x1
        self.x1 = self.x0
            
        self.F2 = self.F1
        self.F1 = self.F0

        b1 = (2*self.dt2 + self.dt1)/2./self.dt2
        b2 = -self.dt1/2/self.dt2

        self.x0 = self.x1 + self.dt1*(b1*self.F1 + b2*self.F2)
        self.F0 = self.F(self.x0, self.t0)
    
class AdamsBashforth3:
    """
    Esquema de Adams-Bashforth de orden 3, de paso variable.
    Nota:
        No calcula el paso para integrar hasta el instante
        t, solo da un paso, por lo que es responsabilidad
        de otra funcion elegir el paso.
    """

    def __init__(self):
        pass
    
    def init(self, F, ti, xi):
        """
        Coloca las condiciones iniciales:
            
            x(ti) = xi
        
        y define F(x, t), tal que:
            
            dx
            -- = F(x, t)
            dt
        """
        self.F, self.ti, self.xi = F, ti, xi
        # x0, t0 guardan los ultimos resultados, de momento
        # esos son las condiciones iniciales
        self.x0 = xi
        self.t0 = ti
        self.F0 = self.F(xi, ti)

        self.ci = 0
    
    def __call__(self, t):
        """
        Integra hasta el instante de tiempo t
        El resultado e instante final se guardan en:
            x0
            t0 
        Si el esquemas es multipaso de orden p, los anteriores
        pasos se encuenran en:
            x1, x2, ..., xp
            t1, t2, ..., tp
        """
        dtmin = 1e-8
        
        # Si pedimos un tiempo equivocado, ignoramos
        if t <= self.ti or ( self.ci != 0 and t<self.t0 ):
            return 
        
        # Hemos arrancado el esquema?
        if self.ci == 0:
            self.ci = 1
            
            self.x1 = self.x0
            self.t1 = self.t0
            self.t0 = self.t1 + dtmin

            self.dt1 = t - self.t0
            self.F1 = self.F0
            
            # Un Euler para la primera aproximacion
            self.x0 = self.x1 + self.dt1*self.F1
            self.F0 = self.F(self.x0, self.t0)
        
        elif self.ci == 1:
            self.ci = 2

            self.t2 = self.t1
            self.t1 = self.t0
            
            self.dt2 = self.dt1
            self.dt1 = t - self.t1
            
            self.x2 = self.x1
            self.x1 = self.x0 
            
            self.F2 = self.F1
            self.F1 = self.F0

            b1 = (2*self.dt2 + self.dt1)/2./self.dt2
            b2 = -self.dt1/2/self.dt2

            # Un AdamsBashForth2
            self.x0 = self.x1 + self.dt1*(b1*self.F1 + b2*self.F2)
            self.t0 = t
            self.F0 = self.F(self.x0, self.t0)

        else:
            self.t3 = self.t2
            self.t2 = self.t1
            self.t1 = self.t0

            self.dt3 = self.dt2
            self.dt2 = self.dt1
            self.dt1 = t - self.t1

            self.x3 = self.x2
            self.x2 = self.x1
            self.x1 = self.x0
            
            self.F3 = self.F2
            self.F2 = self.F1
            self.F1 = self.F0

            dt1 = self.dt1
            dt2 = self.dt2
            dt3 = self.dt3

            b1 =  1. + dt1*(2*dt1 + 6*dt2 + 3*dt3)/(6*dt2*(dt2 + dt3))
            b2 = -dt1*(2*dt1 + 3*dt2 + 3*dt3)/(6*dt2*dt3)
            b3 =  dt1*(2*dt1 + 3*dt2)/(6*dt3*(dt2 + dt3))

            self.x0 = self.x1 + self.dt1*(b1*self.F1 + b2*self.F2 
                    + b3*self.F3)
            self.t0 = t 
            self.F0 = self.F(self.x0, self.t0)
            
class ABM2:
    """
    Predictor-Corrector basado en un Adams-Bashforth 2, y un Adams-Moulton
    1 de paso variable.
    Nota:
        No calcula el paso para integrar hasta el instante
        t, solo da un paso, por lo que es responsabilidad
        de otra funcion elegir el paso.
    """

    def __init__(self):
        pass
    
    def init(self, F, ti, xi):
        """
        Coloca las condiciones iniciales:
            
            x(ti) = xi
        
        y define F(x, t), tal que:
            
            dx
            -- = F(x, t)
            dt
        """
        self.F, self.ti, self.xi = F, ti, xi
        # x0, t0 guardan los ultimos resultados, de momento
        # esos son las condiciones iniciales
        self.x0 = xi
        self.t0 = ti
        self.F0 = self.F(xi, ti)

        self.ci = 0
    
    def __call__(self, t):
        """
        Integra hasta el instante de tiempo t
        El resultado e instante final se guardan en:
            x0
            t0 
        Si el esquemas es multipaso de orden p, los anteriores
        pasos se encuenran en:
            x1, x2, ..., xp
            t1, t2, ..., tp
        """
        # Si pedimos un tiempo equivocado, ignoramos
        if t <= self.ti or ( self.ci != 0 and t<self.t0 ):
            return 
        
        # Hemos arrancado el esquema?
        if self.ci == 0:
            self.ci = 1
            
            self.x1 = self.x0
            self.t1 = self.t0

            self.dt1 = t - self.t1
            self.F1 = self.F0
            
            # Un Euler para la primera aproximacion
            self.x0 = self.x1 + self.dt1*self.F1
            self.t0 = t
            self.F0 = self.F(self.x0, self.t0)
        else:
            self.t2 = self.t1
            self.t1 = self.t0
            
            self.dt2 = self.dt1
            self.dt1 = t - self.t1
            
            self.x2 = self.x1
            self.x1 = self.x0
            
            self.F2 = self.F1
            self.F1 = self.F0

            b1 = (2*self.dt2 + self.dt1)/2./self.dt2
            b2 = -self.dt1/2/self.dt2

            # Punto estimado por el el predictor
            xp = self.x1 + self.dt1*(b1*self.F1 + b2*self.F2)
            Fp = self.F(xp, self.t0)
            
            # Paso del corrector
            self.x0 = self.x1 + self.dt1*0.5*(Fp + self.F1)
            self.t0 = t
            self.F0 = self.F(self.x0, self.t0)

