# -*- coding: latin-1 -*-
from numarray import *
from historial import *
from math import *

class Controles:
    def __init__(self, modelo):
        self.modelo = modelo
        # para actualizar todos los relojes al mismo tiempo, utilizamos
        # esta version absurda de reloj
        def reloj_tonto():
            return self.modelo.t
        
        # Colectivo
        self.et0p  = SerieTemporal(1000, reloj_tonto) 
        # Ciclico lateral
        self.et1cp = SerieTemporal(1000, reloj_tonto) 
        # Ciclico longitudinal
        self.et1sp = SerieTemporal(1000, reloj_tonto) 
        # Colectivo de cola (pedal)
        self.etpp  = SerieTemporal(1000, reloj_tonto) 
    
    def init(self):
        def func(x):
            if not hasattr(x, '__call__'):
                return lambda y: x
            else:
                return x
        self.th_0 = func(self.th_0)
        self.fi_0 = func(self.fi_0)
        self.ch_0 = func(self.ch_0)

        self.c30_bak = self.c1[3,0]
        self.c00_bak = self.c1[0,0]

    def calcula_controles(self):
        t = self.modelo.t
        # Controles
        et0p  = self.et0p(t)
        et1sp = self.et1sp(t)
        et1cp = self.et1cp(t)
        etpp  = self.etpp(t)

        #L = self.modelo.LBI
        #Ca = -(L[0,0]*L[1,1] - L[0,1]*L[1,0])/sqrt(L[1,0]**2 + L[1,1]**2)
        #if Ca>1:
        #    Ca = 1.0
        #elif Ca<-1:
        #    Ca = -1.0
        #a = acos(Ca)
        #if L[0,2]<0:
        #    a = -a
        #print a

        et = dot(self.S0, [ et0p, et1sp, et1cp, etpp ]) +\
                dot(self.S1, [ self.modelo.DTfi, 
                               self.modelo.DTth, 
                               self.modelo.r ]) +\
                dot(self.S2, [ self.modelo.th - self.th_0(et1sp), 
                               self.modelo.fi - self.fi_0(et1cp), 
                               self.modelo.ch - self.ch_0(etpp)   ])

        # Transforma la posicion del control en un porcentaje
        # Estas formulas no son validas para el colectivo que 
        # unicamente se mueve en el rango de 0 a 1
        def per(x):
            return (x + 1.)/2.
        
        # Transorma un porcentaje en posicion de un control
        def iper(x):
            return 2*x - 1
        
        p1sp = per(et1sp) 
        p1cp = per(et1cp) 
        ppp  = per(etpp)
        p1s  = per(et[1])
        p1c  = per(et[2])
        pp   = per(et[3])
        
        limite = self.L
        
        if p1s>p1sp + limite:
            et[1] = iper(p1sp + limite)
        elif p1s<p1sp - limite:
            et[1] = iper(p1sp - limite)
        
        if p1c>p1cp + limite:
            et[2] = iper(p1cp + limite)
        elif p1c<p1cp - limite:
            et[2] = iper(p1cp - limite)

        if pp>ppp + limite:
            et[3] = iper(ppp + limite)
        elif pp<ppp - limite:
            et[3] = iper(ppp - limite)

        # Pasos
        if self.modelo.motor.nmotores<=0:
            self.c1[3,0] = 0.0
            self.c1[0,0] = self.c00_bak*1
        else:
            self.c1[0,0] = self.c00_bak
            self.c1[3,0] = self.c30_bak


        th = self.c0 + dot(self.c1, et)
        self.th0, self.th1s, self.th1c, self.th0T = th

        self.l_et0p  = et0p
        self.l_et1sp = et1sp
        self.l_et1cp = et1cp
        self.l_etpp  = etpp
