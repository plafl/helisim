# -*- coding: latin-1 -*-
from numarray import *

class Motor:
    def __init__(self, modelo):
        self.modelo = modelo
    
    def init(self):
        self.start = False
        self.K3_bak = self.K3
    
    def calcula_K3(self):
        if self.start:
            self.K3 = 1e3
        else:
            self.K3 = self.K3_bak
        if self.modelo.Om>0.2*self.Omi:
            self.start = False
        
    
    def calcula_Ir(self):
        """ 
        Inercia de giro, referida al rotor, del sistema transmision + rotor. 
        """

        # Inercia del sistema de transmision (Irt) mas inercia de las palas 
        # (Nb*Ibe)
        self.Ir = self.Irt[self.nmotores] +\
                self.modelo.rotor.Nb*self.modelo.rotor.Ibe

    def calcula_tae3(self):
        """ 
        Constante de tiempos. 
        """
        Q =  self.Qef/self.Qmc
        if Q<=1:
            x1 = self.Qid/self.Qmc
            self.tae3 = self.a3 + (self.b3 - self.a3)/(1 - x1)*(Q - x1)
        else:
            x2 = self.Qto/self.Qmc
            self.tae3 = self.b3 + (self.c3 - self.b3)/(x2 - 1)*(Q - 1)
    
    def calcula_tae2(self):
        """ 
        Constante de tiempos. 
        """
        Q =  self.Qef/self.Qmc
        if Q<=1:
            x1 = self.Qid/self.Qmc
            self.tae2 = self.a2 + (self.b2 - self.a2)/(1 - x1)*(Q - x1)
        else:
            x2 = self.Qto/self.Qmc
            self.tae2 = self.b2 + (self.c2 - self.b2)/(x2 - 1)*(Q - 1)

    def calcula_tae1(self):
        """ 
        Constante de tiempos. 
        """
        pass

    def calcula_ecs_motor(self):
        # Corregimos velocidad angular
        if self.modelo.Om<0:
            self.modelo.integrador.x1[10] = self.modelo.Om = 0
        
        # Evidentemente para Om pequena las formulas siguientes no
        # son validas, pero tiramos adelante igual
        Om = self.modelo.Om
        if Om<1.0:
            Om = 1.0
        k = (self.modelo.ro/1.225)**(self.x)
        self.Wto = self.Wto_0*k
        self.Wmc = self.Wmc_0*k
        # Par maximo
        self.Qto = self.Wto/Om
        # Para maximo continuo
        self.Qmc = self.Wmc/Om
        # Par minimo
        self.Qid = self.Wid/Om
        
        # Corregimos par
        if self.modelo.Q1>self.Qto:
            self.modelo.integrador.x0[12] = self.modelo.Q1 = self.Qto

        # Par efectivo para calculo de coeficientes
        # Qid <= Qef <= Qto
        if self.modelo.Q1<self.Qid:
            self.Qef = self.Qid
        else:
            self.Qef = self.modelo.Q1
    
        # Calculamos constantes de tiempo
        self.calcula_tae1()
        self.calcula_tae2()
        self.calcula_tae3()

        if Om>0.8*self.Omi:
            self.start = False
        
        # Rigidez del motor
        self.calcula_K3()
        
        # Inercia total del sistema
        self.calcula_Ir()
        
        # Mucho ojo con esto, porque si se modifica x hay que modificar 
        # tambien esto
        ecs_1 = zeros(shape = (3, 16), type = 'Float64')
        ecs_1[0, 5 ] = -1.0
        ecs_1[0, 10] =  1.0
        ecs_1[0] += 1./self.Ir*self.modelo.rotor.QR[1]
        ecs_1[1, 10:13] = [0.0, 0.0, 1.0]
        ecs_1[2, 10:13] = [self.K3*self.tae2, self.tae1*self.tae3, 
                self.tae1 + self.tae3]

        # Para de rotor
        QR =   self.modelo.rotor.QR[0]
        # Para de cola
        QT =   self.modelo.rotor_cola.QT
        # Par total
        Qtot = (1. + self.P)*(QR + self.modelo.rotor_cola.gT*QT)
        if self.modelo.Om>0:
            Qtot += self.Qid

        # Diferencia entre par de motor y total
        Qdif = 1./self.Ir*(self.nmotores*self.modelo.Q1 - Qtot)

        ecs_0 = array(shape = (3,), type = 'Float64')
        ecs_0[0] =  Qdif
        ecs_0[1] =  self.modelo.DTQ1
        ecs_0[2] = -self.modelo.Q1 + self.K3*( self.Omi - Om )
        
        self.ecs_motor = ecs_0, ecs_1
