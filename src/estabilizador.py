from matematicas import *
from numarray import *

class EstabilizadorHorizontal:
    def __init__(self, modelo):
        self.modelo = modelo

    def init(self):
        x = self.modelo.rotor.hR - self.modelo.rotor_cola.hT
        # angulos entre los cuales la estela afecta al estabilizador
        # xi1<x<xi2
        self.xi1 = ang(x, self.modelo.rotor.xcg + self.ltp -
                self.modelo.rotor.R)
        self.xi2 = ang(x, self.modelo.rotor.xcg + self.ltp +
                self.modelo.rotor.R)
    
    def calcula_FMtp(self):
        """
        Calcula las fuerzas y momentos del estabilizador horizontal.

        NECESITA:
            modelo.Vrw_b
            rotor.la0
            rotor.xi

        PROVEE:
            estabilizador_horizontal.FMtp
        """
        
        # Velocidad relativa al viento de la cola, en ejes cuerpo
        urwtp_b = self.modelo.Vrw_b[0] - self.modelo.q*self.htp
        vrwtp_b = ( self.modelo.Vrw_b[1] 
                - self.modelo.r*(self.ltp + self.modelo.rotor.xcg) 
                + self.modelo.p*self.htp  )
        wrwtp_b = self.modelo.Vrw_b[2] + self.modelo.q*\
                (self.ltp + self.modelo.rotor.xcg) 
        
        # Angulo de ataque
        al = self.altp0 + ang(urwtp_b, wrwtp_b)
        Sal, Cal = sin(al), cos(al)
        # Angulo de resbalamiento
        be = ang(urwtp_b, vrwtp_b)
        # Coeficientes en ejes viento
        cL, cD = self.aero.coefs(al, be)
        # Coeficientes en ejes cuerpo
        cZ = -cD*Sal - cL*Cal
        cX =  cL*Sal - cD*Cal
        # Presion dinamica
        pd = 0.5*self.modelo.ro*(urwtp_b**2 + wrwtp_b**2)
        Ztp = pd*self.Stp*cZ
        Xtp = pd*self.Stp*cX

        self.FMtp = (
                array([Xtp, 0.0, Ztp, 0.0, 
                    (self.ltp + self.modelo.rotor.xcg)*Ztp -
                    self.htp*Xtp, 0.0],  type='Float64'), 
                zeros(shape=(6, 16), type='Float64')
        )

        # logging
        self.al = al
        self.be = be
        self.cL, self.cD = cL, cD
        ( self.Xtp_b, self.Ytp_b, self.Ztp_b, 
                self.Ltp_b, self.Mtp_b, self.Ntp_b ) =\
                self.FMtp[0]

class EstabilizadorVertical:
    def __init__(self, modelo):
        self.modelo = modelo

    def calcula_FMfn(self):
        """
        Calcula las fuerzas y momentos del estabilizador vertical.

        NECESITA:
            modelo.Vrw_b
        PROVEE:
            estabilizador_vertical.FMfn
        """

        # Velocidad relativa al viento del estab vert en ejes cuerpo
        vrwfn_b = self.modelo.Vrw_b[1] - self.modelo.r*\
                (self.lfn + self.modelo.rotor.xcg) + self.hfn*self.modelo.p
        urwfn_b = self.modelo.Vrw_b[0] - self.modelo.q*self.hfn
        wrwfn_b = self.modelo.Vrw_b[2] + self.modelo.q*\
                (self.lfn + self.modelo.rotor.xcg) 
        # Angulo de ataque del estabilizador
        be = self.befn0 + ang(urwfn_b, vrwfn_b)
        Sbe, Cbe = sin(be), cos(be)
        # Angulo de resbalamiento
        al = ang(urwfn_b, wrwfn_b)
        # Coeficientes en ejes viento
        cL, cD = self.aero.coefs(be, al)
        # Coeficientes en ejes cuerpo
        cY = -cD*Sbe - cL*Cbe
        cX =  cL*Sbe - cD*Cbe
        # Fuerzas
        pd = 0.5*self.modelo.ro*(urwfn_b**2 + vrwfn_b**2)
        Xfn = pd*self.Sfn*cX
        Yfn = pd*self.Sfn*cY
        Zfn = 0.0

        Lfn, Mfn, Nfn = cross(
                [-(self.lfn + self.modelo.rotor.xcg), 0.0, -self.hfn],
                [  Xfn,                               Yfn,  Zfn     ])

        #Nfn = -Yfn*(self.lfn + self.modelo.rotor.xcg)
        self.FMfn = ( 
                array([Xfn, Yfn, Zfn, Lfn, Mfn, Nfn], type='Float64'),
                zeros(shape=(6, 16), type='Float64')
        )

        # logging
        self.be = be
        self.al = al
        self.cL = cL
        self.cD = cD
        self.Xfn_b = Xfn
        self.Yfn_b = Yfn
        self.Zfn_b = Zfn
        self.Lfn_b = Lfn
        self.Mfn_b = Mfn
        self.Nfn_b = Nfn

