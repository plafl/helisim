from numarray import *
import numarray.linear_algebra as la
from matematicas import *

class RotorCola:
    def __init__(self, modelo):
        self.params = {}
        self.modelo = modelo

    def init(self):
        """ 
        Inicializa algunos parametros de la cola. 
        """
        
        self.params['bat'] = {'ila0T': -1.0, 'Nmax': 100, 'eps': 1e-4}
        SK, CK = sin(self.K), cos(self.K)
        # Matriz de cambio de base de ejes cuerpo a ejes cola
        self.LBT = array([
            [ 1.0, 0.0, 0.0 ],
            [ 0.0, SK, -CK  ],
            [ 0.0, CK,  SK  ]], type='Float64')
        self.LTB = transpose(self.LBT)
        self.sT = 4*self.cT/pi/self.RT
    
    def calcula_FMT(self):
        """ 
        Calcula las fuerzas y momentos del rotor de cola
        
        NECESITA:
            modelo.Vrw_b
            controles.th0T
        PROVEE:
            rotor_cola.FMT
            rotor_cola.QT

        """
        ########################################################################
        #   CAMBIO A EJES VIENTO DEL ROTOR                                     #
        ########################################################################
        
        # Velocidad de giro del rotor de cola
        OmT = self.gT*self.modelo.Om
        if OmT<1.0:
            self.FMT = (
                    zeros(shape=(6, ),   type='Float64'),
                    zeros(shape=(6, 16), type='Float64')
            )
            self.QT = 0.0
            return 
            
        # Velocidad de la cola, relativa al viento en ejes cuerpo
        VTrw_b =  self.modelo.Vrw_b +\
                cross(
                    [   self.modelo.p, self.modelo.q,  self.modelo.r ], 
                    [ -(self.modelo.rotor.xcg + self.lT), 0.0, -self.hT ]
                )

        # Velocidad de la cola, relativa al viento en ejes cola
        uTrw_T, vTrw_T, wTrw_T = dot(self.LTB, VTrw_b)
        # Velocidad de la cola, relativa al viento en ejes viento
        uTrw_w = sqrt( uTrw_T**2 + vTrw_T**2 )
        wTrw_w = wTrw_T
        # Angulo ejes viento-ejes cola
        if uTrw_w!=0:
            CchwT = uTrw_T/uTrw_w
            SchwT = vTrw_T/uTrw_w
        else:
            if vTrw_T==0:
                CchwT = 1.0
                SchwT = 0.0
            elif vTrw_T<0:
                CchwT =  0.0
                SchwT = -1.0
            else:
                CchwT = 0.0
                SchwT = 1.0
        # Parametros de avance
        nuT =  uTrw_w/( OmT*self.RT )
        nuzT = wTrw_w/( OmT*self.RT )
        
        
        ########################################################################
        #   CALCULO DEL BATIMIENTO                                             #
        ########################################################################
        
        # Numero de Locke
        gaT = self.modelo.ro*self.cT*self.a0T*self.RT**4/self.IbeT
        # Frecuencia natural de batimiento al cuadrado
        labeT2 = 1. + self.KbeT/( self.IbeT*(OmT**2) )
 
        # Lado izquierdo de la ecuacion de batimiento
        A_bebeT = array(shape=(3, 3), type='Float64')
        
        A_bebeT[0, 0] = (8.*labeT2/gaT) - self.k3*(1. + nuT**2)
        A_bebeT[0, 1] =  0. 
        A_bebeT[0, 2] = -self.k3*4./3.*nuT 
        
        A_bebeT[1, 0] =  4./3.*nuT
        A_bebeT[1, 1] =  8.*(labeT2 - 1.)/gaT - self.k3*(1. + 0.5*nuT**2)
        A_bebeT[1, 2] =  1 + 0.5*nuT**2 
        
        A_bebeT[2, 0] = -self.k3*8./3.*nuT
        A_bebeT[2, 1] = -1. + 0.5*nuT**2
        A_bebeT[2, 2] =  8.*(labeT2 - 1.)/gaT - self.k3*(1. + 1.5*nuT**2)

        # Matriz de influencia del colectivo de cola sobre el batimiento
        A_beth0T = array([ 
            1. + nuT**2, 
            0., 
            8./3.*nuT ], type='Float64')

        # Matriz de influencia de la velocidad normal al rotor sobre el 
        # batimiento
        A_bela0T = array([
            4./3., 
            0., 
            2*nuT ], type='Float64')
        
        # Matriz de influencia de la torsion de cola sobre el batimiento
        A_bethtT = array([ 
            4*(1./5. + nuT**2/6.), 
            0., 
            2.*nuT ], type='Float64')

        # Despejamos el batimiento en funcion de la velocidad inducida
        invA_bebeT = la.inverse(A_bebeT)
        beth0T = dot(invA_bebeT, A_beth0T)
        bela0T = dot(invA_bebeT, A_bela0T)
        bethtT = dot(invA_bebeT, A_bethtT)
        
        # Estas funciones tambien seran utiles mas adelante, por eso no
        # se declaran dentro de h
        def f_be0T(la0T):
            return beth0T[0]*self.modelo.controles.th0T +\
                    bela0T[0]*(nuzT - la0T) +\
                    bethtT[0]*self.modelo.rotor_cola.tht
        
        def f_be1swT(la0T):
            return beth0T[2]*self.modelo.controles.th0T +\
                    bela0T[2]*(nuzT - la0T) +\
                    bethtT[2]*self.modelo.rotor_cola.tht
        
        # No interviene en el calculo de la0T pero si en el de las fuerzas
        def f_be1cwT(la0T):
            return beth0T[1]*self.modelo.controles.th0T +\
                    bela0T[1]*(nuzT - la0T) +\
                    bethtT[1]*self.modelo.rotor_cola.tht

        # Angulos de paso efectivos
        def f_th0T(la0T):
            return self.modelo.controles.th0T + self.k3*f_be0T(la0T)
        
        def f_th1swT(la0T):
            return self.k3*f_be1swT(la0T)

        # Coeficiente de sustentacion de la cola
        def f_cTT(la0T):
            return self.a0T*self.sT/2.*(f_th0T(la0T)*(1./3. + 0.5*nuT**2) +\
                    0.5*nuT*f_th1swT(la0T) +\
                    0.5*(nuzT - la0T))
        
        # Constantes de la TCMM
        ki  = (9./5.)**0.25
        knu = (5./4.)**0.25
        # Cuadrado de la velocidad inducida en el rotor
        def f_VT2(la0T):
            return (nuT/knu)**2 + (1./knu**2 - 1./ki**2)*(nuzT**2) +\
                    ((nuzT - la0T)/ki)**2

        def f_h(la0T):
            #return la0T/ki - f_cTT(la0T)/(2*sqrt(f_VT2(la0T)))
            return la0T/ki*2*sqrt(f_VT2(la0T)) - f_cTT(la0T)

        la0T = solveSecante(f_h, -0.2, 0.2, 100, 1e-12)
        # Calculamos otras magnitudes
        # Coeficiente de sustentacion
        cTT = f_cTT(la0T)
        # Coeficiente medio de resistencia
        deT = self.de0T + self.de2T*(cTT**2)
        # Coeficiente de par
        cQT = -(nuzT - la0T)*cTT + deT*self.sT/8.*(1. + 3*nuT**2)
        be0T = f_be0T(la0T)
        be1swT = f_be1swT(la0T)
        be1cwT = f_be1cwT(la0T)
        # Pasamos los angulos de batimiento a ejes rotor de cola
        be1cT = be1cwT*CchwT + be1swT*SchwT
        be1sT = be1swT*CchwT - be1cwT*SchwT
        
        ########################################################################
        #   CALCULO DE LAS FUERZAS DEL ROTOR DE COLA EN EJES CUERPO            #
        ########################################################################
        # Factor empirico de bloqueo del rotor cola debido al estabilizador 
        # vertical
        FT = 1. - 3./4*self.modelo.estabilizador_vertical.Sfn/( pi*self.RT**2 )
        # Dimension de fuerza para el rotor de cola
        dim = self.modelo.ro*( OmT*self.RT )**2*( pi*self.RT**2 )
        # Fuerzas y momentos con dimesiones
        TT = dim*FT*cTT
        QT = dim*self.RT*cQT
        # Fuerzas y momentos en ejes cola
        XT_T =  TT*be1cT
        YT_T = -TT*be1sT
        ZT_T = -TT
        LT_T =  0.0
        MT_T =  0.0
        NT_T =  QT
        # Fuerzas y momentos en ejes cuerpo
        XT_b, YT_b, ZT_b = dot(self.LBT, [XT_T, YT_T, ZT_T])
        LT_b, MT_b, NT_b = dot(self.LBT, [LT_T, MT_T, NT_T])
        # Fuerzas y momentos en ejes cuerpo respecto al centro de masas
        l = self.lT + self.modelo.rotor.xcg
        XT = XT_b
        YT = YT_b
        ZT = ZT_b
        LT = LT_b + YT_b*self.hT
        MT = MT_b - XT_b*self.hT + ZT_b*l
        NT = NT_b - YT_b*l
        
        FMT_0 = array([XT, YT, ZT, LT, MT, NT], type='Float64')
        FMT_1 = zeros(shape = (6, 16), type = 'Float64')
        
        ########################################################################
        #   HACEMOS NOTAR LOS CAMBIOS                                          #
        ########################################################################
        
        #----------------------------------------------------------------------#
        self.FMT = FMT_0, FMT_1
        #----------------------------------------------------------------------#    
        
        # Otros subsistemas necesitan esto
        # self.QT = -M_b[1]
        self.QT = QT
        
        # Logging
        self.deT = deT
        self.la0T = la0T
        self.cTT = cTT
        self.be1swT = be1swT
        self.be1cwT = be1cwT
        self.OmT = OmT
        self.XT_b, self.YT_b, self.ZT_b = XT, YT, ZT
        self.LT_b, self.MT_b, self.NT_b = LT, MT, NT
        self.nuT, self.nuzT = nuT, nuzT

