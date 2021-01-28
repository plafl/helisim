# -*- coding: latin-1 -*-
from numarray import *
import numarray.linear_algebra as la
from matematicas import *


class Rotor:
    def __init__(self, modelo):
        self.modelo = modelo
        self.params = {}
    
    def init(self):
        """ 
        Inicializa algunas variables que permanecen constantes y declara 
        algunos parametros. 
        """
        Cgas = cos(self.gas)
        Sgas = sin(self.gas)
        # Matriz de cambio de base ejes rotor-cuerpo
        self.LHB = array([
            [  Cgas, 0.0, Sgas ], 
            [  0.0,  1.0, 0.0  ],
            [ -Sgas, 0.0, Cgas ]], type = 'Float64' )

        self.LBH = transpose(self.LHB)
        
        # Solidez del rotor
        self.s = self.Nb*self.c/(pi*self.R) 
        
        # Coeficientes de la TCMM
        self.ki =  (9./5.)**0.25
        self.knu = (5./4.)**0.25

        # Parametros de la solución numérica de la velocidad inducida
        self.params['la0'] = {'ila0': -1.0, 'Nmax': 50, 'eps': 1e-6}


    def calcula_FMR(self):
        """ 
        Calcula las fuerzas y momentos del rotor
        
        NECESITA:
            modelo.Vrw_b
            controles.
                  th0
                  th1c
                  th1s
        
        PROVEE:
            rotor.FMR
            rotor.la0
            rotor.xi
            rotor.QR
        """
        # Unos alias para mejorar la legibilidad
        a0, s = self.a0, self.s
        Om = self.modelo.Om 
        
        # Si el rotor no gira nada de esto tiene sentido y ademas no se puede
        # adimensionalizar con Om
        if Om < 1:
            self.FMR = (
                    zeros(shape=(6, ),   type='Float64'),
                    zeros(shape=(6, 16), type='Float64') )
            self.la0 = 0.0
            self.xi = 0.0
            self.QR = ( 
                    0.0,
                    zeros(shape=(16, ), type='Float64') )
            return 
        ########################################################################
        #   CAMBIO A EJES VIENTO DEL ROTOR                                     #
        ########################################################################
        
        # Velocidad angular en ejes cuerpo
        W_b = array([ self.modelo.p, self.modelo.q, self.modelo.r ], 
                type='Float64')
        
        # Velocidad del rotor relativa al viento, en ejes cuerpo
        Vhrw_b = self.modelo.Vrw_b +\
                cross(W_b, [-self.xcg, 0.0, -self.hR])

        # en ejes rotor
        uhrw_h, vhrw_h, whrw_h = dot(self.LHB, Vhrw_b)
        
        # Velocidad angular en ejes rotor
        p_h, q_h, r_h = dot(self.LHB, W_b)
        
        # Velocidad relativa del rotor respecto al viento, en ejes viento
        uhrw_w = sqrt(uhrw_h**2 + vhrw_h**2)
        whrw_w = whrw_h
        
        # Angulo que forman los ejes viento con los ejes rotor. 
        # Si no hay velocidad se supone 0.
        if uhrw_w!=0:
            Cchw = uhrw_h/uhrw_w
            Schw = vhrw_h/uhrw_w
        else:
            # Por convenio chw = 0 si u = v = 0
            if vhrw_h>0:
                Cchw = 0.0
                Schw = 1.0
            elif vhrw_h<0:
                Cchw =  0.0
                Schw = -1.0
            else:
                Cchw = 1.0
                Schw = 0.0
        
        # Matriz cambio de base de ejes viento a ejes cuerpo
        # nos es util para pasar las fuerzas y momentos a ejes
        # cuerpo
        LHW = array([
            [ Cchw, -Schw, 0.0 ],
            [ Schw,  Cchw, 0.0 ],
            [ 0.0,   0.0,  1.0 ]], type = 'Float64')
        LBW = dot(self.LBH, LHW)
        
        # Cambiamos a ejes viento los angulos de paso del rotor
        th0 = self.modelo.controles.th0
        tht = self.tht
        th1c_w = self.modelo.controles.th1c*Cchw -\
                self.modelo.controles.th1s*Schw
        th1s_w = self.modelo.controles.th1c*Schw +\
                self.modelo.controles.th1s*Cchw

        # Parametros de avance 
        nu  = uhrw_w/( Om*self.R )
        nuz = whrw_w/( Om*self.R )
        
        # Velocidad angular adimensionalizada en ejes viento
        p_w = (  Cchw*p_h + Schw*q_h )/Om
        q_w = ( -Schw*p_h + Cchw*q_h )/Om

        # Las derivadas respecto al tiempo del angulo chw y de las
        # velocidad son matrices, que al multiplicar por la derivada del
        # vector de estado devuelven el valor que corresponda.

        # Derivada respecto al tiempo del angulo chw, en funcion del
        # vector de estado
        m = uhrw_h**2 + vhrw_h**2
        # Por convenio DTchw = 0 si u=v=0
        DTchw = zeros( shape = (16,), type = 'Float64' )
        if m!=0:
            DTchw[:2] = 1./m*dot( [-vhrw_h, uhrw_h ], self.LHB[:2, :2] )
        # Derivada respecto al tiempo de W_h, en funcion del vector de
        # estado
        DTW_h = zeros(shape = (3, 16), type = 'Float64')
        DTW_h[:, 3:6] = self.LHB
        
        # Derivada respecto al azimuth de p_w y q_w (adimensionalizadas), 
        # en funcion del  vector de estado
        DChip_w = 1./(Om**2)*( DTchw*( -Schw*p_h + Cchw*q_h ) +\
                Cchw*DTW_h[0] + Schw*DTW_h[1])
        DChiq_w = 1./(Om**2)*( DTchw*( -Cchw*p_h - Schw*q_h ) -\
                Schw*DTW_h[0] + Cchw*DTW_h[1])

        ########################################################################
        #   CALCULO DE LA VELOCIDAD INDUCIDA UNIFORME                          #
        ########################################################################
        # Coeficientes de la TCMM
        ki, knu = self.ki, self.knu
        # z = altura del rotor sobre el suelo
        # En el instante inicial con helicoptero tocando suelo z_i = 0
        # h0 = distancia rotor a punto mas bajo del helicoptero
        z = self.modelo.alt - self.modelo.hG + self.h0

        # Resolvemos h(la0) = 0 mediante el metodo de la secante
        def h(la0):
            # Efecto suelo
            def K0G(la0):
                # Si no el helicoptero no se estrella nunca!!!
                if z>self.R/3. and la0!=0.0:
                    return 1. - 1./(16.*(z/self.R)**2)/(1 + (nu/la0)**2)
                else:
                    return 1.0
            # Coeficiente de sustentacion del rotor
            def cT(la0):
                return a0*s/2.*(th0*(1./3. + nu**2/2.) +\
                        nu/2.*(th1s_w + p_w/2.) + K0G(la0)*(nuz - la0)/2. +\
                        1./4.*(1. + nu**2)*tht)
            # Velocidad en el rotor al cuadrado
            def VT2(la0):
                return (nu/knu)**2 + (1./knu**2 - 1./ki**2)*nuz**2  +\
                        ((nuz - la0)/ki)**2
            # Ecuacion que combina la TEP con la TCMM
            # return la0/ki - cT(la0)/(2*sqrt(VT2(la0)))
            return 2*sqrt(VT2(la0))*la0/ki - cT(la0)
        # la0NG no incluye efecto suelo
        la0NG = solveSecante(h, -0.2, 0.2, 100, 1e-12)
        # Si la iteracion no ha convergido asignamos un valor razonable
        # para que no se note :)
        if la0NG == None:
            la0NG = 0.06

        ########################################################################
        #   CALCULO DEL BATIMIENTO                                             #
        ########################################################################
        # VT = Velocidad del aire en el rotor
        VT2 = (nu/knu)**2 + (1/knu**2 - 1/ki**2)*(nuz**2) +\
                ((nuz - la0NG)/ki)**2
        VT = sqrt(VT2)
        # Velocidad del aire en el rotor corregida para terminos lineales
        barV = VT*(1 - la0NG*(nuz - la0NG)/ki**2/VT2)
        # Efecto suelo para la componente no lineal de la velocidad inducida
        if z>self.R/10:
            K0G = 1 - 1./(16.*(z/self.R)**2)/(1 + (nu/la0NG)**2)
        else:
            K0G = 1.0
        # Velocidad con efecto suelo
        la0 = K0G*la0NG + nuz*(1-K0G)
        # Angulo de la estela
        xi = ang(abs(la0*Om*self.R - whrw_w), uhrw_w)
        # Numero de Locke, varia con la densidad
        ga = self.modelo.ro*self.c*a0*self.R**4/self.Ibe
        # Frecuencia natural de batimiento, al cuadrado
        # Varia con la velocidad de rotacion del rotor
        labe2 = 1. + self.Kbe/( self.Ibe*Om**2)

        # Matrices de las 3 ecuaciones de batimiento
        # Lado izquierdo
        A_bebe = array([
            [ 8*labe2/ga,  0,              0                ],
            [ 4./3.*nu,    8*(labe2-1)/ga, 1 + 0.5*nu**2    ],
            [ 0,          -1. + 0.5*nu**2, 8*(labe2 - 1)/ga ]], type='Float64')
        
        A_bela = array([
            [ 2./3.*nu, 0  ],
            [ 0.,       1. ],
            [ 1.,       0. ]], type='Float64')
        
        # Lado derecho
        # Influencia del paso sobre el batimiento
        A_beth = array([
            [ 1. + nu**2, 4./5. + 2./3.*nu**2, 4./3.*nu,       0.             ],
            [ 0.,         0.,                  0.,             1. + 0.5*nu**2 ],
            [ 8./3.*nu,   2*nu,                1. + 1.5*nu**2, 0              ]
            ], type='Float64')
        
        # Influencia de la aceleracion angular sobre el batimiento
        A_beDChiom = 8./ga*array([
            [0., 0.],
            [0., 1.],
            [1., 0.]], type='Float64')

        # Influencia de la velocidad angular sobre el batimiento
        A_beom = array([
            [ 2./3.*nu, 0      ],
            [ 16./ga,   1.     ],
            [ 1.,      -16./ga ]], type='Float64')
        
        # Influencia de la velocidad inducida media sobre el batimiento
        A_bela0 = array([
            4./3.,
            0.,
            2*nu], type='Float64')

        # Matrices de las 2 ecuaciones de velocidad inducida
        # Matriz L de Pitt-Peters
        X = tan(abs(xi/2.))
        K = array([
            [ 0,                 2*(1. + X**2)/barV, 0.                  ],
            [ 15.*pi/(64.*VT)*X, 0.,                 2.*(1. - X**2)/barV ]], 
            type='Float64')
        
        # Matriz de fuerzas y momentos aerodinamicos
        # N = { F0/2, M1s, M1c}
        Nth = array([
            [ 2./3. + nu**2, 0.5*(1. + nu**2), nu,                    0. ],
            [ 2./3.*nu,      0.5*nu,           0.25*(1. + 1.5*nu**2), 0  ],
            [ 0.,            0., 0.,           0.25*(1. + 0.5*nu**2)     ]], 
            type='Float64')

        Nbe = array([
            [  0.,    0.,                     0.                    ],
            [  0.,    0.25*(1. - 0.5*nu**2),  0.                    ],
            [ -nu/3., 0,                     -0.25*(1. + 0.5*nu**2) ]], 
            type='Float64')

        Nla = -0.5*array([
            [ nu,  0. ],
            [ 0.5, 0. ],
            [ 0.,  0.5]], type='Float64')

        Nom = -Nla

        Nla0 = array([ 
            1., 
            0.5*nu, 
            0. ], type='Float64')
        
        # Lado izquierdo de las ecuaciones de la velocidad inducida
        m = a0*s/4.
        A_labe = -m*dot(K, Nbe)
        A_lala = identity(2) - m*dot(K, Nla)
        # Matriz de influencia de controles, velocidad angular y velocidad
        # inducida uniforme sobre la velocidad inducida
        A_lath, A_laom, A_lala0 = [m*dot(K, x) for x in [Nth, Nom, Nla0]]
        # Matriz de influencia de la aceleracion angular sobre la velocidad 
        # inducida
        A_laDChiom = zeros(shape=(2,2), type='Float64')
        
        # Ensamblamos el sistema
        # Lado izquierdo de las ecuaciones
        lhs = array(shape=(5, 5), type='Float64')
        lhs[0:3, 0:3], lhs[0:3, 3:5] = A_bebe, A_bela
        lhs[3:5, 0:3], lhs[3:5, 3:5] = A_labe, A_lala

        # Lado derecho de las ecuaciones
        # Terminos en funcion del vector de estado
        rhs_0 = dot(concatenate([A_beth, A_lath]), [th0, tht, th1s_w, th1c_w]) \
                + dot(concatenate([A_beom, A_laom]), [p_w, q_w]) \
                + (nuz - la0)*concatenate([A_bela0, A_lala0]) 

        # Terminos dependientes linealmente de la derivada del vector de
        # estado
        rhs_1 = dot(concatenate([A_beDChiom, A_laDChiom]), 
                array([DChip_w, DChiq_w]))
        
        # Y resolvemos el batimiento
        lhs_inv = la.inverse(lhs)
        # {be0, be1c_w, be1s_w, la1s_w, la1c_w} = bat[0] + bat[1]DTx
        bat = [ dot(lhs_inv, rhs_0), dot(lhs_inv, rhs_1) ]
        
        ########################################################################
        #   CALCULO DE LAS FUERZAS DEL ROTOR                                   #
        ########################################################################
        # Alias
        # Para los resultados del batimiento 
        be0, be1c_w, be1s_w, la1s_w, la1c_w = bat[0]
        # Angulo de ataque efectivo en el perfil de la pala
        al1s_w = p_w - la1s_w + be1c_w + th1s_w
        al1c_w = q_w - la1c_w - be1s_w + th1c_w
        
        # F1 = Fuerza aerodinámica normal al rotor
        # F1 = int(UT^2*th + UP*UT, r, 0, 1)
        # F1 = F10 + F11c*cos(ch) + F11s*sin(ch) + ...
        F10 = th0*( 1./3. + 0.5*nu**2 ) + 0.5*nu*( th1s_w + p_w/2. ) +\
                0.5*( nuz - la0 ) + 0.25*( 1. + nu**2)*tht
        # Coeficiente de resistencia medio
        de = self.de0 + self.de2*( 1./2.*a0*s*F10 )**2
        # Coeficiente de sustentacion
        cT = F10/2*a0*s
        # Primeros armonicos 
        F11s = al1s_w/3. + nu*( th0 + nuz - la0 + 2./3.*tht )
        F11c = al1c_w/3. - 0.5*nu*be0
        # Segundos armonicos
        F12s =  0.5*nu*( 0.5*al1c_w + 0.5*( th1c_w - be1s_w ) - nu*be0 )
        F12c = -0.5*nu*( 0.5*al1s_w + 0.5*( th1s_w + be1c_w ) +\
                nu*( th0 + 0.5*tht ) )
        
        # F2 = Componente tangencial de la fuerza aerodinamica
        # F2 = int(UP*UT*th + UP^2 - de*UT^2/a0, r, 0, 1)
        # F2 = F20 + F21c*cos(ch) + F21s*sin(ch)
        # Solo necesitamos los primeros armonicos
        F21s = 0.5*(nu**2)*be0*be1s_w  + ( nuz - la0 - 0.25*nu*be1c_w ) *\
                   ( al1s_w - th1s_w ) - 0.25*nu*be1s_w*( al1c_w - th1c_w ) +\
               th0*(
                 (al1s_w - th1s_w)/3. + 
                  nu*( nuz - la0 ) -
                  0.25*(nu**2)*be1c_w
               ) +\
               tht*( 
                     (al1s_w - th1s_w)/4. + 
                  0.5*nu*(nuz - la0 - 0.25*nu*be1c_w) 
               )+ \
               th1s_w*( 
                    0.5*( nuz - la0 ) + 
                    nu*( 3./8.*( p_w - la1s_w ) + 
                    0.25*be1c_w) 
                  ) +\
               th1c_w*0.25*nu*( 
                        0.5*(q_w - la1c_w) - 
                    be1s_w - 
                    nu*be0 
                      )-\
               de*nu/a0
        
        F21c = ( al1c_w - th1c_w - 2*be0*nu )*( nuz - la0 - 0.75*nu*be1c_w )-\
               0.25*nu*be1s_w*( al1s_w - th1s_w ) + \
               th0*( 
                 (al1c_w - th1c_w)/3. - 
                  0.5*nu*( be0 + 0.5*nu*be1s_w ) 
               ) +\
               tht*( 
                     0.25*(al1c_w - th1c_w) - 
                 nu*(be0/3. + 0.125*nu*be1s_w) 
               ) +\
               th1c_w*( 
                        0.5*( nuz - la0 ) + 
                0.25*nu*( 0.5*( p_w - la1s_w ) - 
                be1c_w ) 
                  ) +\
               th1s_w*0.25*nu*( 
                        0.5*( q_w - la1c_w ) -
                    be1s_w - 
                    nu*be0 
                      )

        # Coeficientes de fuerzas
        # _Cxw = 2*Cxw/(a0*s)
        # _Cyw = 2*Cyw/(a0*s)
        # _Czw = 2*Czw/(a0*s)
        _Cxw = ( F10/2. + F12c/4.)*be1c_w + F11c/2.*be0 + F12s/4.*be1s_w +\
                F21s/2.
        _Cyw = (-F10/2. + F12c/4.)*be1s_w - F11s/2.*be0 - F12s/4.*be1c_w +\
                F21c/2.
        _Czw = -F10
        # Dimensiones de las fuerzas * a0s/2
        dim_F = (a0*s)/2.*self.modelo.ro*((Om*self.R)**2)*pi*(self.R**2)
        # Componentes en ejes viento de las fuerzas sobre el rotor
        # Parte estacionaria
        F_w_0 = array([ 
                dim_F*_Cxw,
                dim_F*_Cyw,
                dim_F*_Czw ], type='Float64')
        
        # Para las fuerzas no hay terminos en aceleraciones angulares
        # Parte no estacionaria
        F_w_1 = zeros(shape=(3, 16), type='Float64')
        
        ########################################################################
        #   CALCULO DE LOS MOMENTOS DEL ROTOR                                  #
        ########################################################################
        # Momentos en ejes viento, lado derecho
        M_w_0 = zeros( shape = (3, ), type = 'Float64' )
        # Parece ser que esta formula esta mal en el Padfield, mas abajo la 
        # formula coregida, segun viene en el  Bramwell
        # M_w_0[2] = self.R*( ( nuz - la0 )*F_w_0[2] + nu*F_w_0[0] + \
        # dim_F*de/( 4.*a0 )*( 1 + 7./3.*nu**2 ) )
        M_w_0[2] =  self.R*( ( nuz - la0 )*F_w_0[2] + nu*F_w_0[0] +\
            dim_F*de/( 4.*self.a0 )*( 1 + 3.*nu**2 ) )
        # Debido a la inclinacion del disco el par de motor contribuye un 
        # 10% aproximadamente
        M_w_0[0] = -self.Nb/2.*self.Kbe*be1s_w - M_w_0[2]/2.*be1c_w
        M_w_0[1] = -self.Nb/2.*self.Kbe*be1c_w + M_w_0[2]/2.*be1s_w

        # Los momentos tienen terminos en DTx
        M_w_1 = array( shape = (3, 16), type = 'Float64' )
        M_w_1[0, :] = -self.Nb/2.*self.Kbe*bat[1][2]
        M_w_1[1, :] = -self.Nb/2.*self.Kbe*bat[1][1]
        M_w_1[2, :] =  zeros(shape=(16,), type='Float64')
        M_w_1[2, 10] = self.modelo.motor.Irt[self.modelo.motor.n] +\
            self.Nb*self.Ibe

        ########################################################################
        #   CALCULO DE FM EN EJES CUERPO                                       #
        ########################################################################
        # Pasamos de ejes viento a ejes cuerpo
        # Estas son las fuerzas sobre la cabeza del rotor, pero ya en ejes 
        # cuerpo
        F_b_0 = dot(LBW, F_w_0)
        F_b_1 = dot(LBW, F_w_1)
        M_b_0 = dot(LBW, M_w_0)
        M_b_1 = dot(LBW, M_w_1)
        # Utilizamos notacion matricial para el producto vectorial para
        # poder pasar los terminos que dependen de la derivada del
        # vector de estado (que son matrices)
        r_cross = array([
                [ 0.0,      self.hR,  0.0      ],
                [-self.hR,  0.0,      self.xcg ],
                [ 0.0,     -self.xcg, 0.0      ]], type='Float64')
        
        # Trasladamos al centro de masas
        F_0 = F_b_0
        F_1 = F_b_1
        M_0 = M_b_0 + dot(r_cross, F_b_0)
        M_1 = M_b_1 + dot(r_cross, F_b_1)
        
        # FM: Fuerzas y momentos conjuntos en un solo vector
        FM_0 = array(shape = (6,), type = 'Float64')
        FM_0[0:3] = F_0
        FM_0[3:6] = M_0
        
        FM_1 = array(shape = (6, 16), type = 'Float64')
        FM_1[0:3, :] = F_1
        FM_1[3:6, :] = M_1
        
        ########################################################################
        #   HACEMOS NOTAR LOS CAMBIOS                                          #
        ########################################################################
        
        #----------------------------------------------------------------------#
        self.FMR = FM_0, FM_1
        #----------------------------------------------------------------------#    
        
        # Otras magnitudes usadas en otros subsistemas
        # Angulo de la estela, para el estabilizador horizontal, respecto al 
        # eje z_b
        self.xi = xi + self.gas
        # La velocidad inducida, para el estabilizador horizontal y en
        # un futuro tambien podria ser para el fuselaje
        self.la0 = la0
        # El par del rotor, para el motor
        self.QR = M_w_0[2], M_w_1[2, :]
        
        # Logging de magnitudes intermedias
        self.Xr_w, self.Yr_w, self.Zr_w = F_w_0
        self.Lr_w, self.Mr_w, self.Nr_w = M_w_0
        self.Xr_b, self.Yr_b, self.Zr_b = F_b_0
        self.Lr_b, self.Mr_b, self.Nr_b = M_0
        self.Xr_i, self.Yr_i, self.Zr_i = dot(self.modelo.LIB, F_b_0)
        self.Lr_i, self.Mr_i, self.Nr_i = dot(self.modelo.LIB, M_b_0)

        self.nu, self.nuz = nu, nuz
        self.Cchw, self.Schw = Cchw, Schw
        self.cT, self.de = cT, de
        self.K0G, self.VT = K0G, VT
        self.be0, self.be1c_w, self.be1s_w = bat[0][:3]
        self.la1s_w, self.la1c_w = bat[0][3:]
        self.DTchw, self.DChip_w, self.DChiq_w = DChip_w, DChiq_w, DTchw
        self.th1c_w, self.th1s_w = th1c_w, th1s_w
        self.A_bebe, self.A_bela, self.A_labe, self.A_lala =\
                A_bebe, A_bela, A_labe, A_lala
        self.A_beth, self.A_lath, self.A_beom, self.A_laom =\
                A_beth, A_lath, A_beom, A_laom
        self.A_bela0, self.A_lala0, self.A_beDChiom, self.A_laDChiom =\
                A_bela0, A_lala0, A_beDChiom, A_laDChiom
        self.A = lhs
        self.A_th = concatenate([A_beth, A_lath])
        self.A_om = concatenate([A_beom, A_laom])
        self.A_la0 = concatenate([A_bela0, A_lala0])
        self.A_DChiom = concatenate([A_beDChiom, A_laDChiom])
        self.barV = barV
        self.labe2, self.ga = labe2, ga
        self.p_w, self.q_w = p_w, q_w
        
