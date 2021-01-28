# -*- coding: latin-1 -*-
"""
Implementa la clase Modelo, encargada del calculo de las fuerzas del helicóptero
"""

from pylab import *
from numarray import *
from matematicas import *
import numarray.linear_algebra as la

from rotor import Rotor
from rotor_cola import RotorCola
from motor import Motor
from fuselaje import Fuselaje
from estabilizador import EstabilizadorHorizontal, EstabilizadorVertical
from controles import Controles
from geodesia import Geoide
from colision import Tren
from historial import *

class Modelo:
    """
    La clase principal de la simulacion:

    1º Contiene a todos los subsistemas:
        rotor
        rotor_cola
        fuselaje
        estabilizador_horizontal
        estabilizador_vertical
        motor
        controles
    
    2º Calcula las ecuaciones del tipo
        ecs_modelo[1]*(dx/dt) = ecs_modelo[0]

    3º Para llamar a calcula_ecs_modelo hay que asegurarse de tener
    los valores correctos en el vector de estado:
        Velocidad en ejes cuerpo
            self.u
            self.v
            self.w
        Velocidad angular
            self.p
            self.q
            self.r
        Cuaternio
            self.q0
            self.q1
            self.q2
            self.q3
        Controles:
            self.controles.th0
            self.controles.th1c
            self.controles.th1s
        Ambiente:
            self.ro
            self.hG
            self.Vw_i
        Posicion:
            self.x_i
            self.y_i
            self.z_i
        Motor/Rotor:
            self.Om
            self.Q1
            self.DTQ1
    """
    
    def __init__(self, reloj=None):
        """
        Crea los subsistemas y un reloj para controlar el tiempo
        """

        self.rotor = Rotor(self)
        self.rotor_cola = RotorCola(self)
        self.motor = Motor(self)
        self.fuselaje = Fuselaje(self)
        self.estabilizador_horizontal = EstabilizadorHorizontal(self)
        self.estabilizador_vertical = EstabilizadorVertical(self)
        self.controles = Controles(self)
        self.geoide = Geoide()
        self.tren = Tren(self)
        
        if reloj == None:
            self.reloj = Reloj()
        else:
            self.reloj = reloj
        
        # control de tiempo
        self.t = 0.0
        # altura del suelo
        self.hG = 0.0
        # Inicializamos con un valor cualquiera
        self.geoide.ejesLocales(0.0, 0.0, 0.0)
        
        self.trim = False
        
        # Cada variable de los subsistemas tiene una cadena de texto que el
        # simbolo mediante el cual se le hace referencia inequivocamente.
        self.__logstr = {
                "Om":       ("modelo", "Om"),
                "u":        ("modelo", "u"),
                "v":        ("modelo", "v"),
                "w":        ("modelo", "w"),
                "p":        ("modelo", "p"),
                "q":        ("modelo", "q"),
                "r":        ("modelo", "r"),
                "h":        ("modelo", "h"),
                "x_i":      ("modelo", "x_i"),
                "y_i":      ("modelo", "y_i"),
                "z_i":      ("modelo", "z_i"),
                "q0":       ("modelo", "q0"),
                "q1":       ("modelo", "q1"),
                "q2":       ("modelo", "q2"),
                "q3":       ("modelo", "q3"),
                "th":       ("modelo", "th"),
                "fi":       ("modelo", "fi"),
                "ch":       ("modelo", "ch"),
                "alt":      ("modelo", "alt"),
                "lon":      ("modelo", "lon"),
                "lat":      ("modelo", "lat"),
                
                "la0":      ("rotor", "la0"),
                "la1c_w":   ("rotor", "la1c_w"),
                "la1s_w":   ("rotor", "la1s_w"),
                "XR_b":     ("rotor", "Xr_b"),
                "YR_b":     ("rotor", "Yr_b"),
                "ZR_b":     ("rotor", "Zr_b"),
                "LR_b":     ("rotor", "Lr_b"),
                "MR_b":     ("rotor", "Mr_b"),
                "NR_b":     ("rotor", "Nr_b"),
                "cT":       ("rotor", "cT"),
                "th1c_w":   ("rotor", "th1c_w"),
                "th1s_w":   ("rotor", "th1s_w"),
                "be0":      ("rotor", "be0"),
                "be1c_w":   ("rotor", "be1c_w"),
                "be1s_w":   ("rotor", "be1s_w"),
                
                "Xf_b":     ("fuselaje", "Xf_b"),
                "Yf_b":     ("fuselaje", "Yf_b"),
                "Zf_b":     ("fuselaje", "Zf_b"),
                "Lf_b":     ("fuselaje", "Lf_b"),
                "Mf_b":     ("fuselaje", "Mf_b"),
                "Nf_b":     ("fuselaje", "Nf_b"),
                "alf":      ("fuselaje", "al"),
                "bef":      ("fuselaje", "be"),

                "th0":      ("controles", "th0"),
                "th1c":     ("controles", "th1c"),
                "th1s":     ("controles", "th1s"),
                "th0T":     ("controles", "th0T"),
                "et0p":     ("controles", "l_et0p"),
                "et1sp":    ("controles", "l_et1sp"),
                "et1cp":    ("controles", "l_et1cp"),
                "etpp":     ("controles", "l_etpp"),

                "Q1":       ("modelo", "Q1"),
                "DTQ1":     ("modelo", "DTQ1"),
                
                "QT":       ("rotor_cola", "QT"),
                "la0T":     ("rotor_cola", "la0T"),
                "cTT":      ("rotor_cola", "cTT"),
                "XT_b":     ("rotor_cola", "XT_b"),
                "YT_b":     ("rotor_cola", "YT_b"),
                "ZT_b":     ("rotor_cola", "ZT_b"),
                "LT_b":     ("rotor_cola", "LT_b"),
                "MT_b":     ("rotor_cola", "MT_b"),
                "NT_b":     ("rotor_cola", "NT_b"),
                
                "Xfn_b":    ("estabilizador_vertical", "Xfn_b"),
                "Yfn_b":    ("estabilizador_vertical", "Yfn_b"),
                "Zfn_b":    ("estabilizador_vertical", "Zfn_b"),
                "Lfn_b":    ("estabilizador_vertical", "Lfn_b"),
                "Mfn_b":    ("estabilizador_vertical", "Mfn_b"),
                "Nfn_b":    ("estabilizador_vertical", "Nfn_b"),
                "alfn":     ("estabilizador_vertical", "al"),
                "befn":     ("estabilizador_vertical", "be"),
                
                "Xtp_b":    ("estabilizador_horizontal", "Xtp_b"),
                "Ytp_b":    ("estabilizador_horizontal", "Ytp_b"),
                "Ztp_b":    ("estabilizador_horizontal", "Ztp_b"),
                "Ltp_b":    ("estabilizador_horizontal", "Ltp_b"),
                "Mtp_b":    ("estabilizador_horizontal", "Mtp_b"),
                "Ntp_b":    ("estabilizador_horizontal", "Ntp_b"),
                "altp":     ("estabilizador_horizontal", "al"),
                "betp":     ("estabilizador_horizontal", "be"),
        }
        


    ########################################################################### 
    #   INICIALIZACION DE VARIABLES                                           #
    ########################################################################### 
    def init(self):
        """
        Inicializa los subsistemas
        """

        self.rotor.init()
        self.rotor_cola.init()
        self.motor.init()
        self.fuselaje.init()
        self.estabilizador_horizontal.init()
        self.tren.init()
        self.controles.init()

    ###########################################################################
    #   UTILIDADES                                                            # 
    ###########################################################################
    def euler(self, ch, th, fi):
        """
        Dados los angulos de Euler actualiza el cuaternio de rotacion
        """
        
        ( self.q0, 
          self.q1, 
          self.q2, 
          self.q3 )= Quaternion.euler(
                  ch, th, fi)*Quaternion.rot(pi, 0.0, 1.0, 0.0)

    
    ########################################################################### 
    #   PASO DE TIEMPO                                                        #
    ########################################################################### 
    def paso(self, x, t):
        """ 
        Devuelve el vector f(x, t) tal que: dx/dt = f(x,t). 
        Esta rutina es necesaria para el integrador numerico.
        """
        
        self.t = t
        ( self.u, 
          self.v, 
          self.w, 
          self.p, 
          self.q, 
          self.r, 
          self.q0, 
          self.q1, 
          self.q2,
          self.q3, 
          self.Om, 
          self.DTQ1, 
          self.Q1, 
          self.x_i, 
          self.y_i, 
          self.z_i ) = x
        
        self.calcula_ecs_modelo()
        ecs_0, ecs_1 = self.ecs_modelo
        self.f = dot(
            la.inverse(ecs_1), 
            ecs_0
        )
            
        return self.f
                
    def avanza(self, t):
        """
        Avanza el estado del modelo hasta el instante de tiempo t
        """
        if self.tren.colision == 0:
            # Guardamos el estado para el tren antes de dar un paso de integracion
            self.tren.x[0:6] = (self.u, self.v, self.w, self.p, self.q, self.r)
            
            self.tren.x[6] = self.p0
            self.tren.x[7] = self.p1
            self.tren.x[8] = self.p2
            self.tren.x[9] = self.p3

            self.tren.x[10] = self.x_l
            self.tren.x[11] = self.y_l
            self.tren.x[12] = self.z_l - self.hG

            self.tren.t = self.t
        
        self.tren.FC = [0.0, 0.0, 0.0] 
        self.tren.MC = [0.0, 0.0, 0.0]
        self.tren.FB = self.F0 + self.M*self.G
        self.tren.MB = self.M0

        # Damos un paso de integracion
        self.integrador(t)
        # Damos un paso de colision
        xc0, yc0, zc0 = self.tren.x[10:13]
        self.tren.step(t)
        xc1, yc1, zc1 = self.tren.x[10:13]
        # tren.step afecta a:
        # tren.x, tren.t, tren.F, tren.colision

        
        if self.tren.colision:
            q0 = self.tren.x[6]
            q1, q2, q3 = dot(self.LIL, self.tren.x[7:10])
            
            # !! Parece ser que se acumula un error numerico
            # El codigo es formalmente correcto
#            xc, yc, zc = dot(self.LCL, 
#                    [self.tren.x[10], 
#                        self.tren.x[11], 
#                        self.tren.x[12] + self.hG]) + self.Ol
#            xi, yi, zi = self.geoide.cart2loc(xc, yc, zc)
            
            # Para evitar el error numerico calulamos los desplazamientos 
            # en ejes inerciales
            dxi, dyi, dzi = dot(self.LIL, [xc1 - xc0, yc1 - yc0, zc1 - zc0])
            xi = self.x_i + dxi
            yi = self.y_i + dyi
            zi = self.z_i + dzi
            
            # No es correcto, pero por lo menos no se acumula
            # error
#            xi, yi, zi = dot(self.LIL, 
#                    [self.tren.x[10], 
#                        self.tren.x[11], 
#                        self.tren.x[12] + self.hG])

            x = array(shape=(16,), type='Float64')
            x[0:6  ] = self.tren.x[0:6]
            x[6:10 ] = (q0, q1, q2, q3)
            x[13:16] = (xi, yi, zi)
            x[10   ] = self.integrador.x0[10]
            x[11   ] = self.integrador.x0[11]
            x[12   ] = self.integrador.x0[12]

            ( self.u, self.v, self.w, 
                    self.p, self.q, self.r,
                    self.q0, self.q1, self.q2, self.q3,
                    self.Om, self.DTQ1, self.Q1,
                    self.x_i, self.y_i, self.z_i ) = x

            self.t = self.tren.t
            
            self.calcula_preparativos()
            self.calcula_FM()

            F = array(shape=(16, ), type='Float64')
            F[0:6  ] = self.tren.F[0:6]
            F[6    ] = self.tren.F[6]
            F[7:10 ] = dot(self.LIL, self.tren.F[7:10 ])
            F[13:16] = dot(self.LIL, self.tren.F[10:13])
            F[10   ] = self.integrador.F0[10]
            F[11   ] = self.integrador.F0[11]
            F[12   ] = self.integrador.F0[12]
            
#           self.integrador.init(self.paso, self.tren.t, x)
            self.integrador.x0 = x
            self.integrador.F0 = F
            self.integrador.t0 = self.tren.t
            self.integrador.ci = 0
        
            
        # Comprobamos si se producen colisiones
#        self.colisiones()
        
        self.t = t
        
    def reset(self):
        """
        Coloca el modelo en unas condiciones estandard.
        """
        self.u = 0.0
        self.v = 0.0
        self.w = 0.0
        self.p = 0.0
        self.q = 0.0
        self.r = 0.0

        self.Om = 0.0
        self.Q1 = 0.0
        self.DQ1 = 0.0

        self.x_i = 0.0
        self.y_i = 0.0
        self.z_i = self.fuselaje.ht

        self.nmotores = 0.0

    def calcula_x(self):
        """
        Actualiza el vector de estado.
        """

        self.x = array([
            self.u, 
            self.v, 
            self.w,
            self.p, 
            self.q, 
            self.r,
            self.q0, 
            self.q1, 
            self.q2, 
            self.q3,
            self.Om, 
            self.DTQ1, 
            self.Q1,
            self.x_i, 
            self.y_i, 
            self.z_i], type='Float64')

    def calcula_LBI(self):
        """
        Calcula la matriz de cambio de ejes inercia a ejes cuerpo y 
        viceversa.
        """

        # La matriz de cambio de base es la transpuesta de la matriz de
        # rotacion
        self.LIB = Quaternion(
                self.q0, self.q1, self.q2, self.q3).toMatrix()
        self.LBI = transpose(self.LIB)

    def calcula_G(self):
        """
        Calcula la aceleracion de la gravedad en ejes cuerpo.
        """

        self.Gx, self.Gy, self.Gz = self.G = dot(self.LBI, [0.0, 0.0, -9.81])
    
    def calcula_Vw_i(self):
        """
        Calcula la velocidad del viento en ejes inerciales a partir
        de modulo y azimuth
        """
        self.Vw_i = dot(self.LIL, 
                [ self.viento_vel*cos(self.viento_dir), 
                  self.viento_vel*sin(self.viento_dir),
                  0 ])

    def calcula_Vrw_b(self):
        """
        Velocidad relativa al viento en ejes cuerpo.
        """
        self.Vrw_b = [self.u, self.v, self.w] - dot(self.LBI, self.Vw_i)
    

    def calcula_ro(self):
        """
        Calcula la densidad a considerando al aire como gas perfecto.
        """
        self.T = self.T0 - 0.0065*self.alt
        self.P = self.P0*(self.T/self.T0)**5.2586
        self.ro = self.P/(self.T*286.9)

    def calcula_velocidades(self):
        """
        Diversas velocidades para pasarlas a FlightGear.
        """
        p = self.p
        q = self.q
        r = self.r

        Sfi, Cfi = sin(self.fi), cos(self.fi)
        Sth, Cth = sin(self.th), cos(self.th)

        if Cth!=0:
            self.DTfi = p + (q*Sfi + r*Cfi)*Sth/Cth
        else:
            self.DTfi = 0.0

        self.DTth = q*Cfi - r*Sfi
        if Sth!=0:
            self.DTch = (q*Sfi + r*Cfi)*Cth/Sth
        else:
            self.DTch = 0.0
        
    def colisiones(self):
        """ 
        Esta funcion se tiene que llamar en caso de que se estrelle el 
        helicoptero, o se presente un valor absurdo, ya que coloca
        al helicoptero en el suelo con valores razonables.
        """

        x = self.integrador.x1
        u_i, v_i, w_i = dot(self.LIB, [x[0], x[1], x[2]])
        if self.alt - self.fuselaje.ht - self.hG < 0 and w_i<0.0:
            # No penetrabilidad y no deslizamiento
            x[0] = x[1] = x[2] = 0.0
            x[3] = x[4] = x[5] = 0.0
            
            # Colocamos el helicoptero en posicion horizontal
            E = Quaternion(x[6], x[7], x[8], x[9]).toEuler()
            x[6], x[7], x[8], x[9] = Quaternion.euler(E[0], 0.0, pi)
            
            # Y corregimos el error de penetracion
            self.alt = self.hG + self.fuselaje.ht
            x[15] = self.geoide.geod2loc(
                    self.lat,
                    self.lon,
                    self.alt )[2] - 0.01

            # reseteamos el integrador
            self.integrador.init(F = self.paso,
                    ti = self.t,
                    xi = x)
    
    def calcula_FM(self):
        """
        Calcula las fuerzas y momentos totales que actuan sobre el C.M del 
        helicoptero sumando todas las fuerzas y momentos de los subsistemas.
        """

        # A lo mejor queremos controles fijos, por ejemplo
        # los calculados por el trimado
        if not self.trim:
            self.controles.calcula_controles()
        
        # Calculamos las fuerzas de los diferentes subsistemas
        self.rotor.calcula_FMR()
        self.rotor_cola.calcula_FMT()
        self.fuselaje.calcula_FMf()
        self.estabilizador_horizontal.calcula_FMtp()
        self.estabilizador_vertical.calcula_FMfn()
        
        FMR = self.rotor.FMR
        if self.motor.nmotores <= 0:
            FMR[0][5] = 0.0
            FMR[1][5] = zeros(shape=(16, ), type='Float64')
        
        FM0 = ( FMR[0] +
                self.rotor_cola.FMT[0] +
                self.fuselaje.FMf[0] +
                self.estabilizador_vertical.FMfn[0] +
                self.estabilizador_horizontal.FMtp[0]
        )
        
        FM1 = ( FMR[1] +
                self.rotor_cola.FMT[1] +
                self.fuselaje.FMf[1] +
                self.estabilizador_vertical.FMfn[1] +
                self.estabilizador_horizontal.FMtp[1]
        )

        self.F0 = FM0[0:3]
        self.M0 = FM0[3:6]
        self.F1 = FM1[0:3]
        self.M1 = FM1[3:6]

    def calcula_geo(self):
        self.lat, self.lon, self.alt = \
                self.geoide.loc2geod(
                        self.x_i,
                        self.y_i,
                        self.z_i)

    def calcula_locales(self):
        Slat, Clat = sin(self.lat), cos(self.lat)
        Slon, Clon = sin(self.lon), cos(self.lon)
        # matriz ejes locales-cartesianos geocentricos
        LLC = array( [
            [  Slat*Clon,  Slat*Slon, -Clat ],
            [ -Slon,       Clon,       0    ],
            [  Clat*Clon,  Clat*Slon,  Slat ]] ,type='Float64')

        self.LLC = LLC
        # matriz ejes cartesianos geocentricos-locales
        self.LCL = transpose(LLC)
        
        # matriz ejes locales-ejes inerciales
        LLI = dot(LLC, self.geoide.Lcl)

        # origen en geocentricas de los ejes locales
        self.Ol = self.geoide.geod2cart(self.lat, self.lon, 0.0)
        
        r_c = self.geoide.loc2cart(self.x_i, self.y_i, self.z_i)
        # Coordenadas locales
        ( self.x_l, 
                self.y_l, 
                self.z_l ) = dot(LLC, r_c - self.Ol)

        # Cuaternio en ejes locales
        self.p0 = self.q0
        ( self.p1, 
                self.p2,
                self.p3 ) = dot(LLI, [self.q1, self.q2, self.q3])

        self.LLI = LLI
        self.LIL = transpose(LLI)
        
    def calcula_euler(self):
        Q = Quaternion(self.p0, self.p1, self.p2, self.p3 )
        E = Q.toEuler()
        self.ch =  pi - E[0]
        self.th = -E[1]
        self.fi =  pi + E[2]
        
        # eps y 2pi + eps son equivalentes, pero no para
        # el sistema de control, asi que metemos los valores
        # en el rango -pi, pi
        if self.th > pi:
            self.th -= 2*pi
        elif self.th<-pi:
            self.th += 2*pi
        
        if self.fi > pi:
            self.fi -= 2*pi
        elif self.fi<-pi:
            self.fi += 2*pi

    def calcula_preparativos(self):
        """
        Suponiendo conocidas las magnitudes del vector de estado
        calcula el resto de magnitudes que faltan para ya poder pasar
        al calculo de las fuerzas y moentos y las ecuaciones
        del motor y de solido rigido.
        """

        # Posicion geodesica
        self.calcula_geo()
       
        # Ejes locales
        self.calcula_locales()
        
        # Angulos de Euler
        self.calcula_euler() 
        
        # Orientacion de los ejes cuerpo respecto a los inerciales
        self.calcula_LBI()
        
        # Densidad
        self.calcula_ro()
        
        # Velocidad del viento en ejes inerciales
        self.calcula_Vw_i()

        # Velocidad relativa al viento en ejes cuerpo
        self.calcula_Vrw_b()
        
        # Varias velocidades
        self.calcula_velocidades()

        # Gravedad
        self.calcula_G()
        
        # renormalizamos el cuaternio de rotacion
        modq = sqrt(self.q0**2 + self.q1**2 + self.q2**2 + self.q3**2)
        self.q0 /=modq
        self.q1 /=modq
        self.q2 /=modq
        self.q3 /=modq

    def calcula_ecs_modelo(self):
        """
        Calcula las variables ecs_0 y ecs_1 tal que la ecuacion a integrar
        para el vector de estado es:
            ecs_1*(dx/dt) = ecs_0

        En el proceso calcula tambien multitud de variables intermedias
        """

        if not self.tren.colision:
            self.calcula_preparativos()        
            # Calculamos fuerzas y momentos de subsistemas
            self.calcula_FM()
        
        # Calculamos ecuaciones del motor
        self.motor.calcula_ecs_motor()
        
        # alias para facilitar la legibilidad del codigo
        Ixx = self.Ixx
        Iyy = self.Iyy
        Izz = self.Izz
        Ixz = self.Ixz
        
        M = self.M
        
        u = self.u
        v = self.v
        w = self.w
        p = self.p
        q = self.q
        r = self.r

        q0 = self.q0
        q1 = self.q1 
        q2 = self.q2
        q3 = self.q3
            
        ####################################################################### 
        #   TERMINOS DEPENDIENTES LINEALES                                    #
        ####################################################################### 
        ecs_1 = zeros(shape=(16, 16), type='Float64')
        
        # Ecuaciones de velocidad lineal en ejes cuerpo
        ecs_1[0:3, 0:3] = [
                [ 1.0, 0.0, 0.0 ],
                [ 0.0, 1.0, 0.0 ],
                [ 0.0, 0.0, 1.0 ]]
        # Ecuaciones de velocidad angular en ejes cuerpo
        ecs_1[3:6, 3:6] = [
                [  Ixx, 0.0, -Ixz ],
                [  0.0, Iyy,  0.0 ],
                [ -Ixz, 0.0,  Izz ]]
        # Relaciones cinematicas del cuaternio cambio ejes cuerpo-inerciales
        ecs_1[6:10, 6:10] = [
                [ 1.0, 0.0, 0.0, 0.0 ],
                [ 0.0, 1.0, 0.0, 0.0 ],
                [ 0.0, 0.0, 1.0, 0.0 ],
                [ 0.0, 0.0, 0.0, 1.0 ]]
        
        # Ecuaciones de la posicion del avion
        ecs_1[13:16, 13:16] = [
                [ 1.0, 0.0, 0.0 ],
                [ 0.0, 1.0, 0.0 ],
                [ 0.0, 0.0, 1.0 ]]

        # Añadimos los terminos dependientes a las ecuaciones
        ecs_1[0:3, :]   -= self.F1/M
        ecs_1[3:6, :]   -= self.M1
        ecs_1[10:13, :] += self.motor.ecs_motor[1]
        
        ####################################################################### 
        #   TERMINOS INDEPENDIENTES  NO LINEALES                              #
        ####################################################################### 
        ecs_0 = zeros( shape = (16,), type = 'Float64')
        
        # Aceleraciones y momentos de inercia
        ecs_0[0] = -( w*q - v*r ) 
        ecs_0[1] = -( u*r - w*p ) 
        ecs_0[2] = -( v*p - u*q ) 
        ecs_0[3] =  ( Iyy - Izz )*q*r + Ixz*p*q
        ecs_0[4] =  ( Izz - Ixx )*r*p + Ixz*( r**2 - p**2 )
        ecs_0[5] =  ( Ixx - Iyy )*p*q + Ixz*q*r 
        
        # Gravedad y fuerzas de rotor, cola, etc...
        ecs_0[0:3] += self.G + self.F0/M
        ecs_0[3:6] += self.M0
        
        # Lado derecho de las ecuaciones del quaternio
        # Mucho ojo que hay que pasar la velocidad a angular a ejes inerciales
        p_i, q_i, r_i = dot(self.LIB, [p, q, r])
        
        ecs_0[6] = -0.5*(  q1*p_i + q2*q_i + q3*r_i ) 
        ecs_0[7] = -0.5*( -q0*p_i - q3*q_i + q2*r_i ) 
        ecs_0[8] = -0.5*(  q3*p_i - q0*q_i - q1*r_i ) 
        ecs_0[9] = -0.5*( -q2*p_i + q1*q_i - q0*r_i ) 
        
        # Lado derecho de las ecuaciones del motor
        ecs_0[10:13] = self.motor.ecs_motor[0]
        
        # Lado derecho de las ecuaciones cinematicas
        ecs_0[13:16] = dot(self.LIB, [u,v,w])
        
        #######################################################################
        #           CONTACTO CON EL SUELO                                     #
        #######################################################################
        # ESTO NECESITA MUCHISIMA MEJORA
        # si tocamos el suelo, y vamos hacia el
#        if (self.alt - self.fuselaje.ht - self.hG) <= 0.0 and w>=0:
#            # No penetrabilidad
#            if ecs_0[2]>0:
#                ecs_0[2] = 0.0                      # Fuerzas verticales
#            # No deslizamiento
#            ecs_0[0] = ecs_0[1] = 0.0               # Fuerzas tangenciales
#            ecs_0[3] = ecs_0[4] = ecs_0[5] = 0.0    # Momentos
#
#            ecs_1[0:6] = zeros(shape=(6, len(self.x)), type='Float64')
#            ecs_1[0:6, 0:6] = identity(6)
        
        #------------------------------------------------------------------#    
        self.ecs_modelo = ecs_0, ecs_1
        #------------------------------------------------------------------#    
        
        # logging
        ( self.ecs_u, 
                self.ecs_v, 
                self.ecs_w, 
                self.ecs_p, 
                self.ecs_q, 
                self.ecs_r ) = ecs_0[:6]

        ( self.ecs_q0, 
                self.ecs_q1, 
                self.ecs_q2, 
                self.ecs_q3 ) = ecs_0[6:10]

        self.ecs_Om, self.ecs_DTQ1, self.ecs_Q1 = ecs_0[10:13]
        
    
    
    ###########################################################################
    #           ACCESO MEDIANTE SIMBOLOS                                      #
    ###########################################################################
    def __path(self, s):
        return self.__logstr[s]

    def get(self, s):
        p = self.__path(s)
        try:
            if p[0] == 'modelo':
                x = getattr(self, p[1])
            else:
                x = getattr(getattr(self, p[0]), p[1])
        except AttributeError:
            x = 0.0
        return x
    
    def set(self, s, v):
        p = self.__path(s)
        if p[0] == 'modelo':
            setattr(self, p[1], v)
        else:
            setattr(getattr(self, p[0]), p[1], v)

    ###########################################################################
    #           SERIES TEMPORALES                                             #
    ###########################################################################
    def series(self, T, et0p, et1sp, et1cp, etpp, dt=0.03, *args):
        """
        La siguiente funcion calcula la respuesta de las magnitudes pedidas
        en args (simbolos) desde el instante actual (self.t) hasta el instante
        T, dando un paso de tiempo dt. et0p, et1sp, et1cp, etpp son funciones
        del tiempo o valores constantes. Devuelve los resultados en un 
        diccionario que contiene por lo menos las llaves 't' con los pasos 
        de tiempo dados, 'et0p', 'et1sp', 'et1cp' y 'etpp' con los valores
        en cada instante de tiempo y ademas una llave por cada simbolo que se 
        pidiese calcular.

        Un ejemplo de una posible sesion en python seria:
        
        In [1]: from Lynx import modelo
        In [2]: modelo.init() 
        In [3]: tr = modelo.trimado(0, 0, 0, 0.0, 100., 1.225)
                  ...    (Informacion soltada por el proceso de trimado)
                  ...
                  ...
        In [4]: modelo.t = 0.0
        In [5]: s = modelo.series(0.2, 0, 0, 0, 0, 0.03, 'zi')
        In [6]: s['t']
        Out[6]: array([ 0.03,  0.06,  0.09,  0.12,  0.15,  0.18,  0.21])
        In [7]: s['zi']
        Out[7]:
        [100.0,
         100.0,
         99.986765709951754,
         99.964728999868001,
         99.933898156358978,
         99.894286155391967,
         99.845908960219973]
        """
        
        # Hay que calcular los controles
        self.trim = False

        def reloj_tonto():
            return self.t
        
        # Ya tenemos el instante actual
        # El instante final es T, incluido
        DT = arange(self.t + dt, T + dt, dt)
        
        # Diccionario con los resultados
        r = {
            'et0p':  [],
            'et1sp': [],
            'et1cp': [],
            'etpp':  []
        }
        for k in args:
            r[k] = [] 
        
        # Si alguno de los controles es una constante lo 
        # transformamos en una funcion constante
        def func(x):
            if not hasattr(x, '__call__'):
                return lambda y: x
            else:
                return x
        
        f_et0p  = func(et0p )
        f_etpp  = func(etpp )
        f_et1sp = func(et1sp)
        f_et1cp = func(et1cp)

        # Preparamos el integrador
        self.controles.et0p.append (f_et0p  (self.t))
        self.controles.et1sp.append(f_et1sp (self.t))
        self.controles.et1cp.append(f_et1cp (self.t))
        self.controles.etpp.append (f_etpp  (self.t))

        self.integrador.init(
                    F  = self.paso,
                    ti = self.t,
                    xi = array([
                        self.u, 
                        self.v, 
                        self.w, 
                        self.p, 
                        self.q, 
                        self.r, 
                        self.q0, 
                        self.q1, 
                        self.q2, 
                        self.q3, 
                        self.Om, 
                        self.DTQ1, 
                        self.Q1, 
                        self.x_i, 
                        self.y_i, 
                        self.z_i]
                    )
        )
        
        
        for t in DT:
            # Controles en el instante actual
            self.controles.et0p.append (f_et0p  (t), t)
            self.controles.et1sp.append(f_et1sp (t), t)
            self.controles.et1cp.append(f_et1cp (t), t)
            self.controles.etpp.append (f_etpp  (t), t)
            
            # Avanzamos al siguiente instante
            self.avanza(t)
            self.t = t
            
            # Guardamos resultados
            for k in r.keys():
                r[k].append(self.get(k))
        r['t'] = DT
        return r

    def trimado(self, V_t, ga_t, be_t, Oma_t, z_t, ro_t, metodo='Fijo'):
        """
        Calcula el trimado del helicoptero para los parametros:
            V_t:   Velocidad del helicoptero.
            ga_t:  Angulo de subida
            be_t:  Angulo de resbalamiento.
            Oma_t: Velocidad de giro
            z_t:   Altura 
            ro_t:  Densidad del aire
            metodo: 'Fijo' o 'Newton', modifica el método numérico. 
                     Se recomienda MUCHO utilizar Fijo porque Newton no 
                     funciona bien todavia.
        """

        # motores encendidos!!!
        self.motor.nmotores = 2
        self.motor.start = False
        ######################################################################
        #   ALIAS                                                            #
        ######################################################################
        # Para no tener self-verborrea se acortan los siguientes nombres
        M   = self.M
        Ixx = self.Ixx
        Iyy = self.Iyy
        Izz = self.Izz
        Ixz = self.Ixz

        tht = self.rotor.tht
        a0 =  self.rotor.a0
        s =   self.rotor.s
        
        # Para calcular los mandos del piloto
        c0 = self.controles.c0
        c1 = self.controles.c1
        S0 = self.controles.S0
        S1 = self.controles.S1
        S2 = self.controles.S2
        
        # Constantes de la TCMM
        ki =  (9./5.)**0.25
        knu = (5./4.)**0.25
        g =  9.81
        
        # Precalculamos senos y cosenos que permanecen constantes
        Cbe, Sbe = cos(be_t), sin(be_t)
        Cga, Sga = cos(ga_t), sin(ga_t)
        CK = cos(self.rotor_cola.K)
        SK = sin(self.rotor_cola.K)


        
        ########################################################################
        #   ESTIMACION VALORES INICIALES                                       #
        ########################################################################
        # Coeficientes del fuselaje, practicamente no varian
        cDf, cLf, cYf, clf, cmf, cnf = self.fuselaje.aero.coefs(0.0, 0.0)
        
        # Angulos de euler:
        # th compensa la resistencia del fuselaje
        # fi compensa la aceleracion centrifuga
        th = -0.5*ro_t*V_t**2*self.fuselaje.Sp*cDf/(M*g)
        fi =  atan(Oma_t*V_t/9.81)
        
        # Hay que iniciar una serie de variables, aunque sea a cero
        # Fuerzas y momentos del rotor de cola
        XT = YT = ZT = 0.0
        LT = MT = NT = 0.0
        QT = 0.0
        
        # Fuerzas y momentos del rotor
        XR = YR = LR = MR = NR = 0.0
        LRh_b = NRh_b = YRh_b = 0.0
        LRh_w = NRh_w = QRh_w= YRh_w = 0.0
        
        # Velocidad inducida
        la1c_w = la1s_w = 0.0
        la0 = la0NG = 0.006
        
        # Batimiento del rotor
        be0 = be1c_w = be1s_w = 0.0
        # Batimiento del rotor de cola
        be1cT = be1sT = be1sT_w = 0.0
        
        # Velocidad de giro del rotor
        Om = self.motor.Omi
        # Coeficiente de traccion
        cT = M*g/(ro_t*(Om*self.rotor.R)**2*pi*self.rotor.R**2)
        
        # Pasos
        th0 = th1s = th1c = th0T = 0.0
        th1c_w = th1s_w = 0.0
        # Controles que necesita hacer el piloto
        et0p = et1cp = et1sp = etpp = 0.0   
        
        ########################################################################
        #   FUNCIONES AUXILIARES                                               #
        ########################################################################
        # Esta funcion, que se llama en dos puntos del programa calcula el 
        # vector be0, la1s_w, la1c_w, th0, th1s_w, th1c_w
        def f_controles():
            # Velocidades en el rotor
            C = (nu/knu)**2 + (1/knu**2 - 1/ki**2)*nuz**2 +\
                    ((nuz - la0NG)/ki)**2
            VT = sqrt(C)
            barV = VT*(1 - la0NG*(nuz - la0NG)/ki**2/C)
    
            # Constantes del batimiento
            labe2 = 1. + self.rotor.Kbe/( self.rotor.Ibe*Om**2)
            ga = ro_t*self.rotor.c*a0*self.rotor.R**4/self.rotor.Ibe
        
            #
            # Batimiento y controles del rotor principal    
            #

            A_bebe = array([
                [ 8*labe2/ga, 0,              0                 ],
                [ 4./3.*nu,   8*(labe2-1)/ga, 1 + 0.5*nu**2     ],
                [ 0,         -1. + 0.5*nu**2, 8*(labe2 - 1)/ga  ]], 
                type='Float64')
            
            A_bela = array([
                [ 2./3.*nu, 0. ],
                [ 0.,       1. ],
                [ 1.,       0. ]], type='Float64')
            
            # Lado derecho
            # Influencia del paso sobre el batimiento
            A_beth = array([
                [ 
                    1. + nu**2, 
                    4./5. + 2./3.*nu**2, 
                    4./3.*nu,      
                    0.             
                ],
                [ 
                    0.,         
                    0.,                  
                    0.,             
                    1. + 0.5*nu**2 
                ],
                [ 
                    8./3.*nu,   
                    2*nu,                
                    1. + 1.5*nu**2, 
                    0.              
                ]], type='Float64')
            
            # Influencia de la aceleracion angular sobre el batimiento
            A_beDChiom = 8./ga*array([
                [ 0., 0. ],
                [ 0., 1. ],
                [ 1., 0. ]], type='Float64')
    
            # Influencia de la velocidad angular sobre el batimiento
            A_beom = array([
                [ 2./3.*nu,  0       ],
                [ 16./ga,    1.      ],
                [ 1.,       -16./ga  ]], type='Float64')
            
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
                    [ 15.*pi/(64.*VT)*X, 0.,                 2.*(1. - X**2)/barV ]]
                    , type='Float64')
            
            # Matriz de fuerzas y momentos aerodinamicos
            # N = { F0/2, M1s, M1c }
            Nth = array([
                [ 
                    2./3. + nu**2, 
                    0.5*(1. + nu**2), 
                    nu,                    
                    0.                    
                ],
                [ 
                    2./3.*nu,      
                    0.5*nu,           
                    0.25*(1. + 1.5*nu**2), 
                    0                     
                ],
                [ 
                    0.,            
                    0.,               
                    0.,                    
                    0.25*(1. + 0.5*nu**2) 
                ]], type='Float64')
    
            Nbe = array([
                [ 0.,   0.,                    0.                        ],
                [ 0.,   0.25*(1. - 0.5*nu**2), 0.                        ],
                [-nu/3.,                       0, -0.25*(1. + 0.5*nu**2) ]], 
                type='Float64')
    
            Nla = -0.5*array([
                [ nu,  0.  ],
                [ 0.5, 0.  ],
                [ 0.,  0.5 ]], type='Float64')
    
            Nom = -Nla
    
            Nla0 = array([
                1., 
                0.5*nu, 
                0.], type='Float64')
            
            # Lado izquierdo de las ecuaciones de la velocidad inducida
            m = a0*s/4.
            A_labe = -m*dot(K, Nbe)
            A_lala =  identity(2) - m*dot(K, Nla)
            # Matriz de influencia de controles, velocidad angular y velocidad
            # inducida uniforme sobre la velocidad inducida
            A_lath =  m*dot(K, Nth)
            A_laom =  m*dot(K, Nom)
            A_lala0 = m*dot(K, Nla0)
            # Matriz de influencia de la aceleracion angular sobre la velocidad
            # inducida
            A_laDChiom = zeros(shape=(2,2), type='Float64')
    
            lhs = array(shape=(6,6), type='Float64')
            lhs[0:3, 0  ] =  A_bebe[:, 0]
            lhs[0:3, 1:3] =  A_bela
            lhs[0:3, 3  ] = -A_beth[:, 0] 
            lhs[0:3, 4:6] = -A_beth[:, 2:4]
            lhs[3:5, 0  ] =  A_labe[:, 0]
            lhs[3:5, 1:3] =  A_lala
            lhs[3:5, 3  ] = -A_lath[:, 0] 
            lhs[3:5, 4:6] = -A_lath[:, 2:4]
            lhs[5,   0:6] =  array([0, 0, 0, 1./3. + nu**2/2., nu/2., 0], 
                    type='Float64')

            rhs = array(shape=(6), type='Float64')
            rhs[0:3] = (
                    -dot(A_bebe[:, 1:3], [be1c_w, be1s_w]) +
                        A_beth[:, 1]*self.rotor.tht +
                        dot(A_beom, [p_w, q_w]) +
                        (nuz-la0)*A_bela0
            )
            rhs[3:5] = (
                    -dot(A_labe[:, 1:3], [be1c_w, be1s_w]) +
                    A_lath[:, 1]*self.rotor.tht +
                    dot(A_laom, [p_w, q_w]) +
                    (nuz-la0)*A_lala0
            )
            rhs[5] = (
                    2*cT/(a0*s) -
                    nu/4.*p_w -
                    0.5*(nuz - la0) -
                    0.25*(1. + nu**2)*self.rotor.tht
            )

            # be0, la1s_w, la1c_w, th0, th1s_w, th1c_w
            return la.solve_linear_equations(lhs, rhs)
                    
        # Esta funcion, que se llama en dos puntos del programa calcula 
        # controles a partir de los pasos
        def f_etp(th0, th1s, th1c, th0T, p, q, r, th, fi, ch):
            et =  dot(la.inverse(c1), 
            array([th0, th1s, th1c, th0T], type='Float64') - c0)
            print et
            # Controles del piloto a partir de los controles
            def f(x):
                th_0 = self.controles.th_0(x[1])
                fi_0 = self.controles.fi_0(x[2])
                ch_0 = self.controles.ch_0(x[3])
                return dot(S0, x) + dot(S1, [p, q, r]) +\
                        dot(S2, [th - th_0, fi - fi_0, ch - ch_0]) - et

            etp = solveNewtonMulti(f, [0., 0., 0., 0.], 100, 1e-5)

            return etp  
        
        #######################################################################
        #   CONSTANTES NUMÉRICAS                                              #
        #######################################################################
        
        # Errores numericos: en ciertos puntos del programa si 
        # abs(x)<eps0-->x = 0
        eps0 = 1e-3
        
        # Precision deseada en Om, fi y th
        eps_Om = 1e-3
        eps_fi = 1e-3
        eps_th = 1e-3
        
        # Numero maximo de iteraciones en los respectivos bucles
        Nmax_Om = 50
        Nmax_fi = 50
        Nmax_th = 200
        
        # Numero minimo de iteraciones en los respectivos bucles
        Nmin_Om = 3
        Nmin_fi = 3
        Nmin_th = 3

        # Amortiguadores para asegurar convergencia, valores
        # experimentales
        k_th = 1.0
        k_fi = 0.5
        k_Om = 1.0
        
        # Bandera para controlar el caso ga=pi/2, que numericamente
        # es problematico, si se detecta  entonces
        # se activa y se tiene en cuenta
        resb = False
            
        
        #######################################################################
        #   COMIENZO DE LA ITERACION                                          #
        #######################################################################

        # COMIENZA BUCLE OM #
        n_Om = 0        # Numero de iteraciones realizadas
        err_Om = 1      # error en la ultime iteracion
        while err_Om>eps_Om or n_Om<Nmin_Om:
            if n_Om>Nmax_Om:
                print "No converge Om, amortiguando el proceso!!"
                k_Om -=0.1
                # reseteamos el bucle
                n_Om = 0
                Om = self.motor.Omi
                if k_Om <0.2:
                    print ( "Muy amortiguado ya, no se puede alcanzar " + 
                            "convergencia" )
                    print "Terminando trimado"
                    break
            
            # COMIENZA BUCLE FI #
            n_fi = 0
            err_fi = 1
            while err_fi>eps_fi or n_fi<Nmin_fi:
                if n_fi>Nmax_fi:
                    print "No converge fi, amortiguando el proceso!!"
                    k_fi -=0.1
                    # reseteamos el bucle
                    n_fi = 0
                    fi = atan(Oma_t*V_t/9.81)
                    if k_fi <0.2:
                        print ( "Muy amortiguado ya, no se puede alcanzar " +
                                "convergencia" )
                        print "Saltando trimado de fi"
                        k_th = 1.0
                        break
                
                # No cambia en todo el bucle TH
                Cfi, Sfi = cos(fi), sin(fi)
    
                # COMIENZA BUCLE TH #
                st = 0
                n_th = 0
                err_th = 1
                while err_th>eps_th or n_th<Nmin_th:
                    if n_th>Nmax_th:
                        print "No converge th, amortiguando el proceso!!"
                        k_th -= 0.2
                        # reseteamos el bucle
                        n_th = 0
                        th = -0.5*ro_t*V_t**2*self.fuselaje.Sp*cDf/(M*g)
                        if k_th <0.2:
                            print ( "Muy amortiguado ya, no se puede alcanzar "
                                    + "convergencia" )
                            print "Saltando trimado de theta"
                            k_th = 1.0
                            break
                    Cth, Sth = cos(th), sin(th)                 
                    ############################################################
                    #   VELOCIDAD Y VELOCIDAD ANGULAR                          #
                    ############################################################
                    # Calculamos la velocidad angular
                    p = -Oma_t*Sth
                    q =  Oma_t*Sfi*Cth
                    r =  Oma_t*Cfi*Cth
                    # Calculamos la velocidad
                    # Resolvemos la ecuacion de segundo orden para el seno 
                    # del angulo de tracking: a = b*cos(xi) + c*sin(xi)
                    a = Sbe - Sga*Sfi*Cth
                    b = Cga*Sth*Sfi
                    c = Cga*Cfi
                    # A*sin(xi)^2 + B*sin(xi) + C = 0
                    A =  c**2 + b**2
                    B = -2*a*c
                    C =  a**2 - b**2
                    # sin(xi) = (-B +- D)/(2A)
                    D =  B**2 - 4*A*C
                    # Manejamos errores numericos que nos echan del dominio 
                    if abs(A)<eps0:
                        # Suponemos que ha pasado ga = pi/2, en cuyo caso 
                        # be no esta definido, en cualquier caso el valor 
                        # de xi no afecta el valor de la velocidad
                        # Obligamos a=0 ajustando be
                        Sbe = Sga*Sfi*Cth
                        Sxi = 0.
                        Cxi = 1.
                        # Solo avisamos la primera vez
                        if resb == False:
                            print "¿ga=pi/2? Dejando libre  al resbalamiento"
                            resb = True
                    else:
                        # Corregimos error numerico, -1e-16 daria error
                        if D<0:
                            if -D<eps0:
                                D = 0.0
                            else:
                                # La ecuacion no tiene solucion obligatoriamente
                                sys.exit("Error calculando el ángulo de " +
                                         "tracking, raíz imaginaria")
                        
                        # Las dos raíces
                        Sxi1 = (-B - sqrt(D))/(2.*A)
                        Sxi2 = (-B + sqrt(D))/(2.*A)
                        # Descartamos una, nos quedamos con la mas pequeña
                        if abs(Sxi1) < abs(Sxi2):
                            Sxi = Sxi1
                        else:
                            Sxi = Sxi2
                        if Sxi>1:
                            if (Sxi-1)<eps0:
                                Sxi = 1.0
                            else:
                                sys.exit("Error calculando el ángulo de " +
                                         "tracking, sin(xi)>1")
                        Cxi2 = 1. - Sxi**2
                        # Si es menor que 0 tiene que ser error numerico seguro
                        if Cxi2 < 0:
                            Cxi2 = 0.0
                        Cxi = sqrt(Cxi2)
                    
                    # Velocidad en ejes cuerpo
                    u = V_t*(Cth*Cga*Cxi - Sth*Sga)
                    v = V_t*Sbe
                    w = V_t*(Cga*Sth*Cfi*Cxi + Cga*Sfi*Sxi + Sga*Cfi*Cth)
    
                    # Para aprovechar funciones ya disponibles
                    # Cuaternio de rotacion de ejes inerciales a cuerpo
                    _Q = Quaternion.rotacion([
                        (pi, 0, 1, 0), # Ejes inerciales a auxiliares
                        (th, 0, 1, 0), # Angulos de euler, ch = 0
                        (fi, 1, 0, 0)])
                    
                    # Solo estas variables influyen en el calculo de las 
                    # fuerzas aerodinamicas
                    self.u = u
                    self.v = v
                    self.w = w

                    self.p = p
                    self.q = q
                    self.r = r

                    self.q0 = _Q.q0
                    self.q1 = _Q.q1
                    self.q2 = _Q.q2
                    self.q3 = _Q.q3

                    self.Om = Om

                    self.rotor.la0 = la0
                    
                    ############################################################
                    #   CALCULO DE EJES VIENTO                                 #
                    ############################################################
                    # Calculamos ejes viento del rotor
                    self.calcula_LBI()
                    self.Vw_i = array([0.0, 0.0, 0.0], type='Float64')
                    self.ro = ro_t
                    self.calcula_Vrw_b()
                    
                    # Velocidad angular en ejes cuerpo
                    W_b = array([p, q, r], type='Float64')
                    
                    # Velocidad del rotor relativa al viento, en ejes cuerpo
                    Vhrw_b = self.Vrw_b +\
                            cross(W_b, [-self.rotor.xcg, 0.0, -self.rotor.hR])                  
                    # en ejes rotor
                    uhrw_h, vhrw_h, whrw_h = dot(self.rotor.LHB, Vhrw_b)
                    
                    # Velocidad angular en ejes rotor
                    W_h = dot(self.rotor.LHB, W_b)
                    
                    # Velocidad relativa del rotor respecto al viento, en ejes 
                    # viento
                    uhrw_w = sqrt(uhrw_h**2 + vhrw_h**2)
                    whrw_w = whrw_h
                    
                    # Angulo que forman los ejes viento con los ejes rotor. 
                    # Si no hay velocidad se supone 0
                    if uhrw_w!=0:
                        Cchw = uhrw_h/uhrw_w
                        Schw = vhrw_h/uhrw_w
                    else:
                        if vhrw_h==0:
                            Cchw = 1.0
                            Schw = 0.0
                        elif vhrw_h<0:
                            Cchw =  0.0
                            Schw = -1.0
                        else:
                            Cchw = 0.0
                            Schw = 1.0
                    
                    # Matriz cambio de base de ejes viento a ejes cuerpo
                    # nos es util para pasar las fuerzas y momentos a ejes
                    # cuerpo
                    LHW = array([
                        [Cchw, -Schw, 0.0],
                        [Schw,  Cchw, 0.0],
                        [0.0,   0.0,  1.0]], type = 'Float64')
                    LBW = dot(self.rotor.LBH, LHW)
                    LWB = transpose(LBW)
                
                    # Parametros de avance 
                    nu =  uhrw_w/( Om*self.rotor.R )
                    nuz = whrw_w/( Om*self.rotor.R )
                
                    # Velocidad angular adimensionalizada en ejes viento
                    p_w = (  Cchw*W_h[0] + Schw*W_h[1] )/Om
                    q_w = ( -Schw*W_h[0] + Cchw*W_h[1] )/Om
                    
                    ############################################################
                    #   FUERZAS DE FUSELAJE Y ESTABILIZADORES                  #
                    ############################################################
                    # Vrw_b, [p,q,r] y [q0,q1,q2,q3] ya estan calculados, 
                    # que es lo que cuenta
                    # Angulo de la estela, para el estabilizador horizontal
                    xi = self.rotor.xi = ang(
                            abs(la0*Om*self.rotor.R - whrw_w), uhrw_w)                     

                    # Calculamos las fuerzas aerodinamicas y guardamos los 
                    # resultados, en ejes cuerpo. Aprovechamos funciones ya 
                    # declaradas.
                    
                    # Fuselaje
                    self.fuselaje.calcula_FMf()
                    Xf, Yf, Zf, Lf, Mf, Nf = self.fuselaje.FMf[0]
                    
                    self.estabilizador_horizontal.calcula_FMtp()
                    Xtp, Ytp, Ztp, Ltp, Mtp, Ntp = \
                            self.estabilizador_horizontal.FMtp[0]
                    
                    # Estabilizador vertical
                    self.estabilizador_vertical.calcula_FMfn()
                    Xfn, Yfn, Zfn, Lfn, Mfn, Nfn = \
                            self.estabilizador_vertical.FMfn[0]
                
                
                    ############################################################
                    #   RESOLUCION ECUACIONES LONGITUDINALES                   #
                    ############################################################
                    # Despejamos las fuerzas vertical del rotor
                    ZR = M*(-q*u + p*v) - (Zf + Zfn + Ztp + ZT) - M*g*Cth*Cfi
                    
                    # Despejamos el momento del rotor
                    MR = -(Izz - Ixx)*r*p - Ixz*(r**2 - p**2) - \
                            (Mf + Mfn + Mtp + MT)
                    # Trasladmos las fuerzas y momentos del rotor en ejes cuerpo
                    # en el centro de masas a fuerzas y momentos en ejes viento
                    # en el rotor
                    XRh_b = XR
                    ZRh_b = ZR
                    MRh_b = MR - ZR*self.rotor.xcg + XR*self.rotor.hR
                    
                    # FIX: Si tenemos angulo de resbalamiento habria que 
                    # investigar como de necesario es realizar un buen calculo 
                    # de YRh_b, LRh_b y NRh_b antes de despejar valores en ejes
                    # viento
                    ZRh_w = dot(LWB[2], [XRh_b, YRh_b, ZRh_b])
                    MRh_w = dot(LWB[1], [LRh_b, MRh_b, NRh_b])
                    
                    # Dimension de fuerzas para el rotor
                    dimFR = ro_t*(Om*self.rotor.R)**2*pi*self.rotor.R**2
                    
                    # Calculamos cT
                    cT = -ZRh_w/dimFR
                    
                    # Coeficiente de rozamiento de los perfiles
                    de = self.rotor.de0 + self.rotor.de2*cT**2
        
                    # Calculamos el batimiento be1c_w usando la ecuacion de 
                    # momento de cabeceo
                    # FIX: Como de bueno necesitamos a Qrh_w?
                    be1c_w = -(MRh_w - QRh_w/2.*be1s_w)*2/(
                            self.rotor.Nb*self.rotor.Kbe)
        
                    # Estimamos un nuevo valor para XR, aprovechando el nuevo 
                    # valor de be1c_w
                    # El problema es que necesitamos los valores de los 
                    # controles para poder realizar una buena estimación
                    be0, la1s_w, la1c_w, th0, th1s_w, th1c_w = f_controles()
                    
                    # Utilizamos funciones para calcular los terminos de 
                    # F1 y F2 porque mas adelante los calcularemos de nuevo
                
                    #Angulo de ataque efectivo en el perfil de la pala
                    al1s_w = p_w - la1s_w + be1c_w + th1s_w
                    al1c_w = q_w - la1c_w - be1s_w + th1c_w
        
                    # Primeros armonicos F1 
                    def f_F11s():
                        return al1s_w/3. + nu*( th0 + nuz - la0 + 2./3.*tht )
                    
                    def f_F11c():
                        return al1c_w/3. - 0.5*nu*be0
        
                    # Segundos armonicos F1
                    def f_F12s():
                        return  0.5*nu*( 0.5*al1c_w + 0.5*( th1c_w - be1s_w ) 
                                - nu*be0 )
                    def f_F12c():
                        return -0.5*nu*( 0.5*al1s_w + 0.5*( th1s_w + be1c_w ) 
                                + nu*( th0 + 0.5*tht ) )
        
                    # Primeros armonicos F2
                    def f_F21s():
                        return 0.5*(nu**2)*be0*be1s_w  + \
                                ( nuz - la0 - 0.25*nu*be1c_w ) *\
                                   ( al1s_w - th1s_w ) -\
                                   0.25*nu*be1s_w*( al1c_w - th1c_w ) +\
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
                    def f_F21c():
                        return ( al1c_w - th1c_w - 2*be0*nu )*( 
                                nuz - la0 - 0.75*nu*be1c_w )-\
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
                    # _Cxw = 2*Cxw/(a0*s)
                    # _Cyw = 2*Cyw/(a0*s)
                    # Comunes _Cxw y _Cyw
                    F12c = f_F12c()
                    F12s = f_F12s()
                    
                    _Cxw = ( cT/(a0*s) + F12c/4. )*be1c_w + f_F11c()/2.*be0 +\
                            F12s/4.*be1s_w + f_F21s()/2.
                    _Cyw = (-cT/(a0*s) + F12c/4. )*be1s_w - f_F11s()/2.*be0 -\
                            F12s/4.*be1c_w + f_F21c()/2.
                    
                    XRh_w = dimFR*a0*s/2.*_Cxw
                    YRh_w = dimFR*a0*s/2.*_Cyw
                    XR = XRh_b = dot(LBW[0], [XRh_w, YRh_w, ZRh_w])
                    
                    # Calculamos el angulo theta a partir de la ecuacion 
                    # horizontal y vertical
                    # Fuerzas de la gravedad
                    XG = M*(-r*v + q*w) - (Xf + Xfn + Xtp + XT + XR)
                    ZG = M*(-q*u + p*v) - (Zf + Zfn + Ztp + ZT + ZR)    
                    
                    # Iteramos th
                    if metodo == 'Newton':
                        f_th = atan(-Cfi*XG/ZG)
                        epsNew = 0.001
                        if st==0:
                            f_th0 = f_th
                            th = th - epsNew
                            st += 1
                        elif st==1:
                            f_th1 = f_th
                            th = th + epsNew
                            st += 1
                        elif st==2:
                            f_th2 = f_th
                            th_0 = th
                            th = th_0 - k_th*(th_0 - f_th0)/(
                                    1 - (f_th2 - f_th1)/(2*epsNew))
                            st = 0
                            err_th = abs(th - th_0)
                            n_th += 1
                    elif metodo == 'Fijo':
                        th_0 = th
                        th = th_0 + k_th*(atan(-Cfi*XG/ZG) - th_0)
                        err_th = abs(th - th_0)
                        n_th += 1
                    else:
                        sys.exit("Metodo de iteracion desconocido: %s" % metodo)
                    
                    print ( ( "\t\tth = %f, be1c_w = %f, cT = %f, altp0 = %f, "+
                           "err_th = %f") % 
                            (th, be1c_w, cT, 
                                self.estabilizador_horizontal.altp0, err_th))
                    
                # FIN BUCLE TH #
                Cth, Sth = cos(th), sin(th)
            
                ################################################################
                #   VELOCIDAD INDUCIDA                                         #
                ################################################################
                # Calculamos la velocidad inducida en el rotor
                # ki, knu = 1.0, 1.0
                def f(la0):
                    def VT(la0):
                        return (nu/knu)**2 + (1./knu**2 - 1./ki**2)*nuz**2  +\
                                ((nuz - la0)/ki)**2
                    # return la0/ki - cT/(2*sqrt(VT(la0)))
                    return 2*sqrt(VT(la0))*la0/ki - cT

                la0NG = solveSecante(f, 0.1, 0.5, 100, 1e-12)

                z = z_t + self.rotor.h0
                if z>self.rotor.R/10.:
                    K0G = 1 - 1./(16.*(z/self.rotor.R)**2)/(1 + (nu/la0NG)**2)
                else:
                    K0G = 1.0
                # Velocidad con efecto suelo
                la0 = K0G*la0NG + nuz*(1-K0G)
                
                ################################################################
                #   TRACCION DE COLA                                           #
                ################################################################
                # Calculamos el momento del rotor
                NRh_w = QRh_w = self.rotor.R*(-(nuz - la0)*(-ZRh_w) +\
                        nu*XRh_w + dimFR*de*s/8.*(1. + 3.*nu**2))
                NRh_b = dot(LBW[2], [LRh_w, MRh_w, NRh_w])
                
                # FIX: de nuevo hay que investigar la influencia de YRh_b, 
                # podria dar valores de error altos para helicoptero como el 
                # Lynx o el BlackHawk con xcg grande
                QR = NR = NRh_b - YRh_b*self.rotor.xcg 
                
                # Despejamos la traccion de cola a partir de la ecuación de 
                # momentos de guiñada
                # FIX: que influencia tiene QT en helicopteros con K!=0 como el
                # Black Hawk?
                NT = -(Ixx - Iyy)*p*q + Ixz*q*r - (Nf + Nfn + Ntp + NR)
                TT =  (-NT + SK*QT)/((self.rotor_cola.lT +\
                        self.rotor.xcg)*(CK - be1sT*SK))
                
                # Fuerzas y momentos en ejes cola
                XT_T =  TT*be1cT
                YT_T = -TT*be1sT
                ZT_T = -TT
                LT_T =  0.0
                MT_T =  0.0
                NT_T =  QT
                
                # Fuerzas y momentos en ejes cuerpo
                XT_b, YT_b, ZT_b = dot(self.rotor_cola.LBT, [XT_T, YT_T, ZT_T])
                LT_b, MT_b, NT_b = dot(self.rotor_cola.LBT, [LT_T, MT_T, NT_T])
                
                # Fuerzas y momentos en ejes cuerpo respecto al centro de masas
                l = self.rotor_cola.lT + self.rotor.xcg
                XT, YT, ZT = XT_b, YT_b, ZT_b
                LT = LT_b + YT_b*self.rotor_cola.hT
                MT = MT_b - XT_b*self.rotor_cola.hT + ZT_b*l
                NT = NT_b - YT_b*l
        
                # Despejamos el batimiento lateral a partir del momento de 
                # balance
                LR = (Ixx - Iyy)*q*r - Ixz*p*q - (Lf + Lfn + Ltp + LT)
                LRh_b = LR - YR*self.rotor.hR
                LRh_w = dot(LWB[0], [LRh_b, MRh_b, NRh_b])
                be1s_w = -(LRh_w + QRh_w/2.*be1c_w)*2/self.rotor.Nb\
                        /self.rotor.Kbe
        
                # Recalculamos la fuerza lateral del rotor con la mejor 
                # aproximacion de be1s_w
                # _Cyw = 2*Cyw/(a0*s)
                _Cyw = (-cT/(a0*s) + f_F12c()/4. )*be1s_w - f_F11s()/2.*be0 -\
                        f_F12s()/4.*be1c_w + f_F21c()/2.
                YRh_w = dimFR*a0*s/2.*_Cyw  
                YR = YRh_b = dot(LBW[1], [XRh_w, YRh_w, ZRh_w])
                
                # Despejamos fi de la ecuacion de fuerza lateral
                if Cth == 0:
                    sys.exit("Error: th = pi/2")
                Sfi = 1/(M*g*Cth)*( M*(u*r - w*p) - (Yf + Yfn + Ytp + YT + YR) )
                if abs(Sfi) > 1:
                    sys.exit("Error: sin(fi)>1")
                Cfi = 1 - Sfi**2
                if Cfi<0:
                    Cfi = 0.0
                fi_0 = fi
                fi = fi_0 + k_fi*(asin(Sfi) - fi_0)
                err_fi = abs(fi - fi_0)
                n_fi += 1
                # print "YR = %f, YT = %f, LR = %f, LT = %f" % (YR, YT, LR, LT)
                print "\tfi= %f, la0 = %f, be1s_w = %f, err_fi = %f" \
                        % (fi, la0, be1s_w, err_fi)

                # FIN BUCLE FI #
            
            ####################################################################
            #   VELOCIDAD INDUCIDA DE COLA                                     #
            ####################################################################
            # Velocidad de giro del rotor de cola
            OmT = self.rotor_cola.gT*Om
            
            # Velocidad de la cola, relativa al viento en ejes cuerpo
            VTrw_b = self.Vrw_b +\
                    cross(  W_b, 
                            [-(self.rotor.xcg + self.rotor_cola.lT), 
                                0.0, -self.rotor_cola.hT])
            
            # Velocidad de la cola, relativa al viento en ejes cola
            uTrw_T, vTrw_T, wTrw_T = dot(self.rotor_cola.LTB, VTrw_b)
            
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
            nuT =  uTrw_w/( OmT*self.rotor_cola.RT )
            nuzT = wTrw_w/( OmT*self.rotor_cola.RT )
            
            # Coeficiente de traccion de la cola
            dimFT = ro_t*(OmT*self.rotor_cola.RT)**2*pi*self.rotor_cola.RT**2
            FT = 1. - 3./4*self.estabilizador_vertical.Sfn/( 
                    pi*self.rotor_cola.RT**2 )
            cTT = TT/(dimFT*FT)
            
            # Resistencia
            deT = self.rotor_cola.de0T + self.rotor_cola.de2T*cTT**2
            
            def fT(la0T):
                def VTT(la0):
                    return (nuT/knu)**2 + (1./knu**2 - 1./ki**2)*nuzT**2  +\
                            ((nuzT - la0T)/ki)**2
                return 2*sqrt(VTT(la0T))*la0T/ki - cTT

            la0T = solveSecante(fT, 0.1, 0.5, 100, 1e-12)

            ####################################################################
            #   ECUACION DEL MOTOR                                             #
            ####################################################################
            # Par de la cola
            QT = self.rotor_cola.RT*( -(nuzT - la0T)*TT/FT +\
                    dimFT*deT*self.rotor_cola.sT/8.*( 1. + 3*nuT**2) )
            
            # Par del motor
            Q1 = (1. + self.motor.P)*\
                    (QR + self.rotor_cola.gT*QT)/self.motor.nmotores
            # Revoluciones del rotor
            self.motor.calcula_K3()
            Om_0 = Om
            self.Om = Om = Om_0 +\
                    k_Om*(self.motor.Omi - Q1/self.motor.K3 - Om_0)
            err_Om = abs(Om - Om_0)
            n_Om += 1
            
            ####################################################################
            #   CONTROLES                                                      #
            ####################################################################
            # Resolvemos simulataneamente los controles 
            # th0, th1c_w, th1s_w, th0T y los batimientos
            # t velocidades inducidas que todavia no se han calculado: 
            # be0, be1cT_w, be1sT_w, la1cw_w,
            # la1s_w. Tambien calculamos las  fuerzas y momentos que faltan del
            # rotor de cola

            be0, la1s_w, la1c_w, th0, th1s_w, th1c_w = f_controles()

            th1c =  Cchw*th1c_w + Schw*th1s_w
            th1s = -Schw*th1c_w + Cchw*th1s_w
            
            be1c =  Cchw*be1c_w + Schw*be1s_w
            be1s = -Schw*be1c_w + Cchw*be1s_w
        
            #
            # Batimiento y controles del rotor de cola
            #
            
            # Numero de Locke
            gaT = ro_t*self.rotor_cola.cT*self.rotor_cola.a0T*\
                    self.rotor_cola.RT**4/self.rotor_cola.IbeT
            # Frecuencia natural de batimiento al cuadrado
            labeT2 = 1. + self.rotor_cola.KbeT/( self.rotor_cola.IbeT*OmT**2 )
            
            # Lado izquierdo de la ecuacion de batimiento
            A_bebeT = array(shape=(3, 3), type='Float64')
        
            A_bebeT[0, 0] = (8.*labeT2/gaT) - self.rotor_cola.k3*(1. + nuT**2)
            A_bebeT[0, 1] =  0. 
            A_bebeT[0, 2] = -self.rotor_cola.k3*4./3.*nuT 
        
            A_bebeT[1, 0] =  4./3.*nuT
            A_bebeT[1, 1] =  8.*(labeT2 - 1.)/gaT - \
                    self.rotor_cola.k3*(1. + 0.5*nuT**2)
            A_bebeT[1, 2] =  1 + 0.5*nuT**2 
        
            A_bebeT[2, 0] = -self.rotor_cola.k3*8./3.*nuT
            A_bebeT[2, 1] = -1. + 0.5*nuT**2
            A_bebeT[2, 2] =  8.*(labeT2 - 1.)/gaT -\
                    self.rotor_cola.k3*(1. + 1.5*nuT**2)

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


            lhs = array(shape=(4,4), type='Float64')
            lhs[0:3, 0:3] =  A_bebeT
            lhs[0:3, 3  ] = -A_beth0T
            lhs[3, 0] = self.rotor_cola.k3*(1./3. + nuT**2/2.)
            lhs[3, 1] = 0.
            lhs[3, 2] = nuT/2.*self.rotor_cola.k3
            lhs[3, 3] = 1./3. + nuT**2/2.

            rhs = array(shape=(4), type='Float64')
            rhs[0:3] = (nuzT-la0T)*A_bela0T + self.rotor_cola.tht*A_bethtT
            rhs[3] = 2*cTT/(self.rotor_cola.a0T*self.rotor_cola.sT) -\
                    0.5*(nuzT - la0T)

            be0T, be1cT_w, be1sT_w, th0T = la.solve_linear_equations(lhs, rhs)

            # Batimientos en ejes cola
            be1cT = be1cT_w*CchwT + be1sT_w*SchwT
            be1sT = be1sT_w*CchwT - be1cT_w*SchwT
            
            # Calculamos las fuerzas y momentos del rotor de cola
            XT_T =  TT*be1cT
            YT_T = -TT*be1sT
            ZT_T = -TT
            LT_T =  0.0
            MT_T =  0.0
            NT_T =  QT

            # Fuerzas y momentos en ejes cuerpo
            XT_b, YT_b, ZT_b = dot(self.rotor_cola.LBT, [XT_T, YT_T, ZT_T])
            LT_b, MT_b, NT_b = dot(self.rotor_cola.LBT, [LT_T, MT_T, NT_T])
            
            # Fuerzas y momentos en ejes cuerpo respecto al centro de masas
            l = self.rotor_cola.lT + self.rotor.xcg
            XT, YT, ZT = XT_b, YT_b, ZT_b
            LT = LT_b + YT_b*self.rotor_cola.hT
            MT = MT_b - XT_b*self.rotor_cola.hT + ZT_b*l
            NT = NT_b - YT_b*l

            # Controles del piloto
            et0p, et1sp, et1cp, etp =\
                    f_etp(th0, th1s, th1c, th0T, p, q, r, th, fi, 0.0)
                
            print "be0 = %f, la1s_w = %f, la1c_w = %f" %\
                    (be0, la1s_w, la1c_w)
            print "be0T = %f, be1cT_w = %f, be1sT_w = %f" \
                    % (be0T, be1cT_w, be1sT_w)
            print "cTT = %f, la0T = %f, Om = %f, err_Om = %f" \
                    % (cTT, la0T, Om, err_Om)         
            print ( "********************************************************" +
                    "************************" )
            print "th0 = %f, th1s = %f, th1c = %f, th0T = %f"  \
                    % (th0, th1s, th1c, th0T)
            print ( "********************************************************" +
                    "************************" )
        
        # Colocamos al helicopter en estado de trimado
        _Q = Quaternion.rotacion([
                (pi, 0, 1, 0), # Ejes inerciales a auxiliares
                (th, 0, 1, 0), # Angulos de euler, ch = 0
                (fi, 1, 0, 0)])
        
        self.u, self.v, self.w = u, v, w
        self.p, self.q, self.r = p, q, r
        self.q0, self.q1, self.q2, self.q3 = _Q.q0, _Q.q1, _Q.q2, _Q.q3
        self.Om, self.Q1, self.DTQ1 = Om, Q1, 0.0
        self.z_i = z_t
        
        self.controles.th0 =  th0
        self.controles.th1c = th1c
        self.controles.th1s = th1s
        self.controles.th0T = th0T

        # Varios valores intermedios calculados

        return {
                'th':       th,
                'fi':       fi,
                'Om':       Om,
                'be1c_w':   be1c_w,
                'be1s_w':   be1s_w,
                'be0':      be0,
                'la0':      la0,
                'la1c_w':   la1c_w,
                'la1s_w':   la1s_w,
                'be1c':     be1c,
                'be1s':     be1s,
                'cT':       cT,
                'cTT':      cTT,
                'be1cT':    be1cT,
                'be1sT':    be1sT,
                'be1cT_w':  be1cT_w,
                'be1sT_w':  be1sT_w,
                'la0T':     la0T,
                'u':        u,
                'v':        v,
                'w':        w,
                'p':        p,
                'q':        q,
                'r':        r,
                'th0':      th0,
                'th1c':     th1c,
                'th1s':     th1s,
                'th0T':     th0T,
                'Xf':       Xf,
                'Yf':       Yf,
                'Zf':       Zf,
                'Lf':       Lf,
                'Mf':       Mf,
                'Nf':       Nf,
                'Xtp':      Xtp,
                'Ytp':      Ytp,
                'Ztp':      Ztp,
                'Ltp':      Ltp,
                'Mtp':      Mtp,
                'Ntp':      Ntp,
                'Xfn':      Xfn,
                'Yfn':      Yfn,
                'Zfn':      Zfn,
                'Lfn':      Lfn,
                'Mfn':      Mfn,
                'Nfn':      Nfn,
                'XT':       XT,
                'YT':       YT,
                'ZT':       ZT,
                'LT':       LT,
                'MT':       MT,
                'NT':       NT,
                'XR':       XR,
                'YR':       YR,
                'ZR':       ZR,
                'LR':       LR,
                'MR':       MR,
                'NR':       NR,
                'q0':       _Q.q0,
                'q1':       _Q.q1,
                'q2':       _Q.q2,
                'q3':       _Q.q3,
                'DTQ1':     0.0,
                'Q1':       Q1,
                'QT':       QT,
                'nu':       nu,
                'nuz':      nuz,
                'altp0':    self.estabilizador_horizontal.altp0,
                'et0p':     et0p,
                'et1sp':    et1sp,
                'et1cp':    et1cp,
                'etpp':     etp,
                'altp':     self.estabilizador_horizontal.al,
                'cLtp':     self.estabilizador_horizontal.cL,
                'cDtp':     self.estabilizador_horizontal.cD,
                'befn':     self.estabilizador_vertical.be,
                'cLfn':     self.estabilizador_vertical.cL,
                'cDfn':     self.estabilizador_vertical.cD,
                'alf':      self.fuselaje.al,
                'cLf':      self.fuselaje.cL,
                'cDf':      self.fuselaje.cD,
                'cYf':      self.fuselaje.cY,
                'clf':      self.fuselaje.cl,
                'cmf':      self.fuselaje.cm,
                'cnf':      self.fuselaje.cn,
            }


    
