# -*- coding: latin-1 -*-
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#           NO TOCAR                                                          #
from modelo import Modelo
from matematicas import *
from numarray import *
from util import *
from integrador import *

import logging

modelo = Modelo()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

###############################################################################
# A PARTIR DE ESTE PUNTO INTRODUCIR LOS PARAMETROS CORRECTOS DEL HELICOPTERO  #
###############################################################################

###############################################################################
# MISCELANEO                                                                  #
###############################################################################
modelo.nombre = "UH60"
modelo.integrador = AdamsBashforth2()

###############################################################################
# DATOS MASICOS                                                               #
###############################################################################
modelo.M = 7438.9       # Masa del helicoptero  
modelo.Ixx = 7631.9     # Momentos de inercia
modelo.Iyy = 54232.7    #     [  Ixx   0   -Ixz ]
modelo.Izz = 50436.4    # I = [   0   Iyy    0  ]
modelo.Ixz = 2264.1     #     [ -Ixz   0    Izz ]

###############################################################################
# PARAMETROS DEL ROTOR                                                        #
###############################################################################
modelo.rotor.gas = 0.05236  # Inclinacion del rotor, positivo hacia adelante
modelo.rotor.hR = 1.722     # Distancia vertical del CM al rotor, positivo 
                            # hacia arriba
modelo.rotor.xcg = -0.274   # Distancia horizontal, positivo hacia atras
modelo.rotor.h0 = 2.946     # Distancia del rotor al suelo

modelo.rotor.Nb = 4         # Numero de palas
modelo.rotor.R = 8.178      # Radio del rotor
modelo.rotor.c = 0.527      # Cuerda de las palas
modelo.rotor.Ibe = 2050.81  # Momento de inercia de la pala
modelo.rotor.Kbe = 107120   # Constante elastica del muelle equivalente
modelo.rotor.tht = -0.209   # Torsion lineal de la pala

modelo.rotor.a0 = 5.62      # Pendiente de sustenacion
modelo.rotor.de0 = 0.0087   # Coeficiente de resistencia parasita
modelo.rotor.de1 = -0.2448
modelo.rotor.de2 = 9.5      # Coeficiente de resistencia inducida

###############################################################################
# PARAMETROS DEL ROTOR DE COLA                                                #
###############################################################################
modelo.rotor_cola.hT = 1.968    # Distancia vertical de la cola
modelo.rotor_cola.lT = 9.723    # Distancia horizontal de la cola
modelo.rotor_cola.K = rad(20)   # Inclinacion respecto al plano vertical 

modelo.rotor_cola.RT = 1.676    # Radio del rotor de cola
modelo.rotor_cola.cT = 0.2469   # Cuerda de la cola
modelo.rotor_cola.IbeT = 4.203  # Momento de inercia de la pala de cola
modelo.rotor_cola.KbeT = 0      # Muelle equivalente, valor desconocido
modelo.rotor_cola.k3 = -rad(35) # Acoplamiento paso-batimiento = tan(de3)
modelo.rotor_cola.tht = -0.209  # Torsion lineal
modelo.rotor_cola.gT = 4.616    # Reduccion transmision cola-rotor principal

modelo.rotor_cola.a0T = 5.73    # Pendiente de sustentacion
modelo.rotor_cola.de0T = 0.009  # Coeficiente de resistencia parasita
modelo.rotor_cola.de1T = 0.0    #
modelo.rotor_cola.de2T = 9.356  # Coeficiente de resistencia inducida

###############################################################################
# PARAMETROS DEL MOTOR                                                        #
###############################################################################
modelo.motor.Omi = 27.0         # Velocidad de rotacion del rotor para par motor cero

modelo.motor.Wto = 1240000  # Potencia de despegue de un motor
modelo.motor.Wmc = 1240000  # Potencia maxima continua de un motor
modelo.motor.Wid = 1000     # Potencia consumida para motor sin carga
modelo.motor.Irt = [        # Inercias referidas al rotor del sistema  
                            # de la transmision
             3085.069,      # Ningun motor activo
             3344.716,      # 1 motor activo
             3344.716]      # 2 motores activo

# Primera constante de tiempos 
modelo.motor.tae1 = 0.1
# Segunda constante de tiempos 
modelo.motor.a2 = 0.1
modelo.motor.b2 = 0.1
modelo.motor.c2 = 0.1

# Tercera constante de tiempos 
modelo.motor.a3 = 0.1
modelo.motor.b3 = 0.1
modelo.motor.c3 = 0.1

modelo.motor.K3 = 65e3          # "Rigidez" del sistema motor

modelo.motor.n = 2              # numero de motores
modelo.motor.P  = 0.1           # Perdidas de transmision conjuntas de rotor 
                                # y cola
###############################################################################
# PARAMETROS DEL FUSELAJE                                                     #
###############################################################################
modelo.fuselaje.Ss = 41     # Superficie de adimensionalizacion para 
                            # fuerzas/momentos laterales
modelo.fuselaje.Sp = 6      # Superficie de adimensionalizacion para 
                            # fuerzas/momentos longitudinales
modelo.fuselaje.lf = 15     # Longitud de adimensionalizacion
modelo.fuselaje.xca = 0.094 # Posicion del centro aerodinamico respecto al 
                            # Centro de Masas en ejes cuerpo
modelo.fuselaje.zca = 0.335 
modelo.fuselaje.ht = 1.67
# Para el Black Hawk lo mas comodo es definir la funcion analiticamente
def aero(al, be):
    if al<-pi/2:
        al += pi
    elif al>pi/2:
        al -= pi
    if be<-pi/2:
        be += pi
    elif be>pi/2:
        be -= pi

    Sal, S2al, Sbe, S2be, S4be = map(sin, [al, 2*al, be, 2*be, 4*be])
    Cal, C2al, Cbe, C2be, C4be = map(cos, [al, 2*al, be, 2*be, 4*be])

    cD = 1.3944*(Sal)**2 - 0.6435*Cal + 2.4806 +\
         0.0456*C4be - 1.597*C2be - 8.2893e-9*(be**4)
    cL = 0.4546*Sal + 0.6730*S2al - 1.2680*(Sal)**2 - 1.3029*Cal + 1.3215 +\
         0.00127*be - 0.0465*S4be + 0.0005*(be**2) 
    cY = -0.0802*Sbe -0.1627*S2be + 0.0182*S4be
    
    if be<=rad(10):
        cl = 0.0
    elif rad(10)<be and be<=rad(25):
        cl = sign(be)*0.02098*(Cbe**4) - 0.0197
    else:
        cl = 0.0283*Sbe + sign(be)*(0.0022*C4be - 0.0134*(Cbe)**3 +\
             0.0338*(Cbe)**4 - 0.0308)
    
    cm = 0.0007 + 0.2291*S2al + 0.1343*(Sal)**2 + 0.1095*Cal -\
         0.1606*(Cbe)**3 + 0.0176
    
    if be<=rad(20):
        cn = 0.0128*S2be - 0.0195*S4be - 8.4339e-5
    else:
        cn = -0.0101*S2be + sign(be)*(0.0309*(C4be**4) - 0.0197)
    
    return cD, cL, cY, cl, cm, cn

modelo.fuselaje.aero = FuselajeB(aero)

###############################################################################
# PARAMETROS DEL ESTABILIZADOR HORIZONTAL                                     #
################################################################################ Distancia horizontal del estabilizador horizontal
modelo.estabilizador_horizontal.ltp = 8.920     
# Área del estabilizador horizontal
modelo.estabilizador_horizontal.Stp = 4.18      
# Distancia vertical del estabilizador horizontal
modelo.estabilizador_horizontal.htp = -0.081        
# Ángulo de ataque del estabilizador horizontal
modelo.estabilizador_horizontal.altp0 = 0.0     
modelo.estabilizador_horizontal.aero = perfilA(     
        AR  =  4.6,         # Alargamiento
        a = 2.72,           # Pendiente de sustentacion
        cLmax = 1.03        # Coeficiente de sustentación máximo
)

###############################################################################
# PARAMETROS DEL ESTABILIZADOR VERTICAL                                       #
###############################################################################
# Distancia horizontal del estabilizador vertical
modelo.estabilizador_vertical.lfn = 8.783           
# Distancia vertical del estabilizador vertical
modelo.estabilizador_vertical.hfn = 0.655           
# Área del estabilizador vertical
modelo.estabilizador_vertical.Sfn = 3.0         
# Ángulo de ataque del estabilizador vertical
modelo.estabilizador_vertical.befn0 = 0.0           

modelo.estabilizador_vertical.aero = perfilA(
        AR =   1.92,        # Alargamiento
        cLmax = 0.89,       # Coeficiente de sustentación máximo
        De = 0.7156         # Flecha, rad
)
###############################################################################
# CONTROLES                                       #
###############################################################################
# { th0  }          { et0  }
# { th1s } =  c0 + c1*  { et1s }
# { th1c }              { et1c }
# { th0T }              { etp  }

modelo.controles.c0 = array([ 
             0.1728,    
            -0.0008,    
             0.0000, 
             0.2617, 
            ], type='Float64')


modelo.controles.c1 = array([
            [ 1.0990, 0.0000,  0.0000,  0.0000 ],
            [ 0.2343, 1.8614,  0.0000,  0.0000 ],
            [ 0.0000, 0.0000,  1.0990,  0.0000 ],
            [ 0.0000, 0.0000,  0.0000,  3.8060 ],
            ], type='Float64')

# { et0  }       { et0p  }       { p }      { th }
# { et1s } = S0* { et1sp } + S1* { q } + S2*{ fi }
# { et1c }       { et1cp }       { r }      { ch }
# { etp  }       { etpp  }       

modelo.controles.S0 = array([
    [ 1.0000, 0.0000, 0.0000, 0.0000 ],
    [ 0.0000, 1.0000, 0.0000, 0.0000 ],
    [ 0.0000, 0.0000, 1.0000, 0.0000 ],
    [ 0.0000, 0.0000, 0.0000, 1.0000 ]], type='Float64')

modelo.controles.S1 = array([
    [ 0.0000,  0.0000, 0.0000 ],
    [ 0.0000, -0.6000, 0.0000 ],
    [ 0.3000,  0.0000, 0.0000 ],
    [ 0.0000,  0.0000, 0.6000 ]], type='Float64')

modelo.controles.S2 = array([
    [ 0.0000,  0.0000, 0.0000 ],
    [ 0.0100,  0.0000, 0.0000 ],
    [ 0.0000, -0.0100, 0.0000 ],
    [ 0.0000,  0.0000, 0.0000 ]], type='Float64')


modelo.controles.L = 0.2
#----------------------------------------------------------------------------#
# NO IMPLEMENTADO                                                            #
#----------------------------------------------------------------------------#
# Definimos el comportamiento del estabilizador como una funcion
# que depende de u, q y et0p
#def estabilizador(u, q, et0p):
#    def clamp(x, a, b):
#        if x<a:
#            return a
#        if x>b:
#            return b
#        else:
#            return x
#
#    qB = deg(q)     # qB en deg/sec
#    uB = u/0.514444     # uB en kt
#    
#    # dC entre 0 y 100%, 
#    et0p_min = 0.
#    et0p_max = 0.254
#    dC = (et0p - et0p_min)/(et0p_max - et0p_min)*100    
#    dC = 0
#    
#    u2 = uB - 30
#    u3 = clamp(u2, 0, 50)
#    u1 = uB - 80
#    u4 = clamp(u1, 0, 67)
#    u5 = -0.1119*u4
#
#    q1 = 0.16*qB
#    q2 = q1 + 42
#
#    d1 = -dC
#    d2 = clamp(d1 + 50, 0, 50)
#    d3 = clamp(d1 + 70, 0, 20)
#    d4 = d2*0.003286
#    d5 = d3*0.003286
#    d6 = d4 + d5 + 0.67
#    d7 = -d6
#    d8 = d7*u3
#
#    a1 = q2 + u5 + d8
#    a2 = clamp(a1, -8, 39)
#    return 0.01745*a2
#
#
#modelo.controles.estabilizador_horizontal = estabilizador
#----------------------------------------------------------------------------#
# FIN NO IMPLEMENTADO                                                        #
#----------------------------------------------------------------------------#


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#   FIN DE LA CONFIGURACION                                                   #
#   El siguiente codigo no hace falta tocarlo, se utiliza solo                #
#   para funciones de debug                                                   #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

###############################################################################
# TESTING                                                                     #
###############################################################################

coefs = modelo.fuselaje.aero.coefs
al = arange(-pi, pi, 0.01)
c = [ [] for x in range(6) ]
for x in al:
    cD, cL, cY, cl, cm, cn = coefs(x, 0.0)
    c[0].append(cD)
    c[1].append(cL)
    c[2].append(cY)
    c[3].append(cl)
    c[4].append(cm)
    c[5].append(cn)

# Velocidad
modelo.u = 0.0
modelo.v = 0.0
modelo.w = 0.0
modelo.p = 0.0
modelo.q = 0.0
modelo.r = 0.0

# Velocidad del viento
modelo.Vw_i = array([0, 0, 0], type='Float64')

# Controles
modelo.controles.th0 = 0.25038
modelo.controles.th1c = 0.0055
modelo.controles.th1s = 0.0218
modelo.controles.th0T = 0.24795
modelo.controles.th1sT = 0.0

# Motores activos
modelo.motor.nmotores = 2

# Rotor/Motor
modelo.Om = 35
modelo.Q1 = 0.0
modelo.DTQ1 = 0.0

# Densidad del aire
modelo.ro = 1.225

# Altura
modelo.h = 1.0

# Orientacion
th = -0.011458
fi = -0.0516
Om = 35.134769

modelo.euler(th, fi, 0.0)

# Posicion
modelo.x_i = 0.0
modelo.y_i = 0.0
modelo.z_i = 100.0

# Para debug
if __name__ == "__main__":
    import time
    modelo.init()
    tr = modelo.trimado(10.0, 0.0, 0.0, 0.0, 1000., 1.225)

    #modelo.trim = True
    #modelo.calcula_ecs_modelo()

    modelo.hG = 0
    t0 = time.time()
    T = 10.
    s1 = modelo.series(T, 0, 0, 0, 0, 0.03, 
            'u', 'v', 'w', 'p', 'q', 'r', 'la0','be0','be1c_w', 'be1s_w',
            'la1c_w', 'la1s_w', 'xi')
    t1 = time.time()
    print "tiempo integrado = %f, tiempo real = %f, speedup=%f" %\
            (T, t1 - t0, T/(t1 -t0)) 
    
    tr = modelo.trimado(10.0, 0.0, 0.0, 0.0, 1000., 1.225)
    t0 = time.time()
    T = 10.
    modelo.t = 0.0
    s2 = modelo.series(T, 0, 0, 0, 0, 0.01, 
            'u', 'v', 'w', 'p', 'q', 'r', 'la0','be0','be1c_w', 'be1s_w',
            'la1c_w', 'la1s_w', 'xi')
    t1 = time.time()

    from pylab import *
    plot(s1['t'], s1['u'])
    plot(s2['t'], s2['u'])
    xlabel('t')
    ylabel('u')
    show()
    plot(s1['t'], s1['v'])
    plot(s2['t'], s2['v'])
    show()
    plot(s1['t'], s1['w'])
    plot(s2['t'], s2['w'])
    show()
    plot(s1['t'], s1['p'])
    plot(s2['t'], s2['p'])
    show()
    plot(s1['t'], s1['q'])
    plot(s2['t'], s2['q'])
    show()
    plot(s1['t'], s1['r'])
    plot(s2['t'], s2['r'])
    show()
    plot(s1['t'], s1['be1c_w'])
    plot(s2['t'], s2['be1c_w'])
    show()
    plot(s1['t'], s1['la1s_w'])
    plot(s2['t'], s2['la1s_w'])
    show()

