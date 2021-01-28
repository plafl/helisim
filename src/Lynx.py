# -*- coding: latin-1 -*-
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#           NO TOCAR                          #
from modelo import Modelo
from matematicas import *
from numarray import *
from aero import *
from input import *
from integrador import *

import logging

modelo = Modelo()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

################################################################################
# A PARTIR DE ESTE PUNTO INTRODUCIR LOS PARAMETROS CORRECTOS DEL HELICOPTERO   #
################################################################################

################################################################################
# MISCELANEO                                                                   #
################################################################################
modelo.nombre = "Lynx"
modelo.integrador = AdamsBashforth2()

################################################################################
# DATOS MASICOS                                                                #
################################################################################
modelo.M = 4313.7       # Masa del helicoptero  
modelo.Ixx = 2767.1     # Momentos de inercia
modelo.Iyy = 13904.5    #     [  Ixx   0   -Ixz ]
modelo.Izz = 12208.8    # I = [   0   Iyy    0  ]
modelo.Ixz = 2034.8     #     [ -Ixz   0    Izz ]

################################################################################
# PARAMETROS DEL ROTOR                                                         #
################################################################################
modelo.rotor.gas = rad(4)   # Inclinacion del rotor, positivo hacia adelante
modelo.rotor.hR = 1.274     # Distancia vertical del CM al rotor, positivo 
                            # hacia arriba
modelo.rotor.xcg = -0.0198  # Distancia horizontal, positivo hacia atras
modelo.rotor.h0 = 2.946     # Distancia del rotor al suelo

modelo.rotor.Nb = 4         # Numero de palas
modelo.rotor.R = 6.4        # Radio del rotor
modelo.rotor.c = 0.391      # Cuerda de las palas
modelo.rotor.Ibe = 678.14   # Momento de inercia de la pala
modelo.rotor.Kbe = 166352   # Constante elastica del muelle equivalente
modelo.rotor.tht = -0.14    # Torsion lineal de la pala

modelo.rotor.a0 = 6.0       # Pendiente de sustenacion
modelo.rotor.de0 = 0.009    # Coeficiente de resistencia parasita
modelo.rotor.de1 = 0.0
modelo.rotor.de2 = 37.983   # Coeficiente de resistencia inducida

################################################################################
# PARAMETROS DEL ROTOR DE COLA                                                 #
################################################################################
modelo.rotor_cola.hT = 1.146        # Distancia vertical de la cola
modelo.rotor_cola.lT = 7.66         # Distancia horizontal de la cola
modelo.rotor_cola.K = 0.0           # Inclinacion respecto al plano vertical 

modelo.rotor_cola.RT = 1.106        # Radio del rotor de cola
modelo.rotor_cola.cT = 0.18001      # Cuerda de la cola
modelo.rotor_cola.IbeT = 1.08926    # Momento de inercia de la pala de cola
modelo.rotor_cola.KbeT = 2511.24    # Muelle equivalente
modelo.rotor_cola.k3 = -1.          # Acoplamiento paso-batimiento = tan(de3)
modelo.rotor_cola.tht = 0.0         # Torsion lineal
modelo.rotor_cola.gT = 5.8          # Reduccion transmision cola-rotor principal

modelo.rotor_cola.a0T = 6.0         # Pendiente de sustentacion
modelo.rotor_cola.de0T = 0.008      # Coeficiente de resistencia parasita
modelo.rotor_cola.de1T = 0.0        #
modelo.rotor_cola.de2T = 5.334      # Coeficiente de resistencia inducida

################################################################################
# PARAMETROS DEL MOTOR                                                         #
################################################################################
modelo.motor.Omi = 35.63        # Velocidad de rotacion del rotor para par motor
                                # cero
modelo.motor.Wto_0 = 746000     # Potencia de despegue de un motor a nivel del mar
modelo.motor.Wmc_0 = 664000     # Potencia maxima continua de un motor a nivel del mar
modelo.motor.x     = 0.85       # Variacion de la potencia disponible
modelo.motor.Wid   = 1000       # Potencia consumida para motor sin carga
                    
modelo.motor.Irt = [            # Inercias referidas al rotor del sistema 
                                # de la transmision
                3085.069,       # Ningun motor activo
                3344.716,       # 1 motor activo
                3344.716]       # 2 motores activo

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

################################################################################
# PARAMETROS DEL FUSELAJE                                                      #
################################################################################
modelo.fuselaje.Ss = 32.        # Superficie de adimensionalizacion para 
                                # fuerzas/momentos laterales
modelo.fuselaje.Sp = 24.        # Superficie de adimensionalizacion para 
                                # fuerzas/momentos longitudinales
modelo.fuselaje.lf = 12.        # Longitud de adimensionalizacion
modelo.fuselaje.xca = 0.139     # Posicion del centro aerodinamico respecto al 
                                # Centro de Masas en ejes cuerpo
modelo.fuselaje.zca = 0.190     
                                # Distancia tren de aterrizaje-centro de masas, 
                                # positivo hacia abajo
modelo.fuselaje.ht = 0.762

modelo.fuselaje.aero = FuselajeA(  # Modelo aerodinamico de fuselaje
    unidades = 'deg',       # Unidades:  grados:   deg
                            #            radianes: rad  
    cDa90 = 0.3480,         # Coeficiente de resistencia para al = 90º
    cDb90 = 0.4784,         # Coeficiente de resistencia para be = 90º
    cm90 =  0.03,           # Coeficiente de cabeceo para al = 90º
    cl90 =  0.01,           # Coeficiente de balance para be = 90º
    cn90 =  0.1,            # Coeficiente de guinada para be = 90º
    
    # cDa, cL, cm definidos entre al1, al2
    al1 = -21,
    al2 =  21,
    # cDb, cY, cl, cn definidos entre be1, be2
    be1 = -21,
    be2 =  21,
    
    # Coeficientes para angulos pequenos
    # Coeficiente de resistencia en funcion de angulo de ataque 
    cDa = Lista(
        -21, 0.06661,
        -18, 0.06010,
        -15, 0.05531,
        -12, 0.04991,
         -9, 0.04707,
         -6, 0.04561,
         -3, 0.04421,
          0, 0.04233,
          3, 0.04245,
          6, 0.04372,
          9, 0.04486,
         12, 0.04685,
         15, 0.04926,
         18, 0.05170,
         21, 0.05613 
    ),
    
    # Coeficiente de resistencia en funcion de angulo de resbalamiento
    cDb = Lista(
        -21, 0.05532,
        -18, 0.04096,
        -15, 0.02863,
        -12, 0.01842,
         -9, 0.01040,
         -6, 0.00464,
         -3, 0.00116,
          0, 0.00000,
          3, 0.00116,
          6, 0.00464,
          9, 0.01040,
         12, 0.01842,
         15, 0.02863,
         18, 0.04096,
         21, 0.05532 
    ),

    # Coeficiente de cabeceo en funcion de angulo de ataque
    cm = Lista(
        -21, -0.02415,
        -18, -0.02210,
        -15, -0.02015,
        -12, -0.01803,
         -9, -0.01588,
         -6, -0.01351,
         -3, -0.01074,
          0, -0.00747,
          3, -0.00421,
          6, -0.00085,
          9,  0.00237,
         12,  0.00573,
         15,  0.00875,
         18,  0.01197,
         21,  0.01433 
    ),
    
    # Coeficiente de sustentacion en funcion de angulo de ataque
    cL = Lista(
        -21, -0.02415,
        -18, -0.02210,
        -15, -0.02015,
        -12, -0.01803,
         -9, -0.01588,
         -6, -0.01351,
         -3, -0.01074,
          0, -0.00747,
          3, -0.00421,
          6, -0.00085,
          9,  0.00237,
         12,  0.00573,
         15,  0.00875,
         18,  0.01197,
         21,  0.01433 
    ),
        
    # Coeficiente de fuerza lateral en funcion de angulo de resbalamiento
    cY = Lista( 
        -21,  0.17547,
        -18,  0.15138,
        -15,  0.12677,
        -12,  0.10178,
         -9,  0.07653,
         -6,  0.05110,
         -3,  0.02557,
          0,  0.00000,
          3, -0.02557,
          6, -0.05110,
          9, -0.07653,
         12, -0.10178,
         15, -0.12677,
         18, -0.15138,
         21, -0.17547 
    ),
    
    # Coeficiente de balance en funcion de angulo de resbalamiento
    cl = Lista(
        -21, -0.00697,
        -18, -0.00472,
        -15, -0.00270,
        -12, -0.00097,
         -9,  0.00000,
         -6,  0.00000,
         -3,  0.00000,
          0,  0.00000,
          3,  0.00000,
          6,  0.00000,
          9,  0.00000,
         12,  0.00097,
         15,  0.00270,
         18,  0.00472,
         21,  0.00697 
    ),
    
    # Coeficiente de guinada en funcion de angulo de resbalamiento
    cn = Lista(
        -21, -0.01704,
        -18, -0.01461,
        -15, -0.01217,
        -12, -0.00974,
         -9, -0.00730,
         -6, -0.00487,
         -3, -0.00243,
          0,  0.00000,
          3,  0.00243,
          6,  0.00487,
          9,  0.00730,
         12,  0.00974,
         15,  0.01217,
         18,  0.01461,
         21,  0.01704
    ),
)


################################################################################
# PARAMETROS DEL ESTABILIZADOR HORIZONTAL                                      #
################################################################################
# Distancia horizontal del estabilizador horizontal
modelo.estabilizador_horizontal.ltp = 7.66  
# Distancia vertical del estabilizador horizontal
modelo.estabilizador_horizontal.htp = 1.146        
# Area del estabilizador horizontal
modelo.estabilizador_horizontal.Stp = 1.197        
# Angulo de ataque del estabilizador horizontal
modelo.estabilizador_horizontal.altp0 = rad(-1.0)      
# Modelo de estabilizador horizontal
modelo.estabilizador_horizontal.aero = perfilA(         
        AR  =  2.7,         # Alargamiento
        a   =  2.3663,      # Pendiente de sustentacion
        cLmax = 0.6,        # Coeficiente de sustentacion maximo
        de0 =  1.065e-3,    # Coeficientes de resistencia:
        de1 = -8.4703e-2,   # de = de0 + de1*al + de2*al^2
        de2 =  1.46981
)

################################################################################
# PARAMETROS DEL ESTABILIZADOR VERTICAL                                        #
################################################################################
# Distancia horizontal del estabilizador vertical
modelo.estabilizador_vertical.lfn = 7.48           
# Distancia vertical del estabilizador vertical
modelo.estabilizador_vertical.hfn = 0.57           
# Area del estabilizador vertical
modelo.estabilizador_vertical.Sfn = 1.107          
# Angulo de ataque del estabilizador vertical
modelo.estabilizador_vertical.befn0 = 0.0           

modelo.estabilizador_vertical.aero = perfilA( # Modelo de estab. horizontal
        AR =   2.7,         # Alargamiento
        a =    2.5,         # Pendiente de sustentacion
        cL0 =  0.,          # Sustentacion para angulo de ataque nulo
        cLmax = 0.9,        # Coeficiente de sustentacion maximo
        de0 =  1.065e-3,    # Coeficientes de resistencia:
        de1 = -8.4703e-2,   # de = de0 + de1*al + de2*al^2
        de2 =  1.46981
)

################################################################################
# CONTROLES                                                                    #
################################################################################
# { th0  }              { et0  }
# { th1s } =  c0 + c1*  { et1s }
# { th1c }              { et1c }
# { th0T }              { etp  }

modelo.controles.c0 = array([ 
              0.059,    
              0.000,    
              0.000, 
             -0.061, 
            ], type='Float64')


modelo.controles.c1 = array([
            [ 0.300, 0.000,  0.000,  0.000 ],
            [ 0.000, 0.100,  0.000,  0.000 ],
            [ 0.000, 0.000,  0.100,  0.000 ],
            [ 0.442, 0.000,  0.000,  0.200 ],
            ], type='Float64')

# { et0  }       { et0p  }       { p }      { th - th_0}
# { et1s } = S0* { et1sp } + S1* { q } + S2*{ fi - fi_0}
# { et1c }       { et1cp }       { r }      { ch - ch_0}
# { etp  }       { etpp  }       

modelo.controles.S0 = array([
    [ 1.0000, 0.0000, 0.0000, 0.0000 ],
    [ 0.0000, 0.0000, 0.0000, 0.0000 ],
    [ 0.0000, 0.0000, 0.0000, 0.0000 ],
    [ 0.0000, 0.0000, 0.0000, 1.0000 ]], type='Float64')

modelo.controles.S1 = array([
    [  0.0000,   0.0000, 0.0000 ],
    [  0.0000,  -2.0000, 0.0000 ],
    [ +2.0000,   0.0000, 0.0000 ],
    [  0.0000,   0.0000, 2.0000 ]], type='Float64')

modelo.controles.S2 = array([
    [  0.0000,  0.0000, 0.0000 ],
    [ -2.5000,  0.0000, 0.0000 ],
    [  0.0000, +2.5000, 0.0000 ],
    [  0.0000,  0.0000, 0.0000 ]], type='Float64')

def th_0(et1sp):
    r = 0.0603 + pi*et1sp
    if r>pi:
        r -= 2*pi
    elif r<-pi:
        r += 2*pi
    return r

#modelo.controles.th_0 =  th_0
modelo.controles.th_0 =  lambda et1sp:  0.0603 + 1.7*et1sp
modelo.controles.fi_0 =  lambda et1cp: -0.0476 - 1.7*et1cp
modelo.controles.th_0 =  lambda et1sp:  1.0*et1sp
modelo.controles.fi_0 =  lambda et1cp: -1.0*et1cp
modelo.controles.ch_0 =  0.0

# Limite para el control automático
modelo.controles.L = 1.0

################################################################################
# PARAMETROS DEL TREN                                                          #
################################################################################
# x: posicion en ejes cuerpo del punto del tren
# n: coeficiente de rozamiento
# a, b: constante de rigidez y amortiguamiento para la fuerza normal
# c, d: constante de rigidez y amortiguamiento para la primera fuerza tangente
# e, f: constante de rigidez y amortiguamiento para la segunda fuerza tangente
modelo.tren.punto(
    x = [ 1.867, 1.01, 0.762],
    n = 1.0,
    a = 5e5, b = 1e5, c = 5e5, d = 1e5, e = 5e5, f = 1e5 )

modelo.tren.punto(
    x = [-1.147, 1.01, 0.762],
    n = 1.0,
    a = 5e5, b = 1e5, c = 5e5, d = 1e5, e = 5e5, f = 1e5 )

modelo.tren.punto(
    x = [ 1.867, -1.01, 0.762],
    n = 1.0,
    a = 5e5, b = 1e5, c = 5e5, d = 1e5, e = 5e5, f = 1e5 )

modelo.tren.punto(
    x = [-1.147, -1.01, 0.762],
    n = 1.0,
    a = 5e5, b = 1e5, c = 5e5, d = 1e5, e = 5e5, f = 1e5 )

modelo.tren.T0 = 1e-5

modelo.P0 = 101325.
modelo.T0 = 283.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#   FIN DE LA CONFIGURACION                                                    #
#   El siguiente codigo no hace falta tocarlo, se utiliza solo                 #
#   para funciones de debug                                                    #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

################################################################################
# TESTING                                                                      #
################################################################################
if __name__ == "__main__":
    import time
    modelo.init()

    # Propiedades menos vector de estado y controles
    modelo.hG = 0
    modelo.viento_vel = 0
    modelo.viento_dir = 0
    modelo.P = 103100
    modelo.T = 283
    modelo.motor.nmotores = 2
    
    # Vector de estado
    modelo.t = 0.0
    modelo.x_i = modelo.y_i = 0.0
    modelo.z_i = 1000
    tr0 = modelo.trimado(10.0, 0.0, 0.0, 0.0, 1000, 1.225)
    
    t0 = time.time()
    T = 10.
    
    # modelo.series se ocupa de los controles
    s1 = modelo.series(T, 0, 0, 0, 0, 0.02, 
            'u', 'v', 'w', 'p', 'q', 'r', 'la0','be0','be1c_w', 'be1s_w',
            'la1c_w', 'la1s_w', 'x_i', 'th', 'fi', 'ch')
    t1 = time.time()
    print "tiempo integrado = %f, tiempo real = %f, speedup=%f" %\
            (T, t1 - t0, T/(t1 -t0)) 
    
    # De nuevo vector de estado
    tr1 = modelo.trimado(10.0, 0.0, 0.0, 0.0, 1000., 1.225)
    modelo.t = 0.0
    modelo.x_i = modelo.y_i = 0.0
    modelo.z_i = 1000
    
    t0 = time.time()
    T = 10.
    s2 = modelo.series(T, 0, 0, 0, 0, 0.01, 
            'u', 'v', 'w', 'p', 'q', 'r', 'la0','be0','be1c_w', 'be1s_w',
            'la1c_w', 'la1s_w', 'x_i', 'th', 'fi', 'ch')
    t1 = time.time()

    from pylab import *
    plot(s1['t'], s1['u'])
    plot(s2['t'], s2['u'])
    xlabel('t')
    ylabel('u')
    show()
    
    plot(s1['t'], s1['v'])
    plot(s2['t'], s2['v'])
    xlabel('t')
    ylabel('v')
    show()
    
    plot(s1['t'], s1['w'])
    plot(s2['t'], s2['w'])
    xlabel('t')
    ylabel('w')
    show()
    
    plot(s1['t'], s1['p'])
    plot(s2['t'], s2['p'])
    xlabel('t')
    ylabel('p')
    show()
    
    plot(s1['t'], s1['q'])
    plot(s2['t'], s2['q'])
    xlabel('t')
    ylabel('q')
    show()
    
    plot(s1['t'], s1['r'])
    plot(s2['t'], s2['r'])
    xlabel('t')
    ylabel('r')
    show()
    
    plot(s1['t'], s1['th'])
    plot(s2['t'], s2['th'])
    xlabel('t')
    ylabel(r'$\theta$')
    show()
    
    plot(s1['t'], s1['fi'])
    plot(s2['t'], s2['fi'])
    xlabel('t')
    ylabel(r'$\phi$')
    show()

