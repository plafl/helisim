# -*- coding: latin-1 -*-
""" 
Protocol nativo tal como viene especificado en el codigo de FlightGear 

Para mas detalles el protocolo viene descrito en:
    FlightGear/src/Newtork/net_ctrls.hxx
    FlightGear/src/Newtork/net_fdm.hxx

El protocolo de comandos aunque se encuentra en la clase HTTPClient definida e 
implementada en:
    FlightGear/src/FDM/ExternalNet/ExternalNet.hxx
"""

import struct
import re

# formato de de llegada de los controles
ctrlFmt = ">I 4x 9d 2I I 16I 4x 12d 4I 4d 8I 24I I 8I 5I I 4x 5d I I 4d 3d 2d \
2d I 2I 100x" 

# tamano de llegada de controles
ctrlBufSize = struct.calcsize(ctrlFmt)

# formato de salida del estado del FDM
fdmFmt = ">2I 3d 6f 11f 3f 2f 5I 36f I 4f 4I 9f I i f 10f"

# tamano de salida
fdmBufSize = struct.calcsize(fdmFmt)

# formato de linea de comandos
reCmd = re.compile('GET /(?P<prop>\w+-?\w*)\?value=(?P<value>.*) HTTP/1.0')

# hasta que haya generadores de python con parametros de llamada utilizamos esta
# clase

class SequenceGet:
    """ 
    Obtiene los siguientes n elementos de una secuencia 
    """
    
    def __init__(self, tuple):
        self.t = tuple
        self.c = 0
    
    def reset(self):
        self.c = 0
    
    def __call__(self, n):
        if n == 1:
            r = self.t[self.c]
        else:
            r = self.t[self.c:self.c + n]
        self.c += n
        return r
    
def ctrlUnpack(str):
    """ 
    Convierte la cadena de llegada de controles en un diccionario 
    """
    
    # desenpaqueta la cadena en una tupla
    t = struct.unpack(ctrlFmt, str)
    
    # transformamos la tupla en un diccionario
    get = SequenceGet(t)
    
    return {
        'version':              get(1),
        
        # control aerodinamico
        'aileron':              get(1),
        'elevator':             get(1),
        'rudder':               get(1),
        'aileron_trim':         get(1),
        'elevator_trim':        get(1),
        'rudder_trim':          get(1),
        'flaps':                get(1),
        'spoilers':             get(1),
        'speedbrake':           get(1),
        
        # fallos de controles aerodinamicos   
        'flaps_power':          get(1),
        'flap_motor_ok':        get(1),
        
        # controles de motor
        'num_engines':          get(1),
        'master_bat':           get(4),
        'master_alt':           get(4),
        'magnetos':             get(4),
        'starter_power':        get(4),
        'throttle':             get(4),
        'mixture':              get(4),
        'condition':            get(4),
        'fuel_pump_power':      get(4),
        'prop_advance':         get(4),
        'feed_tank':            get(4),
        'reverse':              get(4),
        
        # fallos de motor
        'engine_ok':            get(4),
        'mag_left_ok':          get(4),
        'mag_right_ok':         get(4),
        'spark_plugs_ok':       get(4),
        'oil_press_status':     get(4),
        'fuel_pump_ok':         get(4),

        # combustible   
        'num_tanks':            get(1),
        'fuel_selector':        get(8),
        'xfer_pump':            get(5),

        'cross_feed':           get(1),
        
        # frenos
        'brake_left':           get(1),
        'brake_right':          get(1),
        'copilot_brake_left':   get(1),
        'copilot_brake_right':  get(1),
        'brake_parking':        get(1),

        # tren de aterrizaje
        'gear_handle':          get(1),
        
        # interruptores
        'master_avionics':      get(1),

        # navegacion y comunicaciones
        'comm_1':               get(1),
        'comm_2':               get(1),
        'nav_1':                get(1),
        'nav_2':                get(1),
        
        # viento y turbulencia
        'wind_speed':           get(1),
        'wind_dir_deg':         get(1),
        'turbulence_norm':      get(1),

        # temperatura y presion
        'temp_c':               get(1),
        'press_inhg':           get(1),
    
        # otra informacion de entorno
        'hground':              get(1),
        'magvar':               get(1),

        # hielo
        'icing':                get(1),

        # control de la simulacion
        'speedup':              get(1),
        'freeze':               get(1)
    }


# propiedades del FDM en el orden correcto para enpaquetar
fdmProps = [
        'version',
        'padding',
        
        # Posiciones
        'longitude', 
        'latitude', 
        'altitude', 
        'agl', 
        'phi', 
        'theta', 
        'psi', 
        'alpha', 
        'beta', 
        
        # Velocidades
        'phidot', 
        'thetadot', 
        'psidot', 
        'vcas', 
        'climb_rate', 
        'v_north',
        'v_east', 
        'v_down', 
        'v_wind_body_north', 
        'v_wind_body_east',
        'v_wind_body_down',

        # Aceleraciones
        'A_X_pilot', 
        'A_Y_pilot', 
        'A_Z_pilot',

        # Entrada en perdida
        'stall_warning', 
        'slip_deg',

        # Estatus de motor
        'num_engines', 
        'eng_state', 
        'rpm', 
        'fuel_flow', 
        'fuel_px',
        'egt', 
        'cht', 
        'mp_osi', 
        'tit', 
        'oil_temp', 
        'oil_px', 

        # Consumibles
        'num_tanks', 
        'fuel_quantity',

        # Estado del tren de aterrizaje
        'num_wheels', 
        'wow', 
        'gear_pos', 
        'gear_steer',
        'gear_compression',

        # Entorno
        'cur_time', 
        'warp', 
        'visibility',

        # Posicion de las superficies de control
        'elevator', 
        'elevator_trim_tab', 
        'left_flap', 
        'right_flap',
        'left_aileron', 
        'right_aileron', 
        'rudder', 
        'nose_wheel',
        'speedbrake', 
        'spoilers' 
]


def fdmPack(dict):
    """ 
    Transforma el diccionario de propiedades en una cadena lista para
    transmitir 
    """

    # Chequea si la variable es un diccionario o lista
    # FIX: deberiamos chequear si implementa un protocolo de contenedor
    def isContainer(con):
        if type(con) in map(type, [[],{}]):
            return True
        else:
            return False
    
    t = []
    for prop in fdmProps:
        val = dict[prop]
        if isContainer(val):
            t.extend(val)
        else:
            t.append(val)
    
    return struct.pack(fdmFmt, *t)

