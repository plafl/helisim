# -*- coding: latin-1 -*-
""" 
Gestion de la conexion externa con FlightGear 
"""

import logging
import socket
import time
import os

import protocolo
from historial import *
from matematicas import *
from geodesia import Geoide


def ft2m(feet):
    """
    Convierte de pies a metros.
    """

    return 0.3048*feet

class FlightGear:
    def __init__(self, modelo, logfile='vuelo.dat', props=[], host='localhost',
            port1=5001, port2=5002, port3=5003):
        """ 
        Gestiona la conexion con FlightGear. 
            modelo: helicoptero que se quiera, una instancia de clase modelo.
            logfile: archivo donde se guarden las variables de logging
            props: lista con las variables a guardar
            host: direccion del ordenador que corre FlightGear.
            port1: puerto de controles.
            port2: puerto de FDM.
            port3: puerto de comandos.
        """

        # LOGGING DE ERRORES
        self.logger = logging.getLogger('FlightGear')
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("# %(levelname)s: %(message)s"))
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.INFO)
        
        # CARGAMOS EL MODELO
        self.logger.info('Leyendo modelo de helicoptero %s' % modelo.nombre)
        self.modelo = modelo
        self.modelo.init()

        # COMUNICACIONES
        self.logger.info('Iniciando conexion de entrada de controles')
        self.sockCtrls = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.sockFDM = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.sockCmd = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.cmds_recibidos = False
        
        try:
            self.sockCtrls.setblocking(False)
            self.sockCtrls.bind((host, port1))
        except socket.error, msg:
            self.logger.warning("No se pudo abrir socket de \
            entrada de controles: %s" % msg)
        try:
            self.sockFDM.setblocking(False)
            self.sockFDM.connect((host, port2))
        except socket.error, msg:
            self.logger.warning("No se pudo abrir socket de salida\
            de estado de FDM: %s" % msg)
        try:
            self.sockCmd.setblocking(False)
            self.sockCmd.bind((host, port3))
            self.sockCmd.listen(1)
            # Todavia hay que aceptar las conexiones
        except socket.error, msg:
            self.logger.warning("No se pudo abrir socket de \
            comandos: %s" % msg) 
    
        self.props = Historial(10000)
        self.props.append({})
        self.cmds = Historial(10000)
        self.cmds.append({})

        # LOGGING
        self.logProps = props
        if logfile != None:
            self.logger.info('Abriendo archivo log: %s' % logfile)
            try:
                self.archivoLog = open(logfile, 'w')
            except:
                self.archivoLog = None
                self.logger.warning('No se pudo abrir archivo para logging')
        else:
            self.archivoLog = None
        self.logger.info('Entrando en el bucle de la simulacion. ')
        self.logger.info('Presione Ctrl-C para cancelar la simulacion. ')
    
    def close(self):
        """ 
        Cierra todos los sockets y archivos. 
        """
        
        self.logger.info("Cerrando conexiones...")
        self.sockCtrls.close()
        self.sockFDM.close()
        self.sockCmd.close()
        
        self.logger.info("Cerrando archivos...")
        if self.archivoLog != None:
            self.archivoLog.close()

    def sendFDM(self):
        """ 
        Manda el diccionario de las propiedades. 
        """

        try:
            self.sockFDM.send(protocolo.fdmPack(self.props[0]['fdm']))
        except socket.error, value:
            self.logger.warning("No se pudo mandar estado FDM: %s" % value)
    
    def recvCtrls(self):
        """ 
        Recibe el diccionario de controles. 
        """

        try:
            str = self.sockCtrls.recv(protocolo.ctrlBufSize)
            self.props[0]['ctrls'] = protocolo.ctrlUnpack(str)
        except:
            raise
    
    def recvCmd(self):
        """ 
        Recibe los comandos de inicializacion. Cuando termina de
        recibir todos ajusta los valores en el modelo. 
        """
        
        try:
            conn, addr = self.sockCmd.accept()
            conn.setblocking(True)
        except:
            raise
        else:
            msg = ""
            chunk = conn.recv(1)
            while chunk != "":
                msg += chunk
                chunk = conn.recv(1)
            match = protocolo.reCmd.match(msg)
            # Debemos responder a la inicializacion de las siguientes 
            # propiedades (por orden):
            # longitude-deg, latitude-deg, altitude-ft, ground-m, speed-kts,
            # heading-deg, reset(=ground)
            if match != None:
                def longitud(strVal):
                    self.cmds[0]['longitud'] = rad(float(strVal))
                
                def latitud(strVal):
                    self.cmds[0]['latitud'] = rad(float(strVal))

                def altitud(strVal):
                    self.cmds[0]['altitud'] = 0.0
                    # Por algun motivo FG no nos pasa la altitud correcta
                    # asi que de momento esto lo dejamos en paz, mas adelante
                    # se corrige haciendo altitud = altitud del suelo

                def altura(strVal):
                    # Realmente quiere decir altitud del suelo
                    self.cmds[0]['altura'] = float(strVal)

                def orientacion(strVal):
                    self.cmds[0]['azimuth'] = rad(float(strVal))

                # Todos los valores recibidios
                def reset(strVal):
                    self.logger.info("Posicion inicial recibida:")
                    self.logger.info(
                            "\tlongitud = %f rad" % self.cmds[0]['longitud'])
                    self.logger.info(
                            "\tlatitud  = %f rad" % self.cmds[0]['latitud'] )
                    self.logger.info(
                            "\taltitud  = %f m"   % self.cmds[0]['altura'] )
                    self.logger.info(
                            "\tazimuth  = %f rad" % self.cmds[0]['azimuth'] )
                    self.logger.info(
                            "\taltura del suelo  = %f m" % 
                            self.cmds[0]['altura'])
                    
                    self.cmds_recibidos = True
                    
                prop, val = match.group('prop'), match.group('value')
                try:
                    {
                    'longitude-deg':    longitud, 
                    'latitude-deg':     latitud,
                    'altitude-ft':      altitud,
                    'ground-m':         altura,
                    'heading-deg':      orientacion,
                    'reset':            reset,
                    # No se encuentra implementado
                    'speed-kts':        lambda x: None,
                    }[prop](val)
                except KeyError:
                    self.logger.warning(
                            "Comando inesperado: %s = %s" % (prop, val))
            else:
                self.logger.warning("Error de protocolo")
    

    def loop(self):
        """ 
        Recibe los controles, calcula la respuesta y manda el
        resultado. 
        """
        self.props.append({'ctrls': None, 'fdm':None, 'time': time.time()})
        
        if not self.cmds_recibidos:
            try:    
                self.recvCmd()
            except socket.error:    
                pass
        
        try:
            self.recvCtrls()
        except socket.error:
            pass
        else:
            ####################################################################
            #   CON LA INFORMACION QUE NOS LLEGA ACTUALIZAMOS EL HELICOPTERO   #
            ####################################################################
            
            # CONTROL DE TIEMPOS
            # si no hemos parado la simulacion nos aseguramos de que el reloj 
            # esta corriendo, si ya esta corriendo la siguiente llamada no 
            # afecta el funcionamiento
            if self.props[0]['ctrls']['freeze'] == 0:
                self.modelo.reloj.arranca()
            else:
                self.modelo.reloj.para()
                return 
            t = self.modelo.reloj()
            # hemos recibido los controles y el reloj esta corriendo, ajustamos
            # el modelo. 
            self.modelo.controles.et0p.append(
                    1-self.props[0]['ctrls']['throttle'][1], t)
            self.modelo.controles.et1cp.append(
                    -self.props[0]['ctrls']['aileron'], t)  
            self.modelo.controles.et1sp.append(
                    -self.props[0]['ctrls']['elevator'], t) 
            self.modelo.controles.etpp.append(
                    -self.props[0]['ctrls']['rudder'], t)

            # atmosfera
            # Parece un bug de FG, tiramos con atmosfera ISA
#            self.modelo.T = self.props[0]['ctrls']['temp_c'] + 273.15
#            self.modelo.P = self.props[0]['ctrls']['press_inhg']*3369.39

            
            self.modelo.viento_vel = self.props[0]['ctrls']['wind_speed']*0.814444
            self.modelo.viento_dir = rad(self.props[0]['ctrls']['wind_dir_deg'])

            # Encendemos o apagamos el motor
            if self.props[0]['ctrls']['magnetos'][0] > 0:
                self.modelo.motor.start = True
                self.modelo.motor.nmotores = 2
            else:
                self.modelo.motor.start = False
                self.modelo.motor.nmotores = 0
            
            # Altura del suelo sobre el geoide de referencia
            self.modelo.hG = self.props[0]['ctrls']['hground']

            if self.cmds_recibidos:
                # Vector de estado, condicion inicial
                self.modelo.t   = 0.0
                self.modelo.hG  = self.cmds[0]['altura']
                    
                # Preparamos el integrador
                # Vector de estado
                u = 0.0
                v = 0.0
                w = 0.0
                    
                p = 0.0
                q = 0.0
                r = 0.0
                    
                ( q0, q1, q2, q3 ) = Quaternion.rotacion([
                    (pi - self.cmds[0]['azimuth'], 0, 0, 1),
                    (pi, 1, 0, 0) ])
                    
                Om = 0.0
                DTQ = 0.0
                Q = 0.0
                    
                x_i = 0.0
                y_i = 0.0
                z_i = self.modelo.hG + self.modelo.fuselaje.ht  
                    
                # Motores apagados
                self.modelo.motor.nmotores = 0

                self.modelo.tren.reset()
                 
                self.modelo.integrador.init(
                    F = self.modelo.paso,
                    ti = 0,
                    xi = array(
                        [ u, v, w, 
                          p, q, r, 
                          q0, q1, q2, q3, 
                          Om, DTQ, Q, 
                          x_i, y_i, z_i ]))
                    
                # Posicion inicial y ejes inerciales
                self.modelo.geoide.ejesLocales(
                    self.cmds[0]['latitud'], 
                    self.cmds[0]['longitud'],
                    0.0
                )
                # Reseteamos cmds
                self.cmds.append({})
                self.cmds_recibidos = False

            
            ####################################################################
            #   CALCULAMOS EL RESULTADO                                        #
            ####################################################################          
            # integramos las ecuaciones hasta el instante actual
            self.modelo.avanza(self.modelo.reloj())

            ####################################################################
            #   EJECUTAMOS EL LOGGING DE PROPIEDADES                           #
            ####################################################################
            if self.archivoLog != None:
                if self.logProps != []:
                    str = "%.4f "
                    val = [self.modelo.t]
                for P in self.logProps:
                    str += "%.4f "
                    try:
                        v = self.modelo.get(P)
                    except KeyError:
                        v = 0.0
                    val.append(v)
                   
                str += "\n"
                self.archivoLog.write(str % tuple(val))
                # Si no hacemos un flush no funciona para graficos interactivos
                self.archivoLog.flush()
            
            ####################################################################
            #   MANDAMOS LA INFORMACION CALCULADA                              #
            ####################################################################
            
            # preparamos la informacion del modelo.
            # asignemos o no valores hay que rellenar correctamente el 
            # diccionario
            fdm = {
                'version':      24,
                'padding':      0,
                # Posiciones
                # Coordenadas geodesicas en radianes
                'longitude':            0.0,    
                'latitude':             0.0,
                'altitude':             0.0,
                'agl':                  1.0,    # Above Ground Level
                'phi':                  0.0,    # Euler en radianes
                'theta':                0.0,
                'psi':                  0.0,
                'alpha':                1.0,    # Ataque (rad)
                'beta':                 2.0,    # Resbalamiento (rad)
                # Velocidades
                'phidot':               3.0,
                'thetadot':             4.0,
                'psidot':               5.0,
                'vcas':                 6.0, 
                'climb_rate':           7.0,
                'v_north':              8.0,
                'v_east':               9.0,
                'v_down':               10.0, 
                'v_wind_body_north':    11.0,
                'v_wind_body_east':     12.0,
                'v_wind_body_down':     13.0,
                # Aceleraciones
                'A_X_pilot':            14.0,
                'A_Y_pilot':            15.0,
                'A_Z_pilot':            16.0,
                # Perdida
                'stall_warning':        0.0,
                'slip_deg':             0.0,
                # Motores
                'num_engines':          2,
                'eng_state':            [0, 0, 0, 0],
                'rpm':                  [0.0, 0.0, 0.0, 0.0],
                'fuel_flow':            [0.0, 0.0, 0.0, 0.0],
                'fuel_px':              [0.0, 0.0, 0.0, 0.0],
                'egt':                  [0.0, 0.0, 0.0, 0.0],
                'cht':                  [0.0, 0.0, 0.0, 0.0],
                'mp_osi':               [0.0, 0.0, 0.0, 0.0],
                'tit':                  [0.0, 0.0, 0.0, 0.0],
                'oil_temp':             [0.0, 0.0, 0.0, 0.0],
                'oil_px':               [0.0, 0.0, 0.0, 0.0],
                # Consumibles
                'num_tanks':            0,
                'fuel_quantity':        [0.0, 0.0, 0.0, 0.0],
                # Tren de aterrizaje
                'num_wheels':           0,
                'wow':                  [0, 0, 0],
                'gear_pos':             [0.0, 0.0, 0.0],
                'gear_steer':           [0.0, 0.0, 0.0],
                'gear_compression':     [0.0, 0.0, 0.0],
                # Entorno
                'cur_time':             0,
                'warp':                 0,
                'visibility':           1000.0,
                # Controles
                'elevator':             0.0,
                'elevator_trim_tab':    0.0,
                'left_flap':            0.0,
                'right_flap':           0.0,
                'left_aileron':         0.0,
                'right_aileron':        0.0,
                'rudder':               0.0,
                'nose_wheel':           0.0,
                'speedbrake':           0.0,
                'spoilers':             0.0
            }
            # version
            fdm['version'] = 24
            
            fdm['latitude']  = self.modelo.lat
            fdm['longitude'] = self.modelo.lon
            fdm['altitude']  = self.modelo.alt
            
            fdm['psi']   = self.modelo.ch
            fdm['theta'] = self.modelo.th
            fdm['phi']   = self.modelo.fi

            fdm['agl'] = self.modelo.alt - self.modelo.hG
           
            fdm['phidot']   = self.modelo.DTfi
            fdm['thetadot'] = self.modelo.DTth
            fdm['psidot']   = self.modelo.DTch

            fdm['A_X_pilot'] = self.modelo.ecs_modelo[0][0]
            fdm['A_Y_pilot'] = self.modelo.ecs_modelo[0][1]
            fdm['A_Z_pilot'] = self.modelo.ecs_modelo[0][2]

            fdm['climb_rate'] = self.modelo.ecs_modelo[0][15]/60./0.3048
            
            fdm['beta']  = self.modelo.fuselaje.be
            fdm['alpha'] = self.modelo.fuselaje.al

            fdm['vcas'] = sqrt(self.modelo.u**2 + self.modelo.v**2 +\
                    self.modelo.w**2)/0.3048
            fdm['rpm'][0] = self.modelo.Om
            # Guardamos resultado
            self.props[0]['fdm'] = fdm
            
            # mandamos la informacion del modelo
            self.sendFDM()
        
