#!/usr/bin/env python
# -*- coding:latin 1 -*-
""" 
Contiene el bucle principal 
"""

from flightgear import FlightGear
from optparse import OptionParser
import sys
import imp

# Parseamos las opciones de la linea de comandos
parser = OptionParser(
        usage="%prog --modelo=M, [--logfile=F, --logprop=a, [b,c...]] " +
        "[--host=H], [--puertos=p0,p1,p2]", 
        prog="icaro", 
        version="1.0")

parser.add_option("--modelo",  
        help = "Nombre del archivo que contiene el modelo, sin la " + 
               "extension .py", 
        type="string", 
        dest="modelo")

parser.add_option("--logfile", 
        help = "Fichero opcional de salida de datos para logging",
        type="string", 
        dest="logfile")

parser.add_option("--logprop", 
        help = "Propiedades a escribir en el fichero de logging. \n" +
               "Consultar documentacion para ver la lista de posibilidades",
        type="string", 
        dest="logprop")

parser.add_option("--host", 
        help = "Direccion de la maquina ejecutando FlightGear, si no se " + 
        "especifica se asume localhost",
        type="string", 
        dest="host")

parser.add_option("--puertos", 
        help = "Puertos de conexion con FlightGear, respectivamente son " +
        "de controles, del FDM y de comandos. Por defecto son los " + 
        "puertos 5001, 5002 y 5003 respectivamente",
        type="string", 
        dest="puertos")

opt, args = parser.parse_args()
if opt.modelo == None:
    parser.error("--modelo requerido")

try:
	(file_modelo, path_modelo, desc_modelo) = imp.find_module(opt.modelo)
except ImportError:
	sys.exit("Error: no se encontro el modelo")

modulo = imp.load_module(opt.modelo, file_modelo, path_modelo, desc_modelo)

if opt.logprop!=None:
	logprop = opt.logprop.split(",")
else:
	logprop = []

if opt.host!=None:
    host = opt.host
else:
    host = "localhost"

if opt.puertos!=None:
    try:
        p1, p2, p3 = map(float, opt.puertos.split(","))
    except:
        parser.error("No se especificaron correctamente los puertos")
else:
    p1 = 5001
    p2 = 5002
    p3 = 5003

# Creamos el vuelo
vuelo = FlightGear(modulo.modelo, logfile=opt.logfile, props=logprop,
        host=host, port1=p1, port2=p2, port3=p3)

# Entramos en el bucle de la simulacion
try:	
	while True:
		vuelo.loop()
except KeyboardInterrupt:
	print "# Deteniendo la simulacion..."
	vuelo.close()
	file_modelo.close()

