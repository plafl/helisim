#!/usr/bin/env python

import sys
import pickle
from optparse import OptionParser
from matematicas import InterpoladorLineal

###############################################################################
# PARSEAMOS OPCIONES                                                          #
###############################################################################
parser = OptionParser(
        prog="serializa", 
        version="1.0",
        usage="%prog -i entrada -o salida")

parser.add_option("-i",  type="string", dest="input")
parser.add_option("-o",  type="string", dest="output")

opt, args = parser.parse_args()

# Chequeamos que tenemos todos los datos y son validos
if opt.input == None:
    parser.error("Se requiere fichero de entrada")

if opt.output == None:
    parser.error("Se requiere fichero de salida")

try:
    fi = open(opt.input, 'r')
except IOError, msg:
    sys.exit("No se pudo abrir fichero de entrada: %s" % msg)

try:
    fo = open(opt.output, 'w')
except IOError, msg:
    sys.exit("No se pudo abrir fichero de salida: %s" % msg)

###############################################################################
# LEEMOS FICHERO DE ENTRADA                                                   #
###############################################################################
li = fi.readline()
dat = li.split(" ")
# Cada columna del fichero es una lista contenida en D, que tiene
# tantas listas como columnas el fichero, evidentemente
# El ultimo elemento es el caracter de nueva linea, lo eliminamos
D = [[] for x in dat[:-1]]

while li != "":
    dat = li.split(" ")
    for i in range(len(dat)-1):
        try:
            D[i].append(float(dat[i]))
        except ValueError:
            fi.close()
            fo.close()
            sys.exit("Error en fichero de entrada, no se pudo convertir " + 
                     "a valor numerico")
    li = fi.readline()
fi.close()

###############################################################################
# ESCRIBIMOS FICHERO DE SALIDA                                                #
###############################################################################
t = D[0]
for x in D[1:]:
    a = InterpoladorLineal(t, x)
    pickle.dump(a, fo)
fo.close()
