#!/usr/bin/env python
from Lynx import modelo

import pylab
from matematicas import *
pylab.rc('text', usetex=True)

modelo.init() 

V = 0
be = 0
ga = 0
Om = 0
ro = 1.225

# pylab nos crea listas de tipo float64scalar, y yo quiero
# float a secas

def arange(x0, x1, st=1.):
    r = []
    v = x0
    while (x0<x1 and v<=x1 - st) or (x0>x1 and v>=x1 + st): # cosas raras del error numerico :)
        r.append(v)
        v += st
    return r


sel = raw_input("Seleccion figuras: ")

if "1" in sel:
    z = arange(0, 10, 0.1)
    V = arange(0, 30, 5)
    for y in V:
        th0 = []
        for x in z:
            tr = modelo.trimado(y, ga, be, Om, x, ro)
            th0.append(tr['th0'])

        pylab.plot(z, th0)


    pylab.xlabel(r'$z$')
    pylab.ylabel(r'$\theta_0$')

    pylab.show()

if "2" in sel:
    z = arange(0, 10, 0.1)
    V = 0
    th0, th1c, th1s, th0T = [], [], [], [] 
    for x in z:
        tr = modelo.trimado(V, ga, be, Om, x, ro)
        th0.append(tr['th0'])
        th1c.append(tr['th1c'])
        th1s.append(tr['th1s'])
        th0T.append(tr['th0T'])
    
    pylab.plot(z, th0)
    pylab.xlabel(r'$z$')
    pylab.ylabel(r'$\theta_0$')
    pylab.show()

    pylab.plot(th0, th0T)
    pylab.xlabel(r'$\theta_0$')
    pylab.ylabel(r'$\theta_{0_T}$')
    pylab.show()

    pylab.plot(th0, th1c)
    pylab.plot(th0, th1s)
    pylab.xlabel(r'$\theta_0$')
    pylab.show()

if "3" in sel:
    z = 1000
    V = arange(-30, 30, 1)
    ga = pylab.pi/2.
    th0, th0T = [], [] 
    cT, la0, cTT = [], [], []
    for x in V:
        tr = modelo.trimado(-x, ga, be, Om, z, ro)
        th0.append(tr['th0'])
        th0T.append(tr['th0T'])
        la0.append(tr['la0'])
        cT.append(tr['cT'])
        cTT.append(tr['cTT'])
    
    pylab.plot(V, la0)
    pylab.ylabel(r'$\lambda_0$')
    pylab.xlabel(r'$V_a$')
    pylab.show()
    
    pylab.plot(V, th0)
    pylab.xlabel(r'$V_a$')
    pylab.ylabel(r'$\theta_0$')
    pylab.show()
    
    pylab.plot(th0, th0T)
    pylab.xlabel(r'$\theta_0$')
    pylab.ylabel(r'$\theta_{0_T}$')
    pylab.show()
    
if "4" in sel:
    z = 1000
    V = arange(0, 165, 5)
    ga = 0
    be = 0 
    Om = 0
    ro = 1.04760
    
    resultados = {}
    for x in V:
        metodo = 'Fijo'

        tr = modelo.trimado(x*0.514444, ga, be, Om, z, ro, metodo)
        for k,v in tr.items():
            if not k in resultados:
                resultados[k] = [v]
            else:
                resultados[k].append(v)

    pylab.plot(V, resultados['th0'])
    pylab.ylabel(r'$\theta_0$')
    pylab.xlabel(r'$V$')
    pylab.show()
    
    pylab.plot(V, resultados['la0'])
    pylab.plot(V, resultados['la1c_w'])
    pylab.xlabel(r'$V$')
    pylab.show()
    
    pylab.plot(V, map(deg, resultados['th']))
    pylab.ylabel(r'$\theta$')
    pylab.xlabel(r'$V$')
    pylab.show()
    
    pylab.plot(V, resultados['cTT'])
    pylab.plot(V, resultados['fi'])
    pylab.xlabel(r'$V$')
    pylab.show()
    
    pylab.plot(V, resultados['th0T'])
    pylab.xlabel(r'$V$')
    pylab.ylabel(r'$\theta_{0_T}$')
    pylab.show()

    pylab.plot(V, resultados['th1c'])
    pylab.plot(V, resultados['th1s'])
    pylab.xlabel(r'$V$')
    pylab.show()
    
    def per(x, xmin, xmax):
        return (x - xmin)/(xmax - xmin)*100

#   pylab.plot(V, [ per(x, 0, 0.254) for x in resultados['et0p'] ])
#   pylab.ylabel(r'$\eta_{0p}$')
#   pylab.xlabel(r'$V$')
#   pylab.axis([0, 200, 0, 100])
#   pylab.show()
    
#   pylab.plot(V,  [ per(x, -0.127, 0.127) for x in resultados['et1sp'] ]) 
#   pylab.ylabel(r'$\eta_{1sp}$')
#   pylab.xlabel(r'$V$')
#   pylab.axis([0, 200, 0, 100])
#   pylab.show()
#   
#   pylab.plot(V, [ per(x, -0.127, 0.127) for x in resultados['et1cp'] ])
#   pylab.ylabel(r'$\eta_{1cp}$')
#   pylab.xlabel(r'$V$')
#   pylab.axis([0, 200, 0, 100])
#   pylab.show()
    
#   pylab.plot(V, [ per(x, -0.0683, 0.0683) for x in resultados['etpp'] ])
#   pylab.ylabel(r'$\eta_{pp}$')
#   pylab.xlabel(r'$V$')
#   pylab.axis([0, 200, 0, 100])
#   pylab.show()
    
#   pylab.plot(V, resultados['altp'] )
#   pylab.ylabel(r'$\alpha_{tp}$')
#   pylab.xlabel(r'$V$')
#   pylab.show()
#   
#   pylab.plot(V, resultados['cLtp'] )
#   pylab.ylabel(r'$c_{L_{tp}}$')
#   pylab.xlabel(r'$V$')
#   pylab.show()
#   
#   pylab.plot(V, resultados['cDtp'] )
#   pylab.ylabel(r'$c_{D_{tp}}$')
#   pylab.xlabel(r'$V$')
#   pylab.show()
#   
#   pylab.plot(V, resultados['alf'] )
#   pylab.ylabel(r'$\alpha_f$')
#   pylab.xlabel(r'$V$')
#   pylab.show()
#   
#   pylab.plot(V, resultados['cLf'] )
#   pylab.ylabel(r'$c_{L_f}$')
#   pylab.xlabel(r'$V$')
#   pylab.show()
#   
#   pylab.plot(V, resultados['cDf'] )
#   pylab.ylabel(r'$c_{D_f}$')
#   pylab.xlabel(r'$V$')
#   pylab.show()
#   
#   pylab.plot(V, resultados['cmf'] )
#   pylab.ylabel(r'$c_{m_f}$')
#   pylab.xlabel(r'$V$')
#   pylab.show()
#   
#   #pylab.plot(V, resultados['be1c_w'])
#   #pylab.plot(V, resultados['be1s_w'])
#   #pylab.plot(V, resultados['fi'])
#   #pylab.ylabel(r'$\beta_{1cw}$')
#   #pylab.xlabel(r'$V$')
#   #pylab.show()
#   
#   #pylab.plot(V, resultados['th0'])
#   #pylab.xlabel(r'$V$')
#   #pylab.ylabel(r'$\theta_0$')
#   #pylab.show()
#   

if "5" in sel:
    z = 1000
    V = 30
    ga = 0
    be = 0 
    Om = arange(-0.45, 0.55, 0.05)
    ro = 1.225
    
    resultados = {}
    for x in Om:
        metodo = 'Fijo'

        tr = modelo.trimado(V, ga, be, x, z, ro, metodo)
        for k,v in tr.items():
            if not k in resultados:
                resultados[k] = [v]
            else:
                resultados[k].append(v)

    pylab.plot(Om, resultados['fi'])
    pylab.ylabel(r'$\phi$')
    pylab.xlabel(r'$\Omega_a$')
    pylab.show()
    
    pylab.plot(Om, resultados['th'])
    pylab.ylabel(r'$\theta$')
    pylab.xlabel(r'$\Omega_a$')
    pylab.show()
    
    pylab.plot(Om, resultados['th0'])
    pylab.ylabel(r'$\theta_0$')
    pylab.xlabel(r'$\Omega_a$')
    pylab.show()
    
    pylab.plot(Om, resultados['th0T'])
    pylab.ylabel(r'$\theta_{0_T}$')
    pylab.xlabel(r'$\Omega_a$')
    pylab.show()
    
    pylab.plot(Om, resultados['th1c'], label=r'$\theta_{1c}$')
    pylab.plot(Om, resultados['th1s'], label=r'$\theta_{1s}$')
    pylab.xlabel(r'$\Omega_a$')
    pylab.legend()
    pylab.show()
    

