# -*- coding: latin-1 -*-
"""
Define clases necesarias para el calculo de fuerzas y momentos aerodinámicos 
para todo ángulo de ataque y resbalamiento.
"""

from numarray import *
from matematicas import *

class perfilA:
    """
    Modelo de perfil para todo rango de ángulos de ataque según el 
    NASA Memorandum 84281
    """
    def __init__(self, AR, cLmax=None, a=None, cL0=0, de0=None, de1=None, 
            de2=None, De=0):
        """
        a:     pendiente de sustentación entre 0, pi/4
        cLmax: coeficiente de sustentación máximo
        AR:    alargamiento del ala (Aspect Ratio = b^2/S) 
        De:    flecha del ala
        de0
        de1
        de2:   resistencia del perfil cD = de0 + de1*al + de1*al^2
        """
        
        self.AR = AR
        self.De = De

        # Pendiente de sustentacion, para ángulo de resbalamiento nulo
        # y sin flecha. Si no hay datos se estima
        if a == None:
            self.a = 2*pi/(1 + 2/self.AR)
        else:
            self.a = a
            # se supone que si se da la pendiente del ala ya se ha calculado
            # el efecto de la flecha
            self.De = 0
        
        # Si no se espefica se supone entrada en perdida en pi/4
        if cLmax == None:
            self.cLmax = self.a*pi/4.
        else:
            self.cLmax = cLmax
        
        # Angulo de ataque para sustentacion nula
        self.al0 = -cL0/self.a
        
        # Coeficientes de resistencia, midiendo el angulo de ataque a partir
        # de al0
        A = 0.8*pi*AR
        
        # Angulos pequeños
        if de0 == None:
            self.de0 = 0.009 
        else:
            self.de0 = de0 - self.a**2/A*self.al0**2

        if de1 == None:
            self.de1 = 0.0
        else:
            self.de1 = de1 + 2*self.al0*self.a**2/A

        if de2 == None:
            self.de2 = 0.11 
        else:
            self.de2 = de2 - self.a**2/A
        
        def la(al):
            return self.de0 + self.de1*al + self.de2*al**2
        
        def Dla(al):
            return self.de1 + 2*self.de2*al
            
        # Angulos grandes
        self.ep0 = -0.1254
        self.ep1 =  0.09415
        self.ep2 =  0.977525
        
        def ha(al):
            return self.ep0 + self.ep1*al + self.ep2*(sin(al))**2
        
        def Dha(al):
            return self.ep1 + 2*self.ep2*sin(al)*cos(al)
        
        # Interpolacion intermedia
        
        x1 = -0.60
        x2 = -0.55
        x3 = -0.35
        x4 =  0.35
        x5 =  0.55
        x6 =  0.60

        y5, y6 = map(ha, [x5, x6])
        y1, y2 = y6, y5
        y3, y4 = map(la, [x3, x4])

        D6 = Dha(x6)
        D1 = -D6
        D3, D4 = map(Dla, [x3, x4])

        self.ma1 = InterpoladorSpline(
                x = [ x1, x2, x3 ],
                y = [ y1, y2, y3 ],
                t0 = D1,
                t1 = D3 )
        
        self.ma2 = InterpoladorSpline(
                x = [ x4, x5, x6 ],
                y = [ y4, y5, y6 ],
                t0 = D4,
                t1 = D6 )

        self.la = la
        self.ha = ha


    def coefs(self, al, be):
        """ 
        Calcula (cL, cD) para el angulo de ataque al.
        """

        # Corregimos pendiente de sustentacion con el angulo
        # de resbalamiento y flecha del ala
        a = self.a*(cos(be + self.De))**2
        
        # Angulo de entrada en perdida. Hay que calcularlo aqui porque 
        # varia con el angulo de resbalamiento
        als = self.cLmax/a
        
        # Angulo de ataque auxiliar: de als a al1 se produce el cambio
        # brusco de sustentacion y se interpola linealmente
        al1 = 1.2*als

        # Util para funciones definidas por trozos
        def pick(f, D):
            """ 
            Devuelvo el primer valor de D tal que f aplicada a su 
            llave de verdadero.
            """
            for k,v in D.iteritems():
                if f(k):
                    return v
            return None

        # Se nos puede colar algun error numerico:
        if al<-pi:
            al = -pi
        elif al>pi:
            al = pi
        if be<-pi:
            be = -pi
        elif be>pi:
            be = pi
        
        ali = pick(
        lambda (x, y): x<=al and al<=y,
        {
            ( -pi,  -pi/2 ): lambda x: pi + x,
            ( -pi/2, 0    ): lambda x: -x,
            (  0,    pi/2 ): lambda x: x,
            (  pi/2, pi   ): lambda x: pi - x 
        })(al)
        
        # Damos opcion a que la resistencia no sea simetrica
        aliD = pick(
        lambda (x, y): x<=al and al<=y,
        {
            ( -pi,   -pi/2 ): lambda x: pi + x,
            ( -pi/2, -0.60 ): lambda x: -x,
            ( -0.60,  pi/2 ): lambda x: x,
            (  pi/2,  pi   ): lambda x: pi - x 
        })(al)
        
        cL0 = pick(
        lambda (x, y): x <= ali and ali<= y,
        {
            ( 0,   als  ): lambda x: a*x,
            ( als, al1  ): lambda x: self.cLmax - a*(x - als),
            ( al1, pi/2 ): lambda x:  0.8*self.cLmax*
                                      ( 1 - ((x - al1)/(pi/2 - al1))**2 )
        })(ali)

        cL = pick(
        lambda (x, y): x<=al and al<=y,
        {
            ( -pi,  -pi/2 ):  0.8*cL0,
            ( -pi/2, 0    ): -cL0,
            (  0,    pi/2 ):  cL0,
            (  pi/2, pi   ): -0.8*cL0
        })

        cDp = pick(
        lambda (x, y): x<=aliD and aliD<=y,
        {
            (-0.60, -0.35): self.ma1, 
            (-0.35,  0.35): self.la,
            ( 0.35,  0.60): self.ma2,
            ( 0.60,  pi/2): self.ha
        })(aliD)
        
        cD = cDp + cL**2/(0.8*pi*self.AR)
        
        return (cL, cD)

class perfilB:
    """
    Modelo de perfil para todo rango de ángulos de ataque según el 
    NACA TN 3241
    """

    def __init__(self, a, cLmax, cDmax=1.18, d0=0.09, d2=0.21):
        """
        a: Pendiente de sustentación
        cLmax: máximo coeficiente de sustentación
        cDmax: máximo coeficiente de resistencia
        d0: primer término del coeficiente de resistencia
        d2: segundo término del coeficiente de resistencia
        """
        
        self.a = a
        self.cLmax = cLmax
        self.d0 = d0
        self.d2 = d2
        self.cDmax = cDmax

    def coefs(self, al, be):
        cD = self.d0 + self.d2*(al)**2 
        cL = self.a*sin(al)

        if abs(cL)>self.cLmax:
            cL = self.cLmax*sign(cL)
        
        if abs(cD)>self.cDmax:
            cD = self.cDmax*abs(sin(al))

        return cL, cD

class perfilC:
    """
    Modelo de perfil para todo rango de ángulos de ataque basado en
    tabla para todo ángulo de ataque y resbalamiento
    """
    def __init__(self, cL, cD):
        self.cL = cL
        self.cD = cD

    def coefs(self, al, be):
        return cL(al, be), cD(al, be)


class FuselajeA:
    def __init__(self, cDa90, cDb90, cm90, cl90, cn90, 
            cDa, cDb, cL, cY, cl, cm, cn, 
            al1, al2, be1, be2, unidades = 'rad'):
        """
        Valor de los coeficientes para ángulo de ataque 90º:
            cDa90, cm90
        Valor de los coeficientes para ángulo de resbalamiento 90º:
            cDb90, cl90, cn90

        Coeficientes para pequeños ángulos de ataque:
            cDa, cL, cm
        Definidios entre al1, al2

        Coeficientes para pequeños ángulos de resbalamiento:
            cY, cl, cn
        Definidos entre be1, be2

        """
        
        # Angulos de ataque pequeños
        def la(al):
            if unidades == 'deg':
                al = deg(al)

            return array([
                cDa(al), cL(al), cm(al) ])

        # Angulos de resbalamiento pequeños
        def lb(be):
            if unidades == 'deg':
                be = deg(be)
            return array([
                cDb(be), cY(be), cl(be), cn(be) ])

        # Angulos de ataque grandes
        def ha(al):
            Sal, Cal = sin(al), cos(al)
            
            cDa = cDa90*abs(Sal)*(Sal**2)
            cL  = cDa90*abs(Sal)*Sal*Cal
            cm  = cm90*abs(Sal)*Sal
            
            return array([
                cDa, cL, cm ])
        
        # Angulos de resbalamiento grandes
        def hb(be):
            Sbe, Cbe = sin(be), cos(be)
            
            cDb =  cDb90*abs(Sbe)*(Sbe**2)
            cY  = -cDb90*abs(Sbe)*Sbe*Cbe
            cl  =  cl90*abs(Sbe)*Sbe
            cn  =  cn90*abs(Sbe)*Sbe
            
            return array([
                cDb, cY, cl, cn ])
    
        # Para ángulos de ataque intermedios aproximamos
        # entre los dos casos
        if unidades == 'deg':
            al1, al2 = rad(al1), rad(al2)
            be1, be2 = rad(be1), rad(be2)

        xa0 = al1 - rad(50)
        xa1 = al1 - rad(45)
        xa2 = al1
        xa3 = al2
        xa4 = al2 + rad(45)
        xa5 = al2 + rad(50)

        ya0 = ha(xa0)
        ya1 = ha(xa1) 
        ya2 = la(xa2)
        ya3 = la(xa3)
        ya4 = ha(xa4) 
        ya5 = ha(xa5)

        da0 = deriv_1(ha, xa0)
        da2 = deriv_1(la, xa2)
        da3 = deriv_1(la, xa3)
        da5 = deriv_1(ha, xa5)

        ca0 = []
        for i in range(3):
            ca0.append(
                InterpoladorSpline(
                    x = [xa0, xa1, xa2], 
                    y = [ya0[i], ya1[i], ya2[i]],
                    t0 = da0[i],
                    t1 = da2[i])
            )
        ca1 = []
        for i in range(3):
            ca1.append(
                InterpoladorSpline(
                    x = [xa3, xa4, xa5], 
                    y = [ya3[i], ya4[i], ya5[i]],
                    t0 = da3[i],
                    t1 = da5[i])
            )
        
        def ma0(al):
            return array(
                    [ca0[0](al), ca0[1](al), ca0[2](al)], 
                    type='Float64')

        def ma1(al):
            return array(
                    [ca1[0](al), ca1[1](al), ca1[2](al)], 
                    type='Float64')



        # Interpola linealmente entre (x1, y1) y (x2, y2)
        def lin(x1, y1, x2, y2, x):
            return (y1*(x2 - x) + y2*(x - x1))/(x2 - x1)

        # Coeficientes en funcion del angulo de ataque, 
        # valido entre -pi, pi
        def ca(al):
            if xa2<=al and al<=xa3:
                return la(al)
            elif xa0<=al and al<=xa2:
                #return lin(xa1, ya1, xa2, ya2, al)
                return ma0(al)
            elif xa3<=al and al<=xa5:
                #return lin(xa3, ya3, xa4, ya4, al)
                return ma1(al)
            else: # al<x1 or al>x4
                return ha(al)

        xb0 = be1 - rad(50)
        xb1 = be1 - rad(45)
        xb2 = be1
        xb3 = be2
        xb4 = be2 + rad(45)
        xb5 = be2 + rad(50)
        
        yb0 = hb(xb0)
        yb1 = hb(xb1) 
        yb2 = lb(xb2)
        yb3 = lb(xb3)
        yb4 = hb(xb4) 
        yb5 = hb(xb5)
        
        db0 = deriv_1(hb, xb0)
        db2 = deriv_1(lb, xb2)
        db3 = deriv_1(lb, xb3)
        db5 = deriv_1(hb, xb5)

        cb0 = []
        for i in range(4):
            cb0.append(
                InterpoladorSpline(
                    x = [xb0, xb1, xb2], 
                    y = [yb0[i], yb1[i], yb2[i]],
                    t0 = db0[i],
                    t1 = db2[i])
            )
        cb1 = []
        for i in range(4):
            cb1.append(
                InterpoladorSpline(
                    x = [xb3, xb4, xb5], 
                    y = [yb3[i], yb4[i], yb5[i]],
                    t0 = db3[i],
                    t1 = db5[i])
            )
        
        def mb0(be):
            return array(
                    [cb0[0](be), cb0[1](be), cb0[2](be), cb0[3](be)], 
                    type='Float64')

        def mb1(be):
            return array(
                    [cb1[0](be), cb1[1](be), cb1[2](be), cb1[3](be)], 
                    type='Float64')


        
        # Coeficientes en funcion del angulo de resbalamiento, 
        # valido entre -pi, pi
        def cb(be):
            if xb2<=be and be<=xb3:
                return lb(be)
            elif xb0<=be and be<=xb2:
                return mb0(be)
                #return lin(xb1, yb1, xb2, yb2, be)
            elif xb3<=be and be<=xb5:
                return mb1(be)
                #return lin(xb3, yb3, xb4, yb4, be)
            else:
                return hb(be)
            
        self.call_a = ca
        self.call_b = cb

    def coefs(self, al, be):
        ca = self.call_a(al)
        cb = self.call_b(be)

        cD = ca[0] + cb[0]
        cL = ca[1]
        cY = cb[1]
        cl = cb[2]
        cm = ca[2]
        cn = cb[3]
        return cD, cL, cY, cl, cm, cn 

class FuselajeB:
    def __init__(self, f):
        """ 
        Wrapper alrededor de la funcion f, que devuelve el valor de
        los coeficientes aerodinamicos para todo angulo de ataque y
        resbalamiento:
        cD, cL, cY, cl, cm, cn = f(al, be)
        """

        self.f = f
    
    def coefs(self, al, be):
        return self.f(al, be)
