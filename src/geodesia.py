# -*- coding: latin-1 -*-
""" 
Definicion del Geoide y utilidades de conversion de coordenadas 
"""

from numarray import *
from math import *
from matematicas import *

class Geoide:
    """ 
    Contiene la definicion del elipsoide y conversion de coordenadas 
    """
    def __init__(self, a=6378137, invf=298.257223563, epsLat=1e-6):
        """ 
        Inicializa por defecto el geoide segun el WGS84 
            a: Semieje mayor 
            invf: el inverso del achatamiento 
            epsLat: error numerico en las rutinas de cambio de cordenadas 
        """
        self.a = a
        self.f = 1/invf
        self.e2 = self.f*(2 - self.f)
        self.epsLat = epsLat
        
        # la ultima tangente de latitud calculada por cart2geod
        self.Tlat = None
        # la ultima inversa de la tangente de latitud
        self.iTlat = None
        
    def geod2cart(self, lat, lon, alt):
        """ 
        Transforma de coordenadas geodesicas a cartesianas.
        Unidades en radianes y metros.
        """
        
        N = self.a/sqrt(1 - self.e2*sin(lat)**2)
        x = (N + alt)*cos(lat)*cos(lon)
        y = (N + alt)*cos(lat)*sin(lon)
        z = (N*(1 - self.e2) + alt)*sin(lat)
        
        return x, y, z
    
    def cart2geod(self, x, y, z):
        """ 
        Transforma de coordenadas cartesianas a geodesicas.
        Unidades en radianes y metros.
        """

        # error en latitud
        errLat = 100
        # minimo numero de iteraciones
        minN = 5
        # maximo numero de iteraciones
        maxN = 100
        # numero de iteraciones
        n = 0
        
        if (x**2 + y**2)<(self.a/2)**2:
            if self.Tlat != None:
                self.iTlat = 1./self.Tlat
                self.Tlat = None
            elif self.iTlat == None:
                self.iTlat = (1 - self.e2)*sqrt(x**2 + y**2)/z
            # Ya disponemos de una primera aproximacion para
            # la inversa de la tangente de la latitud
            iTlat1 = self.iTlat

            # Constantes del bucle
            c0 = sqrt(x**2 + y**2)/z
            c1 = self.e2*self.a/z
            c2 = (self.f-1)**2
            while n<minN or errLat>self.epsLat:
                iTlat2 = c0/(1. + c1/sqrt(iTlat1**2 + c2))
                errLat = abs(iTlat2 - iTlat1)
                iTlat1 = iTlat2
                n += 1

            lat = pi/2 - atan(iTlat2)

            lon = ang(x, y)
            Slat2 = 1./(1. + iTlat2**2)
            N = self.a/sqrt(1 - self.f*(2-self.f)*Slat2)
            alt = z/sqrt(Slat2) - N*(1. - self.e2)
        else:
            if self.iTlat != None:
                self.Tlat = 1./self.iTlat
                self.iTlat = None
            elif self.Tlat == None:
                self.Tlat = z/((1 - self.e2)*sqrt(x**2 + y**2))
        
            # Valor inicial para la tangente de la latitud
            Tlat1 = self.Tlat
            # constantes 
            c0 = 1./sqrt(x**2 + y**2)
            c1 = self.e2*self.a
            c2 = (self.f-1)**2
            while n<minN or errLat>self.epsLat:
                Tlat2 = c0*(z + c1*Tlat1/sqrt(1 + c2*Tlat1**2))
                errLat = abs(Tlat2 - Tlat1)
                Tlat1 = Tlat2
                n += 1

            # latitud
            lat = atan(Tlat2)
            Slat1 = Slat = sin(lat)
            N2 = self.a/sqrt(1 - self.f*(2-self.f)*(Slat)**2)
            # longitud
            lon = ang(x, y)
            # altitud
            alt = sqrt(x**2 + y**2)/sqrt(1 - Slat1**2) - N2

        return lat, lon, alt
    
    def ejesLocales(self, lat, lon, alt):
        """ 
        Coloca los ejes locales con origen en lat, lon, alt. Con el
        eje x apuntando al sur, el eje y al este y el z en la vertical. 
        """

        Slat, Clat = sin(lat), cos(lat)
        Slon, Clon = sin(lon), cos(lon)

        # matriz ejes locales-cartesianos geocentricos
        self.Llc = array( [
            [  Slat*Clon,  Slat*Slon, -Clat ],
            [ -Slon,       Clon,       0    ],
            [  Clat*Clon,  Clat*Slon,  Slat ]] ,type='Float64')
        # matriz ejes cartesianos geocentricos-locales
        self.Lcl = transpose(self.Llc)
        # origen de los ejes locales
        self.Ol = array(self.geod2cart(lat, lon, alt), type='Float64')
    
        
    def loc2cart(self, dx, dy, dz):
        """ 
        Calcula las coordenadas cartesianas a partir de los 
        desplazamientos respecto a los ejes locales. 
        """
        return self.Ol + dot(self.Lcl, [dx, dy, dz])

    def cart2loc(self, x, y, z):
        """
        Calcula los desplazamientos en ejes locales
        a partir de las coordenadas cartesianas.
        """
        return dot(self.Llc, [x,y,z] - self.Ol)

    def loc2geod(self, dx, dy, dz):
        """ 
        Calcula las coordenadas geodesicas a partir de los 
        desplazamientos respecto a los ejes locales. 
        """

        return self.cart2geod(*self.loc2cart(dx, dy, dz))
    
    def geod2loc(self, lat, lon, alt):
        """
        Calcula las coordenadas locales a partir de las geodesicas
        """

        # Calculamos los cartesianos geocentricos
        cart = self.geod2cart(lat, lon, alt)
        
        # Pasamos a locales
        return dot(self.Llc, cart - self.Ol)

################################################################################
#       TESTING                                                                #
################################################################################
if __name__ == '__main__':
    lat0, lon0, alt0 = 0.0, 0.5, 1000
    wgs84 = Geoide()
    x, y, z = wgs84.geod2cart(lat0, lon0, alt0)
    lat1, lon1, alt1 = wgs84.cart2geod(x, y, z)
    error = abs(lat1 - lat0) + abs(lon1 - lon0) + abs(alt1 - alt0)
    print "error = %.12f" % error
    
    
