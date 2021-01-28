from matematicas import *
from numarray import *

class Fuselaje:
    def __init__(self, modelo):
        self.modelo = modelo

    def init(self):
        pass
    
    def calcula_FMf(self):
        """
        Calcula las fuerzas aerodinamicas sobre el fuselaje, en ejes viento.
        
        NECESITA:
            modelo.Vrw_b
            
        PROVEE:
            fuselaje.FMf
        """

        # Un alias para le velocidad relativa al viento en ejes cuerpo
        urw_b, vrw_b, wrw_b = self.modelo.Vrw_b
        # Velocidad realativa al viento, al cuadrado, en ejes cuerpo
        Vf2rw_b = urw_b**2 + vrw_b**2 + wrw_b**2
        # Calculamos angulo de resbalamiento (positivo si el viento incide 
        # por la derecha)
        be = ang(sqrt(urw_b**2 + wrw_b**2), vrw_b)
        # Calculamos angulo de ataque
        al = ang(urw_b, wrw_b)
        # Presion dinamica
        pd = 1./2.*self.modelo.ro*Vf2rw_b
        # Calculamos las fuerzas interpolando los coeficientes
        cD, cL, cY, cl, cm, cn = self.aero.coefs(al, be)
        
        Df_w = pd*self.Sp*cD
        Yf_w = pd*self.Ss*cY
        Lf_w = pd*self.Sp*cL
        lf_w = pd*self.Ss*self.lf*cl
        mf_w = pd*self.Sp*self.lf*cm
        nf_w = pd*self.Ss*self.lf*cn

        # Pasamos de ejes viento a ejes cuerpo
        Cal, Sal = cos(al), sin(al)
        Cbe, Sbe = cos(be), sin(be)
        # Matriz de cambio de ejes viento a ejes cuerpo
        Lbw = array( [
            [ Cal*Cbe, -Cal*Sbe, -Sal ],
            [ Sbe,      Cbe,      0.0 ],
            [ Sal*Cbe, -Sal*Sbe,  Cal ]
            ], type = 'Float64')
        # Fuerzas y momentos en el centro aerodinamico, en ejes cuerpo
        Xfca_b, Yfca_b, Zfca_b = dot(Lbw, [-Df_w, Yf_w, -Lf_w])
        Lfca_b, Mfca_b, Nfca_b = dot(Lbw, [lf_w, mf_w, nf_w])

        # Sobre el CM
        Xf, Yf, Zf = Xfca_b, Yfca_b, Zfca_b
        Mf = Mfca_b + Xf*self.zca - Zf*self.xca
        Lf = Lfca_b - Yf*self.zca
        Nf = Nfca_b + Yf*self.xca

        # Fuerzas y Momentos del fuselaje
        self.FMf = ( 
                array([Xf, Yf, Zf, Lf, Mf, Nf], type='Float64'), 
                zeros(shape=(6, 16), type='Float64')
        )

        # logging
        self.be = be
        self.al = al
        self.cD, self.cL, self.cY = cD, cL, cY
        self.cl, self.cm, self.cn = cl, cm, cn
        self.Xf_b = Xf
        self.Yf_b = Yf
        self.Zf_b = Zf
        self.Lf_b = Lf
        self.Mf_b = Mf
        self.Nf_b = Nf

