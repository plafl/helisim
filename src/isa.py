# Atmosfera ISA

# Condiciones a nivel del mar
P0 = 101.3e3    # Pa
T0 = 288        # K

# Gravedad
g  = 9.81       # m/s^2

# Ley de temperatura en la troposfera: T = T0 - L*h
L  = 0.0065     # K/m

# Constante del gas
Ra = 287        # J/(Kg*K)

# h en metros, T en kelvins
def T(h):
    return T0 - L*h

# h en metros, P en pascales
def P(h):
    return P0*(T(h)/T0)**(g/L/Ra)


# Conversion
def inHg_to_Pa(x):
    return 3.386389e3*x

def Pa_to_inHg(x):
    return x/3.386389e3
