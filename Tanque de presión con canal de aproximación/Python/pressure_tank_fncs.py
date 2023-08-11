# Necessary modules
import numpy as np
from scipy.optimize import fsolve

def pressure_tank_design(Qd, D, tp, vd, I, n, alpha, Cf, b, Cd, a, s, a1, a4):
    # Start
    S = {}
    g = 9.81

    ## Calado
    f = lambda y: (1 / n) * ((b * y) / (b + 2 * y)) ** (2 / 3) * I ** (1 / 2) * b * y - Qd
    y_values = np.linspace(.1, 2*b, 100)
    idx = np.argmin(abs(f(y_values)))
    y0 = y_values[idx]

    y = fsolve(f, y0)[0]
    S['Calado'] = y

    ## Velocidad
    v_ca = Qd / (b * y)
    S['CA: Velocidad'] = v_ca

    ## Froude
    Fr_ca = v_ca / np.sqrt(g * y)
    S['CA: Froude'] = Fr_ca

    if Fr_ca < 1: S['CA: Tipo de flujo'] = 'Subcrítico'
    elif Fr_ca > 1: S['CA: Tipo de flujo'] = 'Supercrítico'

    ## Area hidraulica del flujo
    Ah = Qd / (Cd * vd)
    S['Area hidráulica'] = Ah

    ## Altura encima de la tubería
    # Numero de Froude y velocidad de la tuberia
    Fr_t = (4 * Qd) / (np.pi * D ** 2) * (1 / np.sqrt(g * D))
    S['T: Froude'] = Fr_t

    if Fr_t < 1: S['T: Tipo de flujo'] = 'Subcrítico'
    elif Fr_t > 1: S['T: Tipo de flujo'] = 'Supercrítico'

    vt = (4 * Qd) / (np.pi * D ** 2)
    S['T: Velocidad'] = vt

    if 2 <= vt <= 6: S['T: Verificación de velocidad'] = 'OK'
    else: S['T: Verificación de velocidad'] = 'Rediseñar'

    # Propuesta Knauss
    a3_kn = D * (2 * Fr_t + 0.5)

    # Propuesta Krochin
    k = 3
    a3_kr = k * (vt ** 2) / (2 * g)
    if a3_kr < 1: a3_k3 = 0

    # Altura final
    a3 = max([a3_kn, a3_kr])
    S['Altura de sumergencia'] = a3

    ## Ancho efectivo
    Hr = a3
    be = Ah / Hr
    S['Ancho efectivo'] = be

    ## Numero de espacios y numero de barras
    N = np.floor(be/a)
    S['Número de espacios'] = N
    Nb = N - 1
    S['Número de barras'] = Nb

    ## Ancho total
    bt = be + Nb * s
    S['Ancho total reja'] = bt

    ## Perdidas en la reja
    kr = Cf * (s / a) ** (4/3) * np.sin(np.deg2rad(alpha))
    h_lr = kr * vd ** 2 / (2 * g)
    S['Pérdidas reja'] = h_lr

    ## Volumen del tanque
    # Krochin
    A = b * y
    DV_k = (.693 * A * v_ca ** 2) / (g * I)

    # Continuidad
    DV_c = Qd * tp

    # Volumen final
    DV = max([DV_k, DV_c])
    S['Volumen del tanque'] = DV

    ## Altura del tanque
    a2 = D
    h_min = a1 + a2 + a3 + a4

    if (h_min % .1) < 1e-3: h = h_min - (h_min % .1)
    else: h = h_min - (h_min % .1) + .1

    S['Altura del tanque'] = h

    ## Area del tanque
    AT = DV/h
    S['Área del tanque'] = AT

    ## Dimensiones del tanque: [B L]. Y borde extra: bs
    Bs = 1
    B = np.round(bt + 2 * Bs, 0)
    Bs = (B - bt) / 2
    S['Borde extra reja'] = Bs

    L = AT / B
    if (L % .1) < 1e-3: L = L - (L % .1)
    else: L = L - (L % .1) + .1
    S['Dimensiones del tanque'] = [B, L]

    return S