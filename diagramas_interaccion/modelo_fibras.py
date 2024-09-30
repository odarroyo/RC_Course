# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 13:03:33 2024

@author: HOME
"""

import numpy as np
import matplotlib.pyplot as plt
#%% Definicion funciones de los materiales

def elastoplastico(es,Es = 210000,fy = 420):
    fs = es*Es
    if fs > fy:
        fs = fy
    elif fs < -fy:
        fs = -fy
    return fs

def noconf(ec, fc):
    Ec = 4300*np.sqrt(fc)
    e0 = 1.8*fc/Ec # Hognestad
    eu = 0.0038
    m = -0.15*fc/(eu - e0)
    if ec < 0:
        f = 0
    elif 0 < ec < e0:
        f = fc*(2*ec/e0 - (ec/e0)**2)
    elif e0 < ec < eu:
        f = fc + m*(ec-e0)
    else:
        f = 0
    return f

# se convierten las funciones a forma vectorial para facilidad de uso

vec_elastoplastico = np.vectorize(elastoplastico) # se convierten a funciones vectorizables
vec_noconf = np.vectorize(noconf) # se convierten a funciones vectorizables
#%% Verificación de los materiales
eps = np.linspace(-0.004,0.004,1000)
esf_c = np.zeros(len(eps))
esf_s = np.zeros(len(eps))
for ind,ep in enumerate(eps):
    esf_c[ind] = noconf(ep,28)
    esf_s[ind] = elastoplastico(ep)

plt.plot(eps,esf_c)
plt.xlabel('Deformación unitaria')
plt.ylabel('Esfuerzo (MPa)')
plt.show()

plt.plot(eps,esf_s)
plt.xlabel('Deformación unitaria')
plt.ylabel('Esfuerzo (MPa)')
plt.show()
    
#%% Modelo de la sección rectangular
ec = 0.003 # valor que se asume de compresión en fibra extrema
b = 0.3 # base
h = 0.4 # altura
nfib = 8 # numero de fibras a utilizar
fy = 420
fc = 28
dA = b*h/nfib # diferencial de área de las fibras
posfibra = np.linspace(0.025,0.375,nfib) # posición del centroide de cada fibra
posacero = np.array([0.05,0.20,0.35]) # posición de las barras de acero
area_acero = np.array([3*2,2*2,3*2])/10000 # áreas de acero
pos_en = np.linspace(0.01,h*1.3,500) # cantidad de posiciones del eje neutro posición del eje neutro
fcon, fac = [],[] # esfuerzos en el concreto y el acero
P = np.zeros(len(pos_en))
M = np.zeros(len(pos_en))
for ind,c in enumerate(pos_en):
    e_s = ec*(c-posacero)/c
    e_c = ec*(c-posfibra)/c
    f_c = vec_noconf(e_c,fc) # aquí se puede modificar si se desea otro modelo del concreto
    f_s = vec_elastoplastico(e_s) # aquí se puede modificar si se desea otro modelo del acero
    fcon.append(f_c)
    fac.append(f_s)
    F_c = f_c*dA*1000
    F_s = f_s*area_acero*1000
    M_c = F_c*(h/2 - posfibra)
    M_s = F_s*(h/2- posacero)
    P[ind] = np.sum(F_c) + np.sum(F_s)
    M[ind] = np.sum(M_c) + np.sum(M_s)

Pt = -np.sum(area_acero*fy)*1000
P0 = fc*(b*h-np.sum(area_acero))*1000+np.sum(area_acero)*fy*1000
P = np.concatenate(([Pt],P,[P0]))
M = np.concatenate(([0],M,[0]))

plt.plot(M,P)
plt.xlabel('Momento (kN-m)')
plt.ylabel('Fuerza axial (kN)')
plt.show()

fcon = np.array(fcon) # esfuerzos del concreto
fac = np.array(fac) # esfuerzos del acero

    
