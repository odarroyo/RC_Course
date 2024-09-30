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

def conf(ecc,fcc):
    Ec = 4300*np.sqrt(fcc)
    e0 = 1.8*fcc/Ec # Hognestad
    eu = 5*e0 # el modelo asume que la deformación última son 5 veces la de fcc
    m = -0.8*fcc/(eu - e0)
    if ecc < 0:
        f = 0
    elif 0 < ecc < e0:
        f = fcc*(2*ecc/e0 - (ecc/e0)**2)
    elif e0 < ecc < eu:
        f = fcc + m*(ecc-e0)
    else:
        f = 0
    return f

# se convierten las funciones a forma vectorial para facilidad de uso

vec_elastoplastico = np.vectorize(elastoplastico) # se convierten a funciones vectorizables
vec_noconf = np.vectorize(noconf) # se convierten a funciones vectorizables
vec_conf = np.vectorize(conf)

#%% Verificación de los materiales
eps = np.linspace(-0.004,0.015,1000)
esf_c = np.zeros(len(eps))
esf_s = np.zeros(len(eps))
esf_cc = np.zeros(len(eps))

for ind,ep in enumerate(eps):
    esf_c[ind] = noconf(ep,28)
    esf_s[ind] = elastoplastico(ep)
    esf_cc[ind] = conf(ep,28*1.3)

plt.plot(eps,esf_s)
plt.xlabel('Deformación unitaria')
plt.ylabel('Esfuerzo (MPa)')
plt.show()

plt.plot(eps,esf_c,eps,esf_cc)
plt.xlabel('Deformación unitaria')
plt.ylabel('Esfuerzo (MPa)')
plt.legend(['no confinado','confinado'])
plt.show()
    
#%% Modelo de la sección 
ec = 0.003 # valor que se asume de compresión en fibra extrema
b = 0.6 # base
h = 0.4 # altura
nfib = 8 # numero de fibras a utilizar
fy = 420
fc = 28
dA = b*h/nfib # diferencial de área de las fibras
posfibra = np.linspace(0.025,0.375,nfib) # posición del centroide de cada fibra
posacero = np.array([0.05,0.2,0.35]) # posición de las barras de acero
area_acero = np.array([4*2.86,2*2.86,4*2.86])/10000 # áreas de acero
pos_en = np.linspace(0.01,h*1.3,500) # cantidad de posiciones del eje neutroposición del eje neutro
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

#%% Modelo considerando el concreto confinado

ec = 0.003 # valor que se asume de compresión en fibra extrema
b = 0.6 # base
h = 0.4 # altura
bcc = 0.5
hcc = 0.3
r = 0.05 # recubrimiento
nfib = 8 # numero de fibras a utilizar
fy = 420
fc = 28
K = 1.3 # K calculado según Mander
posacero = np.array([0.05,0.2,0.35]) # posición de las barras de acero
area_acero = np.array([4*2.86,2*2.86,4*2.86])/10000 # áreas de acero
fcc = fc*K
dh = h/nfib
posfibra = np.linspace(dh/2,h-dh/2,nfib) # posición del centroide de cada fibra
area_fibra = np.concatenate([np.array([b*dh]),np.array([0.1*dh]*(nfib-2)),np.array([b*dh])])
posfibra_c = np.linspace(3*dh/2,h-3*dh/2,nfib-2) # posición del centroide de las fibras confinadas
area_fibra_c = np.array([(b-0.1)*dh]*(nfib-2))
pos_en = np.linspace(0.01,h*1.3,500) # cantidad de posiciones del eje neutroposición del eje neutro
fcon2, fac2, fconc = [],[],[] # esfuerzos en el concreto y el acero
ec2, es2, ecc2 = [],[],[]
Pc = np.zeros(len(pos_en))
Mc = np.zeros(len(pos_en))
for ind,c in enumerate(pos_en):
    e_s = ec*(c-posacero)/c
    e_c = ec*(c-posfibra)/c
    e_cc = ec*(c-posfibra_c)/c
    f_c = vec_noconf(e_c,fc) # aquí se puede modificar si se desea otro modelo del concreto
    f_s = vec_elastoplastico(e_s) # aquí se puede modificar si se desea otro modelo del acero
    f_cc = vec_conf(e_cc,fcc)
    fcon2.append(f_c)
    fac2.append(f_s)
    fconc.append(f_cc)
    ec2.append(e_c)
    es2.append(e_s)
    ecc2.append(e_cc)
    F_c = f_c*area_fibra*1000 # fuerzas en el concreto no confinado
    F_s = f_s*area_acero*1000 # fuerzas en el acero
    F_cc = f_cc*area_fibra_c*1000 # fuerzas en el concreto confinado
    M_c = F_c*(h/2 - posfibra) # momentos del concreto
    M_s = F_s*(h/2- posacero) # momentos del acero
    M_cc = F_cc*(h/2 - posfibra_c) # momnetos concreto confinado
    Pc[ind] = np.sum(F_c) + np.sum(F_s) + np.sum(F_cc)
    Mc[ind] = np.sum(M_c) + np.sum(M_s) + np.sum(M_cc)

Pt = -np.sum(area_acero*fy)*1000
P0 = fc*(b*h-bcc*hcc)*1000+np.sum(area_acero)*fy*1000+fcc*(bcc*hcc-np.sum(area_acero))*1000
Pc = np.concatenate(([Pt],Pc,[P0]))
Mc = np.concatenate(([0],Mc,[0]))

plt.plot(Mc,Pc)
plt.xlabel('Momento (kN-m)')
plt.ylabel('Fuerza axial (kN)')
plt.show()

fcon = np.array(fcon) # esfuerzos del concreto
fac = np.array(fac) # esfuerzos del acero


#%% comparacion
plt.plot(M,P,Mc,Pc)
plt.xlabel('Momento (kN-m)')
plt.ylabel('Fuerza axial (kN)')
plt.legend(['No confinado','Confinado'])
plt.show()

#%% Gráfica de deformaciones unitarias y esfuerzos

ind = 89
titulo = 'c = '+str(np.round(pos_en[ind],3))
plt.plot(ec2[ind],-posfibra,'o-')
plt.title(titulo)
plt.show()

plt.plot(fcon2[ind],-posfibra,'o-')
plt.title(titulo)
plt.show()

plt.plot(fconc[ind],-posfibra_c,'o-')
plt.title(titulo)
plt.show()
