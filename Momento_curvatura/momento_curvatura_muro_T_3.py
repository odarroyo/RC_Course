# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 08:50:13 2024

@author: User
"""

import opseestools.utilidades as ut # requiere version 0.45 o superior de opseestools
import opseestools.analisis as an 
from openseespy.opensees import *
import opsvis as opsv
import matplotlib.pyplot as plt
import numpy as np


wipe() # es una buena práctica comenzar con un wipe 
model('basic','-ndm',2,'-ndf',3)

#%% Definir materiales
# ======================

fc = 28 # fc del concreto
fy = 420 # fy del acero
noconf,conf,acero = ut.col_materials(fc,fy,'DES',tension='no',conftag=1,unctag=2,steeltag=3)

mallas = 4
p1=509.10*1000
p2=691.50*1000
p3=734.44*1000
p4 = p1*0.02
e1=0.00248
e2=0.005
e3=0.012
e4 = 0.015
pinchX=0.34
pinchY= 0.56
damage1= 0.038
damage2= 0.07
beta= 0.086

# ops.uniaxialMaterial('Hysteretic', mallas, p1, e1, p2, e2, p3, e3, -p1, -e1, -p2, -e2, -p3, -e3, pinchX, pinchY, damage1, damage2, beta)
uniaxialMaterial('HystereticSM', mallas, '-posEnv' ,p1, e1, p2, e2, p3, e3, p4, e4, '-negEnv', -p1, -e1, -p2, -e2, -p3, -e3, -p4, -e4)
#%% Definir seccion en T
# ======================

c = 0.05  # recubrimiento 
bf = 4.0  # ancho de la aleta
tf = 0.1  # espesor de la aleta
Ht = 2.5  # altura total de la sección
h1 = Ht - tf
tw = 0.1  # ancho del alma
yhat = ((bf*tf*tf/2)+(h1*tw*(h1/2+tf)))/(bf*tf+h1*tw) # centroide

# zonas de confinamiento
l1 = 0.3

# estas lineas son para calcular la cantidad de fibras que salgan aprox cada 10cm
nFibZw = int(tw/0.1)
nFibZf = int(bf/0.1)
nFibYw = int(h1/0.1)
nFibYf = int(tf/0.1)

nFibConf = int(l1/0.1)+1

# Las areas de las barras
As4 = 0.000127
As5 = 0.0002
As7 = 0.000387
Asfi7 = np.pi*(0.7/100)**2/4

secT = 1 # tag de la sección a crear

# Se empleará la librería de opsvis que permite visualizar la sección
sM = [['section', 'Fiber', secT, '-GJ', 1.0e6],
             ['patch', 'rect', noconf, nFibYf, nFibZf-nFibConf, yhat-tf, -bf/2+l1, yhat, bf/2-l1], # aleta
             ['patch', 'rect', conf, nFibYf, nFibConf, yhat-tf, -bf/2, yhat, -bf/2+l1], # aleta
             ['patch', 'rect', conf, nFibYf, nFibConf, yhat-tf, bf/2-l1, yhat, bf/2], # aleta
             ['patch', 'rect', conf, nFibConf, nFibZw, -Ht+yhat, -tw/2, -Ht+yhat+l1, tw/2], # alma
              ['patch', 'rect', noconf, nFibYw-nFibConf, nFibZw, -Ht+yhat+l1, -tw/2, yhat-tf, tw/2], # alma
             ['layer', 'straight', mallas, int((bf-2*c-2*l1)/0.1), Asfi7, yhat-c, -bf/2+c+l1, yhat-c, bf/2-c-l1], # aleta
             ['layer', 'straight', mallas, nFibConf, As5, yhat-c, -bf/2+c, yhat-c, -bf/2+l1-c], # aleta
             ['layer', 'straight', mallas, nFibConf, As5, yhat-c, bf/2-l1+c, yhat-c, bf/2-c], # aleta
             ['layer', 'straight', mallas, int((Ht-2*c-tf-l1)/0.1), Asfi7, yhat-c, 0, -Ht+yhat+c+l1, 0], # alma
             ['layer', 'straight', mallas, nFibConf, As5, -Ht+yhat+c, 0, -Ht+yhat+l1-c, 0]] # alma]


# Los siguientes comandos dibujan la sección para saber que nos quedó bien
matcolor = ['y', 'lightgrey', 'gold', 'w', 'w', 'w']
# matcolor=['y', 'b', 'r', 'g', 'm', 'k', 'k', 'k']
#plt.figure()
opsv.plot_fiber_section(sM,matcolor=matcolor)
plt.axis('equal')
plt.show()

# El siguiente comando es el que propiamente crea la sección escribiendo el código definido en la lista
opsv.fib_sec_list_to_cmds(s30x30)
pint = 5 # puntos de integración
beamIntegration('Lobatto', secT, secT,pint) # esquema de integración

#%% Momento curvatura
# ======================
P = 0.0 # Carga axial
phi_u = -0.01 # curvatura hasta la que se intentará llegar
Mom,curv = ut.MomentCurvature(secT, P, phi_u)

print('Momento max: ', np.max(np.abs(Mom)))
print('curvatura mom. max: ', curv[np.argmax(np.abs(Mom))])