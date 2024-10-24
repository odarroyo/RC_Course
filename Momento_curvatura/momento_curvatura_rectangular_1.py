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
#%% Definir seccion
# ======================

# Se presenta el ejemplo de la creación de una sección de fibras
# seccion de 30 x 30 cm con 6 # 4

Bcol = 0.3
Hcol = 0.3
c = 0.05  # recubrimiento 

# creación de la sección de fibra
y1col = Hcol/2.0
z1col = Bcol/2.0

y2col = 0.5*(Hcol-2*c)/3.0

nFibZ = 1
nFibZcore= 10
nFib = 20
nFibCover, nFibCore = 3, 16
As4 = 0.000127
As5 = 0.0002
As6 = 0.00028
As7 = 0.000387
As8 = 0.000508

sec30x30 = 1 # tag de la sección a crear

# Se empleará la librería de opsvis que permite visualizar la sección
s30x30 = [['section', 'Fiber', sec30x30, '-GJ', 1.0e6],
             ['patch', 'rect', noconf, nFib, nFibZcore, -y1col, -z1col, y1col, z1col],
             ['layer', 'straight', acero, 3, As4, y1col-c, z1col-c, y1col-c, c-z1col],
             ['layer', 'straight', acero, 3, As4, c-y1col, z1col-c, c-y1col, c-z1col]]

# Los siguientes comandos dibujan la sección para saber que nos quedó bien
matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
opsv.plot_fiber_section(s30x30,matcolor=matcolor)
plt.axis('equal')
plt.show()

# El siguiente comando es el que propiamente crea la sección escribiendo el código definido en la lista
opsv.fib_sec_list_to_cmds(s30x30)


#%% Momento curvatura
# ======================
P = 0.0 # Carga axial
phi_u = 0.3 # curvatura hasta la que se intentará llegar
ut.MomentCurvature(sec30x30, P, phi_u)