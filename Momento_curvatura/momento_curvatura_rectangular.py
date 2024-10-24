# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 08:50:13 2024

@author: User
"""

import opseestools.utilidades as ut # requiere version 0.45 o superior de opseestools
import opseestools.analisis as an 
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np


wipe() # es una buena práctica comenzar con un wipe 
model('basic','-ndm',2,'-ndf',3)

#%% Definir materiales
# ======================

fc = 28 # fc del concreto
fy = 420 # fy del acero
noconf,conf,acero = ut.col_materials(fc,fy,'DES',tension='no')

#%% Definir seccion
# ======================

sec = 100 # tag de la sección
As4 = 1.27e-4 # area de la barra # 4
ut.create_rect_RC_section(sec, 0.3, 0.3, 0.05, conf, noconf, acero, 3, As4, 3, As4)

#%% Momento curvatura
# ======================
P = 0.0 # Carga axial
phi_u = 0.3 # curvatura hasta la que se intentará llegar
ut.MomentCurvature(sec, P, phi_u)