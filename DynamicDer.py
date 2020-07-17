# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 17:23:03 2020

@author: zhuguohua
"""

import sympy as sp
import scipy.linalg as la
import numpy as np


sp.init_printing(use_unicode=True)
# 车辆参数
m,I_z,l_f,l_r,C_alpha_f,C_alpha_r,V_x,R,delta_ff =  sp.symbols('m I_z l_f l_r C_alpha_f C_alpha_r V_x R delta_ff')
k1,k2,k3,k4 =  sp.symbols('k1 k2 k3 k4')
dT =  sp.symbols('dT')
A = sp.Matrix([
        [0.0,1.0,0.0,0.0],
        [0.0,-2.0*(C_alpha_f + C_alpha_r)/(m*V_x),2.0*(C_alpha_f + C_alpha_r)/m,-2.0*(C_alpha_f*l_f - C_alpha_r*l_r)/(m*V_x)],
        [0.0,0.0,0.0,1.0],
        [0.0,-2.0*(C_alpha_f*l_f - C_alpha_r*l_r)/(I_z*V_x),2.0*(C_alpha_f*l_f - C_alpha_r*l_r)/I_z,-2.0*(C_alpha_f*l_f**2 + C_alpha_r*l_r**2)/(I_z*V_x)]
        ])
#print("A:\r\n",sp.latex(A))

B1 = sp.Matrix([
        [0.0],
        [2.0*C_alpha_f/m],
        [0.0],
        [2.0*l_f*C_alpha_f/I_z]
        ])
#print("B1:\r\n",sp.latex(B1))

#Ad = np.zeros((4, 4))

#Ad = sp.zeros((4,4))

Ad = (1.0*sp.eye(4) - dT * 0.5 * A)**-1 *(1.0*sp.eye(4) + dT * 0.5 * A)

print("B1:\r\n",sp.latex(Ad))