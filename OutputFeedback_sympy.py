# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 20:18:22 2019

@author: zhuguohua
"""

import sympy as sp

sp.init_printing(use_unicode=True)
# 车辆参数
m,I_z,l_f,l_r,C_alpha_f,C_alpha_r,V_x,R,delta_ff =  sp.symbols('m I_z l_f l_r C_alpha_f C_alpha_r V_x R delta_ff')
k1,k2,k3,k4 =  sp.symbols('k1 k2 k3 k4')
s = sp.symbols('s')
A = sp.Matrix([
        [0,1,0,0],
        [0,-2*(C_alpha_f + C_alpha_r)/(m*V_x),2*(C_alpha_f + C_alpha_r)/m,-2*(C_alpha_f*l_f - C_alpha_r*l_r)/(m*V_x)],
        [0,0,0,1],
        [0,-2*(C_alpha_f*l_f - C_alpha_r*l_r)/(I_z*V_x),2*(C_alpha_f*l_f - C_alpha_r*l_r)/I_z,-2*(C_alpha_f*l_f**2 + C_alpha_r*l_r**2)/(I_z*V_x)]
        ])
print("A:\r\n",sp.latex(A))

B1 = sp.Matrix([
        [0],
        [2*C_alpha_f/m],
        [0],
        [2*l_f*C_alpha_f/I_z]
        ])
print("B1:\r\n",sp.latex(B1))

B2 = sp.Matrix(
        [[0],
         [-2*(C_alpha_f*l_f - C_alpha_r*l_r)/(m*V_x) - V_x],
         [0],
         [-2*(C_alpha_f*l_f**2 + C_alpha_r*l_r**2)/(I_z*V_x)]
         ])

S = sp.Matrix(
        [[s,0,0,0],[0,s,0,0],[0,0,s,0],[0,0,0,s]]
        )
print("B2:\r\n",sp.latex(B2))

K = sp.Matrix([[k1,k2,k3,k4]])
print("K:\r\n",sp.latex(K))

D = (S - A)**-1
print("D:\r\n",sp.latex(D))
E = B1
X_ss = D*E


X_factor = sp.factor(X_ss)
print("x_ss:\r\n",sp.latex(X_factor))



