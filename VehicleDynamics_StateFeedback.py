# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:03:25 2019

@author: zhuguohua
"""

import numpy as np
import control as ct
import matplotlib.pyplot as plt

# 车辆参数
M   = 1573.0 #(kg) 总质量 mass
I_z = 2873.0 # (kg/m^2) 绕Z轴的转动惯量 
l_f = 1.10  # (m)
l_r = 1.58  # (m)
C_alpha_f = 80000.0 #(N/rad) 前轮总侧偏刚度
C_alpha_r = 80000.0 #(N/rad) 后轮总侧偏刚度
V_x = 30.0 # (m/s)

A = np.array([
        [0.,1.,0.,0.],
        [0.,-2.*(C_alpha_f + C_alpha_r)/(M*V_x),2.*(C_alpha_f + C_alpha_r)/M,-2.*(C_alpha_f*l_f - C_alpha_r*l_r)/(M*V_x)],
        [0.,0.,0.,1.],
        [0.,-2.*(C_alpha_f*l_f - C_alpha_r*l_r)/(I_z*V_x),2.*(C_alpha_f*l_f - C_alpha_r*l_r)/I_z,-2.*(C_alpha_f*l_f**2 + C_alpha_r*l_r**2)/(I_z*V_x)]
        ])

B1 = np.array([
        [0.],
        [2.*C_alpha_f/M],
        [0.],
        [2.*l_f*C_alpha_f/I_z]
        ])

B2 = np.array([[0.],[-2.*(C_alpha_f*l_f - C_alpha_r*l_r)/(M*V_x) - V_x],[0.],[-2.*(C_alpha_f*l_f**2 + C_alpha_r*l_r**2)/(I_z*V_x)]])

C = np.array([[1., 0., 0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])

D = np.array([[0.],[0.],[0.],[0]])

w,v = np.linalg.eig(A) # 求取数组的特征值
print("特征值：",w)
#print(v)
P = np.array([-5.-3.j,-5.+3.j,-7.,-10.])
K = ct.place(A,B1,P)
print("反馈值K:",K)

sys_init = ct.ss(A,B1,C,D)
ct.pzmap(sys_init)
sys_task = ct.ss(A-B1*K,B2,C,D)

ct.pzmap(sys_task)

t = np.linspace(0, 10, 101)
u = np.zeros(len(t))
u[11:101] = 1.72

#s_t,s_yout = ct.step_response(sys_task,t)

f_t,f_yout,f_xout = ct.forced_response(sys_task,t,u)

#plt.close()
plt.figure()
plt.clf()
plt.grid()
plt.subplot(3,1,1)
plt.plot(f_t,u)
plt.subplot(3,1,2)
plt.plot(f_t,f_yout[0])
plt.subplot(3,1,3)
plt.plot(f_t,f_yout[2])
