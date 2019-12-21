# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:03:25 2019

@author: zhuguohua
"""
# 用于分析output 
import numpy as np
import control as ct
import control.matlab as cmt
import matplotlib.pyplot as plt

# 车辆参数
M   = 1573.0 #(kg) 总质量 mass
I_z = 2873.0 # (kg/m^2) 绕Z轴的转动惯量 
l_f = 1.10  # (m)
l_r = 1.58  # (m)
C_alpha_f = 80000.0 #(N/rad) 前轮总侧偏刚度
C_alpha_r = 80000.0 #(N/rad) 后轮总侧偏刚度
V_x = 25.0 # (m/s)
ds = 2.
K = 0.1
Tn = 0.5
Td = 0.1

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

C = np.array([[1., 0., ds,0.]])

D = np.array([[0.]])


# (1)Plant Transfunction
sys_ss_p = ct.ss(A,B1,C,D)
sys_tf_p = ct.ss2tf(sys_ss_p)

plt.close()
mag,phase,omega = ct.bode_plot(sys_ss_p,np.linspace(0.1,10)*2*np.pi,Plot=False,dB=True,Hz=False,margins=True)

ct.pzmap(sys_ss_p,grid=False)

plt.figure()
plt.subplot(2,1,1)
plt.plot(omega,np.log10(mag))
plt.ylabel("mag")
plt.subplot(2,1,2)
plt.ylabel("phase(deg)")
plt.xlabel("frequency(rad/s)")
plt.plot(omega,phase*57.3-360)

# controllor function
sys_tf_c = K*ct.tf([Tn,1.],[Td,1.])

# PC function
sys_tf_pc = ct.series(sys_tf_p,sys_tf_c)

plt.figure()
real, imag, freq  = ct.nyquist_plot(sys_tf_pc,Plot=True)

#plt.figure()
#plt.grid()
#plt.plot(-1,0,'r+')
#plt.plot(real,imag)
#index = 22
#plt.arrow(real[index],imag[index]+0.01,(real[index+1]-real[index])/2,(imag[index+1]-imag[index])/2,head_width=0.02, head_length=0.02)
#plt.plot(real,-imag)
#plt.xlim(-1.5,0.5)
#plt.ylim(-0.1,0.1)

plt.figure()
mag,phase,omega = ct.bode_plot(sys_tf_pc,np.logspace(-2,1)*2*np.pi,dB=True,Hz=False,margins=True)

plt.figure()
plt.subplot(2,1,1)
plt.plot(omega,np.log10(mag))
plt.ylabel("mag")
plt.subplot(2,1,2)
plt.ylabel("phase(deg)")
plt.xlabel("frequency(rad/s)")
plt.plot(omega,phase*57.3)

ct.root_locus(sys_tf_pc)

