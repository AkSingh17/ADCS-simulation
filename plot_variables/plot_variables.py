import matplotlib.pyplot as plt
import numpy as np
import math

a,b,c,d,e=[],[],[],[],[]
a1,b1,c1,d1,e1=[],[],[],[],[]
a2,b2,c2,d2,e2=[],[],[],[],[]
a3,b3,c3,d3,e3=[],[],[],[],[]
a4,b4,c4,d4,e4=[],[],[],[],[]
a5,b5,c5,d5,e5=[],[],[],[],[]
figure, axis = plt.subplots(3, 2)

for line in open('../ADCS-Codes/q_ri.txt','r'):
    values = [float(s) for s in line.split()]
    a.append(values[0])
    b.append(values[1])
    c.append(values[2])
    d.append(values[3])
    e.append(values[4])

axis[0,0].plot(a,b, label='q1')
axis[0,0].plot(a,c, label='q2')
axis[0,0].plot(a,d, label='q3')
axis[0,0].plot(a,e, label='q4')
axis[0,0].legend(loc='lower right', frameon=False)
axis[0,0].set_xlabel('time(s)')
axis[0,0].set_ylabel('q_ri')
axis[0,0].set_title("q_ri vs time")



for line1 in open('../ADCS-Codes/w_bi.txt','r'):
    values = [float(s) for s in line1.split()]
    a1.append(values[0])
    b1.append(values[1])
    c1.append(values[2])
    d1.append(values[3])

axis[0,1].plot(a1,b1, label='wbi_x')
axis[0,1].plot(a1,c1, label='wbi_y')
axis[0,1].plot(a1,d1, label='wbi_z')
axis[0,1].legend(loc='lower right', frameon=False)
axis[0,1].set_xlabel('time(s)')
axis[0,1].set_ylabel('w_bi')
axis[0,1].set_title("w_bi vs time")


for line2 in open('../ADCS-Codes/q_error.txt','r'):
    values = [float(s) for s in line2.split()]
    a2.append(values[0])
    b2.append(values[1])
    c2.append(values[2])
    d2.append(values[3])
    e2.append(values[4])

axis[1,0].plot(a2,b2, label='q1')
axis[1,0].plot(a2,c2, label='q2')
axis[1,0].plot(a2,d2, label='q3')
axis[1,0].plot(a2,e2, label='q4')
axis[1,0].legend(loc='lower right', frameon=False)
axis[1,0].set_xlabel('time(s)')
axis[1,0].set_ylabel('q_error')
axis[1,0].set_title("q_error vs time")


for line3 in open('../ADCS-Codes/control_torque.txt','r'):
    values = [float(s) for s in line3.split()]
    a3.append(values[0])
    b3.append(values[1])
    c3.append(values[2])
    d3.append(values[3])

axis[1,1].plot(a3,b3, label='tc_x')
axis[1,1].plot(a3,c3, label='tc_y')
axis[1,1].plot(a3,d3, label='tc_z')
axis[1,1].legend(loc='lower right', frameon=False)
axis[1,1].set_xlabel('time(s)')
axis[1,1].set_ylabel('control_torque')
axis[1,1].set_title("control_torque vs time")


for line4 in open('../ADCS-Codes/q_bi.txt','r'):
    values = [float(s) for s in line4.split()]
    a4.append(values[0])
    b4.append(values[1])
    c4.append(values[2])
    d4.append(values[3])
    e4.append(values[4])

axis[2,0].plot(a4,b4, label='q1')
axis[2,0].plot(a4,c4, label='q2')
axis[2,0].plot(a4,d4, label='q3')
axis[2,0].plot(a4,e4, label='q4')
axis[2,0].legend(loc='lower right', frameon=False)
axis[2,0].set_xlabel('time(s)')
axis[2,0].set_ylabel('q_bi')
axis[2,0].set_title("q_bi vs time")

for line5 in open('../ADCS-Codes/w_ri.txt','r'):
    values = [float(s) for s in line5.split()]
    a5.append(values[0])
    b5.append(values[1])
    c5.append(values[2])
    d5.append(values[3])

axis[2,1].plot(a5,b5, label='wri_x')
axis[2,1].plot(a5,c5, label='wri_y')
axis[2,1].plot(a5,d5, label='wri_z')
axis[2,1].legend(loc='lower right', frameon=False)
axis[2,1].set_xlabel('time(s)')
axis[2,1].set_ylabel('w_ri')
axis[2,1].set_title("w_ri vs time")
plt.show()