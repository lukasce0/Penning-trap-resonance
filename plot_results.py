import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d



plt.rcParams.update({'font.size': 13})


"""
Note that there are much more compact ways of doing this, for instance putting the
file reading inside a function or puting the rest into loops and avoid so many variables, 
but I consider the plotting in python only to be a formality as we have been doing this 
for multiple years. Therefore I have only little comments and doing everything in time
efficent way for me, and not necessarily elegant.
"""




#Ploting one particle and comparing to analytical solution
positions1 = np.loadtxt('positions81RK4.txt')
t1 = np.loadtxt('time81RK4.txt')

positions2 = np.loadtxt('positions81FE.txt')
t2 = np.loadtxt('time81FE.txt')

r1 = positions1
r2 = positions2

z1 = []
z2 = []

for i in range(0,len(positions1),3):
    z1.append(r1[i+2])
    
for i in range(0,len(positions2),3):
    z2.append(r2[i+2])
    

x = np.linspace(0,50,1000)
z_T = 20 * np.cos(np.sqrt(0.482) * x) #theoretical z


plt.figure(figsize=(8,6))
plt.xlabel(r'Time [$\mu$s]')
plt.ylabel(r'Position along z [$\mu$m]')
plt.ylim(-30,45)
plt.plot(t1, z1, label=f"RK4 solution with {len(t1)-1} steps")
plt.plot(t2, z2, label=f"FE solution with {len(t2)-1} steps")
plt.plot(x, z_T, linewidth=0.7, label=f'Analytical solution')
plt.legend()
plt.savefig("theor_vs_analytic_81.pdf")




#Two particles with and without interactions
positions1 = np.loadtxt('positions82no.txt')
positions2 = np.loadtxt('positions82yes.txt')

r1n = positions1[:,0]
r2n = positions1[:,1]
r1y = positions2[:,0]
r2y = positions2[:,1]

x1n = []
y1n = []
x2n = []
y2n = []

x1y = []
y1y = []
x2y = []
y2y = []



for i in range(0,len(positions1),3):
    x1n.append(r1n[i])
    y1n.append(r1n[i+1])
    x2n.append(r2n[i])
    y2n.append(r2n[i+1])

for i in range(0,len(positions2),3):
    x1y.append(r1y[i])
    y1y.append(r1y[i+1])
    x2y.append(r2y[i])
    y2y.append(r2y[i+1])




plt.figure(figsize=(14,7))
plt.subplot(1,2,1)
plt.title(f"No interactions between particles {len(x1n)-1} steps")
plt.xlabel(r'Position along x [$\mu$m]')
plt.ylabel(r'Position along y [$\mu$m]')
plt.plot(x1n,y1n)
plt.plot(x2n,y2n)


plt.subplot(1,2,2)
plt.title(f"Interactions between particles {len(x1y)-1} steps")
plt.xlabel(r'Position along x [$\mu$m]')
plt.ylabel(r'Position along y [$\mu$m]')
plt.plot(x1y,y1y)
plt.plot(x2y,y2y)

plt.savefig("inter_noninter82.pdf")







#Phase-space two particles

positions1 = np.loadtxt('positions83no.txt')
positions2 = np.loadtxt('positions83yes.txt')
velocities1 = np.loadtxt('velocities83no.txt')
velocities2 = np.loadtxt('velocities83yes.txt')


r1n = positions1[:,0]
r2n = positions1[:,1]
r1y = positions2[:,0]
r2y = positions2[:,1]

v1n = positions1[:,0]
v2n = positions1[:,1]
v1y = positions2[:,0]
v2y = positions2[:,1]



x1n = []
z1n = []
x2n = []
z2n = []

x1y = []
z1y = []
x2y = []
z2y = []

vx1n = []
vz1n = []
vx2n = []
vz2n = []

vx1y = []
vz1y = []
vx2y = []
vz2y = []




for i in range(0,len(positions1),3):
    x1n.append(r1n[i])
    z1n.append(r1n[i+2])
    x2n.append(r2n[i])
    z2n.append(r2n[i+2])
    
    vx1n.append(v1n[i])
    vz1n.append(v1n[i+2])
    vx2n.append(v2n[i])
    vz2n.append(v2n[i+2])
    
    

for i in range(0,len(positions2),3):
    x1y.append(r1y[i])
    z1y.append(r1y[i+2])
    x2y.append(r2y[i])
    z2y.append(r2y[i+2])
    
    vx1y.append(v1n[i])
    vz1y.append(v1n[i+2])
    vx2y.append(v2n[i])
    vz2y.append(v2n[i+2])



plt.figure(figsize=(14,7))
plt.subplot(1,2,1)
plt.title(f"No interactions between particles {len(x1n)-1} steps")
plt.xlabel(r'Position along x [$\mu$m]')
plt.ylabel(r'Velocity along x [m/s]')
plt.plot(x1n,vx1n)
plt.plot(x2n,vx2n)


plt.subplot(1,2,2)
plt.title(f"Interactions between particles {len(x1y)-1} steps")
plt.xlabel(r'Position along x [$\mu$m]')
plt.ylabel(r'Velocity along x [m/s]')
plt.plot(x1y,vx1y)
plt.plot(x2y,vx2y)

plt.savefig("x_phase_space83.pdf")

plt.figure(figsize=(14,7))
plt.subplot(1,2,1)
plt.title(f"No interactions between particles {len(x1n)-1} steps")
plt.xlabel(r'Position along z [$\mu$m]')
plt.ylabel(r'Velocity along z [m/s]')
plt.plot(z1n,vz1n)
plt.plot(z2n,vz2n)


plt.subplot(1,2,2)
plt.title(f"Interactions between particles {len(x1y)-1} steps")
plt.xlabel(r'Position along z [$\mu$m]')
plt.ylabel(r'Velocity along z [m/s]')
plt.plot(z1y,vz1y)
plt.plot(z2y,vz2y)

plt.savefig("z_phase_space83.pdf")






#3D plot two particles
positions1 = np.loadtxt('positions82no.txt')
positions2 = np.loadtxt('positions82yes.txt')

r1n = positions1[:,0]
r2n = positions1[:,1]
r1y = positions2[:,0]
r2y = positions2[:,1]

x1n = []
y1n = []
z1n = []
x2n = []
y2n = []
z2n = []

x1y = []
y1y = []
z1y = []
x2y = []
y2y = []
z2y = []



for i in range(0,len(positions1),3):
    x1n.append(r1n[i])
    y1n.append(r1n[i+1])
    z1n.append(r1n[i+2])
    x2n.append(r2n[i])
    y2n.append(r2n[i+1])
    z2n.append(r2n[i+2])

for i in range(0,len(positions2),3):
    x1y.append(r1y[i])
    y1y.append(r1y[i+1])
    z1y.append(r1y[i+2])
    x2y.append(r2y[i])
    y2y.append(r2y[i+1])
    z2y.append(r2y[i+2])



fig = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.plot3D(x1n, y1n, z1n)
ax.plot3D(x2n, y2n, z2n)
ax.set_title(f"No interactions between particles {len(x1n)-1} steps")
plt.savefig("3d_no_interactions_84.pdf")

fig = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.plot3D(x1y, y1y, z1y)
ax.plot3D(x2y, y2y, z2y)
ax.set_title(f"Interactions between particles {len(x1n)-1} steps")
plt.savefig("3d_interactions_84.pdf")







#The trap resonance
plt.figure(figsize=(8,8))
plt.xlabel(r'Fluctuation frequency $\omega_{V}$')
plt.ylabel(r'Fraction of particles in trap after 500 $\mu$s')
plt.ylim(-0.1,1.3)


filename_list = ["resonance_frequancy0.200000.txt", "resonance_frequancy0.400000.txt", "resonance_frequancy0.700000.txt"]
label_list = ["f = 0.2", "f = 0.4", "f = 0.7"]

for i in range(3):
    positions = np.loadtxt(filename_list[i])

    w_V = positions[:,0]
    particle_fraction = positions[:,1]

    plt.plot(w_V, particle_fraction, label=label_list[i])

plt.legend()
plt.savefig('resonance_all.pdf')
plt.show()





#The trap resonance zoom in
plt.figure(figsize=(14,7))
plt.subplot(1,2,1)
plt.title('Without interactions')
plt.xlabel(r'Fluctuation frequency $\omega_{V}$')
plt.ylabel(r'Fraction of particles in trap after 500 $\mu$s')
plt.ylim(-0.1,1.3)


filename_list = ["resonance_frequancy_zoom_no_0.2.txt", "resonance_frequancy_zoom_no_0.4.txt", "resonance_frequancy_zoom_no_0.7.txt"]
label_list = ["f = 0.2", "f = 0.4", "f = 0.7"]

for i in range(3):
    positions = np.loadtxt(filename_list[i])

    w_V = positions[:,0]
    particle_fraction = positions[:,1]

    plt.plot(w_V, particle_fraction, label=label_list[i])

plt.legend()




plt.subplot(1,2,2)
plt.title('With interactions')
plt.xlabel(r'Fluctuation frequency $\omega_{V}$')
plt.ylabel(r'Fraction of particles in trap after 500 $\mu$s')
plt.ylim(-0.1,1.3)


filename_list = ["resonance_frequancy_zoom_yes_0.2.txt", "resonance_frequancy_zoom_yes_0.4.txt", "resonance_frequancy_zoom_yes_0.7.txt"]
label_list = ["f = 0.2", "f = 0.4", "f = 0.7"]

for i in range(3):
    positions = np.loadtxt(filename_list[i])

    w_V = positions[:,0]
    particle_fraction = positions[:,1]

    plt.plot(w_V, particle_fraction, label=label_list[i])

plt.legend()


plt.savefig('resonance_narrow.pdf')
plt.show()





