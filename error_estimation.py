import numpy as np
import matplotlib.pyplot as plt




"""
Analytic
"""

T = 50.0 #micro s
V0 = 2.41 * 10**6 #u micro m^2/(micro s ^2 e)
B0 = 96.5 #u/(micro s e)
d = 500.0 #micro m

m = 40 #u
q = 1 #e

x0 = 20 #micro m
z0 = 20 #micro m
v0 = 25 #m/s



def analytic(steps):
    
    t = np.linspace(0,50,steps+1) #each simulation have steps+1 points due to initial vlaues
    
    w_z = np.sqrt(2*q*V0/(m*d**2))
    w_0 = q*B0/m
    w_plus = (w_0 + np.sqrt(w_0**2 - 2 * w_z**2))/2
    w_minus = (w_0 - np.sqrt(w_0**2 - 2 * w_z**2))/2
    
    A_plus = (v0 + w_minus * x0)/(w_minus - w_plus)
    A_minus = - (v0 + w_plus * x0)/(w_minus - w_plus)
    
    x = A_plus * np.cos(w_plus * t) + A_minus * np.cos(w_minus * t) #as phases are zero
    y = - A_plus * np.sin(w_plus * t) - A_minus * np.sin(w_minus * t)
    z = z0 * np.cos(w_z * t)
    
    return x, y, z, t




x4, y4, z4, t4 = analytic(4000)
x8, y8, z8, t8 = analytic(8000)
x16, y16, z16, t16 = analytic(16000)
x32, y32, z32, t32 = analytic(32000)

analytic_list = [[x4,y4,z4,t4],[x8,y8,z8,t8],[x16,y16,z16,t16],[x32,y32,z32,t32]]



"""
Numerical
"""
"""
This could be done more compactly using for loops, but since we only have
four different lengdths, we may just as well write it out
"""
positions4RK4 = np.loadtxt('positions85_4000RK4.txt').T
positions8RK4 = np.loadtxt('positions85_8000RK4.txt').T
positions16RK4 = np.loadtxt('positions85_16000RK4.txt').T
positions32RK4 = np.loadtxt('positions85_32000RK4.txt').T

positions4FE = np.loadtxt('positions85_4000FE.txt').T
positions8FE = np.loadtxt('positions85_8000FE.txt').T
positions16FE = np.loadtxt('positions85_16000FE.txt').T
positions32FE = np.loadtxt('positions85_32000FE.txt').T



RK4_list = [positions4RK4, positions8RK4, positions16RK4, positions32RK4]
FE_list = [positions4FE, positions8FE, positions16FE, positions32FE]



plt.figure(figsize=(14,7))
plt.subplot(1,2,1)
plt.title(f"Runge Kutta relative error")
plt.ylabel(r'Relative error []')
plt.xlabel(r'Time [$\mu$s]')

"""
Inneficent as a lot of alues may be deleted when we will not be using them anymore
"""


ri_RK4_list = []

for i in range(4):
    
    rRK4_rel = np.zeros(0) #relative error
    rRK4_abs = np.zeros(0) #absolute error
    k = 0
    
    for j in range(0,len(RK4_list[i]),3):
        rxRK4 = RK4_list[i][j] - analytic_list[i][0][k] #the relative error in x
        ryRK4 = RK4_list[i][j+1] - analytic_list[i][1][k]
        rzRK4 = RK4_list[i][j+2] - analytic_list[i][2][k]
        
        r_norm = np.sqrt(analytic_list[i][0][k]**2 + analytic_list[i][1][k]**2  + analytic_list[i][2][k]**2)
        
        ri = np.sqrt(rxRK4**2 + ryRK4**2  + rzRK4**2) #absolute value of the error vector

        rRK4_rel = np.append(rRK4_rel, ri/r_norm) #relative error
        rRK4_abs = np.append(rRK4_abs, ri) #absolute error
        
        k += 1

    ri_RK4_list.append(rRK4_abs)
    plt.plot(analytic_list[i][3], rRK4_rel, label=f"{(len(RK4_list[i])-3)//3} steps")
    

plt.legend()
plt.subplot(1,2,2)
plt.title(f"Forward Euler relative error")
plt.ylabel(r'Relative error []')
plt.xlabel(r'Time [$\mu$s]')


ri_FE_list = []

for i in range(4):
    
    rFE_rel = np.zeros(0) #relative error
    rFE_abs = np.zeros(0) #absolute error
    k = 0
    
    for j in range(0,len(RK4_list[i]),3):

        rxFE = FE_list[i][j] - analytic_list[i][0][k] #error vectors
        ryFE = FE_list[i][j+1] - analytic_list[i][1][k]
        rzFE = FE_list[i][j+2] - analytic_list[i][2][k]
        
        r_norm = np.sqrt(analytic_list[i][0][k]**2 + analytic_list[i][1][k]**2  + analytic_list[i][2][k]**2)
        
        ri = np.sqrt(rxFE**2 + ryFE**2  + rzFE**2) #absolute value of the error vector
        
        rFE_rel = np.append(rFE_rel, ri/r_norm) #relative error
        rFE_abs = np.append(rFE_abs, ri) #absolute error
        
        k += 1

    ri_FE_list.append(rFE_abs)
    plt.plot(analytic_list[i][3], rFE_rel, label=f"{(len(RK4_list[i])-3)//3} steps")
    
    
plt.legend()
plt.savefig("error_comparason85.pdf")


"""
Error convergance rate for Runge Kutta 4
"""

RK4_max_values = [] #Here we will store highest values
at_step = 0 

step_size_list = []

for element in ri_RK4_list:

    step_size_list.append(analytic_list[at_step][3][1])

    max_val_loc = max(element) #highest value identical to delta max
        
    RK4_max_values.append(max_val_loc) #store the highest value
    
    at_step += 1

    
r_err = 0 #Error convargence rate
for i in range(3):
    r_err += 1.0/3.0 * np.log10(RK4_max_values[i+1]/RK4_max_values[i])\
        /np.log10(step_size_list[i+1]/step_size_list[i])
            
print(r_err)
        





print('\n')

"""
Error convergance rate for Forward Euler
"""

FE_max_values = [] #Here we will store highest values
at_step = 0 

step_size_list = []

for element in ri_FE_list:

    step_size_list.append(analytic_list[at_step][3][1])

    max_val_loc = max(element) #highest value identical to delta max
        
    FE_max_values.append(max_val_loc) #store the highest value
    
    at_step += 1

    
r_err = 0 #Error convargence rate
for i in range(3):
    r_err += 1.0/3.0 * np.log10(FE_max_values[i+1]/FE_max_values[i])\
        /np.log10(step_size_list[i+1]/step_size_list[i])
            
print(r_err)



