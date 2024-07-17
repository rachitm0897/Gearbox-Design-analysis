import math
torque = 100 #N-m
# table = []
Gear_Ratio = 2
k = 1
power = 100 #kW
angular_speed = 8500 #RPM
phi = 20*0.0174533

import math

ratio = Gear_Ratio
phi = 20
def min_number_teeth(m_g,phi):
    N_p = (2*(m_g + math.sqrt(m_g**2 + (1+2*m_g)*(math.sin(phi*0.0174533)**2))))/((1+2*m_g)*(math.sin(phi*0.0174533)**2))
    return N_p
def max_number_teeth(N_p,phi):
    N_g = ((N_p**2)*(math.sin(phi*0.0174533)**2)-4)/(4-2*N_p*math.sin(phi*0.0174533)**2)
    return N_g

allowable_ratio = []
N_p_min = math.floor(min_number_teeth(Gear_Ratio,phi))
N_g_max = math.ceil(max_number_teeth(min_number_teeth(Gear_Ratio,phi),phi))
print(N_p_min,N_g_max)

N = N_p_min
N_g = N_g_max

import numpy as np
import pandas as pd
from scipy.interpolate import griddata

df = pd.read_csv(r"C:\Users\91987\Downloads\Book1.csv")
x = df['N1'].values
y = df['N2'].values
z = df['Y_j'].values

xi, yi = np.meshgrid(np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100))
zi = griddata((x, y), z, (xi, yi), method='linear')

xi_ext, yi_ext = np.meshgrid(np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100))
zi_ext = griddata((x, y), z, (xi_ext, yi_ext), method='linear')


for module in [1,1.5,2,2.5,3,3.5,4,4.5,5]:
    for N_p in range(N,int((N_g/ratio)+1)):
        for times in [3,4,5]:
            S_c = ((2.41*316) +237) # Contact Strength N/mm^2 
            # S_h =   # AGMA safety Factor
            cycle = 10**7
            Z_n = 1.249*(cycle**(-0.0138))  # Stress cycyle factor for contact strength 
            Y_0 = 1 # Temperature factor
            Y_z = 1.25  # Reliability factor
            Z_w = 1 # Hardness ratio factor
            S_f_all = (S_c*Z_n*Z_w)/(Y_0*Y_z) # Sigma (Surface) Allowed

            d_p = module*N_p
            r_p = d_p/2
            W_t = power*60000/(math.pi*d_p*angular_speed)
            K_s = 1
            K_o = 2.25 # overload factor 
            K_h = 1.6 # Load distribution Factor
            Q_v = 9 # AGMA Quality number 
            B = 0.25*((12-Q_v)**(2/3))
            A = 50 + 56*(1-B)
            V = math.pi*d_p*angular_speed*0.4101049 #pitch line velocity m/s
            K_v = ((A+math.sqrt(200*V))/A)**B #Dynamic Factor
            F = times*math.pi*module
            Z_r = 1 # surface condition factor 
            Z_e = 2300 # Elastic coefficient of steel
            Z_i = (math.sin(phi)*math.cos(phi)*ratio)/(2*(ratio+1))

            S_f = Z_e*math.sqrt((W_t*K_o*K_v*K_s*K_h*Z_r)/(d_p*F*Z_i))

            SF = S_f_all/S_f

            #Bending 

            k_b = 1
            Y_j = zi[int(N_p),int(N_p*ratio)]
            if Y_j == "nan":
                Y_j = zi_ext[int(N_p),int(N_p*ratio)]

            S_b = W_t*K_o*K_v*K_s*K_h*k_b/(F*module*Y_j)

            S_t = (0.533*311) +88.3
            Y_n = 6.1514*(cycle**(-0.1192))

            S_b_all = S_t*Y_n /(Y_0*Y_z)
            
            if S_b_all/S_b > 1.5 :
                if SF > 1.5:
                    print(N_p,N_p*ratio,module,times,int(SF),module*N_p,module*N_g,F,S_b_all/S_b)
