# Sound Speed Profile
# Calculates sound speed profile for given input csv file. Designed for OOI data

import csv
import numpy as np
from matplotlib import pyplot as plt
import math as m
import matplotlib

'''
CSV Files located in ./CSVs
Oregon Slope Base:
    oregon_slope_base_deep_profiler_CTD_20151006_20151017.csv - October 2015
    oregon_slope_base_deep_profiler_CTD_20151221_20151221.csv - December 2015
    oregon_slope_base_deep_profiler_CTD_20190714_20190715.csv - July 2019

'''
file_name_dict = {'Oct2015': './CSVs/oregon_slope_base_deep_profiler_CTD_20151006_20151017.csv', \
                  'Dec2015': './CSVs/oregon_slope_base_deep_profiler_CTD_20151221_20151221.csv', \
                  'Jul2019': './CSVs/oregon_slope_base_deep_profiler_CTD_20190714_20190715.csv'}
file_name = file_name_dict['Oct2015']
with open(file_name) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    conductivities_mS = []
    pressures_dbar = []
    temps_C = []
    salinities_ppt = []
    
    for column in readCSV:       
        pressure_dbar = column[25]
        pressures_dbar.append(pressure_dbar)
        
        temp_C = column [32]
        temps_C.append(temp_C)
        
        salinity_ppt = column[19]
        salinities_ppt.append(salinity_ppt)
        
# Convert to Numpy Array    
pressures_dbar = np.array(pressures_dbar)
temps_C = np.array(temps_C)
salinities_ppt = np.array(salinities_ppt)


# Convert to float instead of string
pressures_dbar = pressures_dbar[1:].astype(np.float)
temps_C = temps_C[1:].astype(np.float)
salinities_ppt = salinities_ppt[1:].astype(np.float)

pressures_MPa = pressures_dbar*0.01


# literally just typing out equation

C00 = 1402.388
A02 = 7.166E-5
C01 = 5.03830
A03 = 2.008E-6
C02 = -5.81090E-2
A04 = -3.21E-8
C03 = 3.3432E-4
A10 = 9.4742E-5
C04 = -1.47797E-6
A11 = -1.2583E-5
C05 = 3.1419E-9
A12 = -6.4928E-8
C10 = 0.153563
A13 = 1.0515E-8
C11 = 6.8999E-4
A14 = -2.0142E-10
C12 = -8.1829E-6
A20 = -3.9064E-7
C13 = 1.3632E-7
A21 = 9.1061E-9
C14 = -6.1260E-10
A22 = -1.6009E-10
C20 = 3.1260E-5
A23 = 7.994E-12

C21 = -1.7111E-6
A30 = 1.100E-10
C22 = 2.5986E-8
A31 = 6.651E-12
C23 = -2.5353E-10
A32 = -3.391E-13
C24 = 1.0415E-12
B00 = -1.922E-2
C30 = -9.7729E-9
B01 = -4.42E-5
C31 = 3.8513E-10
B10 = 7.3637E-5
C32 = -2.3654E-12
B11 = 1.7950E-7
A00 = 1.389
D00 = 1.727E-3
A01 = -1.262E-2
D10 = -7.9836E-6 

T = 3
S = 1
P = 700
T = temps_C
S = salinities_ppt
P = pressures_MPa*10

D = D00 + D10*P 
B = B00 + B01*T + (B10 + B11*T)*P 
A = (A00 + A01*T + A02*T**2 + A03*T**3 + A04*T**4) + (A10 + A11*T + A12*T**2 + A13*T**3 + A14*T**4)*P + (A20 + A21*T + A22*T**2 + A23*T**3)*P**2 + (A30 + A31*T + A32*T**2)*P**3
Cw = (C00 + C01*T + C02*T**2 + C03*T**3 + C04*T**4 + C05*T**5) + (C10 + C11*T + C12*T**2 + C13*T**3 + C14*T**4)*P + (C20 +C21*T +C22*T**2 + C23*T**3 + C24*T**4)*P**2 + (C30 + C31*T + C32*T**2)*P**3


# Calculate Speed of Sound
c = Cw + A*S + B*S**(3/2) + D*S**2

# Calculate Depth from pressure
lat = 44.52757 #deg

# Calculate gravity constant for given latitude
g_phi = 9.780319*(1 + 5.2788E-3*(np.sin(np.deg2rad(lat))**2) + 2.36E-5*(np.sin(np.deg2rad(lat))**4))

# Calculate Depth for Pressure array
depth_m = (9.72659e2*pressures_MPa - 2.512e-1*pressures_MPa**2 + 2.279e-4*pressures_MPa**3 - 1.82e-7*pressures_MPa**4)/(g_phi + 1.092e-4*pressures_MPa)

min_idx = np.argmin(c)
x = c[min_idx]
y = depth_m[min_idx]

fig = plt.figure(figsize=(7,5))

fig.suptitle('Depth of Ocean vs. Speed of Sound', fontsize=16, fontweight='bold')

ax = fig.add_subplot(111)
ax.set_title('Oregon Slope Base: October 6-17, 2015',fontsize=12)
ax.plot(c,depth_m)
fig.subplots_adjust(top=0.87)
ax.set_xlabel('Speed of Sound $(m/s)$')
ax.set_ylabel('Ocean Depth (m)')
plt.gca().invert_yaxis()

ax.text(1485, 2900, 'Minimum speed of sound from data available \nat {0:.4f} meters'.format(y))

ax.plot(x,y,'o',markersize=7)

plt.show()
fig.savefig('Speed_of_sound_profile.png', dpi=500)