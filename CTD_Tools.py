#Import Dependancies
import csv
import numpy as np
from matplotlib import pyplot as plt
import math as m
import matplotlib
from datetime import date, timedelta
import datetime
import pandas as pd
import seaborn as sns
import plotly.express as px


def read_CTD_data(file_name, profiler):
    '''
    read_CTD_data(filename, profiler) - reads CTD Data from csvs

    INPUTS:
    file_name (str) : local path to csv data
    profile (str) : determines which headers to use
        'deep' - reading file from deep profiler
        'shallow' - reading file from shallow profiler

    OUTPUT:
    pressures_MPa, temps_C, salinities_ppt, data_time
    '''

    if profiler == 'deep':
        pressures_dbar = np.array(pd.read_csv(file_name,usecols=['pressure']))
        temps_C = np.array(pd.read_csv(file_name,usecols=['temp']))
        salinities_ppt = np.array(pd.read_csv(file_name,usecols=['practical_salinity']))

        # Convert to float instead of string
        pressures_dbar = pressures_dbar[1:].astype(np.float)
        temps_C = temps_C[1:].astype(np.float)
        salinities_ppt = salinities_ppt[1:].astype(np.float)

        pressures_MPa = pressures_dbar*0.01

    if profiler == 'shallow':
        pressures_dbar = np.array(pd.read_csv(file_name,usecols=['seawater_pressure']))
        temps_C = np.array(pd.read_csv(file_name,usecols=['seawater_temperature']))
        salinities_ppt = np.array(pd.read_csv(file_name,usecols=['practical_salinity']))

        # Convert to float instead of string
        pressures_dbar = pressures_dbar[1:].astype(np.float)
        temps_C = temps_C[1:].astype(np.float)
        salinities_ppt = salinities_ppt[1:].astype(np.float)

        pressures_MPa = pressures_dbar*0.01
    
    time = np.array(pd.read_csv(file_name,usecols=['time']))

    data_time = []
    for k in range(time.shape[0]):
        ref_date = datetime.datetime(1900,1,1)      # This is the "days since" part
        time_since_ref = timedelta(seconds=int(time[k]))     # Create a time delta object from the number of days
        data_time.append(ref_date + time_since_ref)      # Add the specified number of days to 1990

    return [pressures_MPa, temps_C, salinities_ppt, data_time]

def calc_ssp(data_list):
    # unpack data_list
    pressures_MPa = data_list[0]
    temps_C = data_list[1]
    salinities_ppt = data_list[2]
    data_time = data_list[3]

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

    return depth_m, c

def plot_ssp(depth_m, c, plot_title, save_file=False, file_name=None, ):

    sns.set()

    fig = plt.figure(figsize=(7,5))

    fig.suptitle('Depth of Ocean vs. Speed of Sound', fontsize=16, fontweight='bold')

    ax = fig.add_subplot(111)
    ax.set_title(plot_title ,fontsize=12)
    ax.plot(c,depth_m,'.')
    fig.subplots_adjust(top=0.87)
    ax.set_xlabel('Speed of Sound $(m/s)$')
    ax.set_ylabel('Ocean Depth (m)')
    plt.gca().invert_yaxis()

    plt.show()

    if save_file:
        if file_name == None:
            print('Must Provide Filename')
        else:
            fig.savefig(file_name, dpi=500)
    return

def manually_cut_to_8_passes(depth_m_shallow, c_shallow):
    # Cut Shallow data to just 8 passes
    depth_m_shallow = depth_m_shallow[1197:35467]
    c_shallow = c_shallow[1197:35467]
    return depth_m_shallow, c_shallow

def combine_ssps(depth_shallow, speed_shallow, depth_deep, speed_deep):
    c_combined = np.concatenate((speed_deep,speed_shallow))
    depth_combined = np.concatenate((depth_deep,depth_shallow))
    return depth_combined, c_combined

def sort_ssp(depth_combined, c_combined):  
    ssp = np.hstack((depth_combined, c_combined))

    #Sort Data By Depth
    ssp = ssp[ssp[:,0].argsort()]

    depth_sort = ssp[:,0]
    speed_sort = ssp[:,1]
    return depth_sort, speed_sort

def manually_cut_problems(depth, speed):
    ssp = np.vstack((depth, speed))
    #Cut out Redundant data points (hard coded)
    ssp_cut = np.delete(ssp, np.s_[26053:36168], 1)

    ssp_cut = np.delete(ssp_cut, np.s_[7072:20193], 1)
    depth_cut = ssp_cut[0,:]
    speed_cut = ssp_cut[1,:]

    return depth_cut, speed_cut

def calc_moving_average(depth, speed, N=700):

    depth_ma = np.convolve(depth, np.ones((N,))/N, mode='valid')
    speed_ma = np.convolve(speed, np.ones((N,))/N, mode='valid')

    return depth_ma, speed_ma

def get_ssp_SlopeBase():

    file_name_deep = 'CSVs/oregon_slope_base_deep_profiler_CTD_20190714_20190715.csv'
    file_name_shallow = 'CSVs/oregon_slope_base_shallow_profiler_CTD_20190714_20190715.csv'

    profiler1 = 'deep'
    profiler2 = 'shallow'

    data_out_deep = read_CTD_data(file_name_deep, profiler1)
    data_out_shallow = read_CTD_data(file_name_shallow, profiler2)

    depth_shallow, speed_shallow = calc_ssp(data_out_shallow)
    depth_deep, speed_deep = calc_ssp(data_out_deep)

    #Manually Cut Shallow
    depth_shallow, speed_shallow = manually_cut_to_8_passes(depth_shallow, speed_shallow)

    # Combine 2 ssps
    depth, speed = combine_ssps(depth_shallow, speed_shallow, depth_deep, speed_deep)

    #Sort 2 ssps
    depth_sort, speed_sort = sort_ssp(depth, speed)

    # Manually remove problem data
    depth_cut, speed_cut = manually_cut_problems(depth_sort, speed_sort)

    # Calculate Moving Average
    depth_ma, speed_ma = calc_moving_average(depth_cut, speed_cut)

    # Plot Sound Speed Profile
    plot_ssp(depth_ma, speed_ma, 'Oregon Slope Base')

    ssp = np.vstack((depth_ma, speed_ma)).T
    return ssp
   
def get_ssp_Offshore():

    file_name_deep = 'CSVs/oregon_Offshore_deep_profiler_CTD_20190714_20190715.csv'
    file_name_shallow = 'CSVs/oregon_Offshore_shallow_profiler_CTD_20190714_20190715.csv'

    profiler1 = 'deep'
    profiler2 = 'shallow'

    data_out_deep = read_CTD_data(file_name_deep, profiler1)
    data_out_shallow = read_CTD_data(file_name_shallow, profiler2)

    depth_shallow, speed_shallow = calc_ssp(data_out_shallow)
    depth_deep, speed_deep = calc_ssp(data_out_deep)

    # Combine 2 ssps
    depth, speed = combine_ssps(depth_shallow, speed_shallow, depth_deep, speed_deep)

    #Sort 2 ssps
    depth_sort, speed_sort = sort_ssp(depth, speed)


    depth_ma, speed_ma = calc_moving_average(depth_sort, speed_sort)


    # Plot Sound Speed Profile
    plot_ssp(depth_ma, speed_ma, 'Oregon Offshore')

    ssp = np.vstack((depth_ma, speed_ma)).T

    return ssp