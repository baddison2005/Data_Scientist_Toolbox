import scipy.stats
import math
import sys
import os
import re
import glob
import string
import fnmatch
import time
import shutil
import random
import subprocess
import numpy as np

directory_input = '/Volumes/Seagate Backup Plus Drive/Cyclops/RM/'
parameter_filename = 'parameters.txt'

#Lets open up the parameter file and determine the location of the files we are interested in.
input_data_filename = directory_input + parameter_filename

param = np.loadtxt(input_data_filename,dtype='str',delimiter=',')
#parameter_input = open(input_data_filename, "r")
input_RM_filename=param[0].strip() #param(150)
input_my_data_filename=param[1].strip()
input_other_data_filename=param[2].strip()
num_my_rv=int(param[3].strip())
num_other_rv=int(param[4].strip())
directory_location=param[5].strip()
Planet_name=param[6].strip()
my_data_file_name=param[7].strip()
output_temp_data=param[8].strip()
other_RV_files=param[9].strip()
Jupiter_Earth_units=param[10].strip()
Phase_shift=param[11].strip()
subtract_JD_mid_time=param[12].strip()
add_RV_zero_offset=param[13].strip()
data_plot_model_interval=float(param[14].strip())
Phase_angle_start=float(param[15].strip())
Ms_solar=float(param[16].strip())
Inc_prior=float(param[17].strip())
Inc_end=float(param[18].strip())
Inc_begin=float(param[19].strip())
Inc_interval=float(param[20].strip())
Inc_1sigerr=float(param[21].strip())
Inc_fit=param[22].strip()
Number_fit=int(param[23].strip())
Number_orbits=float(param[24].strip())
Number_iterations_fit=int(param[25].strip())
Number_iterations_total=int(param[26].strip())
datafilelength=int(param[27].strip())
omega_arg_periastron_prior=float(param[28].strip())
omega_arg_periastron_end=float(param[29].strip())
omega_arg_periastron_begin=float(param[30].strip())
omega_arg_periastron_interval=float(param[31].strip())
omega_arg_periastron_1sigerr=float(param[32].strip())
omega_arg_periastron_fit=param[33].strip()
Ecc_prior=float(param[34].strip())
Ecc_end=float(param[35].strip())
Ecc_begin=float(param[36].strip())
Ecc_interval=float(param[37].strip())
Ecc_1sigerr=float(param[38].strip())
Ecc_fit=param[39].strip()
Mp_prior=float(param[40].strip())
Mp_end=float(param[41].strip())
Mp_begin=float(param[42].strip())
Mp_interval=float(param[43].strip())
Mp_1sigerr=float(param[44].strip())
Mp_fit=param[45].strip()
Rp_prior=float(param[46].strip())
Rp_end=float(param[47].strip())
Rp_begin=float(param[48].strip())
Rp_interval=float(param[49].strip())
Rp_1sigerr=float(param[50].strip())
Rp_fit=param[51].strip()
Rs_solar=float(param[52].strip())
Orbital_period_prior=float(param[53].strip())
Orbital_period_end=float(param[54].strip())
Orbital_period_begin=float(param[55].strip())
Orbital_period_interval=float(param[56].strip())
Orbital_period_1sigerr=float(param[57].strip())
Orbital_period_fit=param[58].strip()
JD_time_mid_transit_prior=float(param[59].strip())
JD_time_mid_transit_end=float(param[60].strip())
JD_time_mid_transit_begin=float(param[61].strip())
JD_time_mid_transit_interval=float(param[62].strip())
JD_time_mid_transit_1sigerr=float(param[63].strip())
JD_time_mid_transit_fit=param[64].strip()
vsini_prior=float(param[65].strip())
vsini_end=float(param[66].strip())
vsini_begin=float(param[67].strip())
vsini_interval=float(param[68].strip())
vsini_1sigerr=float(param[69].strip())
vsini_fit=param[70].strip()
Number_vsini_mc=param[71].strip()
stellar_rotation_angle_prior=float(param[72].strip())
stellar_rotation_angle_end=float(param[73].strip())
stellar_rotation_angle_begin=float(param[74].strip())
stellar_rotation_angle_interval=float(param[75].strip())
stellar_rotation_angle_1sigerr=float(param[76].strip())
stellar_rotation_angle_fit=param[77].strip()
Number_stellar_rotation_angle_mc=param[78].strip()
vmacro=float(param[79].strip())
vmicro=float(param[80].strip())
RV_zero_offset_prior=float(param[81].strip())
RV_zero_offset_end=float(param[82].strip())
RV_zero_offset_begin=float(param[83].strip())
RV_zero_offset_interval=float(param[84].strip())
RV_zero_offset_1sigerr=float(param[85].strip())
RV_zero_offset_fit=param[86].strip()
RV_offset_datasets_prior=float(param[87].strip())
RV_offset_datasets_end=float(param[88].strip())
RV_offset_datasets_begin=float(param[89].strip())
RV_offset_datasets_interval=float(param[90].strip())
RV_offset_datasets_new_1sigerr=float(param[91].strip())
RV_offset_datasets_fit=param[92].strip()
Time_bef_aft=float(param[93].strip())
Time_compare_vel=float(param[94].strip())
M_pixels=int(param[95].strip())
Albedo=float(param[96].strip())
linear_quadratic=param[97].strip()
u=float(param[98].strip())
q_1=float(param[99].strip())
q_2=float(param[100].strip())
Bessel_function_exit=float(param[101].strip())
chi_squared_change=float(param[102].strip())
chi_squared_change_fit=float(param[103].strip())
num_degrees_freedom_fit=int(param[104].strip())
num_degrees_freedom_total=int(param[105].strip())
mc_sample_size=int(param[106].strip())
my_data_plotting_symbols=param[107].strip()
other_data_plotting_symbols=param[108].strip()

param2 = np.loadtxt(output_temp_data + 'net_degrees_freedom.txt',dtype='int')
net_num_degrees_freedom_total = param2

sigma = 1.0

#confidence intervals these sigmas represent:
conf_int = scipy.stats.chi2.cdf(sigma**2,1)

#The delta chi square value.
delta_chi_squared = scipy.stats.chi2.ppf(conf_int, net_num_degrees_freedom_total)

#Create a file and write out the delta chi squared value.
chi_squared_out = open(output_temp_data + 'net_delta_chi_square.txt', 'w')
chi_squared_out.write('%6.3f' % (delta_chi_squared))
chi_squared_out.close()