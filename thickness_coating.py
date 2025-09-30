#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = "Ilaria Carlomagno"
__license__ = "MIT"
__version__ = "2.0"
__email__   = "ilaria.carlomagno@elettra.eu"


import library_thickness as mylib
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import xraylib as xl

CS_E0 = 'CS_' + str(mylib.E0)
CS_Ef = 'CS_' + str(mylib.Ef)

# Substrate (steel) composition: Fe, Cr, Ni
SS = {}
SS['el']   = [26, 24, 28] # atomic numbers
SS['comp'] = [70, 20, 10] # molar ratios

# Samples composition: Zr, Ti, V
samples_el = [40, 22, 23] # atomic numbers (Z)
sample = {  'S1_1B':{},
            'S1_5B':{},
            'S2_3B':{},
            'S5_7F':{},
            'S7_2B':{},
            'S7_5F':{} } 

# Samples molar ratios of Zr, Ti, V
sample['S1_1B']['el'], sample['S1_1B']['comp'] = samples_el, [20, 35, 45]
sample['S1_5B']['el'], sample['S1_5B']['comp'] = samples_el, [16, 35, 49]
sample['S2_3B']['el'], sample['S2_3B']['comp'] = samples_el, [18, 35, 47]
sample['S5_7F']['el'], sample['S5_7F']['comp'] = samples_el, [21, 34, 45]
sample['S7_2B']['el'], sample['S7_2B']['comp'] = samples_el, [23, 33, 44]
sample['S7_5F']['el'], sample['S7_5F']['comp'] = samples_el, [21, 38, 41]


# here are the experimental data (area of Fe Kalpha emission)
exp_data = {}
exp_inc_angles = np.array([5, 10, 20, 30, 45, 60])
exp_data['data_S1_1B'] = np.multiply([4.17, 4.69, 4.50, 4.05, 3.24, 2.35], 1e5)
exp_data['data_S1_5B'] = np.multiply([2.47, 3.38, 3.60, 3.26, 2.60, 1.78], 1e5)
exp_data['data_S2_3B'] = np.multiply([4.09, 4.63, 4.47, 4.01, 3.22, 2.32], 1e5)
exp_data['data_S5_7F'] = np.multiply([3.52, 4.15, 4.11, 3.73, 2.99, 2.13], 1e5)
exp_data['data_S7_2B'] = np.multiply([5.49, 5.56, 5.08, 4.62, 3.64, 2.70], 1e5)
exp_data['data_S7_5F'] = np.multiply([6.17, 5.99, 5.36, 4.74, 3.83, 2.90], 1e5)
   
def lin_fluo_fit_func(inc_ang, thickness):
    # angle should be in degrees || thickness in cm
    # calculating exit angle and converting all angles to radians
    exit_ang = 90 - inc_ang
    exit_ang = exit_ang * np.pi / 180
    inc_ang = inc_ang * np.pi / 180

    return np.log(mylib.expt_constant / ( substr_density*np.sin(inc_ang) ) / ( substr_CS_E0/np.sin(inc_ang) + substr_CS_Ef/np.sin(exit_ang) ) ) - ( coating_CS_E0/np.sin(inc_ang) + coating_CS_Ef/np.sin(exit_ang) )* coating_density * thickness

thicknesses = [0.25, 0.5, 1]                  # here values are in microns!
thicknesses = np.multiply(thicknesses, 1e-4)  # now thickness is in cm
inc_angles = np.linspace(exp_inc_angles.min(), exp_inc_angles.max(), 50)

print('\tInitialising substrate...')
SS = mylib.simulate_material(SS)
substr_CS_E0  = SS[CS_E0]
substr_CS_Ef  = SS[CS_Ef]
substr_density = SS['density']
    
print('\n\tInitialising samples...')
for key in sample.keys():
    print(key)
    sample[key] = mylib.simulate_material(sample[key])
    coating_CS_E0 = sample[key][CS_E0]
    coating_CS_Ef = sample[key][CS_Ef]
    coating_density = sample[key]['density']

    # calculating and plotting theoretical curves
    theory = {}
    for t in thicknesses:
        point = str(t*1e4)
        theory[point] = []
        for angle in inc_angles:    
            theory[point].append(lin_fluo_fit_func(angle, t))
    for t in theory.keys():
        plt.plot(inc_angles, theory[t], label = t + ' um')

    # now fitting
    print('\t...and now fitting the thickness!')
    key_data = 'data_' + key
    x = exp_inc_angles
    y = exp_data[key_data]
    
    # providing starting value for thickness
    fit_t = 0.5*1e-4
    fit_t, cov = curve_fit(lin_fluo_fit_func, x, np.log(y))

    print('- - - > Sample {0}; thickness = {1} microns \n'. format(key, np.round(fit_t*1e4,3)))
    
    fit_curve = []
    for x in inc_angles:
        fit_curve.append(lin_fluo_fit_func(x, fit_t))
    # plot fit
    plt.plot(inc_angles, fit_curve, "r--", label = 'fit t = {0} um'.format(np.round(fit_t*1e4,3)))
 
    # plot exp data    
    plt.plot(exp_inc_angles, np.log(exp_data[key_data]), "bo", label = 'data')

    plt.xlabel('Incidence angles (deg)')
    plt.ylabel('log(I_fluo_theo) (arb. units)')
    plt.yscale('log')
    plt.title('linearised Fe K_a fluorescence and fit - sample {0}'.format(key))    
    plt.grid()
    plt.legend()
    namefile = 'fit_' + key + '.png'
    plt.savefig(namefile)
    plt.close()

print('END')
