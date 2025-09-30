#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = "Ilaria Carlomagno"
__license__ = "MIT"
__email__ = "ilaria.carlomagno@elettra.eu"

import math
import matplotlib.pyplot as plt
import numpy as np
import xraylib as xl

Ef = 6.4 # fluorescence of Fe Kalpha (keV)
E0 = 10  # primary Xray energy (keV)

i0 = 2.25e9

abs_jump_fact = 0.9
emission_rate_Kline = xl.RadRate(26, xl.KL3_LINE)
Fe_fluo_yield = xl.FluorYield(26, xl.K_SHELL)

# from the previous three, we get the excitation factor
excit_fact = abs_jump_fact * emission_rate_Kline * Fe_fluo_yield

CS_Fe_E0 = xl.CS_Photo(26, 10)
det_accept = 3*50 / ( 4* (math.pi) *30**2)
det_efficiency = 0.995

expt_constant = i0 * excit_fact * CS_Fe_E0 * det_accept * det_efficiency


def compound_CS(material, energy):
    print('\tCalculating cross section at {0} keV...'.format(energy))
    elem_list = material['el']
    mass_fractions = material['mass_fr']

    CS_list = []
    ave_CS = 0
    for el in elem_list:
        CS_list.append(xl.CS_Photo(el, energy)) # returns cross section in cm^2/g)
    
    ave_CS = np.average(CS_list, weights = mass_fractions)
    name = 'CS_' + str(energy)
    material[name] = ave_CS
    return material

def compound_density(material):
    print('\tCalculating density...')
    elem_list = material['el']
    molar_ratios = material['comp']
    density = []
    for el in elem_list:
        density.append(xl.ElementDensity(el))

    comp_density = np.multiply(molar_ratios, density)
    comp_density = np.sum(comp_density)/np.sum(molar_ratios)
    material['density'] = comp_density

    return material

def molar_to_mass_ratio(material):
    print('\tCalculating mass fractions...')
    elem_list = material['el']
    molar_ratios = material['comp']
    mass_fractions = []
    
    denominator = np.multiply(elem_list, molar_ratios)
    denominator = np.sum(denominator)

    for el, i in zip(elem_list, range(len(elem_list))):
        mass_fractions.append( (molar_ratios[i] * xl.AtomicWeight(el)) / denominator )

    material['mass_fr'] = mass_fractions
    return material

def simulate_material(material):
    material = molar_to_mass_ratio(material)
    material = compound_CS(material, E0)  # E0 = incident X-ray energy 
    material = compound_CS(material, Ef)  # Ef = energy of the fluo line 
    material = compound_density(material)
    return material
   
def fluo_int_theo(inc_ang, coating, substr, thickness):
    # angle should be in degrees || thickness in cm
    # calculating exit angle and converting all angles to radians
    exit_ang = 90 - inc_ang
    exit_ang = exit_ang * math.pi / 180
    inc_ang = inc_ang * math.pi / 180
    CS_E0 = 'CS_' + str(E0)
    CS_Ef = 'CS_' + str(Ef)

    expon_fact = -( coating[CS_E0]/np.sin(inc_ang) + coating[CS_Ef]/np.sin(exit_ang) )* coating['density']
    
    I_f_theo = expt_constant / ( substr['density']*np.sin(inc_ang) ) / ( substr[CS_E0]/np.sin(inc_ang) + substr[CS_Ef]/np.sin(exit_ang) ) * np.exp(expon_fact*thickness)

    return I_f_theo
