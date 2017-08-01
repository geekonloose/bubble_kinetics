# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 15:19:22 2017
@author: Parth
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 17:10:21 2017
@author: parth
"""

import numpy as np
import math
from matplotlib import pyplot
'''
#%% Variable Initiallisation

T = 800 #absulute temperature in kelvin
H = 13 # pool height in m
Na = 6.022E23 # Avogrado number
P = 1E5 # Absolute pressure in pascal
k = 1.3807E-23 #boltzman constant in J/k
dp = 1E-6 #particle diameters
lembda = 82.06*T/(float(sqrt(2))*3.14*Na*P) #mean free path of a gas molecule
c = 1 + (2*lembda/dp) * (1.257 + 0.4 * math.exp(-0.55*dp/lembda)) #cunningham slip correction

dp = 1E-6 #particle diameter
dg = 0.5E-6 #effective diameter fo a gas molecular
mug = 1.87E-5 #gas viscosity
#%% viscosity of sodium
t  = 373-273
Al = 1.2162E-5
Bl = 949 - 0.223*( t)-1.75E-5*(t)**2
Cl = 0.6976
mul = Bl**(1.0/3) * Al * math.exp(Bl*Cl/(t+273.15)) #liquid viscosity
rhol =  949 - 0.223*t - 1.75E-5 * t**2 #density of liquid kg/m3
rhop = 1 #density of particle kg/m3

Db =[1E-3] #bubble diameter
E = 1 #bubble evventricity
g = 9.8 #

Pi = 2E5 #pressure inside bubble
Po = 1E5 #pressure outside bubble

theta = k*T*c/(3*3.14*mug*dp)
tau = rhop*dp**2 *c /(18*mug)


x = [1,2,3,4,5,6,7,8,9,10]
sigmal = 0
#%% Program Starts here
'''


#%% data for rate of temperature change
number = int(1E4)
T = np.zeros(number)
Db = np.zeros(number)
A = np.zeros(number)
mass_of_vap_sod = 1E-3 # mass of vaporised sodium in kg
Tb = 1154.1 # Boiling point of sodium in kelvine
Ts = 800 # sodium pool temperature
cp = 1284.4 #kJ/kg

A[i] = 3.14*Db[i] #surface area of bubble


#%% Data for change in velocity

v = np.zeros(number) #bubble rising velocity m/s
rhol = np.zeros(number)













