# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 17:10:21 2017

@author: parth
"""

import numpy as np
import math
from matplotlib import pyplot

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


for i in range(1,10):

    print(i)
    Db.append(Db[0]**3 * (P + (g*rhol*H)/1033.3)/(P+(g*rhol*(H-x[i]))/1033.3)) #bubble diameter

    sigmal = ((Pi-Po)*Db[i])/8 #Surface tension of liquid

    M = g*mul**4/(rhol*sigmal**3) # Morton number

    E0 = g*rhol*Db[i]**2/rhol

    Ut = rhol**2 * Db[i]**4/(576*M**0.7 * mul**2) #terminal velosity of bubble

    Re = Ut*rhol*Db[i]/(mul) # reynolds number

    diffusion = 1.8*(8*theta/(Ut*Db[i]**3)) #diffusion rate of aerosol in bubble

    sedimentation = 1.5*g*tau/(Db[i]*Ut) #sedimentation rate of aerosol in bubble

    inertial_impaction = 18*Ut*tau/Db[i]**2  #Inertial Impaction rate of aerosol in bubble


