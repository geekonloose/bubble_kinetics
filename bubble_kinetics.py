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
from matplotlib import pyplot as plt
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

number = int(1E3)
dt = 1E-2

#%% data for rate of temperature change

T = np.zeros(number)
Db = np.zeros(number)
A = np.zeros(number)
mass_of_vap_sod = 1E-3 # mass of vaporised sodium in kg
Tb = 1154.1 # Boiling point of sodium in kelvine
Ts = 800 # sodium pool temperature
cp = 1284.4 #kJ/kg

lg = 4197E3


#%% Data for change in velocity

v = np.zeros(number) #bubble rising velocity m/s
t = np.zeros(number)
y = np.zeros(number)
rhol = np.zeros(number)
rhog = np.zeros(number)


#%% Data for bubble

Pb = np.zeros(number)
k = 1.380648813131313131E-23 #Boltzman constant J/K
g = 9.8 # gravitational acceleration
P0 = 1E5 # Atmosperic pressure 1 bar = 1E5 pascals
H = 13

#%% Input Parameter Needed

N = 0.23E23 #chack
h = 1E1 #check
Pb[0] = 1E5 #check
T[0]= 3000 #check
rhog[:] = 200 #check the value
rhol[0]=949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2
print('sodium density \t{:0.3e}'.format(rhol[0]))
Db[0] =N*k*T[0]/(3.14*Pb[0])
print('bubble diameter \t{:0.3e}'.format(Db[0]*1E3))
#Db[0] = 2
v[0] = 0 #check value
for i in range(number-3):

    t[i] = i*dt
    y[i] = (i)*dt*v[i]

    Pb[i] = P0 + rhog[i]*g*(H-y[i])
    Db[i] = (N*k*T[i]/(3.14*Pb[i]))**(1/3.0)



    T[i+1] = T[i] + dt* (0-h*3.14*Db[i]**2 * (T[i]-Ts))

    rhol[i+1] = 949-0.223*(T[i+1]-273.15)-1.7E-5*(T[i+1]-273.15)**2

    v[i+1] = v[i] + dt*(6/(rhol[i]*3.14*Db[i]**3))*(  (3.14*Db[i]**3 *g*(rhol[i]-rhog[i]))/6 - 12*3.14*Db[i]*v[i])
plt.figure()
plt.plot(t[:-4],Db[:-4],label='diameter',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('diameter (meter)',color='blue',size=14)
plt.legend()
plt.savefig('diameter')
plt.figure()
plt.plot(t[:-4],T[:-4],label = 'Temperature',color='blue')
plt.xlabel('time',color='blue',size=14)
plt.ylabel('Temp (K)',color='blue',size=14)
plt.legend()
plt.savefig('temperature')
plt.figure()
plt.plot(t[:-4],Pb[:-4],label='pressure',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('Pressure (Pa)',color='blue',size=14)
plt.legend()
plt.savefig('pressure')
plt.figure()
plt.plot(t[:-4],v[:-4],label='velocity',color='blue')
plt.xlabel('time(second)',color='blue',size=14)
plt.ylabel('velocity(m/s)',color='blue',size=14)

plt.legend()
plt.savefig('velocity')
#convert meter to cm
Db =(Db*1E2)
y = y*1E2-250
#convert array to list
Db = list(Db)
y = list(y)
import turtle

def draw_circle(turtle, color, size, x, y):
    turtle.penup()
    turtle.color(color)

    turtle.goto(x,y)
    turtle.pendown()

    turtle.circle(size)


turtle.shape("circle")
#turtle.resizemode("user")
turtle.shapesize(1E-2,1E-2)
#tommy.speed(7)
scale = 100
# Draw three circles:
for i in range(int(len(y)/scale)):
    draw_circle(turtle, "green", Db[scale*i], 25, y[scale*i])







