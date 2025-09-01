#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 14:32:43 2023

@author: corinnekomlodi
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

#question 1 algebra
n = 2.56e19
Z = 336.28e5
sig_r = 1.91e-27
sig_b = 1.79e-26
sig_g = 5e-27

tau_r = n*Z*sig_r
tau_g = n*Z*sig_g
tau_b = n*Z*sig_b
# print(np.exp(-tau_r))
# print(np.exp(-tau_g))
# print(np.exp(-tau_b))

#2
print('2)')
#consider a protostellar disk that accretes onto a solar-type protostar with m and r.
def eq(x):
    return (x/1)**-2.5

initialM, err = integrate.quad(eq, .1, 10)
finalM = 0.4/initialM
print('Mdot at 1 Myr: {}'.format(finalM/1e6))
#1.9E-8 Msun/yr


initialM1, err = integrate.quad(eq, 1, 10)
print('Mdot {}'.format((finalM)*initialM1))

sb = 5.67e-5 #g s-3 K-4
G = 6.67e-8 #cm3 g-1 s-2
Mstar = 1 * 2e33 #g
Rstar = 2 * 7e10 #cm
Rsun = 1 * 7e10 #cm
Mdot = 1.9e-8 * (2e33/3.15e7) #g s-1
r = np.array([1 * 1.5e13, 2.7 * 1.5e13])  #cm

photoa = (3*G*Mstar*Mdot)
photob = (8*np.pi*sb*Rstar**3)
PhotoTeff = ((photoa/photob)**(0.25))*((r/Rstar)**(-0.75))
print('Photospheric Disk Temp at 1 and 2.7: {}'.format(PhotoTeff))
print('e and f: {}'.format(PhotoTeff*4))



#3
print('3)')
Rstar = 2 #solar radius
Teff = 5000 #k
sb = 5.67e-5 #W/m2/K4

Luminosity = (4*np.pi*(Rstar*7e10)**2*sb*Teff**4)/3.8e33
print('solar luminosities {}'.format(Luminosity))
#2.3 Solar Luminosities

wavelength = np.arange(1000,10000010, 10) #angstroms
h = 6.67e-27 #cm2 g s-1
c = 3e10 #cm s-2
K = 1.4e-16 #cm2 g s-2 K-1
BTeff = ((2*c**2*h)/((wavelength*1e-8)**5))*(1/(np.exp((h*c)/((wavelength*(1e-8))*K*Teff))-1))

StelLumCGS = 4*np.pi**2*BTeff*(Rstar*7e10)**2
stelLum = StelLumCGS/3.8e33

fig, ax = plt.subplots(figsize = (15,10))
ax.set_yscale('log')
ax.axvline(x = .73, color = 'darkgrey')
ax.axhline(y = 1.6, color = 'darkgrey')
ax.set_xscale('log')
ax.plot(wavelength/10000, (wavelength/10000/10000)*stelLum, color = "teal")

ax.set_xlabel('wavelength($\mu$m)')
ax.set_ylabel('$L_{*,\lambda}$')
ax.text(.8, 1.7, '(0.73, 1.6)', color = 'k', fontsize = 'x-large')
ax.plot(.73, 1.6, 'ko')

r = np.arange(.05,100.05, .05)
Teffr = ((3*G*Mstar*Mdot)/(8*np.pi*sb*Rstar**3))**.25*(r/Rsun)**-.75
Bteffa = ((2*h*c**2)/((wavelength*1e-8)**5))

BTeffr = np.array([])
for i in range(0, len(Teffr)):
    loop = (((wavelength*(1e-8))*K*Teffr[i])-1)
    exponent = (np.exp((h*c)/loop))
    bTeffr = Bteffa*(1/exponent)
    BTeffr = np.append(BTeffr, bTeffr)
    print(i)
    
BTeffr = np.reshape(BTeffr, [2000, 999901])
FxDensity = np.pi*BTeffr


for i in range(0, 2000):
    delr = .05
    a_r = 2*np.pi*r[i]*delr
    Lcgs = a_r*FxDensity[i:]


LsolLum = Lcgs/3.8e33
summation = np.sum(LsolLum, axis = 0)

ax.plot(wavelength/1000, summation*wavelength/10000/10000)
fig, ax = plt.subplots(figsize = (15,10))
ax.set_yscale('log')
ax.axvline(x = .73, color = 'darkgrey')
ax.axhline(y = 1.6, color = 'darkgrey')
ax.set_xscale('log')
ax.plot(wavelength/10000, (wavelength/10000/10000)*stelLum, color = "teal")

ax.set_xlabel('wavelength($\mu$m)')
ax.set_ylabel('$L_{*,\lambda}$')
ax.plot(range(1000, 10000010, 2000), summation)

print('e')
Mdot10 = 3.81e15 #g s-1
Mstar = 2 * 6.5e8 #g

Teffr = ((3*G*Mstar*Mdot)/(8*np.pi*sb*Rstar**3))**.25*(r/Rsun)**-.75
Bteffa = ((2*h*c**2)/((wavelength*1e-8)**5))

BTeffr = np.array([])
for i in range(0, len(Teffr)):
    loop = (((wavelength*(1e-8))*K*Teffr[i])-1)
    exponent = (np.exp((h*c)/loop))
    bTeffr = Bteffa*(1/exponent)
    BTeffr = np.append(BTeffr, bTeffr)
    print(i)
    
BTeffr = np.reshape(BTeffr, [2000, 999901])
FxDensity = np.pi*BTeffr


for i in range(0, 2000):
    delr = .05
    a_r = 2*np.pi*r[i]*delr
    Lcgs = a_r*FxDensity[i:]


LsolLum2 = Lcgs/3.8e33
summation2 = np.sum(LsolLum2, axis = 0)


fig, ax = plt.subplots(figsize = (15,10))
ax.set_yscale('log')
ax.axvline(x = .73, color = 'darkgrey')
ax.axhline(y = 1.6, color = 'darkgrey')
ax.set_xscale('log')
ax.plot(wavelength/10000, (wavelength/10000/10000)*stelLum, color = "teal")

ax.set_xlabel('wavelength($\mu$m)')
ax.set_ylabel('$L_{*,\lambda}$')
ax.plot(wavelength/1000, summation2*wavelength/10000/10000)

print('f)')
#gap between .2 and 5 AU
r = np.append(np.arange(.05,2.05, .05), np.arange(2, 100.05, .05))

Teffr = ((3*G*Mstar*Mdot)/(8*np.pi*sb*Rstar**3))**.25*(r/Rsun)**-.75
Bteffa = ((2*h*c**2)/((wavelength*1e-8)**5))

BTeffr = np.array([])
for i in range(0, len(Teffr)):
    loop = (((wavelength*(1e-8))*K*Teffr[i])-1)
    exponent = (np.exp((h*c)/loop))
    bTeffr = Bteffa*(1/exponent)
    BTeffr = np.append(BTeffr, bTeffr)
    print(i)
    
BTeffr = np.reshape(BTeffr, [2000, 999901])
FxDensity = np.pi*BTeffr


for i in range(0, 2000):
    delr = .05
    a_r = 2*np.pi*r[i]*delr
    Lcgs = a_r*FxDensity[i:]


LsolLum3 = Lcgs/3.8e33
summation3 = np.sum(LsolLum3, axis = 0)


fig, ax = plt.subplots(figsize = (15,10))
ax.set_yscale('log')
ax.axvline(x = .73, color = 'darkgrey')
ax.axhline(y = 1.6, color = 'darkgrey')
ax.set_xscale('log')
ax.plot(wavelength/10000, (wavelength/10000/10000)*stelLum, color = "teal")

ax.set_xlabel('wavelength($\mu$m)')
ax.set_ylabel('$L_{*,\lambda}$')
ax.plot(wavelength/1000, summation3*wavelength/10000/10000)












