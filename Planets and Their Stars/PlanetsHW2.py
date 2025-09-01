#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 19:56:38 2023

@author: corinnekomlodi
"""

import matplotlib.pyplot as plt
import numpy as np

print('2.')
print('a)')
Mstar = 1 * 1.9884e33 #g
Lstar = 2 * 3.828e33 #solLum 
Rstar = 2 * 6.957e10 #cm
Mdot = 2e-8 * (2e33/3.15e7) #g s-1
a = np.array([0.8, 1, 1.2]) * 1.496e+13 #cm
mu = 2
mh =1.67e-24 #g
k = 1.3807e-16 #cm2 g s-2 K-1
G = 6.6743e-8 #cm3 g-1 s-2
h = 6.67e-27 #cm2 g s-1
sb = 5.67e-5 #g s-3 K-4
Rsun = 6.957e10  #cm

Tdisk = (((3*G*Mstar*Mdot)/(8*np.pi*sb*Rstar**3))**.25)*((a/Rstar)**-.75)
SSpeed = ((k*Tdisk)/(mu*mh))**.5 #cm/s
KepFreq = np.sqrt(G*Mstar/a**3) #hz
ScaleH = SSpeed/KepFreq #cm
ScaleHAU = ScaleH/1.496e13 #au

print('Tdisk = {} K'.format(Tdisk))
print('Speed of Sound = {} cm s-1'.format(SSpeed))
print('Keplerian Frequency = {} Hz'.format(KepFreq))
print('Scale height = {} cm'.format(ScaleH))
print('Scale height = {} AU'.format(ScaleHAU))

print('b)')
theta = np.arctan(ScaleHAU[1]/1)
print('center to scale height = {} degrees'.format(np.degrees(theta)))

delta = (ScaleHAU[2]-ScaleHAU[0])
diskflare = np.arctan(delta/.4)
print('Flaring angle of disk = {} degrees'.format(np.degrees((diskflare))))

incangle = diskflare-theta
print('incident angle = {}'.format(np.degrees(incangle)))

F = (Lstar/(4*np.pi*a[1]**2))*np.sin(incangle) 
print('bolometric Flux = {} erg s-1 cm-2'.format(F))

Ftot = (sb*Tdisk**4)+F
Ttot = (Ftot/sb)**.25

print('T total = {}'.format(Ttot))
deltaK = Ttot-Tdisk
percent = (Tdisk/Ttot) * 100

print('delta K for .8, 1, and 1.2 {}'.format(deltaK))
print('percentage change {}'.format(percent))

###################################################################################

print('3.')
r = np.linspace(1.391e+11, 100*1.496e+13, 200) #cm
Mdisk = .1 * 1.989e+33 #grams
#SurfDen = Mdisk/(np.pi*r**2)
SD = Mdisk/(2*np.pi*r*(100*1.496e+13))

fig, ax = plt.subplots(figsize = (15,10))
ax.plot(r/1.496e+13, SD, label = 'Surface Density')
ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_ylabel('Surface Density')
ax.set_xlabel('Radius AU')
ax.legend()

Kr = 5 #cm2 g-1
tau = (Kr*SD)/2

fig, ax = plt.subplots(figsize = (15,10))
ax.plot(r/1.496e+13, tau, label = 'Optical Depth')
ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_ylabel('Optical Depth')
ax.set_xlabel('Radius AU')
ax.legend()

TT = ((3*tau)/4)**.25

fig, ax = plt.subplots(figsize = (15,10))
ax.plot(r/1.496e+13, TT, label = 'Tmid/Teff')
ax.plot(np.where(max(TT)), max(TT))
ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_ylabel('Tmid/Teff')
ax.set_xlabel('Radius AU')
ax.legend()


Teff = ((3*G*Mstar*Mdot)/(8*np.pi*sb*Rstar**3))**.25*(r/Rstar)**-.75
Tmid = TT**.25 * Teff

fig, ax = plt.subplots(figsize = (15,10))
ax.plot(r/1.496e+13, Teff, label = 'Teff')
ax.plot(r/1.496e+13, Tmid, label = 'Tmid')
ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_ylabel('Temperature K')
ax.set_xlabel('Radius AU')
ax.legend()


Tdisk2 = r**(-3/4) #K
Cs = ((k*Tmid)/(mu*mh))**.5
Omega = (G*Mstar/r**3)**.5
scaleheight = Cs/Omega
Vk = (G*Mstar/r)**.5
Mach = Vk/Cs

fig, ax = plt.subplots(figsize = (15,10))
ax.plot(r/1.496e+13, scaleheight, label = 'scale height')
ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_ylabel('height cm')
ax.set_xlabel('Radius AU')
ax.legend()

fig, ax = plt.subplots(figsize = (15,10))
ax.plot(r/1.496e+13, Mach, label = 'Mach #')
#ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_ylabel('Mach #')
ax.set_xlabel('Radius AU')
ax.legend()


z = np.logspace(0,100,10000)*1.496e+13
#z = np.linspace(0.5*1.496e+13, 100*1.496e+13, 100)
densityy = []
for i in range(len(z)):
    rho = []
    for n in range(0, len(r)):
        p0 = (1/np.sqrt(np.pi*2))*(SD[n]/scaleheight[n])
        density = p0*np.exp(-1*(z[i]**2)/(2*(scaleheight[n]**2)))
        rho.append(density)
    densityy.append(rho)

    
#densityy = np.reshape(densityy, [100,2000])
fig, ax = plt.subplots(figsize = (15,10))
#ax.contourf(r/1.496e+13, z, densityy, cmap = 'inferno')
a = ax.contourf(r/1.496e+13, z/1.496e+13, densityy, 20, cmap = 'inferno')
ax.set_ylim(1,5)
ax.set_xlim(0,100)
# a = ax.contourf(r/1.496e+13, -z/1.496e+13, densityy, 20, cmap = 'inferno')
# a = ax.contourf(-r/1.496e+13, z/1.496e+13, densityy, 20, cmap = 'inferno')
# a = ax.contourf(-r/1.496e+13, -z/1.496e+13, densityy, 20, cmap = 'inferno')
#ax.set_yscale('log')
# ax.set_xscale('log')
cbar = fig.colorbar(a)
ax.set_ylabel('Z AU')
ax.set_xlabel('Radius AU')

#ax.set_xscale('log')
    
    
