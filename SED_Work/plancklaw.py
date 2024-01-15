#code for blackbody curve maker

#imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import scipy as sp
import math

#initialize variables
k = 1.38e-23 #boltzmann constant J/K
h = 6.63e-34 #planck's constant J*s
c = 2.99e8   #speed of light, m/s

T  = 8e3 #effective or peak temperature of the blackbody, K
T2 = 9e3 #effective or peak temperature of the blackbody, K
wl = []  #wavelength, m
E  = []  #energy
E2 = []
 
#filling wavelength array
for i in sp.arange(1e-12, 3e-6, 1e-11):
	wl.append(i)

for j in xrange(len(wl)):
	E.append(2*h*c**2 / wl[j]**5 * 1 / (np.exp(h*c / (wl[j] * k * T)) - 1))
	E2.append(2*h*c**2 / wl[j]**5 * 1 / (np.exp(h*c / (wl[j] * k * T2)) - 1))


#plotting intensity as a function of wavelength
plt.plot(np.array(wl)/1e-9, np.array(E)*1e-9)
#plt.plot(np.array(wl)/1e-9, np.array(E2)*1e-9)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (W / m^2 / nm)')
plt.title('Blackbody Spectrum')
#plt.xticks(sp.arange(0, max(wl), 250))
#tick.set_minor_locator()
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.grid()
plt.show()

print np.trapz(E)