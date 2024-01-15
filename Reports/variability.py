#imports
import matplotlib.pyplot as plt
from numpy import log10
import numpy as np
import sys

markSize = 25	#controlling size of points on diagram
c = 2.998e10	#speed of light - useful for calculations
optical_arr 	= ['U', 'B', 'V', 'R', 'I']
optical_lambdas	= [0.36,0.44,0.55,0.70,0.90]
optical_zp	= [1.79e-20,4.063e-20,3.636e-20,3.064e-20,2.416e-20]

mags_Mendigutia = [9.41, 9.16, 8.68, 8.31, 8.01]
red_Mendigutia = [1.207253521,1.026056338,0.775,0.580704225,0.371126761]
mags_me  = [9.41, 9.37, 8.82, 8.31, 8.01]
red_me = [1.545284507,1.313352113,0.992,0.743301408,0.475042254]
mags_Sylvester = [9.35, 9.2, 8.65, 8.36, 7.98]
red_Sylvester = [1.545284507,1.313352113,0.992,0.743301408,0.475042254]

#Johnson Optical; nu * F_zero * magnitude conversion

#Mendigutia
for ii in xrange(0,len(optical_lambdas)):
	if ii != 0:
		plt.scatter(log10(optical_lambdas[ii]), log10(c/optical_lambdas[ii]/1e-4 * optical_zp[ii] * 10**(-1*(mags_Mendigutia[ii] - red_Mendigutia[ii])/2.5) ) + 0.0, color = 'blue', marker = 'o', s = markSize )
	else:
		plt.scatter(log10(optical_lambdas[ii]), log10(c/optical_lambdas[ii]/1e-4 * optical_zp[ii] * 10**(-1*(mags_Mendigutia[ii] - red_Mendigutia[ii])/2.5) ) + 0.0, color = 'blue', marker = 'o', s = markSize, label = 'Mendigutia - 1998/9')

#Me
for ii in xrange(0,len(optical_lambdas)):
	if ii != 0:
		plt.scatter(log10(optical_lambdas[ii]), log10(c/optical_lambdas[ii]/1e-4 * optical_zp[ii] * 10**(-1*(mags_me[ii] - red_me[ii])/2.5) ) + 0.0, color = 'red', marker = 'o', s = markSize )
	else:
		plt.scatter(log10(optical_lambdas[ii]), log10(c/optical_lambdas[ii]/1e-4 * optical_zp[ii] * 10**(-1*(mags_me[ii] - red_me[ii])/2.5) ) + 0.0, color = 'red', marker = 'o', s = markSize, label = 'Me - 1989/93 and 1998/9' )

#Sylvester
for ii in xrange(0,len(optical_lambdas)):
	if ii != 0:
		plt.scatter(log10(optical_lambdas[ii]), log10(c/optical_lambdas[ii]/1e-4 * optical_zp[ii] * 10**(-1*(mags_Sylvester[ii] - red_Sylvester[ii])/2.5) ) + 0.0, color = 'green', marker = 'o', s = markSize )
	else:
		plt.scatter(log10(optical_lambdas[ii]), log10(c/optical_lambdas[ii]/1e-4 * optical_zp[ii] * 10**(-1*(mags_Sylvester[ii] - red_Sylvester[ii])/2.5) ) + 0.0, color = 'green', marker = 'o', s = markSize, label = 'Sylvester - 1993/4' )



plt.title(r'$\rm \, HD\, 142666\, Variability$', fontsize = 20)	#title of plot
plt.minorticks_on()	 	    																							#display plot with minor ticks

plt.tick_params(direction='in', width=1.1, which='minor', labelsize=7)
plt.tick_params(direction='in', width=1.1, which='major', labelsize=20)

plt.ylabel(r'${\rm log} \, {\lambda F_\lambda ({\rm erg \, cm^{-2} \, s^{-1}})}$', fontsize = 20, fontweight = 'bold')	#y-axis title, flux of light
plt.xlabel(r'${\rm log} \, \lambda {(\rm \mu m)}$', fontsize = 20, fontweight = 'bold') 
plt.legend()

# plt.xlim(-0.5, 0)
# plt.ylim(-8.3, -7.6)
plt.tight_layout(h_pad=-30, w_pad=-40)					#if working with ONE SED
plt.show()											#to display when trying to zoom 