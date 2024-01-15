#program that outputs the associated SED for some star, given the associated row and column in the MAIN_DATA.csv file

#imports
import matplotlib.pyplot as plt
from math import pow
from numpy import nan
from numpy import log10
import numpy as np
import pandas as pd
import sys

#initialization
total = 0
c = 2.998e10
IRS_xobj_arr = [3,6,10,11,12,13,14,16,19,20,25,34,44,47,50,57,58,59,60,63,65,67,75,82,85,87,92,94,95,96,97,98,99,100,101,104,105,106,107,108,109,111,112,115,116,120,121,130,133,136,142,145]
optical_arr = ['U', 'B', 'V', 'R', 'I']
JHK_arr 	= ['J', 'H', 'K']
WISE_arr 	= ['W1', 'W2', 'W3', 'W4']
radio_arr 	= ['25', '60', '100', '450', '850', '890', '1300', '2600', '2700', '7000', '13000', '36000']
reddenOPT_arr	= ['A_U', 'A_B', 'A_V', 'A_R', 'A_I']
redden2MASS_arr	= ['A_J', 'A_H', 'A_K']
reddenWISE_arr	= ['A_W1', 'A_W2', 'A_W3', 'A_W4']


#NEED TO CHECK: WHETHER BV DATA IS AVAILABLE FOR REDDENING, WHETHER FLAGS ARE OK, WHETHER DATA IS MISSING, AND WHETHER FITS WITHIN HERNANDEZ BOX
#function to get and initially process data from CSV file using pandas
def dataGetter(x_obj):
	mag_list = []
	with open('MAIN_DATA.csv', 'rb') as csvfile:
		mycsv = pd.read_csv(csvfile)

		#optical magnitudes
    	for i in xrange(0, len(optical_arr)):
    		mag_list.append(float(mycsv[optical_arr[i]][x_obj]) - float(mycsv[reddenOPT_arr[i]][x_obj]))

		#2MASS magnitudes
    	for i in xrange(0, len(JHK_arr)):
    		spot = mycsv['2MASS_FLAG'][x_obj]
    		if spot[i] == ('0' or 'H'):
    			mag_list.append(float(mycsv[JHK_arr[i]][x_obj]) - float(mycsv[redden2MASS_arr[i]][x_obj]))
    		else: 
    			mag_list.append(float(nan))
    	
    	#WISE magnitudes
    	for i in xrange(0, len(WISE_arr)):
    		spot = mycsv['WISE_FLAG'][x_obj]
    		if spot[i] == ('0' or 'H'):
    			mag_list.append(float(mycsv[WISE_arr[i]][x_obj]) - float(mycsv[reddenWISE_arr[i]][x_obj]))
    		else:
    			mag_list.append(float(nan))
    	
    	#AKARI 9 micron flux
    	if mycsv['F9_FLAG'][x_obj] == '3':
    		mag_list.append(float(mycsv['F9'][x_obj]))
    	else:
    		mag_list.append(nan)

    	#AKARI 18 micron flux
    	if mycsv['F18_FLAG'][x_obj] == '3':
    		mag_list.append(float(mycsv['F18'][x_obj]))
    	else:
    		mag_list.append(nan)

    	#radio magnitudes
    	for i in xrange(0, len(radio_arr)):
			spot = mycsv[radio_arr[i]][x_obj]
			mag_list.append(float(spot))

	mag_list.append(mycsv['MAIN_ID'][x_obj])
	mag_list.append(mycsv['SP_TYPE'][x_obj])
	mag_tup = tuple(mag_list)
	return mag_tup


#Need to switch to do 136, 137...also need to switch i range to from 0 to 2
#for z in [136]:
for z in xrange(0, 144, 4):
	for i in xrange(0, 4):
		plt.subplot(2,2,i+1)
		x_obj = z + i
		#using function call to output data needed
		U,B,V,R,I,J,H,K,w1,w2,w3,w4,f9,f18,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12, identifier, sptype = dataGetter(x_obj)
		print z, x_obj, identifier
		radio_dat = [r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12]
		

		#arrays for getting data and overplotting IRS spectrum 
		IRS_lambda 	= []
		IRS_F		= []
		IRS_err		= []

		#condition to check how to read in file because sometimes file type different and number of lines are different to read from
		if (x_obj+2 in IRS_xobj_arr) == True:
			with open('HAEBE Spectra/csvnum' + str(x_obj + 2) + '_spectra.tbl', 'r') as f:
				lines = f.readlines()
				line_data = lines[110:]

				for line in line_data:
					#specifying data in column to put into array
					IRS_lambda.append(float(line[4:14]))
					IRS_F.append(1e-23 * c/1e-4/float(line[4:14]) * float(line[15:28])) 	#converting from Jy
					if x_obj == 12: 
						print line[15:28]
					IRS_err.append(1e-23 * c/1e-4/float(line[4:14]) * float(line[30:40]))

		elif x_obj+2 == 2:
			with open('HAEBE Spectra/csvnum' + str(2) + '_spectra.txt', 'r') as f:
				lines = f.readlines()
				line_data = lines[8:]

				for line in line_data:
					#specifying data in column to put into array
					IRS_lambda.append(float(line[4:14]))
					IRS_F.append(1e-23 * 1e4 * c/float(line[4:14]) * float(line[15:28])) 	#converting from Jy
					IRS_err.append(1e-23 * 1e4 * c/float(line[4:14]) * float(line[30:40]))

		plt.scatter(log10(IRS_lambda), log10(IRS_F), color = 'black', s = 0.2)
		# plt.errorbar(log10(IRS_lambda), log10(IRS_F), xerr = IRS_err, color = 'black', linewidth = 0.1, zorder = 0)


		#creating colors to act as check from Hernandez et al.
		JH = J - H
		HK = H - K
		energ_B = c/0.44e-4 * 4.26e-20 	* 10 ** (-B/2.5)

		#for setting up flux data, dereddened flux conversion from magnitudes, note the conversion factors change
		if ( (HK) < -3.0/2.0*( (JH) - 23.0/15.0) and (HK) > ((JH) - (331.0/560.0))*-28.0/25.0 and (JH) > (14.0/11.0*(HK) - 3.0/5.0) and (JH) < (65.0/53.0*(HK) + 1.0/53.0) ) or \
			((HK) < -0.890065348664074 * ((JH) - 0.57184917296716) and (HK) > -0.861286747696601 * ((JH) + 0.0147349830253168) and \
				(JH) > 1.62175790018995*(HK)-0.424174147193922 and (JH) < 1.31575306101884*(HK)-0.00950529799682048):
					#plotting in flux values
					plt.scatter(log10(0.36),	log10(c/0.36e-4 * 1.79e-20 	* 10 ** (-U/2.5)), color = 'indigo', marker = "o") #Johnson Optical; nu * F_zero * magnitude conversion
					plt.scatter(log10(0.44),  	log10(c/0.44e-4 * 4.063e-20	* 10 ** (-B/2.5)), color = 'indigo', marker = "o")
					plt.scatter(log10(0.55),  	log10(c/0.55e-4 * 3.636e-20	* 10 ** (-V/2.5)), color = 'indigo', marker = "o")
					plt.scatter(log10(0.70),  	log10(c/0.70e-4 * 3.064e-20	* 10 ** (-R/2.5)), color = 'indigo', marker = "o")
					plt.scatter(log10(0.90),  	log10(c/0.90e-4 * 2.416e-20	* 10 ** (-I/2.5)), color = 'indigo', marker = "o")	
					plt.scatter(log10(1.235),  	log10(1e7 * 1.235 *3.129*pow(10,-13) * 10 ** (-(J/2.5))), color = 'blue', marker = "o")	#2MASS; the 1e7 at the front comes from converting zeropoint to ergs, orig is in Watts
					plt.scatter(log10(1.662),  	log10(1e7 * 1.662 *1.133*pow(10,-13) * 10 ** (-(H/2.5))), color = 'blue', marker = "o") #these 2MASS bands go by lambda*F_zero*mag convert using lambda*F_lambda
					plt.scatter(log10(2.159),  	log10(1e7 * 2.159 *4.283*pow(10,-14) * 10 ** (-(K/2.5))), color = 'blue', marker = "o") 
					plt.scatter(log10(9), 		log10(c/9e-4  * 1e-23 * f9), color = 'green', marker = "o") #Akari
					plt.scatter(log10(18), 		log10(c/18e-4 * 1e-23 * f18), color = 'green', marker = "o")
					plt.scatter(log10(3.368),  	log10(1e7 * 3.368 *8.1787*pow(10,-15) * 10 ** (-(w1/2.5))), color = 'orange', marker = "o") #WISE
					plt.scatter(log10(4.618),  	log10(1e7 * 4.618 *2.4150*pow(10,-15) * 10 ** (-(w2/2.5))), color = 'orange', marker = "o") 
					plt.scatter(log10(12.082), 	log10(1e7 * 12.082*6.5151*pow(10,-17) * 10 ** (-(w3/2.5))), color = 'orange', marker = "o") 
					plt.scatter(log10(22.194), 	log10(1e7 * 22.194*5.0901*pow(10,-18) * 10 ** (-(w4/2.5))), color = 'orange', marker = "o")
					for k in xrange(0, len(radio_dat)):
						plt.scatter(log10(float(radio_arr[k])), log10(1e-23 * c/1e-4/float(radio_arr[k]) * radio_dat[k]), color = 'red', marker = "o") #Radio data need to do a conversion from Jansky to nu*F_nu
					
					total += 1
		else:
			print identifier #checking which stars are not in Hernandez box		


		#Plotting photosphere
		if sptype[0] == 'A' or sptype[0] == 'B' or sptype[0] == 'F':
			#reading in photospheric data for stars
			with open('flujosnormJ_stnd_mamajek.txt') as f:
				lines = f.readlines()
				#loop for grabbing data from mamajek using csv sptypes
				for a in xrange(0,27):
					#loop for getting data
					line_sp = lines[a * 9]

					sp = str(line_sp[0:2])
					sp = sp.rstrip('\n')

					#If no spectral type, guess
					if sptype == 'Ae':
						sptype = 'A5'
					elif sptype == ('Be' or '[B]'):
						sptype = 'B5'

					#Using own spectral type, then grabbing from mamajek as appropriate
					if sp == sptype[0:2]:
						lambda_phot = []
						energ_phot 	= []

						flux_ex = lines[a * 9 + 2]
						lam_F_lam_B = float(flux_ex[14:21])

						for b in xrange(1, 9):
							line_data = lines[a * 9 + b]
							#putting together arrays of data for plotting photosphere
							lambda_phot.append(float( line_data[5:11] ))
							energ_phot.append(float( line_data[14:22] ) - lam_F_lam_B + log10(energ_B))

						#adding Rayleigh-Jeans tail
						line_const = 3.0 * log10(float( line_data[5:11] )) + float( line_data[14:21] ) - lam_F_lam_B + log10(energ_B)
						for val in np.arange(1.1 * (lambda_phot[-1]), (50 * lambda_phot[-1]), 0.1):
							lambda_phot.append(val)
							energ_phot.append(line_const - 3.0 * log10(val) )

						#plotting in photosphere
						plt.plot(log10(lambda_phot), energ_phot, color = 'pink')
		

		#Model Herbig Ae/Be Star SED
		with open('A2 HAEBE Model/A2_t3_modA/A2_t3_e0p001/modA.sed.A2.t3.mu0p25.mp1em7.e0p001.alb0.dat') as f:
			#Declaring associated variables for model SED
			lambda_model 	= []
			nu_F_nu 		= []
			
			#Cutting out some extraneous lines from file
			lines = f.readlines()
			line_data = lines[2:] 

			#Reading in data
			for line in line_data:
				lambda_model.append(float(line[4:18]))
				nu_F_nu.append(float(line[36:50]))

			#Plotting
			plt.plot(log10(lambda_model), log10(nu_F_nu), color = 'cyan') 
		
		#formatting graphing
		plt.title(identifier.upper(), fontsize = 9, fontweight = 'bold')										#title of plot

		plt.minorticks_on()	 	    																	#display plot with minor ticks

		plt.tick_params(direction='in', width=1.1, which='minor', labelsize=7)
		plt.tick_params(direction='in', width=1.1, which='major', labelsize=7)

		plt.ylabel(r'Log $\bf{\lambda F_\lambda (\frac{erg}{s*cm^2})}$', fontsize = 9, fontweight = 'bold')	#y-axis title, flux of light
		plt.xlabel(r'Log Wavelength $\bf{(\mu m)}$', fontsize = 9, fontweight = 'bold')         			#x-axis title, lambda or wavelength in microns
	plt.tight_layout()
	plt.savefig('HAEBE SEDs/SED_set' + str(z) + '.jpeg')
	plt.clf()
#plt.show()

print total
sys.exit()