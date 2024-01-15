#Python Script to plot radial profiles based on text files from CASA
#Purpose is to take slice along major, minor axis of a model fits file and compare with observation fits file for image

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

font = {'family' : 'serif',
        'weight' : 'medium',
        'size'   : 13}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2


#Initializing variables
name_start = 'shifted_mod'
name_end = '_radprof.txt'
num = ['480_', '480_tw400_']
#,'274_tw1000_','274_tw900_','274_tw800_'
obs_offset = 0.3
central_offset = 105./0.04/150.1

#Axis
major = 'major'
minor = 'minor'

#Files for observations
obs_data_major = 'shiftednot_obs_major_radprof.txt'
obs_data_minor = 'shiftednot_obs_minor_radprof.txt'


#Importing data with numpy
obs_major_vals = np.genfromtxt(obs_data_major, skip_header=6)	#For obs
obs_minor_vals = np.genfromtxt(obs_data_minor, skip_header=6)   #For obs



#If working with minor axis
plt.scatter((obs_minor_vals[:,0]-central_offset-obs_offset+0.4)*0.04*150.1, obs_minor_vals[:,1], label = 'Observation', color = 'black')

#Collecting model radial profiles
for i in num:
    mod_name = name_start + i + minor + name_end
    model_minor_vals = np.genfromtxt(mod_name, skip_header=6)

    if i.find('tw') == -1:
        plt.plot((model_minor_vals[:,0]-central_offset)*0.04*150.1, model_minor_vals[:,1], label = 'Full Disk', linestyle = '--', color = 'blue')
    else:
        plt.plot((model_minor_vals[:,0]-central_offset)*0.04*150.1, model_minor_vals[:,1], label = 'Gapped Disk', linestyle = '--', color = 'red')

	#plt.plot((model_minor_vals[:,0]-central_offset)*0.04*150.1, model_minor_vals[:,1], label = 'Model', linestyle = '--')


#Formatting
#plt.title('Minor Axis Radial Profile')
plt.plot([0,0],[0,0.05], 'k:')
plt.ylim(-0.001,0.05)
plt.ylabel(r'Brightness (Jy beam$^{-1}$)', fontsize = 16)
plt.xlabel('Radius (AU)', fontsize = 16)
plt.legend(scatterpoints = 1)
plt.savefig('radprof_minor_gap.png', bbox_inches = 'tight')
plt.close()


#If working with major axis
plt.scatter((obs_major_vals[:,0]-central_offset-obs_offset-0.4)*0.04*150.1, obs_major_vals[:,1], label = 'Observation', color = 'black')

#Collecting model radial profiles
for i in num:
        mod_name = name_start + i + major + name_end
        model_major_vals = np.genfromtxt(mod_name, skip_header=6)
        
        if i.find('tw') == -1:
            plt.plot((model_major_vals[:,0]-central_offset)*0.04*150.1, model_major_vals[:,1], label = 'Full Disk', linestyle = '--', color = 'blue')
        else:
            plt.plot((model_major_vals[:,0]-central_offset)*0.04*150.1, model_major_vals[:,1], label = 'Gapped Disk', linestyle = '--', color = 'red')

#Formatting
#plt.title('Major Axis Radial Profile')
plt.plot([0,0],[0,0.05], 'k:')
plt.ylim(-0.001,0.05)
plt.ylabel(r'Brightness (Jy beam$^{-1}$)', fontsize = 16)
plt.xlabel('Radius (AU)', fontsize = 16)
plt.legend(scatterpoints = 1)
plt.savefig('radprof_major_gap.png', bbox_inches = 'tight')
plt.clf() 
