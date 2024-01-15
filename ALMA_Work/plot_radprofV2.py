#Python Script to plot radial profiles based on text files from CASA
#Purpose is to take slice along major, minor axis of a model fits file and compare with observation fits file for image

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

font = {'family' : 'serif',
        'weight' : 'medium',
        'size'   : 25.5}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2

fig = plt.figure(figsize=(18,8))
nrow,ncol=1,2

lwid = 2.1


#Initializing variables
name_start = 'shifted_mod'
name_end = '_radprof.txt'
num = ['480_tw400_', '480_']
#,'274_tw1000_','274_tw900_','274_tw800_'
central_offset = 325./0.04/150.1

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
ax = fig.add_subplot(nrow, ncol, 1)

#Collecting model radial profiles
for i in num:
    mod_name = name_start + i + minor + name_end
    model_minor_vals = np.genfromtxt(mod_name, skip_header=6)

    if i.find('tw') != -1:
        ax.plot((model_minor_vals[:,0]-central_offset+5./10.)*0.04*150.1, model_minor_vals[:,1]*1000, label = 'Gapped Disk', linestyle = '--', color = 'red', linewidth = lwid)
    else:
        ax.plot((model_minor_vals[:,0]-central_offset+5./10.)*0.04*150.1, model_minor_vals[:,1]*1000, label = 'Full Disk', linestyle = '--', color = 'blue', linewidth = lwid)

	#plt.plot((model_minor_vals[:,0]-central_offset)*0.04*150.1, model_minor_vals[:,1], label = 'Model', linestyle = '--')


#ax.scatter((obs_minor_vals[:,0]-central_offset-obs_offset+0.4)*0.04*150.1, obs_minor_vals[:,1]*1000, label = 'Observation', color = 'black')
ax.errorbar((obs_minor_vals[:,0]-central_offset+10./10.)*0.04*150.1, obs_minor_vals[:,1]*1000, yerr=1.7e-4 * 1e3, label = 'Observation', fmt = 'o', color = 'black', capsize = 8, elinewidth = 3.5)


#Formatting
#plt.title('Minor Axis Radial Profile')
ax.plot([0,0],[0,0.05*1000], 'k:')
#ax.text(77.6, 47.6, 'Minor Axis')
ax.text(35+50, 47.1-10, 'Minor Axis')
ax.set_xlim(-250, 250)
ax.set_ylim(-0.005,0.04*1000)
ax.set_ylabel(r'Brightness (mJy beam$^{-1}$)', fontsize = 25.5)
ax.set_xlabel('Radius (au)', fontsize = 25.5)
plt.legend(loc=2, numpoints = 1, fontsize = 19.5, handlelength = 1)
#plt.savefig('radprof_minor_gap.png', bbox_inches = 'tight')
#plt.close()


#If working with major axis
ax = fig.add_subplot(nrow, ncol, 2)

#Collecting model radial profiles
for i in num:
        mod_name = name_start + i + major + name_end
        model_major_vals = np.genfromtxt(mod_name, skip_header=6)
        
        if i.find('tw') != -1:
            ax.plot((model_major_vals[:,0]-central_offset+2./10.)*0.04*150.1, model_major_vals[:,1]*1000, label = 'Gapped Disk', linestyle = '--', color = 'red', linewidth = lwid)
        else:
            ax.plot((model_major_vals[:,0]-central_offset+2./10.)*0.04*150.1, model_major_vals[:,1]*1000, label = 'Full Disk', linestyle = '--', color = 'blue', linewidth = lwid)


#ax.scatter((obs_major_vals[:,0]-central_offset-obs_offset-0.4)*0.04*150.1, obs_major_vals[:,1]*1000, label = 'Observation', color = 'black')
ax.errorbar((obs_major_vals[:,0]-central_offset+5./10.)*0.04*150.1, obs_major_vals[:,1]*1000, yerr=1.7e-4 * 1e3, label = 'Observation', fmt = 'o', color = 'black', capsize = 8, elinewidth = 3.5)

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticklabels(())
ax.set_yticklabels(())

#Formatting
#plt.title('Major Axis Radial Profile')
ax.plot([0,0],[0,0.05*1000], 'k:')
ax.text(35.5+50, 47.1-10, 'Major Axis')
ax.set_xlim(-250, 250)
ax.set_ylim(-0.005,0.04*1000)
#ax.set_ylabel(r'Brightness (Jy beam$^{-1}$)', fontsize = 16)
#ax.set_xlabel('Radius (AU)', fontsize = 16)

#plt.legend(loc=2, numpoints = 1, fontsize = 19)
fig.subplots_adjust(wspace=0.0, hspace=0)
#plt.show()
fig.savefig('radprof_comb.png')
plt.clf() 
