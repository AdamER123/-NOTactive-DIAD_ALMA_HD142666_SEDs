from astropy.io import fits
from matplotlib.patches import Ellipse

import numpy as np
import nc
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import os
import copy
from copy import deepcopy
import time
from scipy.interpolate import griddata


def read_fits(fitsfile,model=False):
    dum = imgObj()
    dum.readFits(fitsfile=fitsfile,model=model)

    return dum



class imgObj():
    def __init__(self):
        self.image       = 0
        self.x           = 0
        self.y           = 0
        self.xx          = 0
        self.yy          = 0
        self.nx          = 0
        self.ny          = 0
        self.sizepix_x   = 0
        self.sizepix_y   = 0
        self.nfreq       = 0
        self.freq        = 0
        self.v           = 0           #km/s
        self.beam        = [0.,0.,0.]  #arcsec
        self.frest       = 0           #Hz
        self.moment0     = 0.
        self.moment1     = 0.
        self.rms         = 0.

# ---------------------------------------Read and write ---------------------------------------
    
    def readFits(self,fitsfile='',model=False,v_def='RADIO'):
        """read fits data and make an image object"""
        
        hdulist        = fits.open(fitsfile)
        self.image     = (hdulist[0].data)

                
        if len(self.image.shape) >3:
            self.image = self.image[0,:,:,:]
            
        self.image = fits_2_np_array(self.image)
            
        prihdr         = hdulist[0].header
        self.nx        = prihdr['NAXIS1']
        self.ny        = prihdr['NAXIS2']
        self.nfreq     = prihdr['NAXIS3']
        self.nv        = prihdr['NAXIS3']
        self.sizepix_x = prihdr['CDELT1']*3600.
        self.sizepix_y = prihdr['CDELT2']*3600. # degree->arc
        df             = prihdr['CDELT3']
        f0             = prihdr['CRVAL3'] #Hz
        
        self.frest     = prihdr['RESTFRQ']#Hz

        if model==False:
            beam_maj       = prihdr['BMAJ']*3600.
            beam_min       = prihdr['BMIN']*3600.
            beam_pa        = prihdr['BPA']
            self.beam      = [beam_maj,beam_min,beam_pa]
        
        self.x         = np.linspace(-0.5*self.nx,0.5*self.nx,self.nx)*self.sizepix_x
        self.y         = np.linspace(-0.5*self.ny,0.5*self.ny,self.ny)*self.sizepix_y
        
        self.freq      = f0+np.arange(self.nfreq)*df

        #
        xx,yy    = np.meshgrid(self.x,self.y, indexing='ij')
        self.xx        = xx
        self.yy        = yy
        
        if v_def=='RADIO': self.v = (self.frest-self.freq)/self.frest*nc.cc*1e-5
        if v_def=='OPTICAL':self.v= (self.frest/self.freq-1.)*nc.cc*1e-5
        
        hdulist.close()

        


    def write_fits(self,fitsname='',header_dic=''):

        # the order of axes on an Numpy array are opposite
        # of the order specified in the FITS file.
        
        data    = np_2_fits_array(self.image)

        hdu     = fits.PrimaryHDU(data)
        hdulist = fits.HDUList([hdu])
        prihdr  = hdulist[0].header

        basic_header = generate_basic_header(self)

        if header_dic!='':
            dum = basic_header.update(header_dic)

        prihdr.update(basic_header)

        safe_write_fits(hdu,fitsname)

        

# ---------------------------------------image statistics --------------------------------------
    def calculate_rms(self,radius=[4.,8],chan_id=[''],moment=10):
        """calculate rms in each channel"""
        
        xx,yy=self.xx,self.yy
        rr   = np.sqrt(xx*xx+yy*yy)
        

        if chan_id[0] !='':
            ind = (rr>=radius[0]) & (rr<radius[1])
            dum = self.image[chan_id]
            dum_std = np.std(dum[ind])
            rms= dum_std
        if moment ==0:
            ind = (rr>=radius[0]) & (rr<radius[1])
            dum = self.moment0
            dum_std = np.std(dum[ind])
            rms= dum_std
            
        return rms
                      
            
        

 
    
#-------------------------- ----------Plotting  --------------------------

    def plot_channel_map(self,fig,chan_range=[-1],nrow=2,ncol=6,vmax='',vmin='',cmap='ocean_r',xlim=2.,\
                         beam_kwg = {'xloc':1.5,'yloc':-1.5,'facecolor':None,'color':'black','label_beam':True},\
                         axiscolor='k',verbose=False,label_v=False,plot_contour=False,\
                         contour_level=np.array([3,6,9]),c_colors='white',**kwg):
        # swap axes for imshow
        
        img_temp = np.swapaxes(self.image,0,1)*1e3

        # if plot_contour==True:
        #     self.calculate_rms()
        #     rms = self.rms
            
        
        plt.subplots_adjust(hspace = .001,wspace=0.001)
        if chan_range[0] == -1:
            chan_range = np.arange(self.nfreq)

        
            
        n_channel = len(chan_range)
        if vmax=='':vmax = self.image[:,:,chan_range[0]:chan_range[-1]+1].max()*1e3
        if vmin=='':vmin = self.image[:,:,chan_range[0]:chan_range[-1]+1].min()*1e3

        if verbose==True:
            print 'n_channel=',n_channel
            print 'vmax,vmin=',vmax,vmin
        
        for ii in range(n_channel):
            tt0 = time.time()
            ax = fig.add_subplot(nrow,ncol,ii+1)
            cs = ax.imshow(img_temp[:,:,chan_range[ii]],extent=(self.x[0],self.x[-1],self.y[0],self.y[-1]),\
                           aspect='equal',vmax=vmax,vmin=vmin,cmap=cmap,origin='lower',**kwg)
            
            if plot_contour==True:
                levels = contour_level
                cs_contour = ax.contour(self.xx,self.yy,self.image[:,:,chan_range[ii]]*1e3,levels,\
                                        origin='lower',colors=c_colors)

                #print contour_level*rms[chan_range[ii]]*1e3
                
            if label_v==True:
                dum_v = self.v[chan_range[ii]]
                plt.text(xlim*(0.5),xlim*0.7,'{:.2f}'.format(dum_v),color=beam_kwg['color'])

            if ii == 0:    
                beam_kwg.update({'ax':ax,'xloc':xlim*0.8,'yloc':xlim*(-0.8)})
                self.plot_beam(**beam_kwg)

            plt.xlim(xlim,-xlim)
            plt.ylim(-xlim,xlim)

            if axiscolor!= 'k':
                ax.spines['bottom'].set_color(axiscolor)
                ax.spines['top'].set_color(axiscolor) 
                ax.spines['right'].set_color(axiscolor)
                ax.spines['left'].set_color(axiscolor)

            if (ii != (nrow-1)*ncol and n_channel>1):
                ax.set_xticklabels(())
                ax.set_yticklabels(())
            else:
                ax.set_xlabel(r'$\Delta\alpha$ [$\prime\prime$]')
                ax.set_ylabel(r'$\Delta\delta$ [$\prime\prime$]')
                

        if n_channel%2!=0 and n_channel>1:
            if n_channel ==1:
                index = 1
            else:
                index = 0
                
            ax = fig.add_subplot(nrow,ncol,index)
            ax.set_xticklabels(())
            ax.set_yticklabels(())
                   
        
        add_color_bar(fig,ax,cs,height_no=nrow)

# --------------------------------

    def plot_one_channel(self,fig,ax,chan_ind=0,vmax='',vmin='',cmap='ocean_r',xlim=2.,\
                         axiscolor='k',verbose=False,label_v=False,plot_contour=False,\
                         contour_level=np.array([3,6,9]),c_colors='white',\
                         units_label='mJy beam$^{-1}$ ',color_bar_width=0.02,\
                         dark_background=False,mom1 = False,label_beam=False,\
                         **kwg):
        # swap axes for imshow
        factor = 1e3
        if mom1 ==True: factor = 1
        img_temp = np.swapaxes(self.image[:,:,chan_ind],0,1)*factor           
        
        if vmax=='':vmax = img_temp.max()
        if vmin=='':vmin = img_temp.min()

        if dark_background==True:

            beam_kwg = {'facecolor':None,'color':'w','label_beam':label_beam}
        else:
            beam_kwg = {'facecolor':None,'color':'k','label_beam':label_beam}

            

        cs = ax.imshow(img_temp,extent=(self.x[0],self.x[-1],self.y[0],self.y[-1]),\
                           aspect='equal',vmax=vmax,vmin=vmin,cmap=cmap,origin='lower',**kwg)
            
        if plot_contour==True:
            levels = contour_level #*1e3*rms
            cs_contour = ax.contour(self.xx,self.yy,self.image[:,:,chan_ind]*1e3,levels,\
                                        origin='lower',colors=c_colors)

                
        if label_v==True:
                dum_v = self.v[chan_ind]
                plt.text(xlim*(0.5),xlim*0.7,'{:.2f}'.format(dum_v),color=beam_kwg['color'])

         
        beam_kwg.update({'ax':ax,'xloc':xlim*0.75,'yloc':xlim*(-0.75)})
        self.plot_beam(**beam_kwg)

        plt.xlim(xlim,-xlim)
        plt.ylim(-xlim,xlim)

        
        ax.set_xlabel(r'$\Delta\alpha$ [$\prime\prime$]')
        ax.set_ylabel(r'$\Delta\delta$ [$\prime\prime$]')

        return cs
                               
        #add_color_bar(fig,ax,cs,height_no=1,loc='top',ylabel=units_label)
    
#--------------------------------------        


    def plot_beam(self,xloc=0.1,yloc=0.1,ax=None,label_beam=True,color='k',**kwarg):

        width  = self.beam[0]
        height = self.beam[1]

        ee=Ellipse(xy=(xloc,yloc),width = width,height=height,angle=90.-self.beam[2],\
                   fill=False,color=color,**kwarg)
        ax.add_artist(ee)

        if label_beam == True:
            dum='{:.2}" X {:.2}"'.format(width, height)
            print xloc-width*1.5,yloc-height*0.4
            ax.text(xloc-width*1.5,yloc-height*0.4,dum,color=color)
        
        

    def plot_moment_map(self,fig,ax,moment=0,color_bar=True,\
                        beam_kwg = {'xloc':1.5,'yloc':-1.5,'facecolor':None,'color':'black'},\
                        xlim=4.,**kwg):
        
        if moment==0: moment_map = self.moment0
        if moment==1: moment_map = self.moment1

        #moment_map = moment_map*1e3
        moment_map = np.swapaxes(moment_map,0,1)*1e3
        
        cs = ax.imshow(moment_map,extent=(self.x[0],self.x[-1],self.y[0],self.y[-1]),\
                           aspect='equal',origin='lower',**kwg)
        
        if color_bar==True:
            add_color_bar(fig,ax,cs,height_no=1.,ylabel='mJy km $^{-1}$ beam$^{-1}$ ')
        
        #beam_kwg.update({'ax':ax})
        
        self.plot_beam(ax=ax,**beam_kwg)
        ax.set_xlim(xlim,-xlim)
        ax.set_ylim(-xlim,xlim)

        return cs

    def plot_scale_bar(self,loc=[0.,0.],scale=1.,color='k',label=None):
        x = np.linspace(0,1,10)*-1.*scale+loc[0]
        y = np.zeros((10))+loc[1]
        plt.plot(x,y,color=color)
        plt.plot(x[0]+np.zeros((2)),[loc[1]-0.1*scale,loc[1]+0.1*scale],color=color)
        plt.plot(x[-1]+np.zeros((2)),[loc[1]-0.1*scale,loc[1]+0.1*scale],color=color)
    
        if label!=None:
            plt.text(x[-1]+1.1*scale,loc[1]+0.25*scale,label,fontsize=15,color=color)

##########################################################################
# --------------------------------------------------------------------------------------------------
# supporting functions


def add_color_bar(fig,ax,cs,xoffset=0.1,yoffset=0.0,height_no=1.,aspect=100.,ylabel='mJy beam$^{-1}$ ',loc='right'):
    
    pos1 = ax.get_position()
    dx = pos1.x1-pos1.x0
    dy = pos1.y1-pos1.y0


    if loc == 'right':
        cbaxes = fig.add_axes([pos1.x1+dx*xoffset, pos1.y0+yoffset,pos1.height*height_no/aspect, pos1.height*height_no])
        cb = plt.colorbar(cs, cax = cbaxes,orientation ='vertical',extend='min')
        #cb.ax.tick_params(axis='y',direction='in',labeltop='on',labelbottom='off')
        if ylabel!='':
            cb.ax.set_ylabel(ylabel,rotation=-90,labelpad=20)
            cb.ax.yaxis.set_label_position('right')
        
    if loc=='top':
        dx = pos1.x1-pos1.x0
        dy = pos1.y1-pos1.y0
        off = 0.02
        width = 0.02
        
        cbaxes = fig.add_axes([pos1.x0+dx*0.025, pos1.y1+dy*0.01-off,dx*0.95,width])
        cb = plt.colorbar(cs, cax = cbaxes,orientation ='horizontal',extend='both')
        
        cb.ax.tick_params(axis='x',direction='in',labeltop='on',labelbottom='off',\
                          labelsize=11,color='k',length=10)
        
        if ylabel!='':
            cb.ax.set_xlabel(ylabel)
            cb.ax.xaxis.set_label_position('top')    


def safe_write_fits(hdu,fname):
    if os.path.exists(fname):
        print fname+' already exists'
        dum = raw_input('Do you want to overwrite it (yes/no)?')
        if (dum.strip()[0]=='y')|(dum.strip()[0]=='Y'):
            os.remove(fname)
            hdu.writeto(fname)
        else:
            print 'No image has been written'
    else:
        hdu.writeto(fname)

        
    
def np_2_fits_array(image,casa=True):
    data = image.transpose()
    
    return data


def fits_2_np_array(image,casa=True):
    """ reshape axis order"""

    data = image.transpose()
            
    return data



def generate_basic_header(img_obj):

    bunit = 'Jy/pixel'

    if img_obj.beam[0] > 0.:bunit ='Jy/beam'

    header = {'BUNIT':bunit,\
              'RESTFRQ':img_obj.frest,\
              'CDELT1': img_obj.sizepix_x/3600.,\
              'CDELT2': img_obj.sizepix_y/3600.,\
              'CTYPE1':'RA---SIN',\
              'CTYPE1':'DEC---SIN',\
              'CRPIX1':(img_obj.nx+1.)/2.,\
              'CRPIX2':(img_obj.ny+1.)/2.,\
              'CRPIX3':1.,\
              'CRVAL3':img_obj.freq[0],\
              'CDELT3':img_obj.freq[1]-img_obj.freq[0],\
              'CUNIT3':'Hz      ' }
    
    if img_obj.beam[0]>0:
        beam_info={'BMAJ':img_obj.beam[0]/3600.,\
                   'BMIN':img_obj.beam[1]/3600.,\
                   'BPA': img_obj.beam[2]}
        dum = header.update(beam_info)
    
    return header



