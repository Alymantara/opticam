import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants as c
import astropy.coordinates as coord
from astropy import wcs

import glob
from pathlib import Path
import pandas as pd
import os
import aplpy

from astropy.time import Time
from statistics import mode
from matplotlib.gridspec import GridSpec
from lmfit import Parameters, fit_report, minimize


class Analysis:
    '''
    Object that analyses the extracted photometric catalogues.

    The `Analysis` class is used to compute the differential photometry 
    over a particular dataset

    Parameters
    ----------
    workdir : str, optional
        Working directory where the data and catalogues are stored

    catalogue : str, optional
        Directory of all the stars catalogue as measured by SExtractor

    name : str, optional
        Name of the target

    Attributes
    ----------

    aper_size : float
        Aperture size used in the 

    airmass : float
        The airmass. This parameter is related to the altitude of the target.

    raw_data : data frame
        Dataframe that contains all the photometry from all stars in the field.

    all_stars : int, arrray
        Unique identifier of each star in the field

    epochs : int, array
        Array of epochs/images

    num_epochs : int
        Total number of epochs/images

    comp_factor : float, array
        Transmission factors at each epoch
    '''

    def __init__(self,workdir=None,catalogue = None,name=None,rule = None):

        if workdir is None: 
            self.workdir = './'
        else:
            self.workdir = workdir
        if catalogue is None: 
            self.catalogue = 'astro'
        else:
            self.catalogue = catalogue

        if name is None: 
            self.name = 'astro'
        else:
            self.name = name
        self.marker = '_'+rule[:-6]
        self.aper_size = 5
        self.raw_data = pd.read_pickle(self.workdir+self.name+self.marker+'_photo.pkl') #.sort_values("MJD")

        self.all_stars = np.unique(self.raw_data.id_apass)


        self.phot_floor = 1.00
        self.phot_factr = 1.00
        self.epochs,self.num_epochs = np.unique(self.raw_data.epoch,return_counts=True)

        self.apertures = np.array([3,5,8,11,13,15,18,21,24,27,30,33])
        print("Num Stars: {}, Num Epochs: {}".format(self.all_stars.size,self.epochs.size))


    def differential_photo(self,target=None,ignore=None,save=True):
        """
        Performs the differential photometry for a specific target.

        Parameters
        ----------
        target : int, optional
            Target's unique number in the reference list

        ignore : int or int,array optional
            Star(s) to be ignored as comparison stars

        save : bool, optional
            Save the corrected photometric data for the target

        """
        mask = np.zeros(self.all_stars.size,dtype=bool)
        self.used = np.zeros(self.all_stars.size,dtype=bool)
        if target is not None:
            mask += self.all_stars == target
        if ignore is not None:
            for i in ignore:
                mask += self.all_stars == i
        mask = ~mask
        

        self.mask = mask
        

        #### Only in those that the target was also detected
        self.epochs_target = np.unique(self.raw_data.epoch[self.raw_data.id_apass == target])

        print('Target epochs: ',self.epochs_target.size)
        self.stds_used = 0
        for jj,i in enumerate(self.all_stars[self.mask]):
            
            ss = (self.raw_data.id_apass == i) 
            #np.equal(self.raw_data.epoch[ss], epochs_target).nonzero()
            #print(i,ss.sum(),self.epochs_target.size)
            if ss.sum() >= self.epochs_target.size:

                xy, x_ind, y_ind = np.intersect1d(self.epochs_target,
                                 self.raw_data.epoch[ss].values, return_indices=True)

                if np.array_equal(self.raw_data.epoch[ss].values[y_ind],self.epochs_target):
                    tmp = np.stack(10**(self.raw_data.mag_aper[ss]/-2.5),axis=0)[:,self.aper_size]
                    if self.stds_used == 0: comp = tmp[y_ind]
                    comp += tmp[y_ind]
                    self.stds_used += 1
                #print(i,np.median(np.stack(lco.mag_aper[ss],axis=0)[:,aper_size]+23))

        self.comp_factor = comp
        
        # Make the photometry of the target
        ss = (self.raw_data.id_apass == target)


        target_mag = np.stack(10**(self.raw_data.mag_aper[ss]/-2.5),axis=0)[:,self.aper_size]
        self.time = self.raw_data.MJD[ss].values
        self.mag = -2.5*np.log10(target_mag/comp)
        self.err = np.stack(self.raw_data.err_aper[ss],axis=0)[:,self.aper_size]

        if save:
            df = pd.DataFrame(data=np.array([self.time,
                                    self.mag,
                                    self.err]).T,columns=["mjd", "mag","err"])
            df.to_csv(self.workdir+self.name+self.marker+'_lc_'+str(target).zfill(2)+'.csv',  index_label=False,index=False)




    def rms_mag(self,target):
        """
        Creates an Magnitude versus RMS of every star in the field. 

        Useful to find other variable stars. Will mark the target used in
        differential_photo.

        Parameters
        ----------
        target : int
            Target's unique number in the reference list
            
        """
        ctr = 0
        fig = plt.figure(figsize=(8,8))

        std_mags = []
        new_errs = []
        rms_mags = []
        #new_errs = np.sum(np.sqrt(1/mct.weight[pp_mct]))/(pp_mct.sum())
        for jj,i in enumerate(self.all_stars):
            ss = (self.raw_data.id_apass == i)

            xy, x_ind, y_ind = np.intersect1d(self.epochs_target,
                             self.raw_data.epoch[ss].values, return_indices=True)
            color='k'
            if (np.array_equal(self.raw_data.epoch[ss].values[y_ind],self.epochs_target)) & (i != target):
                color = 'r'
            
            ctr +=1
            compu = np.stack(10**(self.raw_data.mag_aper[ss]/-2.5),axis=0)[:,self.aper_size]
            compu = compu[y_ind]
            mags = -2.5*np.log10(compu/self.comp_factor[x_ind])
            new_e = np.sum(np.stack(self.raw_data.err_aper[ss],axis=0)[:,self.aper_size])/(ss.sum())
            if color == 'r':
                plt.plot(np.median(mags),np.std(mags),'rs',ls='None',ms=15,mfc='None')
            else:
                plt.plot(np.median(mags),np.std(mags),'k.',ls='None')
            std_mags.append(np.median(mags))
            rms_mags.append(np.std(mags))
            new_errs.append(new_e)

            plt.text(np.median(mags),np.std(mags)*1.5,str(int(i)), 
                 horizontalalignment='center',
                 verticalalignment='center',fontsize=15,color=color)
        self.std_mags = np.array(std_mags)
        rms_mags,std_mags = np.array(rms_mags),np.array(std_mags)
        ss = (rms_mags <1.0) & (std_mags < 5)
        ss[self.all_stars == target] = False

        print(ss.sum(),ss.size)
        std_mags,new_errs,rms_mags,ss = zip(*sorted(zip(std_mags,new_errs,rms_mags,ss)))
        std_mags,new_errs,rms_mags,ss = np.array(std_mags),np.array(new_errs),np.array(rms_mags),np.array(ss)
        

        #######  Fit the CCD model  ######

        def ccd_model(pars, xo, dats=None,sigma=None):
            vals = pars.valuesdict()
            a1 = vals['a1']
            a2 = vals['a2']
            #print(sigma)
            model = np.interp(xo,std_mags,new_errs)*a1 + 10**a2
            if dats is None:
                return model
            if sigma is None:
                return model - dats
            else: return (model - dats)/sigma

        fit_params = Parameters()
        fit_params.add('a1', value=1.0,vary=True,min=0.0)
        fit_params.add('a2', value=-3.0,vary=True,min=-3.5)
                      

        out = minimize(ccd_model, fit_params, args=(std_mags[ss],),
                           kws={'dats': rms_mags[ss]},scale_covar=True,
                           method='nelder')

        print("a1: {:.3f}, a2: {:.5f}".format(out.params['a1'].value,
                                        10**out.params['a2'].value))

        #plt.plot(std_mags[ss],rms_mags[ss],'rx',ls='None')
        plt.plot(np.median(self.mag),np.std(self.mag),
                'bo',ls='None',mfc='None',ms=20,label=self.name)

        self.phot_floor = 10**(out.params['a2'].value)
        self.phot_factr = out.params['a1'].value
        plt.axhline(y=self.phot_floor,ls='--',color='g')
        plt.yscale('log')
        plt.ylim(8e-4,1.1)
        print("{}; SNR ccd: {:8.2f}, Observed: {:8.2f}".format(self.name,
                            1/(self.err.sum()/self.err.size),
                            1/ccd_model(out.params,np.mean(self.mag))))
        print("Mean Mag: {:.3f}".format(np.mean(self.mag)))
        print("RMS     : {:.3f}".format(np.std(self.mag)))
        print("Exposure time: {:.3f} s".format(np.nanmedian(self.raw_data.exptime)))
        
        plt.plot(std_mags,new_errs,'g-',label='SExtractor')
        plt.plot(std_mags,ccd_model(out.params,std_mags),'r',lw=3,alpha=0.5,
            label='Observed')

        plt.ylabel('Light curve RMS')
        plt.xlabel('Median Magnitude')
        lg= plt.legend()
        plt.tight_layout()
        plt.savefig(self.name+self.marker+'_rms_mag.pdf')

        datos = pd.read_csv(self.workdir+self.name+self.marker+'_lc_'+str(target).zfill(2)+'.csv')
        datos['err2'] = np.sqrt((self.err*self.phot_factr)**2 + (self.phot_floor)**2) 
        datos.to_csv(self.workdir+self.name+self.marker+'_lc_'+str(target).zfill(2)+'.csv',
                     index_label=False,index=False)

    def lightcurve(self,comp=None,std = True):
        """
        Plots the lightcurve of the target as well as a comparison star
        either with similar brightness (program finds it automatically)
        or by choosing one yourself (comp parameter).

        Parameters
        ----------
        comp : int, optional
            Comparison star's unique number in the reference list

        std : bool
            Plot comparison star. Default=True


            
        """
        fig = plt.figure(figsize=(14,8))
        gs = GridSpec(6, 1, figure=fig)

        ax1 = fig.add_subplot(gs[:4, 0])
        ax1.plot((self.time-self.time[0])*24*60,self.mag - np.median(self.mag)+\
                np.std(self.mag)*5.0,'.',ls='None',
                 alpha=0.5,color='b',
                 label=self.name+' (using {} comp stars)'.format(self.stds_used))
        ax2 = fig.add_subplot(gs[4:6, 0])
        if std:

            dd = np.abs(np.median(self.mag) - self.std_mags[self.mask])

            ll = dd == np.min(dd)
            print("Using STD star #{:3.0f} in plot".format(self.all_stars[self.mask][ll][0]))
            ss = (self.raw_data.id_apass == self.all_stars[self.mask][ll][0])
            #print(i,ss.sum(),self.epochs_target.size)
            #if ss.sum() >= self.epochs_target.size:
            
            xy, x_ind, y_ind = np.intersect1d(self.epochs_target,
                             self.raw_data.epoch[ss].values, return_indices=True)
            compu = np.stack(10**(self.raw_data.mag_aper[ss]/-2.5),axis=0)[:,self.aper_size]
            compu = compu[y_ind]
            mags = -2.5*np.log10(compu/self.comp_factor[x_ind])

            #print(self.raw_data.MJD[ss][y_ind])
            timer = self.raw_data.MJD[ss].values
            timer = timer[y_ind]
            
            faint = 'fainter'
            if np.median(mags) < np.median(self.mag): faint = 'brighter'
            ax2.plot((timer-self.time[0])*24*60,mags-np.median(mags),'k.',alpha=0.6,
                label='Field Star #{:3.0f}'.format(self.all_stars[self.mask][ll][0])+ \
                r' $\Delta$m='+'{:.3f} mag {}'.format(np.min(dd),faint))
            ax2.axhline(y=0,ls='--',color='r')

            ax2.invert_yaxis()
        ax1.invert_yaxis()
        lg = plt.legend(loc=1)
        ax1.set_title(self.name)
        ax2.set_xlabel('Minutes')
        ax1.set_ylabel(r'$\Delta$m')
        ax2.set_ylabel(r'$\Delta$m')
        plt.tight_layout()
        plt.savefig(self.name+self.marker+'_lc.pdf')

    def ccd_noise(self,image=0,aper=5):
        #aper+=4
        fl2 = self.raw_data.flname.values[image]
        fl1 = self.workdir+self.catalogue+fl2.split('/')[-1][:-5]+'_cat.fits'
        aperture = self.apertures[aper]
        #print(self.apertures[aper-4]/2,self.apertures[aper]/2)
        #print((self.apertures[aper-4]/self.apertures[aper]))
        PIX_EDGE = 30
        texp = fits.getval(fl2,'EXPOSURE')
        gain = fits.getval(fl2,'GAIN')
        darkcurr = fits.getval(fl2,'DARKCURR')
        rnoise = 1.1
        naxis1 = fits.getval(fl2,"NAXIS1")
        naxis2 = fits.getval(fl2,"NAXIS2")
        satlevel = fits.getval(fl2,"SATLEVEL")
        binn = (fits.getval(fl2,"CCDXBIN"))**2

        print("Binning: {}x{}".format(fits.getval(fl2,"CCDXBIN"),
                            fits.getval(fl2,"CCDYBIN")))
        circ2pix = 0.78 # approx from circle to pix
        #binn = 1.0
        data_tmp = fits.getdata(fl1)
        ss = (data_tmp['FLUX_APER'][:,-2] >0)
        ss&= ((data_tmp['X_IMAGE'] > PIX_EDGE ) & \
             (data_tmp['X_IMAGE'] < naxis1 -PIX_EDGE)  & \
             (data_tmp['Y_IMAGE'] > PIX_EDGE ) & \
             (data_tmp['Y_IMAGE'] < naxis2 -PIX_EDGE ))
        snr_test = snr(data_tmp['FLUX_APER'][:,-2]/texp,data_tmp['BACKGROUND']/texp,
                        texp,np.round(np.pi*(aperture/2)**2)*circ2pix,
                        rnoise,gain,dark=darkcurr,binning=binn)

        mean_bkg = np.nanmedian(data_tmp['BACKGROUND'])
        snr_fill = snr(np.logspace(0,6,1000),mean_bkg/texp,texp,
                        np.round(np.pi*(aperture/2.0)**2)*circ2pix,
                        rnoise,gain,dark=darkcurr,binning=binn)
        c1,c2,c3,c4 = snr_all(np.logspace(0,6,1000),mean_bkg/texp,texp,
                        np.round(np.pi*(aperture/2.0)**2)*circ2pix,
                        rnoise,gain,dark=darkcurr,binning=binn) #

        ll = data_tmp['FLUX_APER'][:,-2] >0.0
        mm = -2.5*np.log10(data_tmp['FLUX_APER'][:,-2][ll]/texp)
        plt.figure(figsize=(12,10))
        ax1 = plt.subplot2grid((6, 1), (0, 0), rowspan=4)
        ax1.set_title("OPTICam CCD performance vs observations")
        #plt.plot(-2.5*np.log10(data_tmp['FLUX_APER'][:,-2]/texp),1.0875/snr_test,'.',color='k',label='Theory')
        plt.plot(mm,
                data_tmp['MAGERR_APER'][:,-2][ll],'o',lw=2,ms=10,
                 color='r',label='SExtractor')
        plt.plot(-2.5*np.log10(np.logspace(0,6,1000)),1.0875/snr_fill,
                        'k-',label='Theory')
        plt.plot(-2.5*np.log10(np.logspace(0,6,1000)),1.0875/c1,
                        'm--',label='Shot Noise')
        plt.plot(-2.5*np.log10(np.logspace(0,6,1000)),1.0875/c2,
                        'g--',label='Sky Noise')
        plt.plot(-2.5*np.log10(np.logspace(0,6,1000)),1.0875/c3,
                        'c--',label='Readout Noise')
        plt.plot(-2.5*np.log10(np.logspace(0,6,1000)),1.0875/c4,
                        'r--',label='Dark Current')
        #plt.plot(-2.5*np.log10(np.logspace(0,6,1000)),1.0875/c4,'k--')
        lg = plt.legend()
        plt.yscale('log')
        plt.xlim(np.nanmin(mm)-0.5,np.nanmax(mm)+0.5)
        plt.ylim(1e-3,0.8)
        plt.ylabel('$\sigma_{mag}$')
        ax1.xaxis.set_ticklabels([])
        #plt.axvline(x= -2.5*np.log10(satlevel/texp),ls='--',color='k',alpha=0.4,lw=3)
        #print("Saturation count rate = {:.2f} ct/s".format(-2.5*np.log10(satlevel/texp)))

        ax2 = plt.subplot2grid((6, 1), (4, 0), rowspan=2)
        plt.plot(-2.5*np.log10(data_tmp['FLUX_APER'][:,-2][ll]/texp),
            data_tmp['MAGERR_APER'][:,-2][ll]/(1.0875/snr_test[ll]),'k.')
        plt.axhline(y=1.0,ls='--',color='k')
        plt.axhline(y=1.05,ls='--',color='r')
        plt.axhline(y=0.95,ls='--',color='r')
        plt.text(-13.,1.055,'5%',color='r')
        #plt.text(-13.,1.045,'5%',color='r')
        plt.xlim(np.nanmin(mm)-0.5,np.nanmax(mm)+0.5)
        plt.ylim(0.9,1.1)
        plt.xlabel('-2.5 log(Count rate)')
        plt.ylabel(r'$\sigma_{SExt}$ / $\sigma_{theory}$')
        plt.tight_layout(h_pad=0)
        plt.savefig(self.name+self.marker+'_SNR_andor.pdf')

def snr(rate,bkg,time,npix,rn,gain,dark=0.0,binning=1):
    source = rate * time
    shot_noise = rate * time
    sky_noise = bkg * npix * time*binning
    ro_noise = (rn**2 + (gain/2.0)**2 * npix*binning)
    dark_noise = dark * npix * time*binning
    
    return source / np.sqrt(shot_noise + sky_noise + ro_noise + dark_noise)

def snr_all(rate,bkg,time,npix,rn,gain,dark=0.0,binning=1):
    source = rate * time
    shot_noise = rate * time
    sky_noise = bkg * npix * time * binning
    ro_noise = (rn**2 + (gain/2.0)**2 * npix)*binning
    dark_noise = dark * npix * time * binning
    
    return source / np.sqrt(shot_noise ), source / np.sqrt(sky_noise),source / np.sqrt(ro_noise), source / np.sqrt(dark_noise)




