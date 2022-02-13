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
import astroalign as aa

class Reduction:
    '''
    Object that performs the aperture photometry
    with SExtractor and creates a list of stars.

    The `Reduction` class is used to compute the 
    differential photometry over a particular dataset

    Parameters
    ----------
    workdir : str, optional
        Working directory where the data and catalogues are stored

    savedir : str, optional
        Directory of all the stars catalogue as measured by SExtractor

    name : str, optional
        Name of the target

    rule : str, optional
        File rule to be used when collecting all the fits files. 
        Default '*.fits'

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

    def __init__(self,workdir=None,rawdata = None,savedir=None,
                name=None,rule='*.fits'):

        if workdir is None: 
            self.workdir = './'
        else:
            self.workdir = workdir
        if savedir is None: 
            self.savedir = 'catalogues/'
        else:
            self.savedir = savedir
        if rawdata is None: 
            self.rawdata = 'raw_data/'
        else:
            self.rawdata = rawdata
        if name is None:
            self.name = 'astro'
        else:
            self.name = name
        
        self.rule = rule

        self.flns = self.get_files(self.rule)
        #print("Loading")
        self._ROOT = os.path.abspath(os.path.dirname(__file__))

    def get_files(self,rule):

        print('Looking in: ',self.workdir+self.rawdata+rule)
        self.flns = np.sort(glob.glob(self.workdir+self.rawdata+rule))

        if len(self.flns) == 0: 
            print('WARNING! >> No fits files detected')
        else:
            print('Found {} fits files.'.format(len(self.flns)))

        return self.flns


    def sextractor(self):
        """
        Routine that uses SExtractor to perform
        aperture photometry and create a catalogue of 
        stars for each file.
        """
        sext_def= self._ROOT+'/sextractor_defaults/*'
        os.system('cp '+sext_def+' '+self.workdir)
        if not os.path.isdir(self.workdir+self.savedir):
            os.system('mkdir '+self.workdir+self.savedir)

        for i,fln in enumerate(self.flns[:]):
            exists = os.path.isfile(self.workdir+self.savedir+ \
                    (fln.split(".fits")[0]+"_cat.fits").split("/")[-1])
            if not exists:
                os.system("rm temp_sextractor_file.fits")

                hdul = fits.open(fln)
                hdu1 = fits.PrimaryHDU(data=hdul[0].data)
                hdu1.header = hdul[0].header
                new_hdul = fits.HDUList([hdu1])
                new_hdul.writeto('temp_sextractor_file.fits', overwrite=True)

                hdul.close()
                os.system("chmod 777 temp_sextractor_file.fits")
                try:
                    gain = fits.getval(fln,"GAIN",0)
                except:
                    gain = 1.0

                sex_out = "sex temp_sextractor_file.fits  -c default.sex -CATALOG_NAME "+ \
                          fln.split(".fits")[0]+"_cat.fits"+" -GAIN "+str(gain)
                #print(sex_out)
                os.system(sex_out)
                print("mv "+fln.split(".fits")[0]+"_cat.fits "+\
                            self.workdir+self.savedir+".")
                os.system("mv "+fln.split(".fits")[0]+"_cat.fits "+\
                            self.workdir+self.savedir+".")
                print("{:4.0f} / {:4.0f} -- {}".format(i+1,len(self.flns),fln))
            else:
                print("{:4.0f} / {:4.0f} -- It exists!".format(i+1,len(self.flns)))

    def creat_ref_list(self,number=0):
        '''
        Create reference star list

        number :: Default. First image of the list
        '''
        fln = self.flns[number].split('/')[-1]

        print(fln)
        fl1 = self.workdir+self.savedir+fln.split(".fits")[0]+"_cat.fits"
        fl2 = self.workdir+self.rawdata+fln
        

        data = fits.getdata(fl1)

        print(data.columns)

        fig = plt.figure(figsize=(14,14))
        gc = aplpy.FITSFigure(fl2,hdu=0,figure=fig)
        gc.show_grayscale(pmin=40,pmax=99,stretch="log",invert=True)

        gc.show_circles(data['X_IMAGE'], data['Y_IMAGE'], radius=13,color='g',lw=3)

        for i in range(data['X_IMAGE'].size):
            plt.text(data['X_IMAGE'][i]+5, data['Y_IMAGE'][i]+5,data['NUMBER'][i])
        #gc.show_circles(coo_image.ra.deg[ss][pp], coo_image.dec.deg[ss][pp], 
        plt.show()

        gc.savefig(self.workdir+self.name+'_fov.pdf')

        df = pd.DataFrame(data=np.array([data['NUMBER'],
                                data['X_IMAGE'],
                                data['Y_IMAGE']]).T,columns=["id", "x","y"])
        df.to_csv(self.workdir+self.name+'_ref_stars.csv',  index_label=False,index=False)

        self.ref_stars = df

    def get_position(self,num):

        ss = self.ref_stars.id == num

        self.tar_x = self.ref_stars.x.values[ss][0]
        self.tar_y = self.ref_stars.y.values[ss][0]
        print(self.tar_x,self.tar_y)
        print(self.name+ \
            ' position is: X= {:4.0f}, Y={:4.0f}'.format(self.tar_x,self.tar_y))

    def photometry(self):
        """
        Creates a single output file from all the catalogues. 
        Cross-matches the positions of each catalogue and assigns
        every star its unique identifier.
        """
        self.photo_file = self.name+'_photo'
        apass = pd.read_csv(self.workdir+self.name+'_ref_stars.csv',
            comment="#")
        apass.set_index('id')

        vrb = True #verbose, Default=True
        plots = False
        save_output = True

        save_standards = True
        save_target = True

        flux_id = 'MAG_APER'
        err_id  = 'MAGERR_APER'

        
        PIX_EDGE = 30
        MAG_LIM = 14.0
        MAG_OFFSET = 2.0
        MAG_LIM_DIFF = 16.0
        DIFF_REF_STAR = 1362
        SEEING_LIMIT = 60
        aper_size = 3 # --> 8,11,13,15,18,21,24,27,30,33 in pixels DIAMETER
        colour_correction = True

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        coo_apass = SkyCoord(apass['x']/1000.*u.deg, apass['y']/1000*u.deg)

        num_flns = len(self.flns)
        fac_col = 0.0
        MAG_LIM_INIT = MAG_LIM

        lco_log = Path(self.workdir+self.photo_file)

        if lco_log.exists():
            print('Photometry file already exists')
            lco = pd.read_csv(lco_log)
            check_flag = True
        else:
            #"Object","Filter","MJD","Flux","Error"
            d = {'flname': [], 'Object': [],'Filter': [],'MJD': [],
                 'zp_mag': [],'zp_err': [],
                 'mag_aper': [],'err_aper': [],
                 'zp0_fit': [],'zp0_err': [],
                 'dist':[],'flag':[],'n_std':[],
                 'n_diff':[],
                 'airmass':[],'exptime':[],'fwhm':[],'observatory': [],'seeing': [],
                 'skymag': [],'telid': [],"psf_fwhm": [],
                 'm_fit': [],'b_fit': [],
                 'X': [],'Y': [],
                 'corrects':[]}
            lco = pd.DataFrame(data=d)

            dd = {'flname': [], 'id_apass': [],'Filter': [],'MJD': [],
                 'mag_aper': [],'err_aper': []
                 }
            lco = pd.DataFrame(data=d)
            sta = pd.DataFrame(data=dd)
            df2 = {}
            df3 = {}
            id3 = 0
            check_flag = False
        print("OPTICAM - Light curve generator")

        for i,flname in enumerate(self.flns[:]):
            cat_flname = self.workdir+self.savedir+flname.split('/')[-1][:-5]+'_cat.fits'
            print(flname,cat_flname)
            if check_flag : # & (True in (lco['flname'] == flname).values)
                if vrb: print(flname+" exists")
            else:
                print("Processing {:5.0f} / {:5.0f} : {}".format(i+1,num_flns,
                        flname.split('/')[-1]))

                filt = fits.getval(flname,"FILTER",0)
                obj = fits.getval(flname,"OBJECT",0)
                exptime = fits.getval(flname,"EXPOSURE",0)
                mjd_t = fits.getval(flname,"DATE-OBS",0)
                mjd = Time(mjd_t, format='isot', scale='ut1').mjd

                #camera = fits.getval(flname,"DATE-OBS",0)
                airmass = fits.getval(flname,"AIRMASS",0)
                #fwhm = fits.getval(flname,"L1FWHM",0)
                naxis1 = fits.getval(flname,"NAXIS1",0)
                naxis2 = fits.getval(flname,"NAXIS2",0)
                observatory = fits.getval(flname,"OBSERVAT",0)
                telid = fits.getval(flname,"ORIGIN",0)
                telid = flname.split("-")[0][-5:]
                pixscale = fits.getval(flname,"PIXSIZE",0)
                pixscale = 0.14
                PSF_FWHM = np.median(fits.getdata(cat_flname,0)['FWHM_IMAGE'])
                try:
                    seeing = fits.getval(flname,"L1FWHM",0)
                except:
                    seeing = PSF_FWHM*pixscale
                if seeing == "UNKNOWN": seeing = PSF_FWHM*pixscale
                if vrb: print("Seeing = {:7.3f} arcsec".format(seeing))
                if vrb: print("PSF FWHM = {:7.3f} arcsec".format(PSF_FWHM*pixscale))
                data = fits.getdata(cat_flname,0)


                #### Align images #####
                d_x,d_y = 0.0,0.0

                if i==0:
                    c_ref = np.array([data['X_IMAGE'],data['Y_IMAGE']]).T
                else:
                    c_tar = np.array([data['X_IMAGE'],data['Y_IMAGE']]).T
                    p, (pos_img, pos_img_rot) = aa.find_transform(c_ref, c_tar)
                    d_x,d_y = p.translation[0],p.translation[1]
                    print("Translation: (x, y) = ({:.2f}, {:.2f})".format(*p.translation))

                coo_image = SkyCoord((data['X_IMAGE']-d_x)/1000*u.deg, 
                                     (data['Y_IMAGE']-d_y)/1000*u.deg)

                if vrb: print("Filter: {}".format(filt))

                idx_apass, d2d_apass, d3d_apass = coo_image.match_to_catalog_sky(coo_apass)

                # Make mask due to separation
                ss = (d2d_apass.deg*1000 < 5)

                #if vrb: print("Separation to {}: {:5.3f} pixels".format(self.name,d2d_fair.deg[0]*1000))
                #if vrb: print("Index in the image catalogue for {}: {}".format(self.name,idx_fair))


                std_mag = apass['x'][idx_apass]
                ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                inst_mag = data[flux_id] + 2.5 * np.log10(exptime)
                inst_err = data[err_id]

                pp = np.isfinite(inst_mag[:,aper_size][ss]) & \
                     (data['X_IMAGE'][ss] > PIX_EDGE ) & \
                     (data['X_IMAGE'][ss] < naxis1 -PIX_EDGE)  & \
                     (data['Y_IMAGE'][ss] > PIX_EDGE ) & \
                     (data['Y_IMAGE'][ss] < naxis2 -PIX_EDGE )
                fwhm = np.median(data['FWHM_IMAGE'])

                if vrb: print("Numer of Absolute calibration stars {}".format(pp.sum()))

        #############################################################################################
        ##############################    WRITE IN PANDAS DATAFRAME    ############################
        #############################################################################################
                if ((pp.sum() >= 3)) & save_target:
                   if save_standards:
                        for jj in np.arange(pp.sum()):
                            #df3 = pd.DataFrame([[flname, std_mag[ss][pp].index.values[jj],
                            #        filt,mjd+exptime/86400./2.,
                            #        std_mag[ss][pp].values[jj],
                            #        zp0, zp0_err,
                            #        zp0_min,zp0_min_err,
                            #        inst_mag[ss][pp][jj],inst_err[ss][pp][jj],
                            #        corrects]],
                            #            columns=list(dd.keys()))
                            #sta = pd.concat([sta,df3], ignore_index=True)
                            #apass[apass.id == 1210]
                            df3[id3] = {'flname': flname, 'id_apass':apass.id.iloc[std_mag[ss][pp].index.values[jj]],
                                 'Filter': filt,'MJD': mjd+exptime/86400./2.,
                                 'epoch':i,
                                 'mag_aper': inst_mag[ss][pp][jj],
                                 'err_aper': inst_err[ss][pp][jj],
                                 'exptime': exptime,
                                 'airmass': airmass
                                 }
                            id3 += 1
        #############################################################################################
                if vrb: print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                #stop

        lco = pd.DataFrame.from_dict(df2,"index")
        sta = pd.DataFrame.from_dict(df3,"index")
        ### Now, let's put all in one log.
        #if save_output & save_target:
        #    lco.to_csv(self.photo_file+".csv")
        #    lco.to_pickle(self.photo_file+".pkl")
        if save_output & save_standards:
            sta.to_csv(self.photo_file+".csv")
            sta.to_pickle(self.photo_file+".pkl")

