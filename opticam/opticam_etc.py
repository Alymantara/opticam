# -*- coding: utf-8 -*-
import os
import yaml
from synphot.models import BlackBody1D, ConstFlux1D, PowerLawFlux1D, Empirical1D
from synphot import units
from astropy import units as u
import numpy as np
from synphot import SourceSpectrum, SpectralElement
from astropy.io import ascii
from scipy import interpolate
from sys import exit



class Sky:
    """Object that represents the sky.

    The `Sky` class is used to compute the transmission and emission of the sky. The transmission of the sky is inferred
    from the airmass, and the emission of the sky is based on the lunar phase.

    Parameters
    ----------
    lunar_phase : float, optional
        Floating point that represents the lunar phase where 0 means new moon and 1 is a full moon. Defaults to 0.

    seeing : float, optional
        The seeing. This parameter is used to calculate the background area in the S/N ratio equation.  Defaults to 1

    airmass : float, optional
        The airmass of the target. Max airmass handled is 3 Defaults to 1.

    Attributes
    ----------

    lunar_phase : float
        The phase of the moon. 0 is a new moon and 1 is a full moon.

    airmass : float
        The airmass. This parameter is related to the altitude of the target.

    seeing : float
        The seeing parameter. For large aperature telescopes, this is typically 1 arcsecond. Defaults to 1 arcsecond

    sky_transmission : Interpolated Object
        The transmission of the sky interpolated from 3000A - 30000A.

    sky_emission : Interpolated Object
        The emission of the sky interpolated from 3000A - 30000A

    """

    def __init__(self, lunar_phase=0, seeing=1, airmass=1):

        self.lunar_phase = lunar_phase
        self.airmass = airmass
        self.seeing = seeing

        self.sky_transmission = self.transmission()
        self.sky_emission = self.emission()

    def transmission(self):
        """Determine the transmission of the sky.

        The transmission of the sky is determined by the seeing. The package includes data files which read the
        appropriate transmission file based on the airmass.

        Returns
        -------
        sky_transmission : Interpolated Object
            The transmission of the sky interpolated over a given wavelength range specified in the data files.
        """

        # Find the appropriate airmass file.
        if self.airmass <= 1.25:
            trans_file = 'trans_1.txt'
        elif 1.75 > self.airmass > 1.25:
            trans_file = 'trans_1_5.txt'
        elif 1.75 <= self.airmass < 2.25:
            trans_file = 'trans_2.txt'
        elif self.airmass >= 2.25:
            trans_file = 'trans_2_5.txt'

        # Load the data file
        transmission = np.loadtxt(get_data('Sky/' + trans_file))

        # Interpolate the transmission
        sky_transmission = interpolate.InterpolatedUnivariateSpline(
            transmission[:, 0] * 10, transmission[:, 1])

        # Return the interpolated transmission.
        return sky_transmission

    def emission(self):
        """Determines the emission of the sky.

        The emission of the sky is primarily based on the lunar phase. This method computes the emission (photon flux)
        of the sky per wavelength based on the ``lunar_phase`` parameter.

        Returns
        -------
        sky_emission : Interpolated Object
            The emission of the sky interpolated over a given wavelength range specified in the data files.
        """

        # Find the appropriate date files.
        if self.lunar_phase < 0.25:
            emission_file = 'moon_00.txt'
        elif 0.25 <= self.lunar_phase < 0.75:
            emission_file = 'moon_50.txt'
        elif self.lunar_phase >= 0.75:
            emission_file = 'moon_100.txt'

        # Load the data files
        emission = np.loadtxt(get_data('Sky/' + emission_file))

        # Interpolate
        sky_emission = interpolate.InterpolatedUnivariateSpline(
            emission[:, 0] * 10, (emission[:, 1] * 1E-8))

        # Return the interpolated emission
        return sky_emission


class Target:
    """This object represents the target star which you wish to compute an exposure time for.

    This class is intended to be used for unresolved or point source objects (i.e., stars) and we do not recommend using
    it for extended objects. The class can compute the spectrum of your target by taking the temperature and scaling a
    black body spectrum to match the specified magnitude.

    Parameters
    ----------
    mag : float
        The magnitude of the target object.

    magsystem : str The magnitude system used in the `mag` parameter. The 3 options available are 'VEGAMAG', 'stmag',
    and 'abnu'. 

    filt_range : tuple or str The wavelength range of the filter you wish to observe in. Default is wavelength range
    corresponding to the Johnson V band. Also use any Johnson filter name.

    sed : arr, optional
        Optional ability to enter your own spectral energy distribution of the target object. Defaults to None.

    temp : float, optional
        The temperature (K) of the target object which is used to compute a black body spectrum. Defaults to 5778.

    Attributes
    ----------
    mag : float
        The magnitude of the target object.

    magsystem : str
        The magnitude system used in the `mag` parameter (i.e., VEGAMAG).

    filt_range : tuple
        The wavelength range of the filter you wish to observe in. !

    sed : arr, optional
        The spectral energy distribution of the target object.

    temp : float, optional
        The temperature (K) of the target object which is used to compute a black body spectrum.

    """

    def __init__(self, mag, magsystem='VEGAMAG', filt_range=None, sed=None, temp=None, index=None):

        # Define the magnitude system.
        if filt_range is None:
            filt_range = [5000, 6000]

        if isinstance(filt_range,str):
            #print('johnson_'+filt_range.lower())
            try:
                #print('johnson_'+filt_range.lower())
                b = SpectralElement.from_filter('johnson_'+filt_range.lower())
                #print(b.avgwave().value)
                filt_range = [b.avgwave().value-b.equivwidth().value/2.,
                              b.avgwave().value+b.equivwidth().value/2.]
            except:
                print("Wrong filter. Only Johnson filters supported")
                exit(1)

        if magsystem.lower() == 'vegamag':
            sys = units.VEGAMAG
        elif magsystem.lower() == 'stmag':
            sys = u.STmag
        elif magsystem.lower() == 'abnu':
            sys = u.ABmag

        # Get Vega's spectrum.
        vega = SourceSpectrum.from_vega()

        # Set attributes.
        self.mag = mag
        self.SED = sed
        self.temp = temp
        self.index = index
        self.inputFlux = units.convert_flux(filt_range, mag * sys, units.FLAM, vegaspec=vega)
        self.range = filt_range
        self.mean_range = np.mean(filt_range)
        self.F_lambda = self.starF_lambda()

    def starF_lambda(self):
        """Compute the wavelength flux of the target object.

        This method creates a black body spectrum of temperature ``temp`` and scaled that spectrum to match the flux of
        a ``mag`` magnitude object.

        Returns
        --------
        F_lambda : Interpolated Object
            The spectrum of the star interpolated from 1000 A to 30000 A.
        """

        # The object will default to a BB, if no options are given
        if (self.SED == None) & (self.temp == None) & (self.index == None): temp = 5778

        # Get the black body spectrum of an object at temperature "temp".
        if self.temp != None: sp = SourceSpectrum(BlackBody1D, temperature=self.temp * u.K)

        #sp = SourceSpectrum(ConstFlux1D, amplitude=1)  # PHOTLAM
        #print(self.index)
        # There is something odd about the index
        if self.index != None: sp = SourceSpectrum(PowerLawFlux1D, amplitude=1,
                                x_0=self.mean_range*u.AA, alpha=self.index-1)

        if self.SED != None: 
            wave = self.SED[:,0]  # Angstrom
            flux = self.SED[:,1] * units.FLAM
            sp = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux, keep_neg=True)
        
        # Scale that black body to match the flux of a "mag" magnitude star.
        sp_new = sp / np.mean(sp(self.range * u.AA, flux_unit=units.FLAM) / self.inputFlux)
        x = sp_new(range(1000, 30000) * u.AA, flux_unit=units.FLAM)

        # Interpolate the flux.
        F_lambda = interpolate.InterpolatedUnivariateSpline(range(1000, 30000), x)
        

        # Return the interpolated flux.
        return F_lambda


class Observation:
    """Creates object for an observation given a certain telescope, instrument, sky conditions, and target.

    This object takes in the three classes specified above to compute a variety of things such as the signal/noise, the
    count rate from the source, the count rate of the sky, etc.

    Parameters
    ----------
    target : Object
        The ``opticam.Target`` class.

    sky : Object
        The ``opticam.Sky`` class.

    instrument : Object
        The ``opticam.Instrument`` class.

    Attributes
    ----------
    detector_qe : Interpolated Object
        The quantum efficiency of the detector.

    telescope_area : float
        The light collecting area of the telescope (in cgs).

    source : Interpolated Object
        The flux of the target object interpolated.

    skySED : Interpolated Object
        The emission of the sky interpolated.

    skyTransmission : Interpolated Object
        The transmission of the sky interpolated.

    seeing : float
        The seeing.

    rdnoise : float
        The readout noise of the instrument.

    isImager : bool
        1 if the object is an imager. 0 if it is a spectrograph.

    gain : float
        The gain of the instrument.
    """

    def __init__(self, target, sky, instrument):


        # telescope_transm = telescope.transmission
        self.q_efficiencies = instrument.efficiencies
        self.InstTransmission = instrument.transmissions
        self.cameras = instrument.cameras
        self.telescope_area = (instrument.Telescope_rad ** 2) * np.pi
        self.source = target.F_lambda
        self.skySED = sky.sky_emission
        self.skyTransmission = sky.sky_transmission
        self.seeing = sky.seeing
        self.rdnoise = instrument.readout_noise
        self.isImager = instrument.isImager
        self.elements = instrument.element_num
        self.names = instrument.names
        self.scale = instrument.scale
        self.slit_height = instrument.slit_height
        self.ranges = instrument.range
        self.Npix_lam = instrument.Npix_lam
        self.name = instrument.name

        self.Npix, self.seeing_area = self.Npix()
        self.counts(self.source, instrument)
        self.skycounts(self.skySED, instrument)

    def Npix(self):
        """The number of pixels covered by the source and sky.

        The number of pixels is used to compute the area covered by the sky on the 
        detector as well as the amount of pixels that contributed to the readout noise. 
        This method takes the seeing and the plate scale of the instrument
        to compute Npix.

        Parameters
        ----------
        instrument : object
            The ``opticam.Instrument`` class.

        Returns
        -------
        Npix : float
            The number of pixels.

        """

        # Determine whether the instrument is an imager or a spectrograph.
        seeing_area = []
        Npix = []

        for i, row in enumerate(self.ranges):
            if self.isImager == 1:
                s = np.pi * ((self.seeing / 2) ** 2)  # set blur size for count equation
                Npix.append(s / (self.scale ** 2))
            else:
                s = self.seeing ** 2  # Slit should be about seeing size anyway
                spec_height = self.slit_height
                if self.slit_height == 'slit':  # If its the echell, then just use hieght of slit in pixels,
                    # otherwise set height in slit to be set by seeing
                    spec_height = self.seeing / self.scale
                Npix.append(spec_height * self.Npix_lam[i](range(int(row[0]), int(row[1]))))

            seeing_area.append(s)

        return Npix, seeing_area

    def skycounts(self, sky, instrument):
        """Computes the amount of counts received by the sky.

        Parameters
        ----------
        sky : object
            The ``opticam.Sky`` class.

        instrument : object
            The ``opticam.Instrument`` class.

        """

        sky_prime_dlam = []
        for i, row in enumerate(self.names):
            s_integrand = InterpolationMultiplier(
                [sky, self.q_efficiencies[i], self.skyTransmission, self.InstTransmission[i], self.cameras[i]],
                self.ranges[i])

            sky_prime_dlam.append([self.telescope_area * self.seeing_area[i] * s_integrand[1], s_integrand[0]])

        self.sky_prime_dlam = sky_prime_dlam

    def counts(self, source, instrument):
        """The counts received from the source.

        Parameters
        -----------
        source : Interpolated Object
            The wavelength flux received from the source.

        instrument : object
            The ``opticam.Instrument`` class.
        """

        h = 6.626 * 10 ** (-27)  # ergs*s
        c = 2.9979 * 10 ** (18)  # A/s
        s_prime_dlam = []
        for i, row in enumerate(self.names):
            s_integrand = InterpolationMultiplier(
                [source, self.q_efficiencies[i], self.skyTransmission, self.InstTransmission[i], self.cameras[i]],
                self.ranges[i])

            s_prime_dlam.append([self.telescope_area * (1 / (h * c)) * s_integrand[1] * s_integrand[0], s_integrand[0]])

        self.s_prime_dlam = s_prime_dlam

    def SNfromTime(self, exptime):
        """Computes the signal to noise ratio for a given exposure time.

        Parameters
        ----------
        exptime : float
            The exposure time for which you wish to compute the signal to noise ratio.

        Returns
        --------
        returnList : bytearray
            Array containing signal to noise and filter/dispersion names for each 
            filter/dispersion of instrument

        """
        returnList = []
        self.exptime = exptime
        for i, row in enumerate(self.names):
            if self.isImager == 0:
                SN_d_lam = (self.s_prime_dlam[i][0] * exptime) / np.sqrt(self.s_prime_dlam[i][0] * exptime +
                                                                         self.sky_prime_dlam[i][0] * exptime +
                                                                         (self.Npix[i] * self.rdnoise ** 2))
                returnList.append([np.array(self.s_prime_dlam[i][1]), SN_d_lam, row])
            else:
                s_prime = np.trapz(self.s_prime_dlam[i][0], self.s_prime_dlam[i][1])
                sky_prime = np.trapz(self.sky_prime_dlam[i][0], self.sky_prime_dlam[i][1])
                SN = (s_prime * exptime) / np.sqrt(s_prime * exptime + sky_prime * exptime
                                                   + self.Npix[i] * self.rdnoise ** 2)
                returnList.append([SN, row])

        self.SN = returnList
        return returnList

    def TimefromSN(self, SN):
        """Computes the exposure time need to achieve a desired signal to noise ratio.

        Parameters
        ----------
        SN : float
            The desired signal to noise ratio.

        Returns
        --------
        returnList : bytearray
            Array containing Exposure time and filter/dispersion names for each filter/dispersion of instrument

        """

        returnList = []
        self.SigToNoise = SN
        for i, row in enumerate(self.names):
            if self.isImager == 0:
                t_d_lam = (1. / (2. * self.s_prime_dlam[i][0] ** 2)) * (
                            SN ** 2 * (self.s_prime_dlam[i][0] + self.sky_prime_dlam[i][0]) + np.sqrt(
                    SN ** 4 * (self.s_prime_dlam[i][0] + self.sky_prime_dlam[i][0]) ** 2 + 4. * self.Npix[i] * (
                            self.s_prime_dlam[i][0] * SN * self.rdnoise) ** 2))
                returnList.append([np.array(self.s_prime_dlam[i][1]), t_d_lam, row])
            else:
                s_prime = np.trapz(self.s_prime_dlam[i][0], self.s_prime_dlam[i][1])
                sky_prime = np.trapz(self.sky_prime_dlam[i][0], self.sky_prime_dlam[i][1])
                t = (1. / (2. * s_prime ** 2)) * (SN ** 2 * (s_prime + sky_prime) + np.sqrt(
                    SN ** 4 * (s_prime + sky_prime) ** 2 + 4. * self.Npix[i] * (s_prime * SN * self.rdnoise) ** 2))
                returnList.append([t, row])

        self.Time = returnList
        return returnList


class Instrument:
    """Object that represents the instrument used to observe.

    It is important to note that since this exposure time calculator is 
    designed for the Observatorio Astronomico Nacional, San Pedro Martir, Mexico 
    (SPM) 2.1m telescope, the list of instruments available is exclusive to this telescope. 
    That list is::
        * OPTICAM        (Imager)


    Parameters
    -----------
    Instr_name : (str)
        The name of the instrument used.

    Attributes
    ----------
    efficiency: Interpolated Object
        UnivariateInterpolatedSpline of the instrument efficiency.

    readout_noise : float
        Value of instrument readout noise.

    filter_num : int
        Number of filters for the instrument.

    gain : float
        Gain of the instrument.

    scale : float
        The plate scale of the instrument.
    """

    def __init__(self, Instr_name, Telescope_name='spm2_1m'):       # will add multi telescope support in future

        path = get_data(Telescope_name + '/' + Instr_name + "/" + Instr_name + '_param.yaml')
        with open(r''+ path) as file:
            param = yaml.full_load(file)
        self.isImager = param['isImager']
        self.readout_noise = param['readoutnoise[electrons]']
        self.scale = param['plate_scale[arcsec/pix]']
        self.slit_height = param['Slit_height']
        self.element_num = param['filter/dispersion_Num']
        self.name = Instr_name
        self.Telescope_name = Telescope_name

        #if self.Telescope_name == 'apo3_5m': self.Telescope_rad = 350./2.0
        if self.Telescope_name == 'spm2_1m': self.Telescope_rad = 210./2.0

        names = []
        efficiencies = []
        cameras = []
        transmissions = []
        lambda_range = []
        Npix_lam = []

        for row in param['filters/dispersions']:
            names.append(row[0].split('.data')[0])
            transmission = ascii.read(get_data(Telescope_name + '/' + Instr_name + "/" + row[0]))
            q_efficiency = ascii.read(get_data(Telescope_name + '/' + Instr_name + "/" + row[1]))
            cam_efficiency = ascii.read(get_data(Telescope_name + '/' + Instr_name + "/" + row[2]))

            efficiencies.append(
                interpolate.InterpolatedUnivariateSpline(q_efficiency['col1'] * 10, q_efficiency["col2"] / 100))
            cameras.append(
                interpolate.InterpolatedUnivariateSpline(cam_efficiency['col1'] * 10, cam_efficiency["col2"]))
            ff = np.diff(transmission['col1'].data) > 0
            transmissions.append(
                interpolate.InterpolatedUnivariateSpline(transmission['col1'].data, transmission["col2"].data / 100))
            lambda_range.append([transmission['col1'].min(), transmission['col1'].max()])

            if param['isImager'] == 0:
                dispersion_file = row[0].split('_effic.data')[0] + '_disp.data'
                dispersion = ascii.read(get_data(Telescope_name + '/' + Instr_name + "/" + dispersion_file))
                Npix_lam.append(interpolate.InterpolatedUnivariateSpline(dispersion['col2'],
                                                                         (dispersion['col1'] ** (-1))))

            self.transmissions = transmissions
            self.cameras = cameras
            self.names = names
            self.efficiencies = efficiencies
            self.range = lambda_range
            self.Npix_lam = Npix_lam


def InterpolationMultiplier(functions, interpolation_range):
    """The integrand of the count equation.

    This objects takes in all of the interpolated objects that goes into the count equation (sky transmission,
    sky emission, telescope throughput, instrument effiency, etc.) and multiplies them together. It then Outputs an
    array with the values of the product vs wavelength

    Parameters
    ----------
    functions : arr-like
        List of interpolated objects to go into the count equation.

    interpolation_range : tuple
        The range that wish you to interpolate over.

    Returns
    -------
    interpolation_range, x : tuple
        Tuple where the first element is the interpolation range and the
        second element is the product array from the multiplication.


    """
    r = range(int(interpolation_range[0]), int(interpolation_range[1]))
    for i, f in enumerate(functions):
        if i == 0:
            x = np.ones(len(r))
        x = f(r) * x

    return [r, x]


_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)