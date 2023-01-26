import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]



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
