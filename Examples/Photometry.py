import opticam

datadir = 'bl_cam/'
catdir = 'bl_cam_cat/'
name = 'BL_Cam_i'

extract = False
analyse = True


if extract:
	op = opticam.Reduction(rawdata=datadir,savedir=catdir,
							name=name,rule='C3*.fits')
	op.sextractor()   		# Perform aperture photometry
	op.creat_ref_list()		# Make master star list & FoV image
	op.photometry()			# Cross-match between all images

target = 20
if analyse:
	photo = opticam.Analysis(catalogue=catdir,name=name)
	photo.differential_photo(target=target,ignore=None)
	photo.rms_mag(target=target)
	photo.lightcurve(std=True)
	photo.ccd_noise()


datadir = 'bl_cam/'
catdir = 'bl_cam_cat/'
name = 'BL_Cam_r'

extract = False
analyse = True


if extract:
	op = opticam.Reduction(rawdata=datadir,savedir=catdir,
							name=name,rule='C2*.fits')
	op.sextractor()   		# Perform aperture photometry
	op.creat_ref_list()		# Make master star list & FoV image
	op.photometry()			# Cross-match between all images

target = 21
if analyse:
	photo = opticam.Analysis(catalogue=catdir,name=name)
	photo.differential_photo(target=target,ignore=None)
	photo.rms_mag(target=target)
	photo.lightcurve(std=True)
	photo.ccd_noise()

#### Change here #####
datadir = 'gk_per/'
catdir = 'gk_per_cat/'
name = 'astro'

extract = False
analyse = True


if extract:
	op = opticam.Reduction(rawdata=datadir,savedir=catdir,
							name=name,rule='C3*.fits')
	op.sextractor()   		# Perform aperture photometry
	op.creat_ref_list()		# Make master star list & FoV image
	op.photometry()			# Cross-match between all images

target = 29
if analyse:
	photo = opticam.Analysis(catalogue=catdir,name=name)
	photo.differential_photo(target=target,ignore=None)
	photo.rms_mag(target=target)
	photo.lightcurve(std=True)
	photo.ccd_noise()