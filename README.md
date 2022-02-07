OPTICAM Exposure Time Calculator
======
This package is used to calculate an exposure time based on a desired signal to noise ratio for OPTICAM `<https://www.southampton.ac.uk/opticam/>`. This code has been modified and adapted specifically for OPTICAM: The original code can be downloaded here: `here <https://apoexposuretimecalculator.github.io/APOExptime/>`_.


##  Section 1:  Usage
Import the four distinct objects inside OPTICAM:

```python
from opticam import Sky, Target, Instrument, Observation
```

### Sky
First, let's create the <strong>Sky</strong> conditions, where lunar phase can go from 0 (New moon) to 1 (Full moon). Seeing is in arcseconds:
```python
sky = Sky(lunar_phase=0.5, seeing=1.5)
```

### Target
Then, let's create our <strong>Target</strong> Object. First entry is the target's magnitude, second entry is a string that specifies the magnitude system of the input magnitude [vegamag/stmag/abmag], third entry is the band pass for entered magnitude (only Johnson filters are supported at the moment). Default object is a blackbody spectrum with a solar temperature T=5778 K. You can choose between a blackbody, a powerlaw F<sub>&lambda;</sub>&Proportional;&lambda;<sup>-index</sup>
```python
#BlackBody T=5000 K
star1 = Target(16.5, 'VEGAMAG', 'V', temp=5000)			
# PowerLaw F<sub>&lambda;</sub>, index = -4
star2 = Target(20.1, 'stmag', 'B', index=4)				
```

You can also upload a custom SED
```python
import numpy as np

# custom_sed.txt has two columns: 1) Wavelenght [AA]; 2) Flux [erg/s/cm^2/AA]
custom_sed = np.loadtxt('custom_sed.txt')				
star3 = Target(13.8, 'abmag', 'B', sed=custom_sed)	
```

### Instrument
Now, let's load the instrument. Only Opticam is supported at the moment.
```python
inst = Instrument('Opticam')
```

### Observation
We can now combine all three objects to generate and observation.
```python
obs = Observation(star, sky, inst)
```
Now, we can use this object to generate either the SNR for a given exposure time or viceversa:
```python
snr_1 = obs.SNfromTime(200) #value in seconds
time_1 = obs.TimefromSN(50) #value in S/N ratio

print(snr_1)
print(time_1)
```
This will generate for every filter in OPTICam the desired SNR (assuming a unique exposure time for all filters), or the exposure time required in each filter to achieve the desired SNR.
```
[[35.7873250731623, 'uprime_filter'], [141.40232164373245, 'gprime_filter'], [182.41241192694474, 'rprime_filter'], [183.60018496695776, 'iprime_filter'], [114.7288335187604, 'zprime_filter']]

[[387.1788766085068, 'uprime_filter'], [25.201236621119214, 'gprime_filter'], [15.179251054637014, 'rprime_filter'], [14.98133450353503, 'iprime_filter'], [38.189538299358965, 'zprime_filter']]
```

##  Section 2:  Plotting
A more convenient way of using this calculator is to export the information as a plot. 
```python
from opticam import makeplots

ob3 = Observation(star, sky, inst)
# 10 second integrations
ob3.SNfromTime(10)   					

# 'SN', will create a SNR plot for the 10 sec exposures
dd = makeplots(ob3, 'SN')				
```
<p align="middle">
 <img src="Examples/SN_plot.png" width="350" height="450" />
</p>

```python
ob4 = Observation(star, sky, inst)
# SNR = 50
ob4.TimefromSN(50)

# 'Time', will create a Exposure time plot for SNR=50
dd = makeplots(ob4, 'Time')					
											
```
<p align="middle">
 <img src="Examples/EXP_plot.png" width="350" height="450" />
</p>

Limitations
------------
This package is  only for the 2.1m telescope at the Observatorio Astronomico Nacional.
The current supported instruments are: OPTICam
