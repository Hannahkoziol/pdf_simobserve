#!/usr/bin/env python
# coding: utf-8

# In[2]:


import aplpy
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy import units as u


# In[3]:


#Creating histograms


# In[747]:


sim = fits.open('/Users/hannahkoziol/Desktop/feather10_B1_000_clean16.fits')
fig = plt.figure(figsize=(5, 5))
subplot = aplpy.FITSFigure(sim, figure=fig, convention='calabretta')
subplot.show_colorscale(cmap='magma')
cmin, cmax = 23, 26
#nticks = 10
subplot.add_colorbar()

#opening the simulated observation file 


# In[748]:


sim_array = sim[0].data
flat_sim_array = sim_array.flatten()

#taking the data from the file and making an array
#flattening the array to use in a histogram


# In[749]:


def Sr_conversion(bmaj, bmin):
    return ((np.pi * bmaj *bmin)/(4*np.log(2)))/((180/(np.pi))**2)

#defining a function to convert the units from Jy/beam to Jy/Sr to get column density from flux


# In[750]:


data, header = fits.getdata("/Users/hannahkoziol/Desktop/feather10_B1_000_clean16.fits", header=True)
bmaj_sim = header['BMAJ']
bmin_sim = header['BMIN']

#grabbing bmaj and bmin values directly from fits file header


# In[751]:


conversion = Sr_conversion(bmaj_sim, bmin_sim)
converted_array = sim_array / conversion

#applying the converstion to the array 


# In[752]:


sim_header = header.copy()
sim_header['BUNIT'] = 'Jy/Sr'
fits.writeto('converted.fits', converted_array, sim_header, overwrite = True)

#grab header, copy to edit, write new file with updated header


# In[753]:


flux_sim = fits.open('converted.fits')
flux_data = flux_sim[0].data

#opening new converted file and getting data 


# In[754]:


# Image information and assumptions
distance        = 8178. # distance to GC; (GRAVITY collaboration 2019)
Wave            = (3.0e8/226.0e9) # define the wavelength from the SMA dust continuum frequency. 
Wave0           = 1.3e-3
k0              = 0.899
nu              = 3.e08/Wave
nu0             = 3.e08/Wave0
beta            = 1.75
Kappag2d        = k0*((nu/nu0)**beta)
g2d             = 100.0
Kappa           = Kappag2d / g2d # this kappa is the dust opacity, it's kind of complicated where it comes from, but if you're interested I can find you some stuff about it!
mu              = 2.8 # This mu is the mean atomic weight. 
hplanck = 6.626176e-34 #Joule seconds
clight = 3e8 #m/s
kboltzmann = 1.380649e-23 #Joule/Kelvin
mh = 1.6737236e-27

#defining constants and assumptions


# In[755]:


def planck_wave(Wave, Temp):
    planck_conv_wave = 1.e-26 *clight /Wave**2.0
    planck = ((2.0*hplanck*clight**2.0)/(Wave**5.0))*(1.0/(np.exp((hplanck*clight)/(Wave*kboltzmann*Temp))-1.0))
    planck = planck/planck_conv_wave
    return planck
#the modified blackbody curve for dust emission
#input Wave: the wavelength of emission in meters
#input Temp: the assumed dust temperature in kelvin
#returns flux density in Jy/S


# In[756]:


def column_density(Wave, Temp, Kappa, Flux_Density, mu):
    B = planck_wave(Wave, Temp)
    N = Flux_Density / (mu * (mh*1.e3) * Kappa * B)
    return N
#calculate the column density of a source based on its flux density.
#input Wave: the wavelength of emission in Meters
#input Temp: the dust temperature in kelvin
#input Kappa: the dust opacity
#input Flux_Density: the flux density of the source in Jy/Sr
#input mu: the assumed mean molecular weight.
#returns column density in cm^-2


# In[757]:


column_density_array = column_density(Wave, 25, Kappa, flux_data, mu)
flat_column_density = column_density_array.flatten()

#running the function and making a new array of data


# In[758]:


data, header = fits.getdata("/Users/hannahkoziol/Desktop/N_dens_B1snap_000.fits", header=True)
header2 = header.copy()
fits.writeto('/Users/hannahkoziol/Desktop/input_flux_column_density.fits', column_density_array, header2, overwrite = True)

#creating column density file from the cleaned image 


# In[759]:


fig, ax = plt.subplots()
plt.hist(flat_column_density, bins=50, color='firebrick', alpha=0.5)
plt.xlabel('Column Density $\#/cm^2$')
plt.ylabel('Frequency')
plt.show()

#histogram in linear space


# In[760]:


fig, ax = plt.subplots()
hist2 = plt.hist(flat_column_density, bins=np.logspace(22, 23.4, 30), log=True, color='firebrick', alpha=0.5)
ax.set_xscale('log')
plt.xlabel('$Log_{10}$ Column Density $\#/cm^2$')
plt.ylabel('$Log_{10} Frequency$')
plt.show()

#histogram in log10 space


# In[761]:


original_sim = fits.open('/Users/hannahkoziol/Desktop/N_dens_B1snap_000.fits')
fig = plt.figure(figsize=(5, 5))
subplot = aplpy.FITSFigure(original_sim, figure=fig, convention='calabretta')
subplot.show_colorscale(cmap='magma')
cmin, cmax = 23, 26
#nticks = 10
subplot.add_colorbar()

#opening the original simulation snapshot, before running simobserve and tclean


# In[762]:


original_data = original_sim[0].data * 10 
#multiplying by 10 for every flux image that was multiplied by 10

column_density_array2 = original_data/(1e4)
#factor of 1e4 is corrected


# In[763]:


data, header = fits.getdata("/Users/hannahkoziol/Desktop/N_dens_B1snap_000.fits", header=True)
header2 = header.copy()
fits.writeto('final2.fits', column_density_array2, header2, overwrite = True)
original_converted_sim = fits.open('final2.fits')

#writing a converted fits file for the original snapshot


# In[764]:


flat_column_density2 = column_density_array2.flatten()


# In[765]:


fig, ax = plt.subplots()
plt.hist(flat_column_density2, bins=60, color='firebrick', alpha=0.5)
plt.xlabel('Column Density $\#/cm^2$')
plt.ylabel('Frequency')
plt.show()

#original snapshot histogram in linear space


# In[766]:


fig, ax = plt.subplots()
hist2 = plt.hist(flat_column_density2, bins=np.logspace(22.5, 23.5, 25), log=True, color='firebrick', alpha=0.5)
ax.set_xscale('log')
plt.xlabel('$Log_{10}$ Column Density $\#/cm^2$')
plt.ylabel('$Log_{10} Frequency$')
plt.show()

#original snapshot histogram in log10 space


# In[767]:


fig, ax = plt.subplots()
plt.hist(flat_column_density2, bins=60, label = "Before Simobserve", color='firebrick', alpha=0.5)
plt.hist(flat_column_density, bins=60, label = "After Simobserve", color='cornflowerblue', alpha=0.5)
plt.xlabel('Column Density $\#/cm^2$')
plt.ylabel('Frequency')
plt.legend()
plt.show()

#combined histograms in linear space


# In[768]:


fig, ax = plt.subplots()
plt.hist(flat_column_density2, bins=60, density = True, label = "Before Simobserve", color='firebrick', alpha=0.5)
plt.hist(flat_column_density, bins=60, density = True, label = "After Simobserve", color='cornflowerblue', alpha=0.5)
plt.xlabel('Column Density $\#/cm^2$')
plt.ylabel('Frequency')
plt.legend()
plt.show()

#conbined histograms in linear space with density = True


# In[769]:


fig, ax = plt.subplots()
hist2 = plt.hist(flat_column_density, bins=np.logspace(22, 23.7, 30), log=True, label = "After Simobserve", color='cornflowerblue', alpha=0.5)
hist2 = plt.hist(flat_column_density2, bins=np.logspace(22, 23.7, 30), log=True, label = "Before Simobserve", color='firebrick', alpha=0.5)
ax.set_xscale('log')
plt.xlabel('$Log_{10}$ Column Density $\#/cm^2$')
plt.ylabel('$Log_{10} Frequency$')
plt.legend()
plt.show()

#combined histograms in log10 space


# In[770]:


fig, ax = plt.subplots()
hist2 = plt.hist(flat_column_density, bins=np.logspace(22, 23.5, 30), density=True, log=True, label = "After Simobserve", color='cornflowerblue', alpha=0.5)
hist2 = plt.hist(flat_column_density2, bins=np.logspace(22, 23.5, 30), density=True, log=True, label = "Before Simobserve", color='firebrick', alpha=0.5)
ax.set_xscale('log')
plt.xlabel('$Log_{10}$ Column Density $\#/cm^2$')
plt.ylabel('$Log_{10} Frequency$')
plt.legend()
plt.show()

#combined historgams in log10 space with density = True 


# In[76]:


#Creating amplified flux files


# In[14]:


flux_base = fits.open('/Users/hannahkoziol/Desktop/flux_B1snap_051.fits')
flux_array1 = flux_base[0].data
flux_array2 = flux_array1 * 10

#beginning of creating the amplified flux files
#taking the original flux data and multiplying by a factor of 10


# In[15]:


data, header = fits.getdata("/Users/hannahkoziol/Desktop/flux_B1snap_051.fits", header=True)
header3 = header.copy()
#fits.writeto('/Users/hannahkoziol/Desktop/flux10_B1snap_051.fits', flux_array2, header3, overwrite = True)

#creating the flux*10 file


# In[77]:


#Converting column density snapshots to flux


# In[69]:


column_dens_file = fits.open('/Users/hannahkoziol/Desktop/N_dens_B10snap_000.fits')
N_data = column_dens_file[0].data / (1e4)

#converting the column density snapshots into flux files to run the simobserve/tclean code on 
#grabbing column density data to convert


# In[71]:


pix_scale_cm = 3.08e18*5.0/column_dens_file[0].header['NAXIS1']
pc2cm = 3.086e18
distance = 8178.

#useful conversion constants


# In[73]:


def flux(Wave, Temp, Kappa, N, mu):
    B = planck_wave(Wave, Temp)
    flux_density = (N*(mu * (mh*1.e3) * Kappa * B))
    return flux_density
flux_density_jysr = flux(Wave, 25, Kappa, N_data, mu)
flux_density_jydeg2 = flux_density_jysr / 3282.8
flux_density_array = flux_density_jydeg2*((pix_scale_cm/(pc2cm*distance))/(2*np.pi/360.))**2

#converstion for column density to flux in jy/pixel


# In[74]:


#fits.writeto('flux_B1snap_000.fits', flux_density_array, header3, overwrite = True)

#writing the new flux file


# In[78]:


#Creating convolution fits files, Bolocam 


# In[70]:


original_flux = fits.open('/Users/hannahkoziol/Desktop/flux10_B10snap_066.fits')
original_flux_array = original_flux[0].data

#Grabbing data from the snapshot the Bolocam image is being created from


# In[71]:


gc_dist_cm = 8178 * 3.08567758e18
pix_scale_cm = 3.08e18*5.0/original_flux[0].header['NAXIS1']
#useful constants for calculating the pixel value for the standard deviation

def std(fwhm):
    standard_deviation = fwhm / (2*(np.sqrt(2*np.log(2))))
    return standard_deviation
#creating a function to get the standard deviation from the fwhm

def pixel_std(std, gc_dist, conv):
    length_cm = (std * (np.pi / (180*3600))) * gc_dist
    pixel_length = length_cm / conv
    return pixel_length
#creating a function for the pixel value of the standard deviation


# In[72]:


bolocam_fwhm_arcsec = 33
#bolocam full width half max value in arcsecs

bolocam_std = std(bolocam_fwhm_arcsec)

bolo_pixel_std = pixel_std(bolocam_std, gc_dist_cm, pix_scale_cm)
print(bolo_pixel_std)

#filling in values and using functions with bolocam


# In[73]:


convolved_flux_array = convolve_fft(original_flux_array, Gaussian2DKernel(bolo_pixel_std), psf_pad=False, 
                                    fft_pad=False, allow_huge=True)

#fits.writeto('/Users/hannahkoziol/Desktop/bolocam_flux10_B10_066.fits', convolved_flux_array, header3, 
#             overwrite = True)

#writing the bolocam fits file
print(np.max(convolved_flux_array))


# In[83]:


#converting Jy/beam to Jy


# In[430]:


def pixel_conversion(bmaj, bmin, cdelt1, cdelt2):
    return ((np.pi * bmaj *bmin)/(4*np.log(2)))/(cdelt1*cdelt2)

#returns units of pixel/beam


# In[431]:


conversion2 = - pixel_conversion(1.271880600188E-03, 7.100287410948E-04, -1.361111111111E-04, 1.361111111111E-04)
print(conversion2)

#inserting values and testing conversion


# In[534]:


sim_before = fits.open('/Users/hannahkoziol/Desktop/flux10_B1snap_000_sim_dirty.fits')
sim_array = sim_before[0].data
converted_sim_array = (sim_array / conversion2) * 1000
data, header = fits.getdata("/Users/hannahkoziol/Desktop/flux10_B1snap_000_sim_dirty.fits", header=True)
header4 = header.copy()
header4['BUNIT'] = 'mJy '

#opening the simulation that needs to be converted
#creating new fits file header


# In[535]:


fits.writeto('/Users/hannahkoziol/Desktop/flux10_B1snap_000_sim_dirty.fits', converted_sim_array, header4, overwrite = True)

#writing over the file with the new values


# In[98]:


#Jy to Jy/beam


# In[74]:


bmaj = 33 / 3600
bmin = 33 / 3600

#defining bolocam beam


# In[75]:


def beam_conversion(bmaj, bmin, cdelt1, cdelt2):
    return (cdelt1*cdelt2)/((np.pi * bmaj *bmin)/(4*np.log(2)))

#returns units of Jy/beam


# In[76]:


conversion3 = - beam_conversion(bmaj, bmin, -1.361111111111E-04, 1.361111111111E-04)
print(1/conversion3)

#implementing and testing


# In[77]:


sim_jy = fits.open('/Users/hannahkoziol/Desktop/bolocam_flux10_B10_066.fits')
jy_array = sim_jy[0].data
converted_beam_array = jy_array / conversion3
data, header = fits.getdata("/Users/hannahkoziol/Desktop/bolocam_flux10_B10_066.fits", header=True)
header5 = header.copy()
header5['BMAJ'] = bmaj
header5['BMIN'] = bmin
header5['BUNIT'] = 'Jy/beam'


# In[78]:


#fits.writeto('/Users/hannahkoziol/Desktop/bolocam_flux10_B10_066.fits', converted_beam_array, header5, overwrite = True)


# In[79]:


print(np.max(converted_beam_array))


# In[16]:


#working with residual files


# In[122]:


residual = fits.open('/Users/hannahkoziol/Desktop/flux10_B1snap_000_residual.fits')
residual_array = residual[0].data
residual_mean = np.mean(residual_array)
residual_median = np.median(residual_array)
print(residual_mean)
print(residual_median)
print(np.log10(residual_mean))


# In[75]:


fig, ax = plt.subplots()
hist2 = plt.hist(flat_column_density, bins=np.logspace(22, 23.6, 30), log=True, label = "After Simobserve", color='cornflowerblue', alpha=0.5)
hist2 = plt.hist(flat_column_density2, bins=np.logspace(22, 23.6, 30), log=True, label = "Before Simobserve", color='firebrick', alpha=0.5)
ax.set_xscale('log')
plt.xlabel('$Log_{10}$ Column Density $\#/cm^2$')
plt.ylabel('$Log_{10} Frequency$')
plt.plot((np.log10(residual_mean), np.log10(residual_mean)), (1, 10000), label = 'mean sigma')
plt.legend()
plt.show()

#combined histograms in log10 space


# In[409]:


clean3 = fits.open('/Users/hannahkoziol/Desktop/flux10_B10snap_000_clean3_residual.fits')
fig = plt.figure(figsize=(5, 5))
subplot = aplpy.FITSFigure(clean3, figure=fig, convention='calabretta')
subplot.show_colorscale(cmap='magma')
cmin, cmax = 23, 26
#nticks = 10
subplot.add_colorbar()


# In[410]:


clean6 = fits.open('/Users/hannahkoziol/Desktop/flux10_B10snap_000_clean6_residual.fits')
fig = plt.figure(figsize=(5, 5))
subplot = aplpy.FITSFigure(clean6, figure=fig, convention='calabretta')
subplot.show_colorscale(cmap='magma')
cmin, cmax = 23, 26
#nticks = 10
subplot.add_colorbar()


# In[412]:


clean3_data = clean3[0].data
clean6_data = clean6[0].data


# In[414]:


difference = clean3_data - clean6_data
data, header = fits.getdata("/Users/hannahkoziol/Desktop/flux10_B10snap_000_clean3_residual.fits", header=True)
header2 = header.copy()
fits.writeto('/Users/hannahkoziol/Desktop/residual_difference.fits', difference, header2, overwrite = True)


# In[415]:


clean6 = fits.open('/Users/hannahkoziol/Desktop/residual_difference.fits')
fig = plt.figure(figsize=(5, 5))
subplot = aplpy.FITSFigure(clean6, figure=fig, convention='calabretta')
subplot.show_colorscale(cmap='magma')
cmin, cmax = 23, 26
#nticks = 10
subplot.add_colorbar()


# In[ ]:




