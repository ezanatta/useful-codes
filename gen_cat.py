import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

#this code generates the catalog in the format needed for the fmap procedures.
#basically it uses the fact that 2MASS have a 1:1 pixel scale to convert 
#tracers RA and DEC to something plotable in the fits images from 2MASS.
#
#as inputs you need:
#---> A catalog with the objetcs coordinates: they need to be in arcsecs, 
#so sometimes this code will require some changes on the initial part. 
#Check fig1_1023.py code for further information.
#
#---> Check NED for galaxy center.

RA = np.loadtxt('2768GC-complete_clean_dec.dat', usecols=(1, ))    #add your catalog
DEC = np.loadtxt('2768GC-complete_clean_dec.dat', usecols=(2, ))

DECdeg = DEC*u.degree
DECrad = DECdeg.to(u.rad)
DECrad = DECrad.value

for i in range(0, len(RA)):
	RA[i] = RA[i]*3600
	DEC[i] = DEC[i]*3600

galcenter = SkyCoord('09h11m37.5s', '+60d02m14s')     #real galaxy center

RAg = galcenter.ra.arcsec
DECg = galcenter.dec.arcsec

gal_center_x = 157.167      #galaxy center IN THE FITS IMAGE (pixels)
gal_center_y = 199.519

DRA = []
DDEC = []

for i in range(0, len(RA)):
	DRA.append((RA[i]-RAg)*np.cos(DECdeg[i]))
	DDEC.append((DEC[i]-DECg))

for i in range(0, len(RA)):
	DRA[i] = DRA[i] + gal_center_x
	DDEC[i] = DDEC[i] + gal_center_y

catname = raw_input('Enter desired output catalog name: ')

with open(catname, 'w') as fcat:
	for i in range(0, len(RA)):
		print>>fcat, DRA[i], DDEC[i]




