from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np


#this code plots all tracers, colored, on a fits file provided. 
#                INPUTS and INFO:
#----> Up to three catalogs: stars, PNe and GCs, and a fits file ( one big enough to fit every tracer inside).
#
#----> Most of time the catalogs I used had different formats and units. 
#      The PN.S catalogs are in this fashion: hh:mm:ss. So I need to load them as strings,
#      use SkyCoord to load these strings as RA and DEC and correct units do decimal degrees.
#----->GC catalogs from Vicenzo Pota are already in decimal degrees, so I load them directly with numpy.
#----->Star catalogs will differ (I guess) between each other, but one of this two methods will work for any case.
#-----> Of course you can use Topcat or any other way to correct the coordinates also. 
#       The important thing is to have everything in decimal degrees so astropy.fits can understand it.



RA = np.loadtxt('2768mag.dat', usecols=(2,))           #GC catalog 
DEC = np.loadtxt('2768mag.dat', usecols=(3,))

with open('star2768.dat', 'r') as f:
    RAstar = [line.split()[1] for line in f]               #star catalog
with open('star2768.dat', 'r') as f:    
    DECstar = [line.split()[2] for line in f]

with open('2768cat.dat', 'r') as f:
    RApne = [line.split()[0] for line in f]              #PNe catalog
with open('2768cat.dat', 'r') as f:    
    DECpne = [line.split()[1] for line in f]

for i in range(0, len(RAstar)):
	c = SkyCoord(RAstar[i], DECstar[i])
	RAstar[i], DECstar[i] = c.ra.degree, c.dec.degree

for i in range(0, len(RApne)):
	c = SkyCoord(RApne[i], DECpne[i])
	RApne[i], DECpne[i] = c.ra.degree, c.dec.degree

coords = zip(RA, DEC)
coords_star = zip(RAstar, DECstar)
coordsPNE = zip(RApne, DECpne)

fits_file = 'fits-null.fits'
hdu = fits.open(fits_file)[0]
wcs = WCS(hdu.header)

ra = []
dec = []
ra_star = []
dec_star = []
ra_pne = []
dec_pne = []

pix = wcs.wcs_world2pix(coords, 1)                          #here we transform the RA and DEC coords in pixel coords 
pix_star = wcs.wcs_world2pix(coords_star, 1)                #to overplot on the fits file in the right places. 
pix_pne = wcs.wcs_world2pix(coordsPNE, 1)                   #This is done using the WCS loaded before.
                                                            #so it is important that the fits you used had the WCS.

#if you really need to use a fits file without WCS (like when you plot the f-map), check how it is done in gen_cat.py

for i in range(0,len(pix)):
    ra.append(pix[i][0])
    dec.append(pix[i][1])
for i in range(0,len(pix_star)):
    ra_star.append(pix_star[i][0])
    dec_star.append(pix_star[i][1])
for i in range(0, len(pix_pne)):
    ra_pne.append(pix_pne[i][0])
    dec_pne.append(pix_pne[i][1])    


allgc = plt.subplot(1,1,1, projection=wcs)
plt.imshow(hdu.data, origin='lower',cmap='gray_r', label='All GCs') #you can change the contrast in the cmap options.
plt.plot(ra, dec, marker='o',markerfacecolor='None', linestyle='none', markeredgewidth=3,markeredgecolor='magenta')
plt.plot(ra_star, dec_star, marker='s', color='orange', markeredgewidth=1, markeredgecolor='orange', linestyle='none')
plt.plot(ra_pne, dec_pne, marker='o', markeredgewidth=2,markeredgecolor='green',markerfacecolor='None',linestyle='none')
plt.title('All Tracers')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.savefig('allgcs-vmap.png', dpi=300)
plt.show()
