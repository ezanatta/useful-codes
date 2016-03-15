def binning(RA, DEC, V, galcenter, pa, i, nbins):
     import numpy as np
     from astropy import units as u
     from astropy.coordinates import SkyCoord

     pa_rad = pa.to(u.rad)
     pa_rad = pa_rad.value

     i_rad = i.to(u.rad)
     i_rad = i_rad.value

     n = len(RA)

     ra = np.linspace(0,0,n)
     dec = np.linspace(0,0,n)

     for i in range(0,n):
	  aux = SkyCoord(RA[i], DEC[i])             #reading RA 	and DEC from GC
	  ra[i] = aux.ra.arcsec
	  dec[i] = aux.dec.arcsec
	
     y_rad = galcenter.dec.radian   # galaxy center DEC 		coordinate in radians
        
     xmc = (ra - galcenter.ra.arcsec)
     xm = xmc*(np.cos(y_rad))                       # setting 	the x coordinates of GC from its RA value relative to 	galactic center      
     ym = (dec - galcenter.dec.arcsec)             # setting 	the y coordinates of the GC from its DEC value relative to 	galactic center
	
     xs = xm*np.cos(pa_rad)-ym*np.sin(pa_rad)            # 	rotating coordinates 
     ys = xm*np.sin(pa_rad)+ym*np.cos(pa_rad)    
	
     cos_i = np.cos(i_rad)               # cosine of 	inclination angle, used to shrink coordinates
	
     r=np.sqrt(((xs**2)*cos_i)+((ys**2)/cos_i))   #distributing 	the GC along the radius
	
     ###binning ####
     t = 0
     l = 0
     n = r.size   #lenght of r --- number of objects to iterate on the following loops
     for i in range(0, n):
         for j in range(i+1,n):
             if r[i] > r[j]:
               			t = r[i]
               			r[i] = r[j] #ordering r and V
               			r[j] = t
				l = V[i]
				V[i] = V[j]
				V[j] = l
	
     rgal = np.linspace(0,0,nbins+1)       #rgal contains the limits of each bin
     for h in range(1, nbins):
    		rgal[h] = r[int(n*h/nbins)]     #binning
		rgal[nbins] = r[n-1]
		rgal[0] = r[0]
	
     return r, rgal, V 
 
def density(nbin, r, rgal, ):
    import numpy as np
    import uncertainties.unumpy as unp
    n = len(r)
    NGC = np.zeros(nbin)    #number of GC in each bin
    rmed = [[] for i in range(0,nbin)]
    for h in range(0, nbin):                              
         for i in range(0, n):
             if (rgal[h+1] >= r[i] and r[i] > rgal[h]):
                 NGC[h] = NGC[h] + 1
                 rmed[h].append(r[i])
             
    #median value of the bins
    radius = [[] for i in range(0,nbin)]
    for i in range(0,nbin):
            radius[i] = np.linspace(rgal[i], rgal[i+1], 1000)   

    median = np.linspace(0,0,nbin)
    medianGC = np.linspace(0,0,nbin)
    for i in range(0,nbin):
        median[i] = np.median(radius[i])
        medianGC[i] = np.median(rmed[i])

    area = []
    for h in range(0, nbin):
        area.append((np.pi*rgal[h+1]**2)-(np.pi*rgal[h]**2)) #area per bin
    
    binsize = np.linspace(0,0,nbin)    
    for h in range(0,nbin):
        binsize[h] = rgal[h+1]-rgal[h]  
        binsize = binsize/2
    
    dens2 = np.linspace(0,0,nbin)
    for h in range(0, nbin):
        dens2[h] = NGC[h]/area[h]

    
    #error in density

    poi_err = (np.sqrt(NGC)/area)
    dens2err = unp.uarray(dens2, poi_err)
    denslog2 = -2.5*np.log10(dens2)    
    denslog2err = -2.5*unp.log10(dens2err)  
    poi_err2 = unp.std_devs(denslog2err)
    
    return denslog2, poi_err2, area, NGC, median, rgal, binsize 
    