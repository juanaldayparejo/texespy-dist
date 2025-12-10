import numpy as np
from struct import *
import sys,os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import texespy as texes

##############################################################################################
##############################################################################################
#                                   GEOMETRY CALIBRATION
##############################################################################################
##############################################################################################

def process_file(fitsname,x_ini_offset=0.0,y_ini_offset=0.0,force_slit_pa=None,MakePlot=True):
    """
    FUNCTION NAME : process_file()

    DESCRIPTION : Process the FITS file from a TEXES observation

    INPUTS : 

        fitsname :: Name of the FITS file

    OPTIONAL INPUTS:
    
        x_ini_offset :: Offset to apply to the TEXES image in the x direction (East is negative) (arcsec)
        y_ini_offset :: Offset to apply to the TEXES image in the y direction (East is negative) (arcsec)
        force_slit_pa :: Force the slit position angle to a specific value (deg)
        MakePlot :: Make a plot of the geometry of the observation
            
    OUTPUTS : 
 
        geometry :: Dictionary with the geometry of the observation
        x_vals(nscans,nslitpix) :: Angular position of each pixel across the x direction (E is negative) (arcsec)
        y_vals(nscans,nslitpix) :: Angular position of each pixel across the y direction (N is positive) (arcsec)
        waven(nwave,norder) :: Wavenumber array for each diffraction order (cm-1)
        data(nscans,nslitpix,nwave,norder) :: Spectra in each pixel and diffraction order (erg s-1 cm-2 sr-1 (cm-1)-1)

    CALLING SEQUENCE:

        geometry,x_vals,y_vals,waven,data = process_file(fitsname)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    """
    
    from astropy.io import fits
    
    #Reading the FITS file
    ####################################################################################################
    
    # Open the FITS file
    hdul = fits.open(fitsname)
    
    header = hdul[0].header
    planet = header['OBJECT']
    
    #Reading date and time of the observations
    date_obs = header['DATE-OBS']   #Date of observation (UTC)
    time_obs = header['TIME-OBS']   #Time of observation (UTC)
    
    #Reading the spectral data and main array sizes
    data = hdul[0].data    #(NSCANS,NSLITPIX,NWAVE)
    datasum = np.sum(data,axis=2)
    waveimg = hdul[1].data    
    sky = hdul[2].data     #(NSLITPIX,NWAVE)
    noise = hdul[3].data   #(NSLITPIX,NWAVE)
        
    nwaveor = header['NSPEC']           #Number of wavelengths in each diffraction order
    nscans = data.shape[0]              #Number of scans
    nslitpix = data.shape[1]            #Number of pixel across the slit
    nwavetot = data.shape[2]            #Total number of spectral points
    norders = int(nwavetot / nwaveor)   #Number of diffraction orders
    
    #Calculating wavenumbers
    waven = np.zeros(nwavetot)
    trans_est = np.zeros(nwavetot)
    for i in range(nwavetot):
        waven[i] = waveimg[i][0] 
        trans_est[i] = waveimg[i][3]    #Estimated transmission at each wavelength
    
    #Reading important parameters for the geometry
    ini_slit_pos = header['OFFSET']
    slit_step = header['TELSTEP']
    plate_scale = header['PLTSCL']   #Plate scale
    slit_width = header['SLITWID']
    slit_height = plate_scale * nslitpix
    slit_angle = header['INSTRPA']
    resolving_power = header['RESOLV']
    
    if force_slit_pa is not None:
        slit_angle = force_slit_pa
            
    hdul.close()
    
    #Calculating the geometry of the observation from spice
    ####################################################################################################
    
    #SPICE inputs
    body = planet.upper()
    ref = 'IAU_'+body
    abcorr = 'LT+S'

    #Mauna Kea
    lat_observer = 19.8263
    lon_observer = -155.473
    alt_observer = 4205.

    ra,dec,ra_np,dec_np = texes.geometry.calc_radec(date_obs,time_obs,lat_observer,lon_observer,alt_observer/1.0e3,body=body,ref=ref,abcorr=abcorr)
    pa = texes.geometry.calc_position_angle(date_obs,time_obs,lat_observer,lon_observer,alt_observer/1.0e3,body=body,ref=ref,abcorr=abcorr)
    emiss_ang_observer = texes.geometry.calc_emission_angle_earth(date_obs,time_obs,lat_observer,lon_observer,alt_observer/1.0e3,body=body,ref=ref,abcorr=abcorr)

    angsize = texes.geometry.calc_angsize_body(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)
    lat_subsol,lon_subsol = texes.geometry.calc_subsol_point(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)
    lat_subobs,lon_subobs = texes.geometry.calc_subobs_point(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)
    phase = texes.geometry.calc_phase_angle(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)
    v_doppler = texes.geometry.calc_vdoppler_earth_body(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)
    dist_sun = texes.geometry.calc_sun_body_dist(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)
    dist_earth = texes.geometry.calc_earth_body_dist(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)


    if MakePlot==True:
        
        #Calculating the geometry
        angmax = float(int(angsize)) ; res = 0.1
        
        rotan = -pa
        hx,hy,lat,lon = texes.geometry.xy2latlon(angmax,res,angsize,lat_subobs,lon_subobs,rot_angle=rotan)
        emiss_ang = texes.geometry.latlon2emiss(lat,lon,lat_subobs,lon_subobs)
        sol_ang = texes.geometry.latlon2sza(lat,lon,lat_subsol,lon_subsol)
        azi_ang = texes.geometry.phase2azi(emiss_ang,sol_ang,phase)
        lst = texes.geometry.lon2lst(lon, lon_subsol)
        
        
        #Plotting some parameters
        fig,([ax1,ax2],[ax3,ax4]) = plt.subplots(2,2,figsize=(8,6))

        cmap = 'turbo'

        nlevels = 51
        
        im1 = ax1.contourf(hx,hy,lat,cmap=cmap,levels=np.linspace(-90.,90.,nlevels),vmin=-90.,vmax=90.)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='black')
        ax1.add_artist( Drawing_uncolored_circle )
        texes.geometry.add_compass(ax1)
        ax1.grid()

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im1, cax=cax)
        cbar.set_label('Latitude ($^\circ$)')

        im2 = ax2.contourf(hx,hy,lst,cmap=cmap,levels=np.linspace(0.,24.,nlevels),vmin=0.,vmax=24.)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='black' )
        ax2.add_artist( Drawing_uncolored_circle )
        texes.geometry.add_compass(ax2)
        ax2.grid()

        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im2, cax=cax)
        cbar.set_label('LST (hour)')


        im3 = ax3.contourf(hx,hy,emiss_ang,cmap=cmap+'_r',levels=np.linspace(0.,90.,nlevels),vmin=0.,vmax=90.)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='black' )
        ax3.add_artist( Drawing_uncolored_circle )
        texes.geometry.add_compass(ax3)
        ax3.grid()

        divider = make_axes_locatable(ax3)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im3, cax=cax)
        cbar.set_label('Emission angle ($^\circ$)')

        im4 = ax4.contourf(hx,hy,sol_ang,cmap=cmap+'_r',levels=np.linspace(0.,180.,nlevels),vmin=0.,vmax=180.)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='black' )
        ax4.add_artist( Drawing_uncolored_circle )
        texes.geometry.add_compass(ax4)
        ax4.grid()

        divider = make_axes_locatable(ax4)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im4, cax=cax)
        cbar.set_label('Solar zenith angle ($^\circ$)')

        plt.tight_layout()
                
        print('SPICE PARAMETERS')
        print('****************************************************************************************')
        print('Date of observation = ',date_obs)
        print('Time of observation = ',time_obs)
        print('Sun-planet distance = ',dist_sun,'AU')
        print('Earth-planet distance = ',dist_earth,'AU')
        print('Planet angular size = ',angsize,'arcsec')
        print('Sub-observer latitude = ',lat_subobs,'degrees','   ','Sub-observer longitude = ',lon_subobs,'degrees')
        print('Sub-solar latitude = ',lat_subsol,'degrees','   ','Sub-solar longitude = ',lon_subsol,'degrees')
        print('Doppler velocity between Earth and planet :: ',v_doppler,'km s-1')
        print('Phase angle :: ',phase,'deg')
        print('Right ascension :: ',ra,'degrees','   ','Declination :: ',dec,'deg')
        print('Right ascension of North Pole :: ',ra_np,'degrees','   ','Declination of North Pole :: ',dec_np,'deg')
        print('Position angle :: ',pa,'deg')
        print('Emission angle of observer on Earth :: ',emiss_ang_observer,'deg')
        print('')
        
        print('TEXES PARAMETERS')
        print('****************************************************************************************')
        print('Slit width = ',slit_width,'arcsec')
        print('Plate scale = ',plate_scale,'arcsec/pixel')
        print('Slit height = ',slit_height,'arcsec')
        print('Slit position angle = ',slit_angle,'deg')
        print('Initial offset = ',ini_slit_pos)
        print('Slit step = ',slit_step)
        print('')



    # Creating the geometry dictionary
    geometry = {
        "date_obs": date_obs,  #date of the observation
        "time_obs": time_obs,   #time of the observation
        "ra": ra,    #right ascension of planet centre (deg)
        "dec": dec,  #declination of planet centre (deg)
        "ra_np": ra_np,   #right ascension of north pole (deg)
        "dec_np": dec_np, #declination of north pole (deg)
        "pa": pa, #position angle (deg)
        "emiss_ang_observer": emiss_ang_observer, #emission angle from Earth (deg)
        "angsize": angsize, #Planet angular size (deg)
        "lat_subsol": lat_subsol, #sub-solar latitude (deg)
        "lon_subsol": lon_subsol, #sub-solar longitude (deg)
        "lat_subobs": lat_subobs, #sub-observer latitude (deg)
        "lon_subobs": lon_subobs, #sub-observer longitude (deg)
        "phase": phase, #phase angle (deg)
        "v_doppler": v_doppler, #Doppler velocity (km s-1)
        "dist_sun": dist_sun, #Planet-Sun distance (AU)
        "dist_earth": dist_earth #Planet-Earth distance (AU)
    }
    
    #Computing first guess of the geometry
    ####################################################################################################
    
    x_ini,y_ini = extract_values_telescope(ini_slit_pos)
    x_step,y_step = extract_values_telescope(slit_step)

    y_ini -= slit_height / 2.
    
    x0_scans = np.arange(0,nscans,1) * x_step + x_ini
    y0_scans = np.arange(0,nscans,1) * y_step + y_ini

    x_vals = np.zeros((nscans,nslitpix))
    y_vals = np.zeros((nscans,nslitpix))

    for i in range(nslitpix):

        x_vals[:,i] = x0_scans - (i * plate_scale) * np.sin(slit_angle/180.*np.pi)
        y_vals[:,i] = y0_scans + (i * plate_scale) * np.cos(slit_angle/180.*np.pi)


    #Filtering the values the pixels across the slit with zero radiance
    inotnan = np.where(datasum[int(nscans/2),:]!=0.0)[0]

    x_vals = x_vals[:,inotnan] + x_ini_offset
    y_vals = y_vals[:,inotnan] + y_ini_offset
    data = data[:,inotnan,:]
    sky = sky[inotnan,:]
    noise = noise[inotnan,:]
    
    #Removing the first and last pixel of the slit
    npixslit = x_vals.shape[1]
    x_vals = x_vals[:,1:npixslit-2]
    y_vals = y_vals[:,1:npixslit-2]
    data = data[:,1:npixslit-2,:]
    sky = sky[1:npixslit-2,:]
    noise = noise[1:npixslit-2,:]
    
    data_order = np.zeros((data.shape[0],data.shape[1],nwaveor,norders))
    waven_order = np.zeros((nwaveor,norders))
    sky_order = np.zeros((sky.shape[0],nwaveor,norders))
    noise_order = np.zeros((noise.shape[0],nwaveor,norders))
    trans_est_order = np.zeros((nwaveor,norders))
    ix = 0
    for iorder in range(norders):
        data_order[:,:,:,iorder] = data[:,:,ix:ix+nwaveor]
        waven_order[:,iorder] = waven[ix:ix+nwaveor]
        sky_order[:,:,iorder] = sky[:,ix:ix+nwaveor]
        noise_order[:,:,iorder] = noise[:,ix:ix+nwaveor]
        trans_est_order[:,iorder] = trans_est[ix:ix+nwaveor]
        ix += nwaveor
    
    return geometry,x_vals,y_vals,waven_order,data_order,sky_order,noise_order,trans_est_order,resolving_power
    

##############################################################################################

def calibrate_geometry(geometry,x_vals,y_vals,waven,image,fitting_method='radiance',max_radiance=None,print_info=False,MakePlot=False,max_off=4.,res_off=0.1,convol=True):
    """
    FUNCTION NAME : calibrate_geometry()

    DESCRIPTION : Calibrate the geometry of the observation

    INPUTS : 

        geometry :: Dictionary with the geometrical parameters calculated with SPICE with process_file()
        x_vals(nscans,npixslit) :: Angular position of each pixel across x-axis as specified in FITS files
        y_vals(nscans,npixslit) :: Angular position of each pixel across y-axis as specified in FITS files
        waven :: Wavenumber of the image (cm-1)
        image(nscans,npixslit) :: Image to calibrate

    OPTIONAL INPUTS:
    
        fitting_method :: This string indicates the method that will be used to calibrate the geometry.
        
                          'sza' == The radiance from the observations will be compared to the cos(sza)
                          'uniform' == The radiance from the observations will be compared to a uniform disk 
                          
        max_off :: Maximum offset (arcsec)
        res_off :: Resolution of the offset steps used in the fitting process (arcsec)
            
    OUTPUTS : 
 
        x_ini_offset :: Offset to apply to the TEXES image in the x direction (East is negative) (arcsec)
        y_ini_offset :: Offset to apply to the TEXES image in the y direction (East is negative) (arcsec)

    CALLING SEQUENCE:

        x_ini_offset,y_ini_offset = calibrate_geometry(geometry,x_vals,y_vals,waven,image)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    """
    
    from astropy.convolution import convolve, Gaussian2DKernel
    from scipy.interpolate import RectBivariateSpline
    
    #Calculating the geometry of the observation using the geometry parameters
    ####################################################################################################
    
    date_obs = geometry['date_obs'] ; time_obs = geometry['time_obs']
    angsize = geometry['angsize'] ; phase = geometry['phase'] ; pa = geometry['pa']
    lat_subobs = geometry['lat_subobs'] ; lat_subsol = geometry['lat_subsol'] ; lon_subobs = geometry['lon_subobs'] ; lon_subsol = geometry['lon_subsol']
    
    angmax = float(int(angsize)) ; res = 0.1
    
    rotan = -pa
    hx,hy,lat,lon = texes.geometry.xy2latlon(angmax,res,angsize,lat_subobs,lon_subobs,rot_angle=rotan)
    emiss_ang = texes.geometry.latlon2emiss(lat,lon,lat_subobs,lon_subobs)
    sol_ang = texes.geometry.latlon2sza(lat,lon,lat_subsol,lon_subsol)
    azi_ang = texes.geometry.phase2azi(emiss_ang,sol_ang,phase)
    lst = texes.geometry.lon2lst(lon, lon_subsol)
            
    #Creating the model to fit the geometry
    ##############################################################################################################
            
    if fitting_method == 'sza':
        
        bsurf = np.cos( sol_ang / 180. * np.pi )
        bsurf[np.isnan(bsurf)==True] = 0.0
        
    elif fitting_method == 'uniform':
        
        sza = np.cos( sol_ang / 180. * np.pi )
        bsurf = np.ones(sza.shape)
        bsurf[np.isnan(sza)==True] = 0.0
        
    if convol is True:
        
        #Convolving the image with a 2D gaussian
        fwhm_texes = 1.  #arcsec
        sigma_texes = fwhm_texes / 2.355  # σ = FWHM / 2.355
        
        # Convert sigma from arcsec to pixels
        sigma_pixels = sigma_texes / res
        
        # Create 2D Gaussian kernel
        gaussian_kernel = Gaussian2DKernel(x_stddev=sigma_pixels)
        
        #Convolving image
        bsurf_conv = convolve(bsurf, gaussian_kernel, normalize_kernel=True)
        
    else:
        
        bsurf_conv = bsurf    
        
    if max_radiance is not None:
        bsurf_conv = bsurf_conv / bsurf_conv.max() * max_radiance
    
        
    if MakePlot is True:
        
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,3))
        
        cmap = 'turbo'
        
        im = ax1.contourf(hx,hy,bsurf*1.0e7,cmap=cmap,levels=21)
        ax1.set_xlabel('Angular distance (arcsec)')
        ax1.set_ylabel('Angular distance (arcsec)')
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax1.add_artist( Drawing_uncolored_circle )
        ax1.set_title('Model image')
        ax1.grid()
        ax1.set_facecolor('lightgray')
        
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('Radiance (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$)')
        
        im = ax2.contourf(hx,hy,bsurf_conv*1.0e7,cmap=cmap,levels=21)
        ax2.set_xlabel('Angular distance (arcsec)')
        ax2.set_ylabel('Angular distance (arcsec)')
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax2.add_artist( Drawing_uncolored_circle )
        ax2.set_title('Model image')
        ax2.grid()
        ax2.set_facecolor('lightgray')
        
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('Radiance (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$)')
                
        plt.tight_layout()
        
        
    #Fitting the geometry of the TEXES image
    ##########################################################################################
    
    xoff = np.arange(-max_off,max_off+res_off,res_off)
    yoff = np.arange(-max_off,max_off+res_off,res_off)

    x = hx[0,:]
    y = hy[:,0]

    s = RectBivariateSpline(x,y,bsurf_conv.T)

    rms = np.zeros((len(xoff),len(yoff)))
    for ioff in range(len(xoff)):
        for joff in range(len(yoff)):

            fm = np.zeros(x_vals.shape)
            for i in range(x_vals.shape[0]):
                for j in range(x_vals.shape[1]):
            
                    fm[i,j] = s(x_vals[i,j]+xoff[ioff],y_vals[i,j]+yoff[joff])[0,0]

            if max_radiance is not None:
                rms[ioff,joff] = np.sum( (fm - image)**2. )
            else:
                rms[ioff,joff] = np.sum( (fm/fm.max() - image/image.max())**2. )
    
    imax,jmax = np.unravel_index(rms.argmin(), rms.shape)
    x_ini_offset = xoff[imax] ; y_ini_offset = yoff[jmax]
    
    if MakePlot==True:
        
        fig,ax1 = plt.subplots(1,1,figsize=(3,3))
        ax1.contourf(xoff,yoff,rms.T,levels=101,cmap='turbo')
        imax,jmax = np.unravel_index(rms.argmin(), rms.shape)
        ax1.scatter(xoff[imax],yoff[jmax],c='tab:red')
        ax1.set_xlabel('x offset (arcsec)')
        ax1.set_ylabel('y offset (arcsec)')
        plt.tight_layout()
        

        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(6,3))

        psize = 5.

        ax1.scatter(x_vals+x_ini_offset,y_vals+y_ini_offset,s=psize,c=image,cmap='turbo',vmin=0.0)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax1.add_artist( Drawing_uncolored_circle )
        ax1.set_title('TEXES map')
        ax1.set_xlabel('Angular distance (arcsec)')
        ax1.set_ylabel('Angular distance (arcsec)')
        ax1.grid()
        texes.geometry.add_compass(ax1)

        lim = angsize
        ax1.set_xlim(-lim,lim)
        ax1.set_ylim(-lim,lim)
        ax1.set_facecolor('lightgray')

        fm = np.zeros(x_vals.shape)
        for i in range(x_vals.shape[0]):
            for j in range(x_vals.shape[1]):
                fm[i,j] = s(x_vals[i,j]+x_ini_offset,y_vals[i,j]+y_ini_offset)[0,0]

        ax2.scatter(x_vals+x_ini_offset,y_vals+y_ini_offset,s=psize,c=fm,cmap='turbo',vmin=0.0)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax2.add_artist( Drawing_uncolored_circle )
        ax2.set_title('Modelled map')
        ax2.set_xlabel('Angular distance (arcsec)')
        ax2.set_ylabel('Angular distance (arcsec)')
        ax2.grid()
        texes.geometry.add_compass(ax2)

        lim = angsize
        ax2.set_xlim(-lim,lim)
        ax2.set_ylim(-lim,lim)
        ax2.set_facecolor('lightgray')

        plt.tight_layout()
    
    return x_ini_offset,y_ini_offset
    
##############################################################################################
##############################################################################################
#                                           UTILS
##############################################################################################
##############################################################################################
    
def extract_values_telescope(value_fits):
    '''
    Change the format of N and E co-ordinates to numerical values (North is positive y, East is negative x).
    Supports values in decimal degrees, optionally followed by a quote (") symbol.
    '''
    
    import re

    # Regex now allows optional trailing double-quote (") after numbers
    pattern = r'([-+]?\d+(?:\.\d+)?)"?\s*([ENWS])'

    matches = re.findall(pattern, value_fits)

    x, y = None, None  # Default to None to distinguish missing values

    for value, direction in matches:
        value = float(value)  # Convert extracted string to float

        if direction == 'E':
            x = -value  # East is negative
        elif direction == 'W':
            x = value  # West is positive
        elif direction == 'N':
            y = value  # North is positive
        elif direction == 'S':
            y = -value  # South is negative

    # Return 0 if x or y is missing
    return (x if x is not None else 0, y if y is not None else 0)