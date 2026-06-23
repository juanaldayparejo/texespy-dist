# NAME:
#       geometry.py (TEXES)
#
# DESCRIPTION:
#
#       This library contains functions to calculate the geeometry of the TEXES measurements
#
# CATEGORY:
#
#       ACS
#
# FUNCTIONS:
#
#       SPICE routines
#       ######################################
#       
#       calc_ra_dec_body()
#
#       MAPPING geometry routines
#       ######################################
#       
#       
#
# MODIFICATION HISTORY: Juan Alday 24/04/2023

import numpy as np
from struct import *
import sys,os,errno,shutil
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.font_manager as font_manager
import matplotlib as matplotlib
import matplotlib as mpl
from texespy.paths import spice_mk
import spiceypy as spice

##############################################################################################
##############################################################################################
#                                       SPICE ROUTINES
##############################################################################################
##############################################################################################

#Loading spice kernels
spice.furnsh(spice_mk)


###############################################################################################

def calc_sun_body_dist(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):

    """
    FUNCTION NAME : calc_sun_body_dist()

    DESCRIPTION : Calculate the distance between the Sun and the observed body  

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        dist_sun_mars :: Sun-Mars distance / AU

    CALLING SEQUENCE:

        dist_sun_body = calc_sun_body_dist(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    et0 = spice.str2et(date_obs+' '+time_obs)

    AU = 1.495978707e11   #Astronomical Unit / m

    positions, lightTime = spice.spkpos(body,et0,ref,abcorr,'SUN')
    dist_sun_mars = np.sqrt( positions[0]**2. + positions[1]**2. + positions[2]**2.)   #Distance between Sun and Mars (km)
    dist_sun_mars = dist_sun_mars * 1.0e3 / AU    #Sun-Mars distance (AU)
    
    return dist_sun_mars

###############################################################################################

def calc_earth_body_dist(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):

    """
    FUNCTION NAME : calc_earth_mars_dist()

    DESCRIPTION : Calculate the distance between the Earth and Mars 

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        dist_earth_mars :: Sun-Mars distance / AU

    CALLING SEQUENCE:

        dist_earth_mars = calc_earth_body_dist(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    et0 = spice.str2et(date_obs+' '+time_obs)
    
    AU = 1.495978707e11   #Astronomical Unit / m

    positions, lightTime = spice.spkpos(body,et0,ref,abcorr,'EARTH')
    dist_earth_mars = np.sqrt( positions[0]**2. + positions[1]**2. + positions[2]**2.)   #Distance between Earth and Mars (km)
    dist_earth_mars = dist_earth_mars * 1.0e3 / AU    #Earth-Mars distance (AU)
    
    return dist_earth_mars

###############################################################################################

def calc_angsize_body(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):

    """
    FUNCTION NAME : calc_angsize_body()

    DESCRIPTION : Calculate the angular size of the observed body

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        calc_angsize_body :: Angular size of the observed body / arcsec

    CALLING SEQUENCE:

        angsize_mars = calc_angsize_body(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    #Calculating Mars radius in each direction
    npos, radii = spice.bodvrd(body, "RADII", 3)
    rx_body  =  radii[0] * 1.0e3
    ry_body  =  radii[1] * 1.0e3
    rz_body  =  radii[2] * 1.0e3
    f_mars   =  ( rx_body - rz_body ) / rx_body   #Flattening coefficient
    
    #Calculating the distance between the body and the Earth
    dist_body_earth = calc_earth_body_dist(date_obs,time_obs,body=body,ref=ref,abcorr=abcorr)
    
    AU = 1.495978707e11   #Astronomical Unit / m
    angsize_body_x = np.arctan(2.0*rx_body / (dist_body_earth*AU)) / (2.0*np.pi) * 360. * 3600.  #arcsec
    angsize_body_y = np.arctan(2.0*ry_body / (dist_body_earth*AU)) / (2.0*np.pi) * 360. * 3600.  #arcsec
    angsize_body_z = np.arctan(2.0*rz_body / (dist_body_earth*AU)) / (2.0*np.pi) * 360. * 3600.  #arcsec
    
    angsize_body = np.mean([angsize_body_x,angsize_body_z])
    
    return angsize_body

###############################################################################################

def calc_subobs_point(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_subobs_point()

    DESCRIPTION : Calculate the sub-observer point on the obsered body centered reference frame

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        lat_subobs :: Sub-observer latitude / degrees
        lon_subobs :: Sub-observer longitude / degrees

    CALLING SEQUENCE:

        lat_subobs,lon_subobs = calc_subobs_point(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    et0 = spice.str2et(date_obs+' '+time_obs)

    spoint,trgepc,srfvec = spice.subpnt("INTERCEPT/ELLIPSOID",body,et0,ref,abcorr,'EARTH')
    
    #Calculating Mars' radius and flattening coefficient
    npos, radii = spice.bodvrd(body, "RADII", 3)
    re  =  radii[0]
    rp  =  radii[2]
    f   =  ( re - rp ) / re

    #Convert Cartesian coordinates to spherical coordinates in IAU_MARS
    lon_subobs,lat_subobs,alt = spice.recgeo(spoint,re,f)
    #radius,lon_subobs,lat_subobs = spice.reclat(spoint)

    lat_subobs = lat_subobs / np.pi * 180.
    lon_subobs = lon_subobs / np.pi * 180.
    
    return lat_subobs,lon_subobs


###############################################################################################

def calc_subsol_point(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_subsol_point()

    DESCRIPTION : Calculate the sub-solar point on the obsered body centered reference frame

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        lat_subsol :: Sub-solar latitude / degrees
        lon_subsol :: Sub-solar longitude / degrees

    CALLING SEQUENCE:

        lat_subsol,lon_subsol = calc_subsol_point(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    
    et0 = spice.str2et(date_obs+' '+time_obs)

    spoint,trgepc,srfvec = spice.subslr("INTERCEPT/ELLIPSOID",body,et0,ref,abcorr,'EARTH')

    #Calculating Mars' radius and flattening coefficient
    npos, radii = spice.bodvrd(body, "RADII", 3)
    re  =  radii[0]
    rp  =  radii[2]
    f   =  ( re - rp ) / re

    #Convert Cartesian coordinates to spherical coordinates in IAU_MARS
    lon_subsol,lat_subsol,alt = spice.recgeo(spoint,re,f)
    #radius,lon_subsol,lat_subsolx = spice.reclat(spoint)

    lat_subsol = lat_subsol / np.pi * 180.
    lon_subsol = lon_subsol / np.pi * 180.
    
    return lat_subsol,lon_subsol

###############################################################################################

def calc_phase_angle(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_phase_angle()

    DESCRIPTION : Calculate the phase angle between the Sun-body and Earth-body vectors

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        phase_angle :: Phase angle / degrees

    CALLING SEQUENCE:

        phase_angle = calc_phase_angle(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    et0 = spice.str2et(date_obs+' '+time_obs)

    #Calculating the Sun-Mars and Earth-Mars vectors
    r_body_sun, ltime_sun_body = spice.spkpos(body,et0,ref,abcorr,'SUN')
    r_body_earth, ltime_earth_body = spice.spkpos(body,et0,ref,abcorr,'EARTH')

    #Calculating the phase angle
    phase = spice.vsep(r_body_earth,r_body_sun) / np.pi * 180.
    
    return phase

###############################################################################################

def calc_vdoppler_earth_body(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_vdoppler_earth_body()

    DESCRIPTION : Calculate the Doppler velocity between the Earth and the observed body

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        vdoppler :: Doppler velocity between the Earth and the body / km s-1

    CALLING SEQUENCE:

        vdoppler = calc_vdoppler_earth_body(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    et0 = spice.str2et(date_obs+' '+time_obs)
    
    #position of a target body relative to observer (for correcting the travel time of light)
    positions, tau = spice.spkpos(body,et0,ref,abcorr,'EARTH')


    #Calculating the different positions and velocities
    state_body_earth,ltime_body_earth = spice.spkezr(body,et0,ref,abcorr,'EARTH')  #Position and speed of the observed body relative to the Earth

    c = 3.0e8 #m s-1
    r_body_earth=state_body_earth[0:3]
    v_body_earth=state_body_earth[3:6]*1.0e3               #[m/s]
    vtot_body_earth=np.sqrt(np.sum(v_body_earth**2.0))
    d_body_earth=ltime_body_earth*c/(149597870700.0)       #[AU]
    dir_body_earth=r_body_earth/np.sqrt(np.sum(r_body_earth**2.0))


    #Calculating the doppler velocity between the Earth and the observed body
    v_doppler_earth_body = spice.vdot(v_body_earth,dir_body_earth) / 1.0e3  #km/s
    
    return v_doppler_earth_body


###############################################################################################

def calc_solar_elongation_body(date_obs,time_obs,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_solar_elongation_body()

    DESCRIPTION : Calculate the angle between the Sun and the body as seen from Earth

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        sol_elong :: Angle between the Sun and the planet as seen from Earth (degrees)

    CALLING SEQUENCE:

        sol_elong = calc_solar_elongation_body(date_obs,time_obs)

    MODIFICATION HISTORY : Juan Alday (19/03/2025)

    """
    
    et0 = spice.str2et(date_obs+' '+time_obs)
    
    #Calculating the different positions and velocities
    state_body_earth,ltime_body_earth = spice.spkezr(body,et0,ref,abcorr,'EARTH')  #Position and speed of the observed body relative to the Earth

    #Calculating the different positions and velocities
    state_sun_earth,ltime_sun_earth = spice.spkezr('SUN',et0,ref,abcorr,'EARTH')  #Position and speed of the Sun relative to the Earth

    #Calculating the angle between the body and the Sun
    sol_elong = spice.vsep(state_body_earth[0:3],state_sun_earth[0:3]) / np.pi * 180. #degrees
    
    return sol_elong

###############################################################################################

def calc_radec(date_obs,time_obs,lat,lon,alt,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_geometry_observation_radecgrid()

    DESCRIPTION : Calculate the Right Ascension and Declination at the centre and North Pole of a planet

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC
        lat :: Latitude of the observer / degrees
        lon :: Longitude of the observer / degrees
        alt :: Altitude of the observer / km

    OPTIONAL INPUTS: 
    
        body :: Observed body
        ref :: Observed body's body-fixed frame
        abcorr :: Correction for aberrations
            
    OUTPUTS : 
 
        RA_centre :: Right ascension at centre of planet
        DEC_centre :: Declination at centre of planet
        RA_pole :: Right ascension at north pole
        DEC_pole :: Declination at north pole

    CALLING SEQUENCE:

        RA,DEC,RA_pole,DEC_pole = calc_radec(date_obs,time_obs,lat,lon,alt,body='MARS',ref='IAU_MARS',abcorr='LT+S')

    MODIFICATION HISTORY : Juan Alday (15/02/2025)

    """
    
    et = spice.str2et(date_obs+' '+time_obs)

    #Calculating Mars radius in each direction
    npos, radii = spice.bodvrd(body, "RADII", 3)
    rx_mars  =  radii[0]
    ry_mars  =  radii[1]
    rz_mars  =  radii[2]
    f_mars   =  ( rx_mars - rz_mars ) / rx_mars   #Flattening coefficient

    #Calculating Earth radius in each direction
    npos, radii = spice.bodvrd("EARTH", "RADII", 3)
    rx_earth  =  radii[0]
    ry_earth  =  radii[1]
    rz_earth  =  radii[2]
    f_earth   =  ( rx_earth - rz_earth ) / rx_earth   #Flattening coefficient

    celestial_frame = 'J2000'

    #Calculating the rotation matrix from IAU_MARS to J2000
    iau_mars_matrix = spice.pxfrm2(ref, celestial_frame, et, et)

    #Calculating the rotation matrix from J2000 to IAU_MARS
    j2000_mars_matrix = spice.pxfrm2(celestial_frame, ref, et, et)

    #Calculating the rotation matrix from IAU_EARTH to J2000
    iau_earth_matrix = spice.pxfrm2('IAU_EARTH', celestial_frame, et, et)

    #Calculating the rotation matrix from IAU_EARTH to IAU_MARS
    iau_earth2mars_matrix = spice.pxfrm2('IAU_EARTH', ref, et, et)
    
    #Calculating the rotation matrix from IAU_EARTH to IAU_MARS
    iau_mars2earth_matrix = spice.pxfrm2(ref, 'IAU_EARTH', et, et)
    
    #Calculate the offset between the centres of Mars-Earth and Mars-Sun in IAU_MARS coordinates
    earth_id = spice.bodn2c('EARTH')
    sun_id = spice.bodn2c('SUN')
    body_id = spice.bodn2c(body)
    offset_mars_earth_iaumars, lt = spice.spkez(earth_id,et,ref,abcorr,body_id)           #Position of Earth in IAU_MARS co-ordinates
    offset_earth_mars_j2000, lt = spice.spkez(body_id,et,celestial_frame,abcorr,earth_id) #Position of Mars in J2000 co-ordinates
    offset_earth_mars_iauearth, lt = spice.spkez(body_id,et,'IAU_EARTH',abcorr,earth_id)  #Position of Mars in IAU_EARTH co-ordinates
    sun_iau_mars, lt = spice.spkez(sun_id,et,ref,abcorr,body_id)                          #Position of Sun in IAU_MARS co-ordinates
    
    #Calculating the Cartesian coordinates of the observer on Earth  (IAU_EARTH and the IAU_MARS and J2000)
    observer_iau_earth = spice.georec(lon/180.*np.pi,lat/180.*np.pi, alt, rx_earth, f_earth)
    observer_iau_mars = spice.mxv(iau_earth2mars_matrix, observer_iau_earth) + offset_mars_earth_iaumars[:3]
    observer_j2000 = spice.mxv(iau_earth_matrix, observer_iau_earth)
    
    #Calculating the Cartesian coordinates of the body's North Pole (IAU_MARS and IAU_EARTH and J2000)
    northpole_iau_mars = spice.georec(0./180.*np.pi,90./180.*np.pi, rz_mars, rx_mars, f_mars)
    northpole_iau_earth = spice.mxv(iau_mars2earth_matrix, northpole_iau_mars) + offset_earth_mars_iauearth[:3]
    northpole_j2000 = spice.mxv(iau_mars_matrix, northpole_iau_mars) + offset_earth_mars_j2000[:3]

    #Calculate the right ascension and declination of the centre of the body as seen from the observer
    ######################################################################################################

    #Calculating the observer-location vector in J2000
    s_obsever_location_J2000 = offset_earth_mars_j2000[:3] - observer_j2000[:3]

    #Calculating the spherical co-ordinates of the observer-location vector in J2000 (right ascension and declination)
    radius, ra_body, dec_body = spice.reclat(s_obsever_location_J2000)

    # Convert right ascension and declination to degrees
    ra_mars = np.rad2deg(ra_body)
    dec_mars = np.rad2deg(dec_body)
    dist_mars = radius
    
    #Calculate the right ascension and declination of the North Pole of the body as seen from the observer
    ######################################################################################################
    
    #Calculating the observer-north pole vector in J2000
    s_obsever_northpole_J2000 = northpole_j2000[:3] - observer_j2000[:3]

    #Calculating the spherical co-ordinates of the observer-location vector in J2000 (right ascension and declination)
    radius, ra_body, dec_body = spice.reclat(s_obsever_northpole_J2000)

    # Convert right ascension and declination to degrees
    ra_northpole = np.rad2deg(ra_body)
    dec_northpole = np.rad2deg(dec_body)
    dist_northpole = radius
    
    return ra_mars,dec_mars,ra_northpole,dec_northpole


###############################################################################################

def calc_position_angle(date_obs,time_obs,lat,lon,alt,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_posangle()

    DESCRIPTION : Calculate the position angle of the planet. The position angle is defined as the angle between the 
                  N-S axis of the planet wrt the N-S axis of the celestial cardinal points. The position angle is defined
                positive towards the East (i.e., left, counter clockwise)

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC
        lat :: Latitude of the observer / degrees
        lon :: Longitude of the observer / degrees
        alt :: Altitude of the observer / km

    OPTIONAL INPUTS: 
    
        body :: Observed body
        ref :: Observed body's body-fixed frame
        abcorr :: Correction for aberrations
            
    OUTPUTS : 
 
        pa :: Position angle / degrees

    CALLING SEQUENCE:

        PA = calc_position_angle(date_obs,time_obs,lat,lon,alt,body='MARS',ref='IAU_MARS',abcorr='LT+S')

    MODIFICATION HISTORY : Juan Alday (15/02/2025)

    """
    
    #Calculate the right ascension and declination of the planet centre and north pole
    ra,dec,ra_pole,dec_pole = calc_radec(date_obs,time_obs,lat,lon,alt,body=body,ref=ref,abcorr=abcorr)
    
    ra, dec = np.radians(ra), np.radians(dec)
    ra_pole, dec_pole = np.radians(ra_pole), np.radians(dec_pole)
    
    # Compute the position angle
    numerator = np.cos(dec) * np.sin(ra - ra_pole)
    denominator = np.sin(dec) * np.cos(dec_pole) - np.cos(dec) * np.sin(dec_pole) * np.cos(ra - ra_pole)
    P = np.degrees(np.arctan2(numerator, denominator))  # Convert back to degrees
    
    PA = 180. + P
    
    return PA
    
    
###############################################################################################

def calc_emission_angle_earth(date_obs,time_obs,lat,lon,alt,body='MARS',ref='IAU_MARS',abcorr='LT+S'):
    
    """
    FUNCTION NAME : calc_emission_angle_earth()

    DESCRIPTION : Calculate the emission angle of the observer at the Earth's surface (or atmosphere)

    INPUTS : 

        date_obs :: Observation date / UTC
        time_obs :: Observation time / UTC
        lat :: Latitude of the observer / degrees
        lon :: Longitude of the observer / degrees
        alt :: Altitude of the observer / km

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        emission_angle :: Emission angle wrt to Earth's atmosphere (0 is looking up, 90 is exactly at the limb)

    CALLING SEQUENCE:

        emiss_angle_earth = calc_emission_angle_earth(date_obs,time_obs,lat,lon,alt)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)

    """
    
    et = spice.str2et(date_obs+' '+time_obs)

    #Calculating Earth radius in each direction
    npos, radii = spice.bodvrd("EARTH", "RADII", 3)
    rx_earth  =  radii[0]
    ry_earth  =  radii[1]
    rz_earth  =  radii[2]
    f_earth   =  ( rx_earth - rz_earth ) / rx_earth   #Flattening coefficient

    #Calculate the offset between the centres of Mars-Earth and Mars-Sun in IAU_MARS coordinates
    earth_id = spice.bodn2c('EARTH')
    body_id = spice.bodn2c(body)

    #Calculating the Cartesian coordinates of the observer on Earth  (IAU_EARTH)
    observer_iau_earth = spice.georec(lon/180.*np.pi,lat/180.*np.pi, alt, rx_earth, f_earth)
    
    body_iau_earth, lt = spice.spkez(body_id,et,'IAU_EARTH',abcorr,earth_id)  #Position of body in IAU_EARTH co-ordinates
    
    #Calculating the position of the observed body in IAU_EARTH coordinates
    observer2body_iau_earth = body_iau_earth[0:3] - observer_iau_earth
    
    #Calculating the angle 
    angle = spice.vsep(observer2body_iau_earth,observer_iau_earth+observer_iau_earth) / np.pi * 180.
    
    return angle
    


##############################################################################################
##############################################################################################
#                                   TIME SERIES ROUTINES
##############################################################################################
##############################################################################################

from datetime import date, timedelta, datetime, timezone

irtf_semesters = {
    "2023A": {
        "initial_date": datetime(2023, 2, 1),
        "end_date": datetime(2023, 7, 31)
    },
    "2023B": {
        "initial_date": datetime(2023, 8, 1),
        "end_date": datetime(2024, 1, 31)
    },
    "2024A": {
        "initial_date": datetime(2024, 2, 1),
        "end_date": datetime(2024, 7, 31)
    },
    "2024B": {
        "initial_date": datetime(2024, 8, 1),
        "end_date": datetime(2025, 1, 31)
    },
    "2025A": {
        "initial_date": datetime(2025, 2, 1),
        "end_date": datetime(2025, 7, 31)
    },
    "2025B": {
        "initial_date": datetime(2025, 8, 1),
        "end_date": datetime(2026, 1, 31)
    },
    "2026A": {
        "initial_date": datetime(2026, 2, 1),
        "end_date": datetime(2026, 7, 31)
    },
    "2026B": {
        "initial_date": datetime(2026, 8, 1),
        "end_date": datetime(2027, 1, 31)
    },
    "2027A": {
        "initial_date": datetime(2027, 2, 1),
        "end_date": datetime(2027, 7, 31)
    },
    "2027B": {
        "initial_date": datetime(2027, 8, 1),
        "end_date": datetime(2028, 1, 31)
    },
    "2028A": {
        "initial_date": datetime(2028, 2, 1),
        "end_date": datetime(2028, 7, 31)
    },
    "2028B": {
        "initial_date": datetime(2028, 8, 1),
        "end_date": datetime(2029, 1, 31)
    },
    "2029A": {
        "initial_date": datetime(2029, 2, 1),
        "end_date": datetime(2029, 7, 31)
    },
    "2029B": {
        "initial_date": datetime(2029, 8, 1),
        "end_date": datetime(2030, 1, 31)
    },
    "2030A": {
        "initial_date": datetime(2030, 2, 1),
        "end_date": datetime(2030, 7, 31)
    },
    "2030B": {
        "initial_date": datetime(2030, 8, 1),
        "end_date": datetime(2031, 1, 31)
    }
}

#######################################################################################################

def evolution_geometry_semester(semester,planet="Mars"):
    """
    FUNCTION NAME : evolution_geometry_semester()

    DESCRIPTION : Calculate the geometry of the ground-based observations for a given IRTF semester

    INPUTS : 

        semester :: Semester (e.g., "2026B")

    OPTIONAL INPUTS: 

        body :: Name of the body we want to observe
            
    OUTPUTS : 
 
        dates(ntimes) :: List of dates
        angsize(ntimes) :: Angular size of the disk (arcsec)
        vdoppler(ntimes) :: Doppler velocity between the Earth and the body (km/s)
        sol_elong(ntimes) :: Solar elongation (degrees)
        phase_ang(ntimes) :: Phase angle (degrees)

    CALLING SEQUENCE:

        dates,angsize,vdoppler,sol_elong,phase_ang = evolution_geometry_semester(semester,body="Mars")

    MODIFICATION HISTORY : Juan Alday (15/03/2025)
    """

    # Define date range for the semester
    start_date = irtf_semesters[semester]["initial_date"]
    end_date = irtf_semesters[semester]["end_date"]

    # Generate dates between start and end dates
    dates = [start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)]

    # store the dates between two dates in a list
    datesx = []

    for i in range(len(dates)):
        datesx.append(dates[i].isoformat()[0:10])

    #Perform the calculations
    angsize = []
    vdoppler = []
    sol_elong = []
    phase_ang = []
    for i in range(len(datesx)):
        
        angsize.append(calc_angsize_body(datesx[i],'00:00:00',body=planet.upper(),ref='IAU_'+planet.upper(),abcorr='LT+S'))
        vdoppler.append(calc_vdoppler_earth_body(datesx[i],'00:00:00',body=planet.upper(),ref='IAU_'+planet.upper(),abcorr='LT+S'))
        sol_elong.append(calc_solar_elongation_body(datesx[i],'00:00:00',body=planet.upper(),ref='IAU_'+planet.upper(),abcorr='LT+S'))
        phase_ang.append(calc_phase_angle(datesx[i],'00:00:00',body=planet.upper(),ref='IAU_'+planet.upper(),abcorr='LT+S'))

    return dates,angsize,vdoppler,sol_elong,phase_ang

#######################################################################################################


def plot_geometry_semester(semester,planet="Mars"):
    """
    FUNCTION NAME : plot_geometry_semester()

    DESCRIPTION : Calculate the geometry of the ground-based observations for a given IRTF semester
                    and make a plot

    INPUTS : 

        semester :: Semester (e.g., "2026B")

    OPTIONAL INPUTS: 

        body :: Name of the body we want to observe
            
    OUTPUTS : 
 
        Plot

    CALLING SEQUENCE:

        plot_geometry_semester(semester,body="Mars")

    MODIFICATION HISTORY : Juan Alday (15/03/2025)
    """

    dates,angsize,vdoppler,sol_elong,phase_ang = evolution_geometry_semester(semester,planet=planet)


    fig, ax1 = plt.subplots(1, 1, figsize=(8, 3))

    #Plotting the data
    ##########################################

    ax1.plot(dates, angsize)

    ax1.grid()

    ax2 = ax1.twinx()
    ax2.plot(dates,vdoppler,c='tab:red')

    ax3 = ax1.twinx()
    ax3.plot(dates,sol_elong,c='tab:green')
    ax3.spines["right"].set_position(("outward", 60))  # Shift outward

    ax4 = ax1.twinx()
    ax4.plot(dates,phase_ang,c='tab:purple')
    ax4.spines["right"].set_position(("outward", 120))  # Shift outward

    ax1.set_facecolor('lightgray')

    ax1.set_title(planet+" - Semester "+semester)
    ax1.set_xlabel('Date (YYYY-MM)')
    ax1.set_ylabel('Angular size (arcsec)',c='tab:blue')
    ax2.set_ylabel('Doppler velocity (km s$^{-1}$)',c='tab:red')
    ax3.set_ylabel('Solar elongation ($^\circ$)',c='tab:green')
    ax4.set_ylabel('Phase angle ($^\circ$)',c='tab:purple')

#######################################################################################################

def evolution_geometry_day(date,planet="Mars"):
    """
    FUNCTION NAME : evolution_geometry_day()

    DESCRIPTION : Calculate the geometry of the ground-based observations for a date from IRTF

    INPUTS : 

        date :: "YYYY-MM-DD"

    OPTIONAL INPUTS: 

        planet :: Name of the body we want to observe
            
    OUTPUTS : 
 
        times_hst(ntimes) :: List of HST times (Hawaii Standard Time)
        emiss_ang_earth(ntimes) :: Emission angle on Earth (degrees)
        

    CALLING SEQUENCE:

        times_hst,emiss_ang_earth = evolution_geometry_day(date,planet="Mars")

    MODIFICATION HISTORY : Juan Alday (15/03/2025)
    """

    from zoneinfo import ZoneInfo

    # Input Date
    hst_zone = ZoneInfo("Pacific/Honolulu")

    # 1. Parse the date string
    y, m, d = map(int, date.split("-"))

    # 2. Create NAIVE start and end (No timezone yet)
    # This ensures NumPy generates exactly 00:00 to 23:00
    start_naive = datetime(y, m, d, 0, 0)
    end_naive = datetime(y, m, d, 23, 0)

    # 3. Generate the range (still naive)
    hourly_naive = np.arange(
        start_naive, 
        end_naive + timedelta(hours=1), 
        timedelta(hours=1)
    ).astype(datetime)

    # 4. Localize to HST and then convert to UTC
    # .replace(tzinfo=...) says "This time IS HST"
    # .astimezone(timezone.utc) says "What is this HST time in UTC?"
    hourly_dates_hst = [dt.replace(tzinfo=hst_zone) for dt in hourly_naive]
    hourly_dates_utc = [dt.astimezone(timezone.utc) for dt in hourly_dates_hst]

    # store the dates between two dates in a list
    datesx = []
    timex = []
    for i in range(len(hourly_dates_utc)):
        datesx.append(hourly_dates_utc[i].isoformat()[0:10])
        timex.append(hourly_dates_utc[i].isoformat()[11:19])

    #Co-ordinates of Mauna Kea
    lat_observer = 19.8263
    lon_observer = -155.473
    alt_observer = 4205.

    #Perform calculations
    emiss_ang_earth = []
    for i in range(len(datesx)):    
        emiss_ang_earth.append(calc_emission_angle_earth(datesx[i],timex[i],lat_observer,lon_observer,alt_observer/1.0e3,body=planet.upper(),ref='IAU_'+planet.upper(),abcorr='LT+S'))

    return hourly_dates_hst, emiss_ang_earth

#######################################################################################################

def plot_geometry_day(date,planet="Mars"):
    """
    FUNCTION NAME : plot_geometry_day()

    DESCRIPTION : Calculate the geometry of the ground-based observations for a 
                    date from IRTF and plot it

    INPUTS : 

        date :: "YYYY-MM-DD"

    OPTIONAL INPUTS: 

        planet :: Name of the body we want to observe
            
    OUTPUTS : 
 
        Summary plot

    CALLING SEQUENCE:

        plot_geometry_day(date,planet="Mars")

    MODIFICATION HISTORY : Juan Alday (15/03/2025)
    """

    times_hst,emiss_ang_earth = evolution_geometry_day(date,planet=planet)

    fig, ax1 = plt.subplots(1, 1, figsize=(10, 4))

    ax1.plot(times_hst, emiss_ang_earth,c='tab:blue')
    ax1.set_ylim(0.,180.)

    ax1.fill_between(times_hst,0.,45.,color="tab:green",alpha=0.1)
    ax1.fill_between(times_hst,45.,90.,color="tab:orange",alpha=0.1)
    ax1.fill_between(times_hst,90.,180.,color="tab:red",alpha=0.1)

    ax1.text(times_hst[-2],5.,"Closer to zenith",verticalalignment="bottom",horizontalalignment="right",color="tab:green")
    ax1.text(times_hst[-2],50.,"Closer to horizon",verticalalignment="bottom",horizontalalignment="right",color="tab:orange")
    ax1.text(times_hst[-2],95.,"Not visible",verticalalignment="bottom",horizontalalignment="right",color="tab:red")

    ax1.grid()
    ax1.set_xlim(times_hst[0],times_hst[-1])
    ax1.set_title(planet+" - "+date)
    ax1.set_ylabel('Emission angle ($^\circ$)')
    ax1.set_xlabel('Date (HST)')
    ax1.set_facecolor('lightgray')
    plt.tight_layout()



##############################################################################################
##############################################################################################
#                                     MAPPING GEOMETRY
##############################################################################################
##############################################################################################


def xy2latlon(angmax,res,angsize_mars,lat_subobs,lon_subobs,rot_angle=0.):
    
    '''
    FUNCTION NAME : xy2latlon()

    DESCRIPTION :   Function to create an angular grid of pixels (in arcsec) and calculate the associated
                    latitude and longitude of the body disk observed in the centre of the grid
                    
                    Formulas were taken from https://en.wikipedia.org/wiki/Orthographic_map_projection

    INPUTS : 

        angmax :: Maximum size to include in the image (image goes from -angmax to angmax) (arcsec)
        res :: Resolution of the image (arcsec/pix)
        angsize :: Angular size of the observed body (arcsec)
        lat_subobs :: Sub-observer latitude (degrees)
        lon_subobs :: Sub-observer lontitude (degrees)
        
    OPTIONAL INPUTS:
    
        rot_angle :: Angle to rotate the disk to align it with the observation (0 angles in North and positive angle is to the East)
            
    OUTPUTS : 
 
        hx :: Angular position of each pixel in x direction (arcsec)
        hy :: Angular position of each pixel in y direction (arcsec)
        lat :: Latitude of each point (degrees)
        lon :: Longitude of each point (degrees)

    CALLING SEQUENCE:

        hx,hy,lat,lon = xy2latlon(angmax,res,angsize,lat_subobs,lon_subobs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''

    lat0 = np.radians(lat_subobs)  # Convert degrees to radians
    lon0 = np.radians(lon_subobs)

    # Image pixel coordinates
    px = np.arange(-angmax, angmax + res, res)
    py = np.arange(-angmax, angmax + res, res)
    hx, hy = np.meshgrid(px, py)

    # If rotation is needed, apply a proper 2D rotation matrix
    if rot_angle != 0.0:
        theta = np.radians(rot_angle)  # Convert degrees to radians
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        # Apply rotation matrix
        hxrot = hx * cos_theta - hy * sin_theta
        hyrot = hx * sin_theta + hy * cos_theta
    else:
        hxrot = hx
        hyrot = hy

    # Compute z of sphere hit position, if pixel's ray hits
    rho = np.sqrt(hxrot**2 + hyrot**2)
    c = np.ones(rho.shape)
    c[rho<=angsize_mars/2.] = np.arcsin(rho[rho<=angsize_mars/2.] / (angsize_mars / 2.0))

    # Compute latitude and longitude
    lat = np.arcsin(np.cos(c) * np.sin(lat0) + (hyrot * np.sin(c) * np.cos(lat0)) / rho)
    lon = lon0 + np.arctan2(hxrot * np.sin(c), (rho * np.cos(c) * np.cos(lat0) - hyrot * np.sin(c) * np.sin(lat0)))

    # Handle singularities (when rho == 0)
    ic = np.where(rho == 0.)
    if len(ic[0]) > 0:
        lat[ic] = lat0
        lon[ic] = lon0

    # Convert back to degrees
    lat = np.degrees(lat)
    lon = np.degrees(lon)

    # Adjust longitudes to stay within [-180,180]
    lon[lon <= -180] += 360
    lon[lon > 180] -= 360

    lat[hx**2. + hy**2. > (angsize_mars/2.)**2.] = np.nan
    lon[hx**2. + hy**2. > (angsize_mars/2.)**2.] = np.nan

    return hx, hy, lat, lon


###############################################################################################

def latlon2emiss(lat,lon,lat_subobs,lon_subobs):
    
    '''
    FUNCTION NAME : latlon2emiss()

    DESCRIPTION :   Given the latitude and longitude of a grid of pixels, and the sub-observer latitude and longitude
                    it calculates the emission angle at each pixel
                    
                    Formulas were taken from https://math.stackexchange.com/questions/2688803/angle-between-two-points-on-a-sphere

    INPUTS : 

        lat :: Latitude of each pixel in the image (deg)
        lon :: Longitude of each pixel in the image (deg)
        lat_subobs :: Sub-observer latitude (deg)
        lon_subobs :: Sub-observer longitude (deg)
    
    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        emiss_ang :: Emission angle of each point (degrees)

    CALLING SEQUENCE:

        emiss_ang = latlon2emiss(lat,lon,lat_subobs,lon_subobs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''
    
    #Changing from degrees to radians
    lat = lat / 180. * np.pi
    lon = lon / 180. * np.pi
    lat_subobs = lat_subobs / 180. * np.pi
    lon_subobs = lon_subobs / 180. * np.pi
    
    #Calculating the angle between each point and the sub-observer point
    #emiss_ang = np.arccos( np.sin(lon)*np.sin(lon_subobs) + np.cos(lon)*np.cos(lon_subobs)*np.cos(lat-lat_subobs) ) / np.pi * 180.
    emiss_ang = np.arccos( np.sin(lat)*np.sin(lat_subobs) + np.cos(lat)*np.cos(lat_subobs)*np.cos(lon-lon_subobs) ) / np.pi * 180.
    
    return emiss_ang

###############################################################################################

def latlon2sza(lat,lon,lat_subsol,lon_subsol):
    
    '''
    FUNCTION NAME : latlon2sza()

    DESCRIPTION :   Given the latitude and longitude of a grid of pixels, and the sub-solar latitude and longitude
                    it calculates the incident angle at each pixel
                    
                    Formulas were taken from https://math.stackexchange.com/questions/2688803/angle-between-two-points-on-a-sphere

    INPUTS : 

        lat :: Latitude of each pixel in the image (deg)
        lon :: Longitude of each pixel in the image (deg)
        lat_subsol :: Sub-solar latitude (deg)
        lon_subsol :: Sub-solar longitude (deg)
    
    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        sol_ang :: Incident solar angle of each point (degrees)

    CALLING SEQUENCE:

        sol_ang = latlon2sza(lat,lon,lat_subsol,lon_subsol)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''
    
    #Changing from degrees to radians
    lat = lat / 180. * np.pi
    lon = lon / 180. * np.pi
    lat_subsol = lat_subsol / 180. * np.pi
    lon_subsol = lon_subsol / 180. * np.pi
    
    #Calculating the angle between each point and the sub-observer point
    #emiss_ang = np.arccos( np.sin(lon)*np.sin(lon_subobs) + np.cos(lon)*np.cos(lon_subobs)*np.cos(lat-lat_subobs) ) / np.pi * 180.
    sza = np.arccos( np.sin(lat)*np.sin(lat_subsol) + np.cos(lat)*np.cos(lat_subsol)*np.cos(lon-lon_subsol) ) / np.pi * 180.
    
    return sza


###############################################################################################

def phase2azi(emiss_ang,sza,phase):
    
    '''
    FUNCTION NAME : phase2azi()

    DESCRIPTION :   Given the emission and incident angles of a grid of pixels, together with the 
                    phase angle it calculates the azimuth angle at each pixel. The azimuth angle here
                    follows the convention used in NEMESIS (phi=0 for forward scattering)

    INPUTS : 

        emiss_ang :: Emission angle (deg)
        sza :: Solar zenith angle (deg)
        phase :: Phase angle (deg)
    
    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        azi_ang :: Azimuth angle of each point (degrees)

    CALLING SEQUENCE:

        azi_ang = phase2azi(emiss_ang,sza,phase)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''
    
    #First of all let's calculate the scattering phase angle
    mu = np.cos(emiss_ang/180.*np.pi)   #Cosine of the reflection angle
    mu0 = np.cos(sza/180.*np.pi)    #Coside of the incidence angle
    cg = np.cos(phase/180.*np.pi)    #Cosine of the phase angle
    
    cazi = ((mu * mu0 - cg)/(np.sqrt(1. - mu**2.) * np.sqrt(1.-mu0**2.)))
    
    iin = np.where( (np.isnan(emiss_ang)==False) & (cazi<=-1.0) )
    cazi[iin] = -1.0
    
    iin = np.where( (np.isnan(emiss_ang)==False) & (cazi>=1.0) )
    cazi[iin] = 1.0
    
    iin = np.where( (np.isnan(cazi)==True) )
    cazi[iin] = 0.0
    
    azi_ang = np.arccos(cazi) / np.pi * 180.
    
    return azi_ang


###############################################################################################

def lon2lst(lon, lon_subsol):
    '''
    FUNCTION NAME : lon2lst()

    DESCRIPTION : Given the longitude in a grid of pixels and the sub-solar longitude, it 
                  calculates the Local Solar Time in each pixel.

    INPUTS : 
        lon :: Longitude of each pixel (deg)
        lon_subsol :: Sub-solar longitude (deg)
    
    OUTPUTS : 
        LST :: Local solar time (h)

    CALLING SEQUENCE:
        lst = lon2lst(lon, lon_subsol)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    '''
    
    # Calculating the Universal Local Time (LST at lon=0)
    ULT = 12. - lon_subsol / 180. * 12.
    
    # Compute Local Solar Time
    lst = ULT + lon / 180. * 12.
    
    # Ensure LST remains in the range [0, 24]
    lst = np.mod(lst, 24.)

    return lst




##############################################################################################
##############################################################################################
#                                     MAPPING GEOMETRY
##############################################################################################
##############################################################################################


# Convert RA from degrees to hours, minutes, seconds
def convert_ra_to_hms(ra_deg):
    '''
    Convert the right ascensions to hours, minutes, second
    '''
    ra_hours = ra_deg / 15
    h = int(ra_hours)
    m = int((ra_hours - h) * 60)
    s = ((ra_hours - h) * 60 - m) * 60
    return f"{h}h {m}m {s:.2f}s"

# Convert Dec to degrees, arcminutes, arcseconds
def convert_dec_to_dms(dec_deg):
    '''
    Convert declination for float to degrees,arcminutes,arcseconds
    '''
    d = int(dec_deg)
    m = int(abs(dec_deg - d) * 60)
    s = (abs(dec_deg - d) * 60 - m) * 60
    sign = "+" if dec_deg >= 0 else "-"
    return f"{sign}{abs(d)}° {m}′ {s:.2f}″"

def add_compass(ax,compass_x=0.15,compass_y=0.15,size=0.075):
    '''
    Add a compass with sky cardinal points to a given axis
    '''

    from matplotlib.patches import FancyArrow
    
    ax.annotate("N", xy=(compass_x, compass_y + size), xycoords='axes fraction',
                fontsize=10, ha="center", va="bottom")
    
    ax.annotate("S", xy=(compass_x, compass_y - size), xycoords='axes fraction',
                fontsize=10, ha="center", va="top")
    
    ax.annotate("W", xy=(compass_x + size, compass_y), xycoords='axes fraction',
                fontsize=10, ha="left", va="center")
    
    ax.annotate("E", xy=(compass_x - size, compass_y), xycoords='axes fraction',
                fontsize=10, ha="right", va="center")
    
    # Adding arrows for the compass
    ax.add_patch(FancyArrow(compass_x, compass_y, 0, size*0.7, transform=ax.transAxes, width=0.005))
    ax.add_patch(FancyArrow(compass_x, compass_y, 0, -size*0.7, transform=ax.transAxes, width=0.005))
    ax.add_patch(FancyArrow(compass_x, compass_y, size*0.7, 0, transform=ax.transAxes, width=0.005))
    ax.add_patch(FancyArrow(compass_x, compass_y, -size*0.7, 0, transform=ax.transAxes, width=0.005))

