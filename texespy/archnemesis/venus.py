import numpy as np
from struct import *
import sys,os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import texespy as texes
import archnemesis as ans

def perform_analysis(filename,
                        waven_min,waven_max,delv,resolving_power=90000.,v_doppler=0.,
                        emiss_ang=0.,sol_ang=0.,azi_ang=0.,
                        id_gases=None,iso_gases=None,vmr_gases=None,
                        hitran_file=texes.paths.archnemesis_hitran24,
                        iscat=0,
                        split_gas_contributions=False,
                        include_telluric=False,emiss_ang_earth=180.):
    """
    FUNCTION NAME : perform_analysis()

    DESCRIPTION : Run Venus forward models to investigate the different contributions to the spectrum

    INPUTS : 

        filename :: Name of the file
        waven_min :: Minimum wavenumber of the measurement in cm-1
        waven_max :: Maximum wavenumber of the measurement in cm-1
        delv :: Wavenumber step for the line-by-line opacity calculation in cm-1

    OPTIONAL INPUTS:
    
        resolving_power :: Spectral resolving power of the measurement (default: 90000)
        v_doppler :: Doppler velocity of the measurement in km/s (default: 0.)
        emiss_ang :: Emission angle of the measurement in degrees (default: 0.)
        sol_ang :: Solar zenith angle of the measurement in degrees (default: 0.)
        azi_ang :: Azimuth angle of the measurement in degrees (default: 0.)

        id_gases :: Array containing the gas IDs of the gases to be included in the spectroscopy (default: None, which means that only the main gases in the Venus atmosphere will be included)
        iso_gases :: Array containing the isotopologue IDs of the gases to be included in the spectroscopy (default: None, which means that only the main isotopologues of the main gases will be included)
        vmr_gases :: Array containing the volume mixing ratios of the gases to be included in the spectroscopy (default: None, which means that the volume mixing ratios of the gases in the reference atmosphere will be used)
        hitran_file :: Path to the HITRAN file to be used for the line-by-line opacity calculation (default: texes.paths.archnemesis_hitran24, which is a custom HITRAN24 file containing the main gases in the Venus atmosphere)
            
        iscat :: Scattering flag


    OUTPUTS : 
 
        vconv :: Convolution wavenumbers (cm-1)
        telluric_transmission(nconv) :: Telluric transmission
        venus_spectrum(nconv) :: Venus spectrum (erg s-1 cm-2 sr-1 (cm-1)-1)
        gaslabels(ngas) :: Names of the gases whose absorption is included
        venus_gas_spectrum(nconv,ngas) :: Venus spectrum with each gas individually (erg s-1 cm-2 sr-1 (cm-1)-1)

    CALLING SEQUENCE:

        vconv, telluric_transmission, venus_spectrum, gaslabels, venus_gas_spectrum = 
                                perform_analysis(filename,
                                waven_min,waven_max,delv,resolving_power=90000.,v_doppler=0.,
                                emiss_ang=0.,sol_ang=0.,azi_ang=0.,
                                id_gases=None,iso_gases=None,vmr_gases=None,
                                hitran_file=texes.paths.archnemesis_hitran24)

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    from copy import deepcopy

    
    #Creating archNEMESIS file
    print("Creating archNEMESIS file")
    create_archnemesis_file(filename,
                            waven_min,waven_max,delv,resolving_power=resolving_power,v_doppler=v_doppler,
                            emiss_ang=emiss_ang,sol_ang=emiss_ang,azi_ang=azi_ang,
                            id_gases=id_gases,iso_gases=iso_gases,vmr_gases=vmr_gases,
                            hitran_file=hitran_file,
                            iscat=iscat,
                            include_telluric=include_telluric,emiss_ang_earth=emiss_ang_earth)

    #Reading reference classes
    print("Reading archNEMESIS file")
    Atmosphere,Measurement,Spectroscopy,Scatter,Stellar,Surface,CIA,Layer,Variables,Retrieval,Telluric = ans.Files.read_input_files_hdf5(filename,calc_SE=False)

    #Calculating telluric contamination
    if include_telluric is True:

        print("Calculating telluric contamination")
        #Calculating high-spectral resolution telluric spectrum
        wavecalc_min_tel,wavecalc_max_tel = Measurement.calc_wave_range(apply_doppler=False,IGEOM=0)
        Telluric.Spectroscopy.read_tables(wavemin=wavecalc_min_tel,wavemax=wavecalc_max_tel)

        wave,telluric_transmissionx = Telluric.calc_transmission()

        #Convolving the telluric spectrum with the instrument lineshape
        telluric_transmission = Measurement.lblconv(wave,telluric_transmissionx,IGEOM=0)

        del wave,telluric_transmissionx
    else:
        telluric_transmission = np.zeros(Measurement.NCONV[0])


    #Calculating the Venus spectrum with all gases
    print("Calculating Venus spectrum")
    ForwardModel = ans.ForwardModel_0(Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables,Telluric=None)
    SPECONV = ForwardModel.nemesisfm()
    venus_spectrum = SPECONV[:,0] * 1.0e7 #(erg s-1 cm-2 sr-1 (cm-1)-1)
    del SPECONV

    #Calculating the contribution from different gases
    print("Calculating gas contributions")
    lbltables = Spectroscopy.LOCATION_LD
    pftables = Spectroscopy.LOCATION_PF
    ngas = len(lbltables)

    venus_gas = np.zeros((Measurement.NCONV[0],ngas))
    gasnames = [""] * ngas

    if split_gas_contributions is True:

        for i in range(ngas):

            gasnames[i] = ans.Data.gas_data.id_to_name(int(Spectroscopy.ID[i]),int(Spectroscopy.ISO[i]))

            #Defining only one gas in the spectroscopy
            Spectroscopy1 = deepcopy(Spectroscopy)
            Spectroscopy1.NGAS = 1
            Spectroscopy1.LOCATION_LD = [lbltables[i]]
            Spectroscopy1.LOCATION_PF = [pftables[i]]
            Spectroscopy1.LOCATION_CD = [lbltables[i]]
            Spectroscopy1.ID = [Spectroscopy.ID[i]]
            Spectroscopy1.ISO = [Spectroscopy.ISO[i]]
            Spectroscopy1.LINE_DATA = [Spectroscopy.LINE_DATA[i]]
            Spectroscopy1.LINE_DATA_PARAMS = [Spectroscopy.LINE_DATA_PARAMS[i]]

            #Running the forward model
            ForwardModel = ans.ForwardModel_0(Atmosphere=Atmosphere,Surface=Surface,Measurement=Measurement,Spectroscopy=Spectroscopy1,Stellar=Stellar,Scatter=Scatter,CIA=CIA,Layer=Layer,Variables=Variables)
            SPECONV1 = ForwardModel.nemesisfm()

            #Saving the results
            venus_gas[:,i] = SPECONV1[:,0] * 1.0e7 #(erg s-1 cm-2 sr-1 (cm-1)-1)

    return Measurement.VCONV[:,0],telluric_transmission, venus_spectrum, gasnames, venus_gas

###########################################################################################################################

def create_archnemesis_file(filename,
                            waven_min,waven_max,delv,resolving_power=90000.,v_doppler=0.,
                            emiss_ang=0.,sol_ang=0.,azi_ang=0.,
                            id_gases=None,iso_gases=None,vmr_gases=None,
                            iscat=0,
                            hitran_file=texes.paths.archnemesis_hitran24,
                            tips_file=texes.paths.archnemesis_tips,
                            include_telluric=False,emiss_ang_earth=180.):
    """
    FUNCTION NAME : create_archnemesis_file()

    DESCRIPTION : Create the archNEMESIS input HDF5 file for running a Venus forward model

    INPUTS : 

        filename :: Name of the file
        waven_min :: Minimum wavenumber of the measurement in cm-1
        waven_max :: Maximum wavenumber of the measurement in cm-1
        delv :: Wavenumber step for the line-by-line opacity calculation in cm-1

    OPTIONAL INPUTS:
    
        resolving_power :: Spectral resolving power of the measurement (default: 90000)
        v_doppler :: Doppler velocity of the measurement in km/s (default: 0.)
        emiss_ang :: Emission angle of the measurement in degrees (default: 0.)
        sol_ang :: Solar zenith angle of the measurement in degrees (default: 0.)
        azi_ang :: Azimuth angle of the measurement in degrees (default: 0.)

        id_gases :: Array containing the gas IDs of the gases to be included in the spectroscopy (default: None, which means that only the main gases in the Venus atmosphere will be included)
        iso_gases :: Array containing the isotopologue IDs of the gases to be included in the spectroscopy (default: None, which means that only the main isotopologues of the main gases will be included)
        vmr_gases :: Array containing the volume mixing ratios of the gases to be included in the spectroscopy (default: None, which means that the volume mixing ratios of the gases in the reference atmosphere will be used)
        hitran_file :: Path to the HITRAN file to be used for the line-by-line opacity calculation (default: texes.paths.archnemesis_hitran24, which is a custom HITRAN24 file containing the main gases in the Venus atmosphere)
        tips_file :: Path to the TIPS file to be used for the partition function calculation (default: texes.paths.archnemesis_tips, which is a custom TIPS file containing the partition function data for the main gases in the Venus atmosphere)
    OUTPUTS : 
 
        ArchNEMESIS HDF5 file

    CALLING SEQUENCE:

        create_archnemesis_file(filename,
                                waven_min,waven_max,delv,resolving_power=90000.,v_doppler=0.,
                                emiss_ang=0.,sol_ang=0.,azi_ang=0.,
                                id_gases=None,iso_gases=None,vmr_gases=None,
                                hitran_file=texes.paths.archnemesis_hitran24)

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    #Loading classes
    Atmosphere = load_reference_atmosphere(id_gases=id_gases,iso_gases=iso_gases,vmr_gases=vmr_gases)
    Measurement = create_measurement_class(waven_min,waven_max,emiss_ang,sol_ang,azi_ang,resolving_power=resolving_power,v_doppler=v_doppler)
    Scatter = create_scatter_class(waven_min,waven_max,iscat=iscat)
    Stellar, solfile = create_stellar_class()
    CIA = create_cia_class()
    Layer = create_layer_class(Atmosphere)
    Spectroscopy = create_spectroscopy_class(waven_min,waven_max,delv,id_gases=id_gases,iso_gases=iso_gases,hitran_file=hitran_file,tips_file=tips_file,resolving_power=resolving_power)
    Surface = create_surface_class()
    Retrieval = create_retrieval_class()
    if include_telluric is True:
        Telluric = create_telluric_class(waven_min,waven_max,delv,hitran_file=hitran_file,tips_file=tips_file,resolving_power=resolving_power,emiss_ang_earth=emiss_ang_earth)
        Telluric.write_hdf5(filename)

    #Writing archnemesis HDF5 file
    Surface.write_hdf5(filename)
    Atmosphere.write_hdf5(filename)
    Spectroscopy.write_hdf5(filename)
    Scatter.write_hdf5(filename)
    Layer.write_hdf5(filename)
    Stellar.write_hdf5(filename,solfile=solfile)
    Retrieval.write_input_hdf5(filename)
    Measurement.write_hdf5(filename)
    CIA.write_hdf5(filename)
    

    #Creating dummy .apr file
    fapr = open(filename+'.apr','w')
    fapr.write('#Forward model with texespy \n')
    fapr.write('\t %i \n' % (1))
    fapr.write('%i \t %i \t %i \n' % (2,0,2))
    fapr.write('\t %7.4f \t %7.4f \n' % (1.0,0.1))
    fapr.close()

###########################################################################################################################

def load_reference_atmosphere(id_gases=None,iso_gases=None,vmr_gases=None):
    """
    FUNCTION NAME : load_reference_atmosphere()

    DESCRIPTION : Load the reference atmosphere class for Venus

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Atmosphere :: An instance of the `Atmosphere_0` class containing the reference atmosphere for Venus

    CALLING SEQUENCE:

        Atmosphere = load_reference_atmosphere()

    MODIFICATION HISTORY : Juan Alday (213/03/2025)
    """

    Atmosphere = ans.Atmosphere_0()
    Atmosphere.read_hdf5(ans.archnemesis_path() + '/archnemesis/Data/reference_profiles/venus_vcd')

    #Updating gases if required
    if id_gases is not None:
        if len(id_gases) != len(vmr_gases):
            raise ValueError('The number of gas IDs and VMRs must be the same')
        for i in range(len(id_gases)):

            igas = np.where( (Atmosphere.ID==id_gases[i]) & (Atmosphere.ISO==iso_gases[i]) )[0]
            if len(igas)>0:
                Atmosphere.update_gas(id_gases[i],iso_gases[i],np.ones(Atmosphere.NP)*vmr_gases[i])
            else:
                Atmosphere.add_gas(id_gases[i],iso_gases[i],np.ones(Atmosphere.NP)*vmr_gases[i])

    return Atmosphere

###########################################################################################################################

def create_measurement_class(waven_min,waven_max,emiss_ang,sol_ang,azi_ang,resolving_power=90000.,v_doppler=0.):
    """
    FUNCTION NAME : create_measurement_class()

    DESCRIPTION : Create an instance of the `Measurement_0` class containing the measurement parameters for the Venus case

    INPUTS : 

        waven_min :: Minimum wavenumber of the measurement in cm-1
        waven_max :: Maximum wavenumber of the measurement in cm-1
        emiss_ang :: Emission angle of the measurement in degrees
        sol_ang :: Solar zenith angle of the measurement in degrees
        azi_ang :: Azimuth angle of the measurement in degrees
        
    OPTIONAL INPUTS:
    
        resolving_power :: Spectral resolving power of the measurement (default: 90000)
        v_doppler :: Doppler velocity of the measurement in km/s (default: 0.)
            
    OUTPUTS : 
 
        Measurement :: An instance of the `Measurement_0` class containing the measurement parameters for the Venus case

    CALLING SEQUENCE:

        Measurement = create_measurement_class(waven_min,waven_max,emiss_ang,sol_ang,azi_ang,resolving_power=90000.)

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Measurement = ans.Measurement_0()

    Measurement.ISPACE = 0 # wavenumber grid in cm-1
    Measurement.V_DOPPLER = v_doppler # Doppler velocity in km/s
    wave_mean = 0.5*(waven_min + waven_max)
    Measurement.ISHAPE = 2 # Gaussian instrumental line shape
    Measurement.FWHM = wave_mean / resolving_power # Full width at half maximum of the instrumental line shape in cm-1

    #Gometry parameters
    ngeom = 1
    Measurement.NGEOM = ngeom # number of geometries
    nav = np.ones(ngeom,dtype='int32') * 1
    flat = np.zeros((ngeom,nav.max()))
    flon = np.zeros((ngeom,nav.max()))
    wgeom = np.ones((ngeom,nav.max()))
    eang = np.zeros((ngeom,nav.max())) + emiss_ang
    aang = np.zeros((ngeom,nav.max())) + azi_ang
    sang = np.zeros((ngeom,nav.max())) + sol_ang

    Measurement.NGEOM = ngeom
    Measurement.NAV = nav
    Measurement.edit_FLAT(flat)
    Measurement.edit_FLON(flon)
    Measurement.edit_WGEOM(wgeom)
    Measurement.edit_EMISS_ANG(eang)
    Measurement.edit_AZI_ANG(aang)
    Measurement.edit_SOL_ANG(sang)

    #Spectral parameters
    delv = Measurement.FWHM / 5.
    vconv = np.arange(waven_min,waven_max+delv,delv)
    nconv = len(vconv)
    Measurement.NCONV = np.ones(ngeom,dtype='int32')*nconv

    VCONVx = np.zeros((Measurement.NCONV.max(),Measurement.NGEOM))
    MEAS = np.zeros((Measurement.NCONV.max(),Measurement.NGEOM))
    ERRMEAS = np.zeros((Measurement.NCONV.max(),Measurement.NGEOM))
    for i in range(Measurement.NGEOM):
        VCONVx[0:Measurement.NCONV[i],i] = vconv[0:Measurement.NCONV[i]]
        MEAS[0:Measurement.NCONV[i],i] = np.ones(Measurement.NCONV[i])
        ERRMEAS[0:Measurement.NCONV[i],i] = np.ones(Measurement.NCONV[i])*0.1
    Measurement.edit_VCONV(VCONVx)
    Measurement.edit_MEAS(MEAS)
    Measurement.edit_ERRMEAS(ERRMEAS)

    return Measurement

###########################################################################################################################

def create_scatter_class(waven_min,waven_max,iscat=0,nmu=5,nf=10,nphi=100):
    """
    FUNCTION NAME : create_cloud_model()

    DESCRIPTION : Create an instance of the `Scatter_0` class containing the cloud model for the Venus case

    INPUTS : 

        waven_min :: Minimum wavenumber of the measurement in cm-1
        waven_max :: Maximum wavenumber of the measurement in cm-1

    OPTIONAL INPUTS:
    
        iscat :: Flag indicating the type of scattering calculations
                    iscat = 0 : No scattering
                    iscat = 1 : Multiple scattering
            
    OUTPUTS : 
 
        Scatter :: An instance of the `Scatter_0` class containing the cloud model for the Venus case

    CALLING SEQUENCE:

        Scatter = create_cloud_model(Atmosphere)

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    #Defining cloud modes
    #Mode 1
    r_g1 = 0.3
    sigma_g1 = 0.445

    #Mode 2
    r_g2 = 1.0
    sigma_g2 = 0.255

    #Mode 2'
    r_g2p = 1.40
    sigma_g2p = 0.207

    #Mode 3
    r_g3 = 3.65
    sigma_g3 = 0.245

    #First of all we need to define our wavelength or wavenumber array
    #and tell the class in what units we want these calculations (wavenumber in cm-1 (ISPACE=0) or wavelength in um (ISPACE=1))
    Scatter = ans.Scatter_0()
    Scatter.read_refind(4)  #Reading optical properties of the main condensates in the Venus atmosphere (sulfuric acid aerosols)
    Scatter.ISCAT = iscat
    Scatter.NF = nf
    Scatter.NMU = nmu
    Scatter.NPHI = nphi
    Scatter.ISPACE = 0 #Wavelength

    Scatter.IMIE = 2 #Legendre polynomial expansion of the phase function
    waven = np.arange(int(waven_min)-1.,waven_max+2.,1.)
    NDUST = 4      #Number of aerosol populations that we want to include in our atmosphere
    NWAVE = len(waven)    #Number of spectral points
    NTHETA = 361    #Number of phase angles for defining the phase function
    theta = np.linspace(0.,180.,NTHETA)

    #Now we initialise the arrays that will be filled with the calculations
    Scatter.initialise_arrays(NDUST,NWAVE,NTHETA,NLPOL=150)
    Scatter.WAVE = waven
    Scatter.THETA = theta

    #Mode 1
    print('Calculating Mode 1')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g1,sigma_g1])

    idust = 0    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    #Mode 2
    print('Calculating Mode 2')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g2,sigma_g2])

    idust = 1    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    #Mode 2'
    print('Calculating Mode 2p')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g2p,sigma_g2p])

    idust = 2    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    #Mode 3
    print('Calculating Mode 3')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g3,sigma_g3])

    idust = 3    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    return Scatter

###########################################################################################################################

def create_stellar_class():
    """
    FUNCTION NAME : create_stellar_class()

    DESCRIPTION : Create an instance of the `Stellar_0` class containing the stellar parameters for the Venus case

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Stellar :: An instance of the `Stellar_0` class containing the stellar parameters for the Venus case

    CALLING SEQUENCE:

        Stellar = create_stellar_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Stellar = ans.Stellar_0()

    #Defining the Sun-Mars distance
    Stellar.DIST = 0.72

    #Defining the file containing the solar spectrum
    solfile = 'houghtonsolarwn.dat'

    return Stellar, solfile

###########################################################################################################################

def create_cia_class():
    """
    FUNCTION NAME : create_cia_class()

    DESCRIPTION : Create an instance of the `CIA_0` class containing the CIA parameters for the Venus case

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        CIA :: An instance of the `CIA_0` class containing the CIA parameters for the Venus case

    CALLING SEQUENCE:

        CIA = create_cia_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    CIA = ans.CIA_0(INORMAL=0,CIATABLE='CO2-CO2_HITRAN.h5')

    return CIA

###########################################################################################################################

def create_layer_class(Atmosphere):
    """
    FUNCTION NAME : create_layer_class()

    DESCRIPTION : Create an instance of the `Layer_0` class containing the layer parameters for the Venus case

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Layer :: An instance of the `Layer_0` class containing the layer parameters for the Venus case

    CALLING SEQUENCE:

        Layer = create_layer_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Layer = ans.Layer_0(Atmosphere.RADIUS)

    Layer.NLAY = 101        #Number of layers
    Layer.LAYTYP = 1       #Layering performed with equal changes in log pressure
    Layer.LAYINT = 1       
    Layer.LAYHT = 0. #m, just below the cloud
    Layer.assess()

    return Layer

###########################################################################################################################

def create_spectroscopy_class(waven_min,waven_max,delv,id_gases=None,iso_gases=None,hitran_file=texes.paths.archnemesis_hitran24,tips_file=texes.paths.archnemesis_tips,resolving_power=90000.):
    """
    FUNCTION NAME : create_spectroscopy_class()

    DESCRIPTION : Create an instance of the `Spectroscopy_0` class containing the spectroscopy parameters for the Venus case

    INPUTS : 

        waven_min :: Minimum wavenumber of the measurement in cm-1
        waven_max :: Maximum wavenumber of the measurement in cm-1
        delv :: Wavenumber step for the line-by-line opacity calculation in cm-1

    OPTIONAL INPUTS:
    
        id_gases :: Array containing the gas IDs of the gases to be included in the spectroscopy (default: None, which means that only the main gases in the Venus atmosphere will be included)
        iso_gases :: Array containing the isotopologue IDs of the gases to be included in the spectroscopy (default: None, which means that only the main isotopologues of the main gases
        hitran_file :: Path to the HITRAN file to be used for the line-by-line opacity calculation (default: texes.paths.archnemesis_hitran24, which is a custom HITRAN24 file containing the main gases in the Venus atmosphere)
        tips_file :: Path to the TIPS file to be used for the partition function calculation (default: texes.paths.archnemesis_tips, which is a custom TIPS file containing the partition function data for the main gases in the Venus atmosphere)
    
    OUTPUTS : 
 
        Spectroscopy :: An instance of the `Spectroscopy_0` class containing the spectroscopy parameters for the Venus case

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Spectroscopy = ans.Spectroscopy_0()

    Spectroscopy.ILBL = 1  #runtime calculation of the line-by-line opacity

    #ids_default = [1,2,9]
    #iso_default = [0,0,0]

    ids_default = [1,2]
    iso_default = [0,0]

    if id_gases is not None:

        ids = ids_default
        iso = iso_default
        for i in range(len(id_gases)):
            if id_gases[i] not in ids_default:
                ids.append(id_gases[i])
                iso.append(iso_gases[i])
    else:
        ids = ids_default
        iso = iso_default

    Spectroscopy.NGAS = len(ids)
    Spectroscopy.ID = ids
    Spectroscopy.ISO = iso
    Spectroscopy.ISPACE = 0
    Spectroscopy.IPROC = np.zeros(Spectroscopy.NGAS,dtype='int32')
    Spectroscopy.LOCATION_LD = [hitran_file] * Spectroscopy.NGAS
    Spectroscopy.LOCATION_PF = [tips_file] * Spectroscopy.NGAS
    Spectroscopy.LOCATION_CD = [hitran_file] * Spectroscopy.NGAS

    Spectroscopy.LINE_DATA_PARAMS = [ans.MolLineDataParams()] * Spectroscopy.NGAS

    #Calculating the spectral grid
    fwhm = np.mean([waven_min,waven_max]) / resolving_power
    waven_minx = waven_min - 5. * fwhm
    waven_maxx = waven_max + 5. * fwhm
    wavemin = np.floor(waven_minx/delv)*delv
    wavemax = np.ceil(waven_maxx/delv)*delv
    nwave = int(np.round((wavemax - wavemin) / delv))
    wave = np.linspace( wavemin , wavemax , nwave )
    Spectroscopy.NWAVE = nwave
    Spectroscopy.WAVE = wave

    return Spectroscopy

###########################################################################################################################

def create_surface_class(tsurf=737.):
    """
    FUNCTION NAME : create_surface_class()

    DESCRIPTION : Create an instance of the `Surface_0` class containing the surface parameters for the Venus case

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        tsurf :: Surface temperature in K (default: 737 K)
            
    OUTPUTS : 
 
        Layer :: An instance of the `Layer_0` class containing the layer parameters for the Venus case

    CALLING SEQUENCE:

        Layer = create_layer_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Surface = ans.Surface_0()
    Surface.LATITUDE = 0.
    Surface.LONGITUDE = 0.
    Surface.TSURF = tsurf
    Surface.LOWBC = 0  #Thermal emission only, no reflection
    Surface.ISPACE = 0 #Wavenumber in cm-1
    Surface.NEM = 2
    Surface.VEM = np.array([0.,100000.]) #Wavenumber grid for the surface emissivity
    Surface.EMISSIVITY = np.ones(Surface.NEM) * 1.
    Surface.assess()

    return Surface

###########################################################################################################################

def create_retrieval_class():
    """
    FUNCTION NAME : create_retrieval_class()

    DESCRIPTION : Create an instance of the `Retrieval_0` class containing the retrieval parameters for the Venus case

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Retrieval :: An instance of the `Retrieval_0` class containing the retrieval parameters for the Venus case

    CALLING SEQUENCE:

        Retrieval = create_retrieval_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Retrieval = ans.OptimalEstimation_0(IRET=0)

    Retrieval.NITER = -1       #Number of iterations
    Retrieval.PHILIMIT = 0.1   #Convergence criterion
    Retrieval.NCORES = 1      #Number of available cores

    Retrieval.assess_input()

    return Retrieval

###########################################################################################################################

def create_telluric_class(waven_min,waven_max,delv,hitran_file=texes.paths.archnemesis_hitran24,tips_file=texes.paths.archnemesis_tips,resolving_power=90000.,emiss_ang_earth=180.):
    """
    FUNCTION NAME : create_telluric_class()

    DESCRIPTION : Create an instance of the `Telluric_0` class containing the telluric parameters for the Venus case

    INPUTS : 

        waven_min :: Minimum wavenumber of the measurement in cm-1
        waven_max :: Maximum wavenumber of the measurement in cm-1
        delv :: Wavenumber step for the line-by-line opacity calculation in cm-1

    OPTIONAL INPUTS:
    
        emiss_ang_earth :: Emission angle of the Earth in degrees (default: 180.0)
        hitran_file :: Path to the HITRAN file to be used for the line-by-line opacity calculation (default: texes.paths.archnemesis_hitran24, which is a custom HITRAN24 file containing the main gases in the Venus atmosphere)
        tips_file :: Path to the TIPS file to be used for the partition function calculation (default: texes.paths.archnemesis_tips, which is a custom TIPS file containing the partition function data for the main gases in the Venus atmosphere)
        resolving_power :: Spectral resolving power of the measurement (default: 90000)

    OUTPUTS : 
 
        Telluric :: An instance of the `Telluric_0` class containing the telluric parameters for the Venus case

    CALLING SEQUENCE:

        Telluric = create_telluric_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Telluric = ans.Telluric_0()

    #Atmosphere
    ##################################################################################

    #Defining the inputs
    Telluric.DATE='01-01-2020'         #UTC date of the observation
    Telluric.TIME='00:00:00'           #UTC time of the observation
    Telluric.LATITUDE=19.82067         #Latitude of the observatory (Mauna Kea, Hawaii)
    Telluric.LONGITUDE=-155.46806      #Longitude of the observatory (Mauna Kea, Hawaii)
    Telluric.ALTITUDE=4207.3           #Altitude of the observatory (Mauna Kea, Hawaii)
    Telluric.EMISS_ANG=emiss_ang_earth #Observing angle, looking straight up to the zenith in this case

    #Extracting the atmosphere from the CIRC reference profiles
    Telluric.extract_atmosphere_circ()

    #Spectroscopy
    ##################################################################################

    Telluric.Spectroscopy = ans.Spectroscopy_0(ILBL=1)

    ids = [1,2,3,4,5,6,7]
    iso = [0,0,0,0,0,0,0]

    Telluric.Spectroscopy.NGAS = len(ids)
    Telluric.Spectroscopy.ID = ids
    Telluric.Spectroscopy.ISO = iso
    Telluric.Spectroscopy.ISPACE = 0
    Telluric.Spectroscopy.IPROC = np.zeros(Telluric.Spectroscopy.NGAS,dtype='int32')
    Telluric.Spectroscopy.LOCATION = [hitran_file] * Telluric.Spectroscopy.NGAS
    Telluric.Spectroscopy.LOCATION_PF = [tips_file] * Telluric.Spectroscopy.NGAS

    #Calculating the spectral grid
    fwhm = np.mean([waven_min,waven_max]) / resolving_power
    waven_minx = waven_min - 5. * fwhm
    waven_maxx = waven_max + 5. * fwhm
    wavemin = np.floor(waven_minx/delv)*delv
    wavemax = np.ceil(waven_maxx/delv)*delv
    nwave = int(np.round((wavemax - wavemin) / delv))
    wave = np.linspace( wavemin , wavemax , nwave )
    Telluric.Spectroscopy.NWAVE = nwave
    Telluric.Spectroscopy.WAVE = wave

    return Telluric

###########################################################################################################################

    