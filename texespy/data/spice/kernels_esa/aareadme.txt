+-----------------------------------------------------------------------------+
|  __   __   __      __   __     __   ___     __   ___  __          __   ___  |
| /__\ /__` '__\    /__` |__) | /  ` |__     /__` |__  |__) \  / | /  ` |__   |
| \__, .__/ \__/    .__/ |    | \__, |___    .__/ |___ |  \  \/  | \__, |___  |
|                                                                             |
+-----------------------------------------------------------------------------+

Welcome to ESA's Solar System Generic SPICE Kernel Dataset
==========================================================

     This directory contains the SPICE kernels Dataset generated for ESA's
     generic Solar System scenario, some of the kernels might be generated
     by NAIF and others by the ESA SPICE Service. Those kernels are not
     tied to a specific mission.

     For more information on SPICE for ESA Planetary Missions, please go to

       http://spice.esac.esa.int

     For information on SPICE, please go to

       https://naif.jpl.nasa.gov

     For any question or suggestion on this site, please contact:

       Marc Costa Sitja,       marc.costa@esa.int
       ESA SPICE Service,      esa_spice@sciops.esa.int


Content of this directory
--------------------------------------------------------

     Each subdirectory within this directory contains all the generic
     kernels of the same type.

        dsk:      Kernels that contain files contain detailed shape models
                  for solar system bodies.

        fk:       Kernels that define reference frames needed for the
                  scenario.

        lsk:      Leapseconds kernel.

        mk:       Meta-kernel files (a.k.a "furnsh" files) that provide
                  lists of kernels suitable for the generic scenario.

        pck:      Kernels that define planetary constants.

        spk:      Orbit kernels, for the ground stations and other solar
                  system bodies.


     Please refer to the aareadme.txt in each of the subdirectories for
     information on the contents for the specific type of kernels and the
     corresponding naming conventions.


Usage of the SPICE Kernels for ESA's generic Solar System scenario
------------------------------------------------------------------

     At least a basic knowledge of the SPICE system is needed in order to
     use these kernels. The SPICE toolkit provides versions in Fortran
     (SPICELIB), C (CSPICE), IDL (icy), Matlab (mice) and Python (SpiceyPy),
     and the user can choose any one that suits him/her.

     The SPICELIB routine FURNSH, CSPICE function furnsh_c, icy routine
     cspice_furnsh and SpiceyPy routine furnsh load a kernel file into the
     kernel pool as shown below.

        CALL FURNSH  ( 'kernel_name' )
        furnsh_c     ( "kernel_name" );
        cspice_furnsh, 'kernel_name'
        cspice_funsh( 'kernel_name' );
        spiceypy.furnsh( 'kernel_name' )

     In the case when two or more files contain data overlapping in time
     for a given object, for binary kernels, the file loaded last takes
     precedence.

     If two (or more) text kernels assign value(s) using the '=' operator
     to identical keywords, the data value(s) associated with the last
     loaded occurrence of the keyword are used -all earlier values have
     been replaced with the last loaded value(s).


-------------------

     This file was last modified on January 14th, 2020 by Marc Costa Sitja.


End of aareadme file.