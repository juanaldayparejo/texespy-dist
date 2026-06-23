KPL/MK

Meta-kernel for Solar System's and Ground Stations' ESA generic kernels
==========================================================================

   This meta-kernel lists the ESA generic kernels for the main Solar System
   bodies and a selection of ground stations SPICE kernels.

   The kernels listed in this meta-kernel and the order in which
   they are listed are picked to provide the best data available and
   the most complete coverage for a generic Solar System scenario.


Usage of the Meta-kernel
-------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make use
   of this kernel must "load" the kernel normally during program
   initialization. Loading the kernel associates the data items with
   their names in a data structure called the "kernel pool". The SPICELIB
   routine FURNSH loads a kernel into the pool.

   The kernels listed below can be obtained from the ESA SPICE FTP server:

      ftp://spiftp.esac.esa.int/data/SPICE/esa_generic_kernels/kernels/


Implementation Notes
-------------------------------------------------------------------------

   It is recommended that users make a local copy of this file and
   modify the value of the PATH_VALUES keyword to point to the actual
   location of the BepiColombo SPICE data set's ``data'' directory on
   their system. Replacing ``/'' with ``\'' and converting line
   terminators to the format native to the user's system may also be
   required if this meta-kernel is to be used on a non-UNIX workstation.


-------------------

   This file was created on January 13, 2020 by Marc Costa Sitja ESA/ESAC.
   The original name of this file was esa_generic_kernels_v01.tm


   \begindata

     PATH_VALUES       = ( '/home/alday/Documents/Projects/texespy-dist/texespy/data/spice/kernels_esa' )

     PATH_SYMBOLS      = ( 'KERNELS' )

     KERNELS_TO_LOAD   = (

                           '$KERNELS/fk/planets/earthfixeditrf93.tf'
                           '$KERNELS/fk/planets/rssd0003.tf'

                           '$KERNELS/fk/stations/estrack_v04.tf'
                           '$KERNELS/fk/stations/earth_topo_050714.tf'
                           '$KERNELS/fk/stations/earthstns_ru_20191031.tf'

                           '$KERNELS/lsk/naif0012.tls'

                           '$KERNELS/pck/de403_masses.tpc'
                           '$KERNELS/pck/gm_de431.tpc'
                           '$KERNELS/pck/pck00010.tpc'

                           '$KERNELS/pck/earth_070425_370426_predict.bpc'
                           '$KERNELS/pck/earth_720101_070426.bpc'

                           '$KERNELS/spk/satellites/jup300_20200101_20400101.bsp'
                           '$KERNELS/spk/satellites/jup310_20200101_20400101.bsp'
                           '$KERNELS/spk/satellites/jup341_20200101_20400101.bsp'
                           '$KERNELS/spk/satellites/noe-5-2017-gal-a-reduced_20200101_20380902.bsp'
                           '$KERNELS/spk/satellites/mar097_20160314_20300101.bsp'

                           '$KERNELS/spk/stations/earthstns_itrf93_050714.bsp'
                           '$KERNELS/spk/stations/estrack_v04.bsp'
                           '$KERNELS/spk/stations/earthstns_ru_20191031.bsp'

                           '$KERNELS/spk/planets/de432s.bsp'
                           '$KERNELS/spk/planets/outerplanets_v0003.bsp'
                           '$KERNELS/spk/planets/outerplanets_v0004.bsp'

                         )

   \begintext


Contact Information
------------------------------------------------------------------------

   If you have any questions regarding this file contact the
   ESA SPICE Service (ESS) at ESAC:

           Marc Costa Sitja
           (+34) 91-8131-457
           esa_spice@sciops.esa.int, marc.costa@esa.int,


End of MK file.
