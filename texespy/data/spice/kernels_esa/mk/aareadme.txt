Solar System Generic MK files
===========================================================================

     This ``aareadme.txt'' file describes the contents of the kernels/mk
     directory of the Solar System Generic SPICE Kernel Dataset.

     It was last modified on August 23rd, 2022 by Alfredo Escalante Lopez, ESAC/ESA.


Brief Summary
--------------------------------------------------------

     This directory contains the Solar System Generic SPICE Meta-Kernel files.
     Meta-kernels (also known as ``FURNSH kernel'') are used
     to name a collection of kernels that are to be loaded into a user's
     application at run-time. Meta-kernels are appropriate to group kernels
     that correspond to a certain study group, phase or kernel state of the
     mission.


File naming conventions
--------------------------------------------------------

   Naming Scheme for Solar System Generic MKs:

     The naming scheme for the Solar System Generic MKs is:

           esa_generic_v[NN].tm

     where

           NNN       a count of the MK version generated for a given TAG
                     (mandatory if TAG is included; e.g 01)


Other directory contents
--------------------------------------------------------

     aareadme.txt         This file.

     former_versions      Directory where versions no longer valid are
                          stored for archive purposes.


Kernel File Details
--------------------------------------------------------

    Name                            Comments
    ---------------------------------------------------------------------

    esa_generic_vNN.tm              Contains the latest available Solar
                                    System generic kernels.


Contact Information
--------------------------------------------------------

   If you have any questions regarding this file contact the
   ESA SPICE Service (ESS) at ESAC:

           Alfredo Escalante Lopez
           (+34) 91-8131-429
           spice@sciops.esa.int


References and required readings
--------------------------------------------------------

     1. ``Kernel Required Reading'', NAIF Document


End of aareadme file.