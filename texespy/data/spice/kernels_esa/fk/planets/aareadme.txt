Solar System Generic Planets FK files
===========================================================================

     This ``aareadme.txt'' file describes the contents of the
     kernels/fk/planets directory of the Solar System Generic SPICE
     data server.

     It was last modified on August 23rd, 2022 by Alfredo Escalante Lopez, ESAC/ESA.


Contact Information
--------------------------------------------------------

     If you have any questions regarding this directory or its
     contents, please contact the ESA SPICE Service at ESAC:

             Alfredo Escalante Lopez
             (+34) 91-8131-429
             spice@sciops.esa.int
      
     
References and required readings
--------------------------------------------------------

     1. ``Frames Required Reading'', NAIF Document

     2. ``Kernel Pool Required Reading'', NAIF

     3. ``C-Kernel Required Reading'', NAIF

          
Brief Summary
--------------------------------------------------------

     This directory contains the SPICE Frames Definition Kernel files for the
     Science frames and Earth ITRF93 default Earth body-fixed frame association.
   

Current FK Kernels Set
--------------------------------------------------------

   earthfixediau.tf

      SPICE FK file that makes the IAU_EARTH frame coincide with the Earth
      fixed reference frame.


   earthfixeditrf93.tf

      SPICE FK file that makes the ITRF93 frame coincide with the Earth
      fixed reference frame.


   rssdNNNN.tf

      SPICE FK file defining a number of cross-mission frames that could be
      used by any of the users of any of the ESA planetary missions and that
      are not ``built'' in the SPICE toolkit.


Other directory contents
--------------------------------------------------------

     aareadme.txt         This file.
   

Kernel File Details
--------------------------------------------------------
 
     The most detailed description of the data in an FK file is provided in
     metadata included inside the descriptive text areas of the file. This
     information can be viewed using any text editor.


End of aareadme file.