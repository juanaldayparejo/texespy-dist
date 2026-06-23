import os

#Path of this package
texespy_path = os.path.dirname(os.path.realpath(__file__))

#Path for general SPICE metakernel
spice_mk = texespy_path+'/data/spice/kernels_esa/mk/esa_generic_v01.tm'

#Path for archNEMESIS HITRAN24 file
archnemesis_hitran24 = '/srv/workspace/data/nemesis/spectroscopy/linedata/hitran24/hitran24.h5'
archnemesis_tips = '/srv/workspace/data/nemesis/spectroscopy/linedata/tips/tips2025.h5'
