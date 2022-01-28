##
## xray_data_test.py - A notebook for playing with X-ray spectra.
##
## 2022.01.21 - liac@umich.edu
##--------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

import pyxsis

## --- Load a Chandra file

TESTFILEDIR = "tgcat/obs_8170_tgid_4247/"

meg_m1 = pyxsis.io.load_chandra_hetg(TESTFILEDIR + "meg_-1.pha.gz",
    rmf=TESTFILEDIR + "meg_-1.rmf.gz",
    arf=TESTFILEDIR + "meg_-1.arf.gz")

pyxsis.binspectrum.group_mincounts(meg_m1, 100)

ax = plt.subplot(111)
pyxsis.plot_unfold(ax, meg_m1, xunit='keV')
plt.semilogy()
plt.xlim(0.1, 2.0)
plt.show()

## --- Load an XMM file

# I have not yet built a loader for XMM data using pyXsis, so this
# will serve as a model for adding to that library
'''
from astropy.io import fits
import astropy.units as u

xmm_hdu = fits.open("xmm-data/XMM_spec_0148220201.FIT")
xmm_hdu.info()

xmm_data = xmm_hdu['SPECTRUM'].data
xmm_hdr = xmm_hdu['SPECTRUM'].header

bin_cen = xmm_hdr['TCRVL1'] # center of the reference bin
bin_cen_i = xmm_hdr['TLMIN1']-1 # reference bin (assuming 1 refers to first bin)
bin_width = xmm_hdr['TCDLT1'] # width of each channel
bin_unit = u.Unit(xmm_hdr['TCUNI1']) # unit for the reference
bin_lo = (np.arange(xmm_hdr['TLMAX1']) * bin_width + bin_cen - bin_width / 2.0) * bin_unit
bin_hi = (np.arange(xmm_hdr['TLMAX1']) * bin_width + bin_cen + bin_width / 2.0) * bin_unit

counts = xmm_data['COUNTS'] * u.Unit('count')
exposure = xmm_hdr['EXPOSURE'] * u.Unit('second')

my_spectrum = pyxsis.XBinSpectrum(bin_lo, bin_hi, counts,
    exposure=exposure, arf=None, rmf='xmm-data/XMM_rmf_0148220201.FIT')
'''

# Test my new pyxsis library

import numpy as np
import matplotlib.pyplot as plt

import pyxsis

rgs1 = pyxsis.io.load_xmm_rgs('xmm-data/XMM_spec_0148220201.FIT',
    rmf='xmm-data/XMM_rmf_0148220201.FIT')

#pyxsis.binspectrum.group_mincounts(rgs1, 30)

ax = plt.subplot(111)
pyxsis.plot_counts(ax, rgs1, xunit='keV')
plt.semilogy()
plt.xlim(0.1, 2.0)
plt.show()
