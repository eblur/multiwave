##
## compare_data.py -- Compare different spectrum data files by opening the FITS
##   files and examining headers, columns, etc
##
## 2022.02.02 - liac@umich.edu
##-----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.units as u

## ----- HETG FILES

TESTFILEDIR = "tgcat/obs_8170_tgid_4247/"

pha_hdu = fits.open(TESTFILEDIR + "meg_-1.pha.gz")
pha_hdu.info()

pha_exposure = pha_hdu['PRIMARY'].header['EXPOSURE']
pha_exposure = pha_hdu['SPECTRUM'].header['EXPOSURE']

pha_data = pha_hdu['SPECTRUM'].data
pha_data.columns

hetg_arf = fits.open(TESTFILEDIR + "meg_-1.arf.gz")
hetg_arf.info()

arf_exposure = hetg_arf['SPECRESP'].header['EXPOSURE']

arf_data = hetg_arf['SPECRESP'].data
arf_data.columns

hetg_rmf = fits.open(TESTFILEDIR + "meg_-1.rmf.gz")
hetg_rmf.info()

rmf_matrix = hetg_rmf['MATRIX'].data
rmf_matrix.columns

rmf_ebounds = hetg_rmf['EBOUNDS'].data
rmf_ebounds.columns

## ------ RGS FILES

xmm_hdu = fits.open("xmm-data/XMM_spec_0148220201.FIT")
xmm_hdu.info()

rgs_data = xmm_hdu['SPECTRUM'].data
rgs_data.columns

rgs_exposure = xmm_hdu['SPECTRUM'].header['EXPOSURE']

rgs_rmf = fits.open('xmm-data/XMM_rmf_0148220201.FIT')
rgs_rmf.info()

rgs_matrix = rgs_rmf['MATRIX'].data
rgs_matrix.columns

rgs_ebounds = rgs_rmf['EBOUNDS'].data
rgs_ebounds.columns
