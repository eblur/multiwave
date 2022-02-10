##
## fit_UV_data.py -- sandbox for evaluating the OI (UV) fit using the COS LSF.
## Apply the LSF to the model that is stored in a separate text file.
##
## 2022.02.04 -- liac@umich.edu
##-----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from cos_LSF import convolve_lsf
from astropy.table import Table
from astropy.io import fits

# if using the spectral information from SPEX saved txt files
from convolve_LSF import SpexDataTxt, SpexModelTxt

# TO DO:
# [ ] Load UV data directly (with specutils?)

# For now: load the SPEX version of the spectrum. pySPEX requires python 3.5 so
# the corresponding verison of specutils is very old

cos_data = SpexDataTxt('spex_output_files/OI_data.txt')

# [ ] Load SPEX model (from SPEX directly)

# Fow now: loading from the SPEX model output text file
cos_model = SpexModelTxt('spex_output_files/OI_spex_model.txt')
cos_model.notice = (cos_model.x > 1290.) & (cos_model.x < 1310.)

# [ ] Do the LSF convolution

cos_model.convolve(cenwave=1291,
	lsf_file='cos_lsf_files/aa_LSFTable_G130M_1291_LP1_cn.dat',
	disptab='cos_lsf_files/05i1639ml_disp.fits',
	detector="FUV")

# [ ] Evaluate the fit statistic between convolved model and data

def lnprob(data, model):
    """
    With Gaussian style data, the log-likelihood function is similar to
    chi-squared. In order to maximize the likelihood, we search for the maximum
    value of negative chi-squared. (See eq 11 of Fitting a Model to data:
    https://arxiv.org/pdf/1008.4686.pdf)

    We have asymmetric errorbars, which is interesting.

    Inputs
    ------
    data : SpexDataTxt object
    model : SpexModelTxt object after convolution has been applied
    """
    if model.y_convolved is None:
        print("Error: Convolution has not been applied to the model")
        return

    # Model is not necessarily on the same grid as the data. Interpolate.
    model_y = np.interp(data.x, model.x_convolved, model.y_convolved)

    # data.yerr has an upper and lower error bar
    # Use the relevant error bar to compute chi-squared
    iupper = (model_y > data.y)
    ilower = (model_y <= data.y)
    sigma  = np.zeros_like(data.x)
    # which element in the array is upper? which is lower?
    sigma[iupper] = data.yerr[1][iupper]
    sigma[ilower] = data.yerr[0][ilower]

    chisqr = ((data.y - model_y) / (2.0 *sigma))**2
    return -0.5 * np.sum(chisqr)

## Test run the statistic!
test = lnprob(cos_data, cos_model)
print(test)
