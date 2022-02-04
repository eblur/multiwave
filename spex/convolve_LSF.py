# This program convolves the model I obtain from spex, when fitting COS data, with the appropriate COS LSF.
# Ioanna Psaradaki
# ipsarad@umich.edu

import numpy as np
import matplotlib.pyplot as plt
from cos_LSF import convolve_lsf
from astropy.table import Table
from astropy.io import fits

# ------------------------------------------------------------------------------------------------ #
# -------------------------------- Load the UV data --------------------------------------------- #
# ------------------------------------------------------------------------------------------------ #
# in spex format

filename0=('cos_data_files/OI_data.txt')
x0,y0=np.loadtxt(filename0, usecols=(0,3), unpack=True)
xei0,xed0=np.loadtxt(filename0, usecols=(1,2), unpack=True)
yei0,yed0=np.loadtxt(filename0, usecols=(4,5), unpack=True)
xfit0,yfit0=np.loadtxt(filename0, usecols=(0,6), unpack=True)

xei0bis=np.zeros(len(x0))
xed0bis=np.zeros(len(x0))

for i in range(len(x0)):

	xed0bis=(-1)*xed0[i]

yei0bis=np.zeros(len(y0))
yed0bis=np.zeros(len(y0))

for i in range(len(y0)):

	yed0bis=(-1)*yed0[i]

xet0=[xei0,(-1)*xed0]
yet0=[yei0,(-1)*yed0]

# ------------------------------------------------------------------------------------------------- #
# ------------------ Load the .txt file that contains the model from spex (pl ty model) ----------- #
# ------------------------------------------------------------------------------------------------- #

x,y=np.loadtxt("cos_data_files/OI_spex_model.txt", usecols=(0,3), unpack=True)

# Fix the range that you want to examine
xnew=x[np.where((x>1290)&(x<1310))]
ynew=y[np.where((x>1290)&(x<1310))]


# ------------------------------------------------------------------------------------------------- #
# ------------------------------------- Convolution ----------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

# Load convolve_lsf from cos_LSF.py to do the convolution
x_convolved, y_convolved = convolve_lsf(wavelength=xnew, spec=ynew, cenwave=1291,
	lsf_file='cos_data_files/aa_LSFTable_G130M_1291_LP1_cn.dat',
	disptab='cos_data_files/05i1639ml_disp.fits',
	detector="FUV")

# ------------------------------------------------------------------------------------------------- #
# --------------------------------------- Plotting ------------------------------------------------ #
# ------------------------------------------------------------------------------------------------- #

plt.plot(x_convolved, y_convolved, color='b', label='real LSF model',linewidth=1.5)
plt.plot(x, y, color='r', label='SPEX model', linewidth=1.5)
plt.errorbar(x0,y0,yet0,xet0,ls='none',capsize=0,color='k', linewidth=1.5)
plt.legend(bbox_to_anchor=(0.88,0.82), loc=4, borderaxespad=0.,fontsize=11,framealpha=1)
plt.xlim(1299.5,1303)
plt.ylim(0,3)
plt.xlabel('Ang')
plt.ylabel('Flux')
plt.title("O I")
plt.show()
