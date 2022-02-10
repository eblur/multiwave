# This program convolves the model I obtain from spex, when fitting COS data, with the appropriate COS LSF.
# Ioanna Psaradaki
# ipsarad@umich.edu

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from cos_LSF import convolve_lsf

# ------------------------------------------------------------------------------------------------ #
# -------------------------------- Load the UV data --------------------------------------------- #
# ------------------------------------------------------------------------------------------------ #
# in spex format

filename0='spex_output_files/OI_data.txt'

class SpexDataTxt(object):
	def __init__(self, filename):
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
		self.x = x0
		self.y = y0
		self.xerr = xet0
		self.yerr = yet0

cos_data = SpexDataTxt(filename0)

# ------------------------------------------------------------------------------------------------- #
# ------------------ Load the .txt file that contains the model from spex (pl ty model) ----------- #
# ------------------------------------------------------------------------------------------------- #

modelfile = "spex_output_files/OI_spex_model.txt"

class SpexModelTxt(object):
	def __init__(self, filename):
		x,y=np.loadtxt(filename, usecols=(0,3), unpack=True)
		self.x = x
		self.y = y
		self.notice = np.ones_like(x, dtype='bool') # an array of True values
		self.x_convolved = None
		self.y_convolved = None
	@property
	def xplot(self):
		return self.x[self.notice]
	@property
	def yplot(self):
		return self.y[self.notice]
	def convolve(self, cenwave, lsf_file, disptab, detector):
		x_convolved, y_convolved = convolve_lsf(wavelength=self.xplot, spec=self.yplot,
			cenwave=cenwave, lsf_file=lsf_file, disptab=disptab, detector=detector)
		self.x_convolved = x_convolved
		self.y_convolved = y_convolved


# Fix the range that you want to examine
model = SpexModelTxt(modelfile)
model.notice = (model.x > 1290.) & (model.x < 1310.)

# ------------------------------------------------------------------------------------------------- #
# ------------------------------------- Convolution ----------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

# Load convolve_lsf from cos_LSF.py to do the convolution
model.convolve(cenwave=1291,
	lsf_file='cos_lsf_files/aa_LSFTable_G130M_1291_LP1_cn.dat',
	disptab='cos_lsf_files/05i1639ml_disp.fits',
	detector="FUV")

# ------------------------------------------------------------------------------------------------- #
# --------------------------------------- Plotting ------------------------------------------------ #
# ------------------------------------------------------------------------------------------------- #
"""
# Make a plot with residuals
fig = plt.figure(figsize=(6,6))
gs  = GridSpec(2, 1, height_ratios=[2,1], hspace=0.0)

ax0 = plt.subplot(gs[0])
ax0.plot(model.x_convolved, model.y_convolved, color='b', label='real LSF model',linewidth=1.5)
ax0.plot(model.xplot, model.yplot, color='r', label='SPEX model', linewidth=1.5)
ax0.errorbar(cos_data.x,cos_data.y,cos_data.yerr,cos_data.xerr,
	ls='none',capsize=0,color='k', linewidth=1.5)
ax0.legend(bbox_to_anchor=(0.88,0.82), loc=4, borderaxespad=0.,fontsize=11,framealpha=1)
ax0.set_ylabel('Flux')
ax0.set_ylim(0,3)
ax0.set_title("O I")

ax1 = plt.subplot(gs[1])
res1 = (cos_data.y - np.interp(cos_data.x, model.x, model.y))
res2 = (cos_data.y - np.interp(cos_data.x, model.x_convolved, model.y_convolved))
ax1.errorbar(cos_data.x, res1, yerr=cos_data.yerr, color='r')
ax1.errorbar(cos_data.x, res2, yerr=cos_data.yerr, color='b')
ax1.axhline(0.0, color='k', ls=':', lw=1)

ax0.set_xlim(1299.5,1303)
ax0.xaxis.set_ticklabels('')
ax1.set_xlim(1299.5,1303)
ax1.set_xlabel('Ang')
ax1.set_ylabel('Residual flux')
plt.show()
"""
