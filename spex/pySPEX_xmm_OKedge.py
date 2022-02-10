#! /opt/anaconda3/envs/spex/bin/python
#
# Fitting xmm data with pyspex --> https://spex-xray.github.io/spex-help/pyspex/commands.html
# Ioanna Psaradaki
# ipsarad@umich.edu

from pyspex import spex
import numpy as np
import matplotlib.pyplot as plt

s=spex.Session()

# Load the data
s.data('spex_input_files/xmm.res','spex_input_files/xmm.spo')

# Bin
s.bin(1,1,19,35,2,'ang') #(inst, reg, elow, ehigh, factor, unit=None)

# Ignore the spectrum, except the oxygen K region
s.ignore(1,1,0,19,'ang')
s.ignore(1,1,35,1000,'ang')

# Set the distance of the source
s.dist(1,2,'kpc')

# Set the abundance reference
a=s.abun('solar')

# Add fitting components
s.com('pow')
s.com('hot')
s.com('hot')
s.com('hot')
s.com('slab')
s.com_rel(1,1, np.array([2,3,4,5]))

# Set the component parameters
s.par(1,1,'norm', 2.5, thawn=True)
s.par(1,1,'gam', 2.57, thawn=True)

#cold ism
s.par(1,2,'nh', 1.9e-3, thawn=True) # 10e28/m**2
s.par(1,2,'t', 5E-4, thawn=True) # keV
s.par(1,2,'26', 0.1, thawn=True)
s.par(1,2,'12', 0.1, thawn=True)
s.par(1,2,'14', 0.1, thawn=True)

#warm ism
s.par(1,3,'nh', 1.8e-5, thawn=True)
s.par(1,3,'t', 3.4e-3, thawn=True)

#hot ism
s.par(1,4,'nh', 6.5e-5, thawn=True)
s.par(1,4,'t', 9.7e-3, thawn=True)

s.par(1,5,'o6', 19.3, thawn=True) # log
s.par(1,5,'o7', 20, thawn=True)
s.par(1,5,'o8', 19.7, thawn=True)

# Calculate the model
s.calc()


# Plot the data and model
s.plot_data(xlog=False, ylog=False, title='SPEX', show=True)


# Plot data with model components
#s.plot_comp(xlog=False, ylog=False, title='SPEX', show=True)


# fit
s.fit_print(False)
s.fit()


# Plot the data and model
s.plot_data(xlog=False, ylog=False, title='SPEX', show=True)


# Show the free parameters fo the fit
print("--------------- FREE PARAMETERS AFTER FITTING --------------")
s.par_show('free')


# Save the parameters of the fit
s.par_write('spex_output_files/xmmparams.com', overwrite=True)


# Statistics
print("------------------- C stat --------------------------------- ")
print(s.fit_cstat()) # (Expected Cstat, Cstat)


# Plot the data and model after fitting
s.plot_data(xlog=False, ylog=False, title='SPEX', show=True)


# Plot the model (after fitting)
s.plot_model(xlog=False, ylog=False, title='SPEX', show=True)
