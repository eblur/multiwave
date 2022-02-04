# Fitting COS data with pyspex --> https://spex-xray.github.io/spex-help/pyspex/commands.html
# Ioanna Psaradaki
# ipsarad@umich.edu

from pyspex import spex
import numpy as np
import matplotlib.pyplot as plt

s=spex.Session()

# Load the data
s.data('data_A_combined.res','data_A_combined.spo')

# Bin
s.bin(1,1,0,1600,2,'ang') #(inst, reg, elow, ehigh, factor, unit=None)

# Ignore the spectrum, except the oxygen region
s.ignore(1,1,0,1300,'ang')
s.ignore(1,1,1303,1800,'ang')

# Set the distance of the source
s.dist(1,2,'kpc')

# Set the abundance reference
a=s.abun('solar')

# Add fitting components
s.com('pow')
s.com('slab')
s.com_rel(1,1, np.array([2]))

# Set the component parameters
s.par(1,1,'norm', 0.5, thawn=True)
s.par(1,1,'gam', 1.2, thawn=True)

s.par(1,2,'o1', 21, thawn=True)
s.par(1,2,'v', 15, thawn=True)
s.par(1,2,'zv', -20, thawn=True)

# Calculate the model
s.calc()

# Plot the data and model
s.plot_data(xlog=False, ylog=False, title='SPEX', show=True)


# Change statistics
s.fit_stat('chi2')

# fit
s.fit_print(True)
s.fit()

# Print statistics
chisq=s.fit_chisq()
print(chisq)

# Plot the data and model
s.plot_data(xlog=False, ylog=False, title='SPEX', show=True)


# Show the free parameters fo the fit
print("--------------- FREE PARAMETERS AFTER FITTING --------------")
s.par_show('free')

# Save the parameters of the fit
s.par_write('cosOIparams.com', overwrite=True)


# Plot the model (after fitting)
s.plot_model(xlog=False, ylog=False, title='SPEX', show=True)
