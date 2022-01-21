##
## xray_data_test.py - A notebook for playing with X-ray spectra.
##
## 2022.01.21 - liac@umich.edu
##--------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

import pyxsis

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
