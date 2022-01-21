##
## test_rmf.py -- Use the 3ML dispersion spectrum tools to examine an RMF file
##
## 2022.01.21 - liac@umich.edu
##-----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

import threeML as tml

from threeML.plugins.DispersionSpectrumLike import DispersionSpectrumLike
from threeML.utils.OGIP.response import OGIPResponse

response = OGIPResponse("../tgcat/obs_8170_tgid_4247/meg_-1.rmf.gz",
                        arf_file="../tgcat/obs_8170_tgid_4247/meg_-1.arf.gz")

source_function = tml.Broken_powerlaw(K=1e-2, alpha=0, beta=-2, xb=2000, piv=200)
background_function = tml.Powerlaw(K=10, index=-1.5, piv=100.0)

dispersion_spectrum_generator = DispersionSpectrumLike.from_function(
    "test",
    source_function=source_function,
    response=response,
    background_function=background_function,
)

tml.threeML_config.response_zero_color = "k"
tml.threemL_config.response_cmap = "magma"

_ = dispersion_spectrum_generator.display_rsp()
