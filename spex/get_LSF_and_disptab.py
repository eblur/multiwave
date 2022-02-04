# This program downloads the LSF and DISPTAB files, that are appropriate for the COS observation I have.
# I need these files to convolve my moddel to the real LSF later on. 
# The LSF changes according to the year of observation.
# More information here ----> https://spacetelescope.github.io/COS-Notebooks/LSF.html#whichL
# What you need to run it is the data that you have downloaded from MAST and the Table file from here --> https://www.stsci.edu/hst/instrumentation/cos/performance/spectral-resolution

# Import for: array manipulation
import numpy as np

# Import for: reading fits files
from astropy.table import Table
from astropy.io import fits

# Import for: line profile functions
from astropy.modeling import functional_models

# Import for: convolutions
from astropy.convolution import convolve

# Import for: downloading the data
#from astroquery.mast import Observations

#Import for: interpolating a discretized function
from scipy.interpolate import interp1d

#Import for: plotting
from matplotlib import pyplot as plt

#Import for: downloading COS' LSF files within python
import urllib 

#Import for: unzipping a compressed file of STIS data
import tarfile

#Import for: working with system paths
from pathlib import Path

# --------------------------------------------------------------------------------------------------------------- #
# --------------------------- THIS FUNCTION WILL DOWNLOAD THE LSF AND DISPTAB FILES ----------------------------- #
# --------------------------------------------------------------------------------------------------------------- #

def fetch_files(det, grating, lpPos, cenwave, disptab):
    """
    Given all the inputs: (detector, grating, LP-POS, cenwave, dispersion table,) this will download both
    the LSF file and Disptab file you should use in the convolution and return their paths.
    Returns:
    LSF_file_name (str): filename of the new downloaded LSF file
    disptab_path (str): path to the new downloaded disptab file
    """
    COS_site_rootname = (
        "https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/cos/"
        "performance/spectral-resolution/_documents/"
    )  # Link to where all the files live - split into 2 lines
    if det == "NUV":  # Only one file for NUV
        LSF_file_name = "nuv_model_lsf.dat"
    elif det == "FUV":  # FUV files follow a naming pattern
        LSF_file_name = f"aa_LSFTable_{grating}_{cenwave}_LP{lpPos}_cn.dat"

    LSF_file_webpath = COS_site_rootname + LSF_file_name  # Where to find file online
    urllib.request.urlretrieve(
        LSF_file_webpath, str(datadir / LSF_file_name)
    )  # Where to save file to locally
    print(f"Downloaded LSF file to {str(datadir/ LSF_file_name)}")

    # And we'll need to get the DISPTAB file as well
    disptab_path = str(datadir / disptab)
    urllib.request.urlretrieve(
        f"https://hst-crds.stsci.edu/unchecked_get/references/hst/{disptab}",
        disptab_path,
    )
    print(f"Downloaded DISPTAB file to {disptab_path}")

    return LSF_file_name, disptab_path 
    
# ------------------------------------------------------------------------------------------------------------------- #
# OPEN THE HEADER OF COS DATAFILE TO SEE THE DISPTAB NAME AND THE DETAILS OF THE OBSERVATION (CENTRAL WAVELENGHT ETC) #
# ------------------------------------------------------------------------------------------------------------------- #

# The file with data that I downloaded from MAST
fuvFile = "data/lb2m02010_x1dsum.fits"

fuvHeader0 = fits.getheader(fuvFile, ext=0)  # Select the primary header
print(f"For the file {fuvFile}, the relevant parameters are: ")
param_dict = {}  # Make a dict to store what you find here

for hdrKeyword in [
    "DETECTOR",
    "OPT_ELEM",
    "LIFE_ADJ",
    "CENWAVE",
    "DISPTAB",
]:  # Print out the relevant values
    try:  # For DISPTAB
        value = fuvHeader0[hdrKeyword].split("$")[1]  # Save the key/value pairs to the dictionary
        param_dict[hdrKeyword] = value                # DISPTAB needs the split here
    except:  # For other params
        value = fuvHeader0[hdrKeyword]  # Save the key/value pairs to the dictionary
        param_dict[hdrKeyword] = value
    print(f"{hdrKeyword} = {value}")  # Print the key/value pairs

# ------------------------------------------------------------------------------------------------------------------- #
# --------------------------- DOWNLOAD THE LSF AND DISPTAB ---------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #

#datadir.mkdir(exist_ok=True), outputdir.mkdir(exist_ok=True), plotsdir.mkdir(exist_ok=True)
cwd = Path(".")
datadir = cwd / "data"
outputdir = cwd / "output"
plotsdir = cwd / "output" / "plots"

# run the function and download the LSF and disptab, appropriate for my dataset
LSF_file_name, disptab_path = fetch_files(*list(param_dict.values()))


