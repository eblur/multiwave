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


def read_lsf(filename):
    # This is the table of all the LSFs: called "lsf"
    # The first column is a list of the wavelengths corresponding to the line profile, so we set our header accordingly
    if "nuv_" in filename:  # If its an NUV file, header starts 1 line later
        ftype = "nuv"

    else:  # assume its an FUV file
        ftype = "fuv"
    hs = 0
    lsf = Table.read(filename, format="ascii", header_start=hs)

    # This is the range of each LSF in pixels (for FUV from -160 to +160, inclusive)
    # middle pixel of the lsf is considered zero ; center is relative zero
    pix = np.arange(len(lsf)) - len(lsf) // 2  # integer division to yield whole pixels

    # the column names returned as integers.
    lsf_wvlns = np.array([int(float(k)) for k in lsf.keys()])

    return lsf, pix, lsf_wvlns


def get_disp_params(disptab, cenwave, segment, x=[]):
    """
    Helper function to redefine_lsf(). Reads through a DISPTAB file and gives relevant\
    dispersion relationship/wavelength solution over input pixels.
    Parameters:
    disptab (str): Path to your DISPTAB file.
    cenwave (str): Cenwave for calculation of dispersion relationship.
    segment (str): FUVA or FUVB?
    x (list): Range in pixels over which to calculate wvln with dispersion relationship (optional).
    Returns:
    disp_coeff (list): Coefficients of the relevant polynomial dispersion relationship
    wavelength (list; if applicable): Wavelengths corresponding to input x pixels 
    """
    with fits.open(disptab) as d:
        wh_disp = np.where(
            (d[1].data["cenwave"] == cenwave)
            & (d[1].data["segment"] == segment)
            & (d[1].data["aperture"] == "PSA")
        )[0]
        disp_coeff = d[1].data[wh_disp]["COEFF"][0] # 0 is needed as this returns nested list [[arr]]
        d_tv03 = d[1].data[wh_disp]["D_TV03"]  # Offset from WCA to PSA in Thermal Vac. 2003 data
        d_orbit = d[1].data[wh_disp]["D"]  # Current offset from WCA to PSA

    delta_d = d_tv03 - d_orbit

    if len(x):  # If given a pixel range, build up a polynomial wvln solution pix -> λ
        wavelength = np.polyval(p=disp_coeff[::-1], x=np.arange(16384))
        return disp_coeff, wavelength
    else:  # If x is empty:
        return disp_coeff    


def redefine_lsf(lsf_file, cenwave, disptab, detector="FUV"):
    """
    Helper function to convolve_lsf(). Converts the LSF kernels in the LSF file from a fn(pixel) -> fn(λ)\
    which can then be used by convolve_lsf() and re-bins the kernels.
    Parameters:
    lsf_file (str): path to your LSF file
    cenwave (str): Cenwave for calculation of dispersion relationship
    disptab (str): path to your DISPTAB file
    detector (str): FUV or NUV?
    Returns:
    new_lsf (numpy.ndarray): Remapped LSF kernels.
    new_w (numpy.ndarray): New LSF kernel's LSF wavelengths.
    step (float): first order coefficient of the FUVA dispersion relationship; proxy for Δλ/Δpixel.
    """

    if detector == "FUV":
        xfull = np.arange(16384)

        # Read in the dispersion relationship here for the segments
        ### FUVA is simple
        disp_coeff_a, wavelength_a = get_disp_params(disptab, cenwave, "FUVA", x=xfull)
        ### FUVB isn't taken for cenwave 1105, nor 800:
        if (cenwave != 1105) & (cenwave != 800):
            disp_coeff_b, wavelength_b = get_disp_params(
                disptab, cenwave, "FUVB", x=xfull)
        elif cenwave == 1105:
            # 1105 doesn't have an FUVB so set it to something arbitrary and clearly not real:
            wavelength_b = [-99.0, 0.0]

        # Get the step size info from the FUVA 1st order dispersion coefficient
        step = disp_coeff_a[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        lsf_array = [np.array(lsf[key]) for key in lsf.keys()]
        if (deltaw < len(pix) * step * 2):  # resamples if the spacing of the original LSF wvlns is too narrow
            # this is all a set up of the bins we want to use
            # The wvln difference between kernels of the new LSF should be about twice their width
            new_deltaw = round(len(pix) * step * 2.0)  
            new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  # nw = number of LSF wavelengths
            new_w = min(w) + np.arange(new_nw) * new_deltaw  # new version of lsf_wvlns

            # populating the lsf with the proper bins
            new_lsf = np.zeros((len(pix), new_nw))  # empty 2-D array to populate
            for i, current_w in enumerate(new_w):
                dist = abs(current_w - w)  # Find closest original LSF wavelength to new LSF wavelength
                lsf_index = np.argmin(dist)
                orig_lsf_wvln_key = lsf.keys()[lsf_index]  # column name corresponding to closest orig LSF wvln
                new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  # assign new LSF wvln the kernel of the closest original lsf wvln
        else:
            new_lsf = lsf
            new_w = w
        return new_lsf, new_w, step

    elif detector == "NUV":
        xfull = np.arange(1024)
        # Read in the dispersion relationship here for the segments
        disp_coeff_a, wavelength_a = get_disp_params(disptab, cenwave, "NUVA", x=xfull)
        disp_coeff_b, wavelength_b = get_disp_params(disptab, cenwave, "NUVB", x=xfull)
        disp_coeff_c, wavelength_c = get_disp_params(disptab, cenwave, "NUVC", x=xfull)

        # Get the step size info from the NUVB 1st order dispersion coefficient
        step = disp_coeff_b[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        lsf_array = [np.array(lsf[key]) for key in lsf.keys()]

        # this section is a set up of the new bins we want to use:
        new_deltaw = round(len(pix) * step * 2.0)  # The wvln difference between kernels of the new LSF should be about twice their width
        new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  # nw = number of LSF wavelengths
        new_w = min(w) + np.arange(new_nw) * new_deltaw  # new version of lsf_wvlns

        # populating the lsf with the proper bins
        new_lsf = np.zeros((len(pix), new_nw))  # empty 2-D array to populate
        for i, current_w in enumerate(new_w):
            dist = abs(current_w - w)  # Find closest original LSF wavelength to new LSF wavelength
            lsf_index = np.argmin(dist)
            orig_lsf_wvln_key = lsf.keys()[lsf_index]  # column name corresponding to closest orig LSF wvln
            new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  # assign new LSF wvln the kernel of the closest original lsf wvln
        return new_lsf, new_w, step

    
def convolve_lsf(wavelength, spec, cenwave, lsf_file, disptab, detector="FUV"):
    """
    Main function; Convolves an input spectrum - i.e. template or STIS spectrum - with the COS LSF.
    Parameters:
    wavelength (list or array): Wavelengths of the spectrum to convolve.
    spec (list or array): Fluxes or intensities of the spectrum to convolve.
    cenwave (str): Cenwave for calculation of dispersion relationship
    lsf_file (str): Path to your LSF file
    disptab (str): Path to your DISPTAB file
    detector (str) : Assumes an FUV detector, but you may specify 'NUV'.
    Returns:
    wave_cos (numpy.ndarray): Wavelengths of convolved spectrum.!Different length from input wvln
    final_spec (numpy.ndarray): New LSF kernel's LSF wavelengths.!Different length from input spec
    """
    # First calls redefine to get right format of LSF kernels
    new_lsf, new_w, step = redefine_lsf(lsf_file, cenwave, disptab, detector=detector)

    # sets up new wavelength scale used in the convolution
    nstep = round((max(wavelength) - min(wavelength)) / step) - 1
    wave_cos = min(wavelength) + np.arange(nstep) * step

    # resampling onto the input spectrum's wavelength scale
    interp_func = interp1d(wavelength, spec)  # builds up interpolated function from input spectrum
    spec_cos = interp_func(wave_cos)  # builds interpolated initial spectrum at COS' wavelength scale for convolution
    final_spec = interp_func(wave_cos)  # Initializes final spectrum to the interpolated input spectrum

    for i, w in enumerate(new_w):  # Loop through the redefined LSF kernels
        # First need to find the boundaries of each kernel's "jurisdiction": where it applies
        # The first and last elements need to be treated separately
        if i == 0:  # First kernel
            diff_wave_left = 500
            diff_wave_right = (new_w[i + 1] - w) / 2.0
        elif i == len(new_w) - 1:  # Last kernel
            diff_wave_right = 500
            diff_wave_left = (w - new_w[i - 1]) / 2.0
        else:  # All other kernels
            diff_wave_left = (w - new_w[i - 1]) / 2.0
            diff_wave_right = (new_w[i + 1] - w) / 2.0

        # splitting up the spectrum into slices around the redefined LSF kernel wvlns
        # will apply the kernel corresponding to that chunk to that region of the spectrum - its "jurisdiction"
        chunk = np.where(
            (wave_cos < w + diff_wave_right) & (wave_cos >= w - diff_wave_left)
        )[0]
        if len(chunk) == 0:
            # off the edge, go to the next chunk
            continue

        current_lsf = new_lsf[:, i]  # selects the current kernel
        
		# Makes sure that the kernel is smaller than the chunk
        if len(chunk) >= len(current_lsf):  
            final_spec[chunk] = convolve(
                spec_cos[chunk],
                current_lsf,  # Applies the actual convolution
                boundary="extend",
                normalize_kernel=True,
            )
        else: 
        	print("The LSF is smaller than the data, cannot convolve")
        	return
		
    return wave_cos, final_spec  # Remember, not the same length as input spectrum data!
    
