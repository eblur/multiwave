# multiwave

My personal Python library for analyzing multiwavelength datasets

## Installed packages for this library

* Python 3.9
* Astropy 5.0
* Specutils 1.5.0
* Numpy
* Scipy
* Matplotlib

### Special (in development) libraries used

* pyXsis (my development version)
* pyatomdb 0.10.8 (`pip install pyatomdb`)

### Special environments

I created a separate environment for using [3ML](https://threeml.readthedocs.io/en/stable/index.html) with Python 3.7, following the instructions from their documentation page.

```
conda create --name threeML -c conda-forge python=3.7 numpy scipy matplotlib
conda activate threeML
conda install -c conda-forge -c threeml astromodels threeml
```

I tried `conda install -c xspecmodels xspec-modelsonly` but there were so many package conflicts that it wouldn't install.
