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

## Special environments

## SPEX

I [installed SPEX for Mac OS](https://spex-xray.github.io/spex-help/getstarted/install.html) using the binary (Administrator) code. I then added the line 
`source /opt/spex/spexdist.sh` 
to my .bash_profile so I can access spex at any time from the command line.

Then I followed the instructions in [Section 1.4.1 on this page](https://spex-xray.github.io/spex-help/getstarted/pyspex.html) to install pySPEX:

```
cp /opt/spex/python/spex.yml .
conda env create -f ~/spex.yml
```

Then I installed my other python packages (above) in the conda spex environment.

## 3ML

I created a separate environment for using [3ML](https://threeml.readthedocs.io/en/stable/index.html) with Python 3.7, following the instructions from their documentation page.

```
conda create --name threeML -c conda-forge python=3.7 numpy scipy matplotlib
conda activate threeML
conda install -c conda-forge -c threeml astromodels threeml
```

I tried `conda install -c xspecmodels xspec-modelsonly` but there were so many package conflicts that it wouldn't install.

*Set up a jupyter notebook kernel*

```
conda activate threeML
python -m ipykernel install --user --name threeML --display-name "threeML"
conda deactivate
```
