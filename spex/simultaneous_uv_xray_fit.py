#! /opt/anaconda3/envs/spex/bin/pythons
##
## simultaneous_uv_xray_fit.py -- sandbox for handling the UV and X-ray data
##    simultaneously
##
## 2022.02.04 -- liac@umich.edu
##-----------------------------------------------------------------------------

from pyspex import spex
import numpy as np
import matplotlib.pyplot as plt

s=spex.Session()

## ---- 1: Set up the X-ray spectrum and model

## ---- 1.1: Load the data
s.data('xmm.res','xmm.spo')

s.bin(1,1,19,35,2,'ang') #(inst, reg, elow, ehigh, factor, unit=None)

# Ignore the spectrum, except the oxygen K region
s.ignore(1,1,0,19,'ang')
s.ignore(1,1,35,1000,'ang')

# Set the distance of the source
s.dist(1,2.0,'kpc') #(inst, value, unit)

## ---- 1.2: Set up the X-ray model
## Before doing this, must run pySPEX_xmm_OKedge.py in order to save the
## parameters of an initial fit

# Set the abundance reference
a=s.abun('solar')

# Add fitting components
s.com('pow')
s.com('hot')
s.com('hot')
s.com('hot')
s.com('slab')
s.com_rel(1,1, np.array([2,3,4,5]))

s.log_exe("xmmparams.com")
s.calc()
s.plot_data(xlog=False, ylog=False, title='O K fit with SPEX', show=True)

## ---- 1.3: Set up a function for evaluating the X-ray model statistic

def eval_xray_stat(ninst=1):
    s.calc()
    cstat = s.fit_cstat() # (cstat, dof)
    return cstat[0]

print("First Cstat value for X-ray spectrum:")
print(eval_xray_stat())

## ---- 2: Load the UV data and setup model with SPEX

## ---- 2.1: Load the data with SPEX for now

s.data('data_A_combined.res','data_A_combined.spo')

# Bin
s.bin(2,1,0,1600,2,'ang') #(inst, reg, elow, ehigh, factor, unit=None)

# Ignore the spectrum, except the oxygen region
s.ignore(2,1,0,1300,'ang')
s.ignore(2,1,1303,1800,'ang')

# Set the distance of the source
s.dist(2,2.0,'kpc')

## ---- 2.2: Set up the model

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
