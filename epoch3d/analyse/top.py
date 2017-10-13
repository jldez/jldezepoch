import numpy as np
import sdf
from pylab import *
import os
import matplotlib.pyplot as plt
import sys
import scipy.constants as c
import math
np.seterr(divide='ignore')
rcParams['ps.useafm'] = True
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
matplotlib.rcParams.update({'font.size': 20})

###############################################################################
dir="../silica_800_3.5_10/pulse1/"
depth0=11
depth=depth0 + 0
f=30
xmax,ymax=1200,1200

data = sdf.read(dir + "%04d.sdf"%f)
#print(data.__dict__.keys())

###############################################################################
key_ez,key_ex,key_ey = "Electric_Field_Ez","Electric_Field_Ex","Electric_Field_Ey"
key_rho,key_temp = "electron_density_Drude","electron_temperature"
key_rho_sfi,key_rho_col = "electron_density_sfi","electron_density_col"
key_mm = "Medium_mask"
key_g,key_g_cb = "dynamic_gamma_drude","dynamic_gamma_drude_cb"
key_ee,key_ep,key_en,key_ei,key_max="g_ee","g_ep","g_en","g_ei","g_max"

c=2.99792458e8
eps0=8.854187817e-12
q=1.602176565e-19
kb=1.38064852e-23
m0=9.10938291e-31
pi=3.14159265359
hbar=1.054571726e-34
atom_density=2.2011e28
bandgap=9.0*q
m_redu=0.5
m_el_cb=1.0
lambd=800.0e-9
omega=2.0*pi*c/lambd
n0=1.45

def get_data(key):
	raw = data.__dict__[key]
	array = raw.data[depth,:,:]
	return array

def fig(dat,norm=1):
	imshow(rot90(dat)/norm,cmap="seismic",extent=[0,xmax,0,ymax])
	colorbar()

mm=get_data(key_mm)
print mm.mean()
fig(mm)

# E_mag=get_data(key_ex)**2.0+get_data(key_ey)**2.0+get_data(key_ez)**2.0
# fig(E_mag)

# rho,Te = get_data(key_rho),get_data(key_temp)
# energy_density = (1.5*kb*Te*rho + bandgap*rho)/(atom_density*q)
# fig(rho,atom_density)

show()
##############################################################################

