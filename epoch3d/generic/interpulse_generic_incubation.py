import numpy as np
import sdf
from pylab import *
import os
import matplotlib.pyplot as plt
import sys
import scipy.constants as c
import math
rcParams['ps.useafm'] = True
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
matplotlib.rcParams.update({'font.size': 20})

###############################################################################
#parametres
zeta = 0.3
eta = 0.3
sigma_max = 10.0
rho_max = 2.2011e28
rho_etoile_0 = 0.5
rho_etoile_min = 0.1

###############################################################################
N=int(sys.argv[1])
f=int(sys.argv[2])
dir=sys.argv[3]
path = dir+"/pulse" + str(N) + "/"
path_next = dir+"/pulse" + str(N+1) + "/"
path_results = dir+"/results/"

print "loading " + path + "%04d.sdf"%f
data = sdf.read(path + "%04d.sdf"%f)
#print(data.__dict__.keys())

###############################################################################
key_ez,key_ex,key_ey = "Electric_Field_Ez","Electric_Field_Ex","Electric_Field_Ey"
key_rho = "electron_density_Drude"
key_mm = "Medium_mask"

def get_data_sdf(key):
	raw = data.__dict__[key]
	array = raw.data[:,:,:]
	return array

def output_data_txt(filename,array):
	file = open(path_next + filename,"w") 
	for i in range(array.shape[0]):
		for j in range(array.shape[1]):
			for k in range(array.shape[2]):
				file.write(str(array[i,j,k])+"\n") 
	file.close() 

def get_data_npy(filename):
	array = np.load(path_results+filename)
	return array

def output_data_npy(filename,array):
	np.save(path_results+filename,array)

# Initialize mm, is and rho_etoile
if (N==1):
	mm_before = get_data_sdf(key_mm)
	output_data_npy("mm/mm0",mm_before)
	output_data_npy("is/is0",mm_before)
	output_data_npy("rho_etoile/rho_etoile0",mm_before*rho_etoile_0)

# Load rho_N, mm_N-1, is_N-1, rho_etoile_N-1
rho = get_data_sdf(key_rho)
mm_before = get_data_npy("mm/mm" + str(N-1) + ".npy")
is_before = get_data_npy("is/is" + str(N-1) + ".npy")
rho_etoile_before = get_data_npy("rho_etoile/rho_etoile" + str(N-1) + ".npy")

# Incubation
is_after = sigma_max-(sigma_max-is_before)*np.exp(-zeta*rho/rho_max)
rho_etoile_after = rho_etoile_min+(rho_etoile_before-rho_etoile_min)*np.exp(-eta*rho/rho_max)

# Dammage
mm_after = copy(mm_before)
for i in range(rho.shape[0]):
		for j in range(rho.shape[1]):
			for k in range(rho.shape[2]):
				if (rho[i,j,k]/rho_max > rho_etoile_before[i,j,k]):
					mm_after[i,j,k] = 0.0

# Output rho_N, mm_N, is_N, rho_etoile_N
output_data_npy("rho/rho" + str(N),rho)
output_data_txt("mm.txt",mm_after)
output_data_npy("mm/mm" + str(N),mm_after)
output_data_txt("is.txt",is_after)
output_data_npy("is/is" + str(N),is_after)
output_data_npy("rho_etoile/rho_etoile" + str(N),rho_etoile_after)