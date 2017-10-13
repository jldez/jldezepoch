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
trap_density = 0.005
xi = 0.05
rho_sat = 2.2011e28
Eb = 1.5
bandGap = 9.0

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
key_rho,key_Te = "electron_density_Drude", "electron_temperature"
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

def mkfig(filename,array,N):
	clf()
	imshow(array[:,1,:],aspect="auto")
	colorbar()
	savefig(path_results+filename+"para/para"+str(N)+".png")
	clf()
	imshow(array[:,:,1],aspect="auto")
	colorbar()
	savefig(path_results+filename+"perp/perp"+str(N)+".png")
	clf()
	imshow(sum(array[:,:,:],axis=0),aspect="auto")
	savefig(path_results+filename+"top/top"+str(N)+".png")

# Initialize mm, ri
if (N==1):
	mm_before = get_data_sdf(key_mm)
	output_data_npy("mm/mm0",mm_before)
	output_data_npy("ri/ri0",mm_before*0.0)

# Load rho_N, mm_N-1, ri_N-1, Te_N-1
rho, Te = get_data_sdf(key_rho), get_data_sdf(key_Te)
mm_before = get_data_npy("mm/mm" + str(N-1) + ".npy")
ri_before = get_data_npy("ri/ri" + str(N-1) + ".npy")

# Incubation
ri_after = (trap_density*rho_sat) * (1.0 - np.exp(-xi*rho/(trap_density*rho_sat)))

# Dammage
mm_after = copy(mm_before)
energy_density = 1.5*Te*1.38064852e-23*rho + rho*bandGap*1.602176565e-19
# for i in range(rho.shape[0]):
# 		for j in range(rho.shape[1]):
# 			for k in range(rho.shape[2]):
# 				if (energy_density[i,j,k] >= Eb*1.602176565e-19*rho_sat):
# 					mm_after[:i,j,k] = 0.0
# 					ri_after[:i,j,k] = 0.0

# Replace the ablated surface max position 30nm 
# from the top boundary of the simulation domain
move_up,n = 0,rho.shape[0]
for i in range(n):
	if (sum(mm_before[i,:,:]) != 0.0 and sum(mm_after[i,:,:]) == 0.0):
		move_up = move_up + 1
	elif (sum(mm_before[i,:,:]) != 0.0 and sum(mm_after[i,:,:]) != 0.0):
		break ;

print "move up by " + str(move_up) + " cells"
mm_moved = copy(mm_after)
ri_moved = copy(ri_after)

if (move_up > 0):
	mm_moved[:-move_up,:,:] = mm_after[move_up:,:,:]
	mm_moved[(n-move_up):,:,:] = [mm_after[n-1,:,:]]*move_up
	ri_moved[:-move_up,:,:] = ri_after[move_up:,:,:]
	ri_moved[(n-move_up):,:,:] = [ri_after[n-1,:,:]]*move_up

# Output rho_N, mm_N, ri_N
output_data_npy("rho/rho" + str(N),rho)
output_data_npy("Te/Te" + str(N),Te)
output_data_txt("mm.txt",mm_moved)
output_data_npy("mm/mm" + str(N),mm_moved)
output_data_txt("ri.txt",ri_moved)
output_data_npy("ri/ri" + str(N),ri_moved)

# figures
mkfig("rho/",rho,N)
mkfig("Te/",Te,N)
mkfig("mm/",mm_moved,N)
mkfig("ri/",ri_moved,N)