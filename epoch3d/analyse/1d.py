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
dir=sys.argv[1]
f=int(sys.argv[2])
path = dir + "%04d.sdf"%f

print "loading " + path
data = sdf.read(path)
#print(data.__dict__.keys())

###############################################################################
key_ez,key_ex,key_ey = "Electric_Field_Ez","Electric_Field_Ex","Electric_Field_Ey"
key_rho = "electron_density_Drude"
key_mm = "Medium_mask"

def get_data(key):
	raw = data.__dict__[key]
	array = raw.data
	return array

def plot_1d(a,l):
	a,amax = a[:,1,1],abs(a)[:,1,1].max()
	print l + " : max value is " + str(amax)
	return plot(a/amax,label=l)

ex,ey,ez = get_data(key_ex),get_data(key_ey),get_data(key_ez)
rho,mm = get_data(key_rho),get_data(key_mm)

arrays,labels = [ex,ey,ez,rho,mm],["ex","ey","ez","rho","mm"]
for a in range(5):
	fig = plot_1d(arrays[a],labels[a])
legend()
ylim(-0.5,1.5)
show()
##############################################################################

