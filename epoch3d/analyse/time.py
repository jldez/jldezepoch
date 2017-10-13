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
dir="../silica_800_2.0_10/pulse1/"
depth0=11
depth=depth0 + 25
Nt=100
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
m_el_cb=0.5
lambd=800.0e-9
omega=2.0*pi*c/lambd
n0=1.45

def get_array(key):
	sol=zeros((Nt))
	for f in range(Nt):
		path = dir + "%04d.sdf"%f
		data = sdf.read(path)
		raw = data.__dict__[key]
		array = raw.data
		sol[f]=array[depth,1,1]
	if (key=="dynamic_gamma_drude" or key=="dynamic_gamma_drude_cb"):
		sol[0]=0
	return sol

def plot1(key,norm=1,lab=" "):
	a=get_array(key)
	max_val=abs(a).max()
	if (max_val==0):
		return 0
	print key + " : " + str(max_val/norm)
	if (norm==1):
		a = a/max_val
	plot(a/norm,label=lab)

def plot2(key,norm=1,lab="",lw=1,ls="-"):
	a=get_array(key)
	max_val=abs(a).max()
	if (max_val==0):
		return 0
	print key + " : " + str(max_val/norm)
	if (norm==1):
		a = a/max_val
	semilogy(a/norm,label=lab,lw=lw,linestyle=ls)

# plasma stuff
plot2(key_rho,2.2011e28,r"$\rho$")
plot2(key_rho_sfi,2.2011e28,r"$\rho_{\mathrm{sfi}}$")
plot2(key_rho_col,2.2011e28,r"$\rho_{\mathrm{col}}$")
ylim(10**-10,10)

# # gamma stuff
# plot2(key_g,1e15,r"$\gamma_{\mathrm{eff}}$",2)
# # plot2(key_g_cb,1e15,r"$\gamma_{\mathrm{cb}}$",2)
# plot2(key_ee,1e15,r"$\gamma_{\mathrm{ee}}$",2,"--")
# plot2(key_ep,1e15,r"$\gamma_{\mathrm{ep}}$",2,"--")
# plot2(key_en,1e15,r"$\gamma_{\mathrm{en}}$",2,"--")
# plot2(key_ei,1e15,r"$\gamma_{\mathrm{ei}}$",2,"--")
# plot2(key_max,1e15,r"$\gamma_{\mathrm{max}}$",2,"--")
# ylim(10**-4,100)

# # mre stuff
# E=(get_array(key_ex)**2.0+get_array(key_ey)**2.0+get_array(key_ez)**2.0)**0.5
# g,g_ei,g_en=get_array(key_g),get_array(key_ei),get_array(key_en)
# pond_energy = (q*abs(E[depth0]).max())**2.0/(4.0*m_redu*m0*omega**2.0)
# critical_energy=(1.0+m_redu)*(bandgap+pond_energy)
# k=np.ceil(critical_energy/(hbar*omega))
# W1pt=q**2.0/((g+1)*m0*m_redu*(1.0+(omega/(g+1.0))**2.0))*E**2.0 \
# 	/(0.69315*critical_energy)/(2.0**(1.0/k)-1.0)
# max_val,norm=abs(W1pt).max(),1e15
# semilogy(W1pt/norm,label=r"$W_{1\mathrm{pt}}$",linestyle="-",lw=1)
# semilogy((g_ei+g_en)/norm,label=r"$\gamma_{\mathrm{en+ei}}$",lw=1)
# g_ib = ((W1pt+1.0)**-1+(g_ei+g_en+1.0)**-1)**-1
# semilogy(g_ib/norm,label=r"$\gamma_{\mathrm{IB}}$",lw=2)
# ylim(1e-4)
# # print "gamma_eff = ", g.max()
# # print "alpha_eff = ", g_ib.max()*(2.0**(1.0/k)-1.0)/(0.5*n0*eps0*c*E.max()**2.0)

# # temperature stuff
# plot1(key_temp,q/kb,r"$T_e$")
# ylim(0,20)

# # energy stuff
# rho=get_array(key_rho_sfi)+get_array(key_rho_col)
# Te=get_array(key_temp)
# energy_density = 1.5*Te*kb*rho + rho*bandgap*q
# max_val,norm=abs(energy_density).max(),atom_density*q
# print "energy_density" + " : " + str(max_val/norm)
# plot(energy_density/norm,label=r"$\mathcal{E}$")

# # field stuff
# plot1(key_ex,1,r"$E_x$")
# plot1(key_ey,1,r"$E_y$")
# plot1(key_ez,1,r"$E_z$")
# ylim(-1,1)

legend(loc=4)
show()
##############################################################################

