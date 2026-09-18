#
#       FieldEvolutionAnimate_EBPolTor.py
#       Produces animation of magnetic and electric field evolution 
#       showing poloidal field as lines and toroidal field as contours.
#
###############################################################################

import re
import numpy as np
import math as mt
import pylab as pylab
from scipy import integrate
from functools import partial
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.tri as tri
from matplotlib import figure
from matplotlib import rc
import matplotlib.animation as animation
import matplotlib as mpl
mpl.rcParams['animation.ffmpeg_path'] = r'/usr/bin/ffmpeg'

from mpl_toolkits.axes_grid1 import make_axes_locatable

import glob
import os.path

import h5py

from mpl_toolkits.axes_grid1 import make_axes_locatable

import glob
import os.path

c = 29979245800 #speed of light in cm/s
t_year = 3600*24*365 #1 year in seconds
t_day = 3600*24 #1 day in seconds
L_km = 1e5 #1 km in cm
yr = 3600*24*365 #1 year in seconds
B_0Local = 1e13 #in G.

OutputFolder = 'EMHD_Sim_Data_LLB'

Nfrac = 5 #use every Nth time step in video

with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', mode='r') as file:

    #T_0 = file.attrs['Advection speed']    
    # print(file['tasks'].keys())
    
    Lx = file.attrs['Lx']
    Ly = file.attrs['Ly']
    B_0 = file.attrs['B_0']
    L_0 = file.attrs['L_0']
    t_0 = file.attrs['t_0']
    
    B_data = file['B'][:]*B_0/B_0Local #in units of B_0 -> B_0Local
    D_data = file['D'][:]*B_0/B_0Local #in units of B_0 -> B_0Local
    x = file['x'][0,:]*L_0/L_km #in units of L_0->km
    y = file['y'][0,:]*L_0/L_km #in units of L_0->km
    t = file['t'][0,:,0] #in units of t_0
    # U_B = file['U_B'][0,:,0] #in units of B_0^2*L_0^2
    # JH = file['JH'][0,:,0] #in units of B_0^2*L_0^2/t_0
    # PF = file['PF'][0,:,0] #in units of B_0^2*L_0^2/t_0
    # DeltaEInt = file['DeltaEInt'][0,:,0] #in units of B_0^2*L_0^2

file.close()

print(f"Number of timesteps in data set: {len(t):.0f}")
    
# fig, ax = plt.subplots(ncols = 1, nrows = 2, figsize=(3.3,6.5), dpi=150)
# fig.subplots_adjust(left=0.2, bottom=0.1, right=0.96, top=0.96, wspace=None, hspace=None)

# #fig.tight_layout()

mpl.rcParams['mathtext.fontset'] = 'cm' #sets font to LaTeX font

deltax_dat = np.empty([0])
for i in range(len(x)):
    if i == 0:
        deltax_dat = np.append( deltax_dat, 2*(x[0]-0) )
    else:
        deltax_dat = np.append( deltax_dat, 2*(x[i]-x[i-1])-deltax_dat[-1] )

deltay_dat = ( y[1] - y[0] )

Bxtemp = B_data[0,0]
Bx = np.empty_like(B_data[0,0,0:-1])
#For Bx, take averages
for j in range(len(Bxtemp[:,0])-1):
    Bx[j] = ( Bxtemp[j] + Bxtemp[j+1] )/2
Bx = np.flip(np.rot90(Bx,axes=(0,1)),axis=0)
By = np.flip(np.rot90(B_data[0,1,0:-1],axes=(0,1)),axis=0)
Bz = np.flip(np.rot90(B_data[0,2,0:-1],axes=(0,1)),axis=0)

ΨB = np.flip( np.rot90( np.empty_like(B_data[0,0,0:-1]) ), axis=0)
for j in range(0,len(Bx[:,0])): #loop over y-direction, which has been rotated into first index
    for i in range(0,len(Bx[0,:])): #loop over x-direction, which has been rotated into second index
        ΨB[j,i] = np.trapezoid(Bx[0:j+1,i],x=y[0:j+1]) - np.trapezoid(By[j,0:i+1],x=x[0:i+1])
ΨB[0,:] = 0

ΨB = ΨB - np.mean(ΨB)
Bmax = np.amax(ΨB) #maximum value of ΨB at t=0. Used to enforce correct periodicity on ΨB throughout simulation
Bmin = np.amin(ΨB) #minimum value of ΨB at t=0. Used to enforce correct periodicity on ΨB throughout simulation

Dxtemp = D_data[0,0]
Dx = np.empty_like(D_data[0,0,0:-1])
#For Dx, take averages
for j in range(len(Dxtemp[:,0])-1):
    Dx[j] = ( Dxtemp[j] + Dxtemp[j+1] )/2
Dx = np.flip(np.rot90(Dx,axes=(0,1)),axis=0)
Dy = np.flip(np.rot90(D_data[0,1,0:-1],axes=(0,1)),axis=0)
Dz = np.flip(np.rot90(D_data[0,2,0:-1],axes=(0,1)),axis=0)

ΨD = np.flip( np.rot90( np.empty_like(D_data[0,0,0:-1]) ), axis=0)
for j in range(0,len(Dx[:,0])): #loop over y-direction, which has been rotated into first index
    for i in range(0,len(Dx[0,:])): #loop over x-direction, which has been rotated into second index
        ΨD[j,i] = np.trapezoid(Dx[0:j+1,i],x=y[0:j+1]) - np.trapezoid(Dy[j,0:i+1],x=x[0:i+1])
ΨD[0,:] = 0

ΨD = ΨD - np.mean(ΨD)
Dmax = np.amax(ΨD) #maximum value of ΨD at t=0. Used to enforce correct periodicity on ΨD throughout simulation
Dmin = np.amin(ΨD) #minimum value of ΨD at t=0. Used to enforce correct periodicity on ΨD throughout simulation


time = t[0].item()*t_0

ngridy = 200
ngridx = 200

# Create grid values first.
if( ngridx ):
    xi = np.linspace(x[0]-deltax_dat[0]/2, x[-2]+deltax_dat[-1]/2, ngridx) #x coordinate in L_0. x[-2]+deltax_dat/2 is the upper edge since the entire final cell is outside the domain
    yi = np.linspace(y[0]-deltay_dat/2, y[-1]+deltay_dat/2, ngridy) #y coordinate in L_0    

# if( ngridx ):
#     xc = xi[0:-1]+(xi[1]-xi[0])/2# np.linspace(x[0]-deltax_dat[0]/2, x[-2]+deltax_dat[-1]/2, ngridx-1) #x coordinate in L_0. x[-2]+deltax_dat/2 is the upper edge since the entire final cell is outside the domain
#     yc = yi[0:-1]+(yi[1]-yi[0])/2 #np.linspace(y[0]-deltay_dat/2, y[-1]+deltay_dat/2, ngridy-1) #y coordinate in L_0    

Yi,Xi = np.meshgrid(yi, xi) #interpolation grid

Bxi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of Bx
Byi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of By
Bzi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of Bz

Dxi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of Dx
Dyi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of Dy
Dzi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of Dz

#Set coordinates for edge grid points to be at edges and not edge cell centers. Makes contour plots look better
y[0] = y[0] - deltay_dat/2
y[-1] = y[-1] + deltay_dat/2

x[0] = x[0] - deltax_dat[0]/2
x[-2] = x[-2] + deltax_dat[-1]/2
#If len(y) is even, make the middle two values almost identical- improves symmetry of plot
# if( len(y) % 2 == 0):
#     y[int(len(y)/2)-1] = y[int(len(y)/2)-1]/1000;
#     y[int(len(y)/2)] = y[int(len(y)/2)]/1000;

ycoords = np.array([])
xcoords = np.array([])
for j in range(len(y)):
    for i in range(len(x)-1):
        ycoords = np.append(ycoords,y[j])
for i in range(len(y)):
   xcoords = np.concatenate((xcoords,x[0:-1]))
triang = tri.Triangulation(ycoords,xcoords)

BzI = Bz.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(BzI,1.))
Bzi_Init = interpolator(Yi,Xi)
Bzi = np.ma.masked_invalid(Bzi_Init)

ΨBI = ΨB.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(ΨBI,1.))
ΨBi_Init = interpolator(Yi,Xi)
ΨBi = np.ma.masked_invalid(ΨBi_Init)

DzI = Dz.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(DzI,1.))
Dzi_Init = interpolator(Yi,Xi)
Dzi = np.ma.masked_invalid(Dzi_Init)

ΨDI = ΨD.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(ΨDI,1.))
ΨDi_Init = interpolator(Yi,Xi)
ΨDi = np.ma.masked_invalid(ΨDi_Init)

fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize=(10,6.5), dpi=200)
fig.subplots_adjust(left=0.075, bottom=0.3, right=1, top=0.95, wspace=None, hspace=None)

ax1 = ax[0]
ax2 = ax[1]

ax1.set_title(r'$B_{\rm pol}$ (lines) and $B_{z,13}$ (color)',fontsize='16') #-B_{\rm pol})/B_0

# Ψcontours = np.linspace( 1.05*np.amin(Ψi), np.amax(Ψi)-0.05*np.amin(Ψi), 20 )
ΨBcontours = np.linspace( np.amin(ΨBi), (np.amax(ΨBi)-np.amin(ΨBi)), 26 )
ΨDcontours = np.linspace( np.amin(ΨDi), (np.amax(ΨDi)-np.amin(ΨDi)), 26 )

if( np.amax(np.abs(Bz)) <= 1e-10 ):
    cntrf1 = ax1.pcolormesh(yi, xi, Bzi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
else:
    cntrf1 = ax1.pcolormesh(yi, xi, Bzi, vmin=np.amin(Bz), vmax=np.amax(Bz), cmap="Spectral_r")
cntr1 = ax1.contour(yi, xi, ΨBi, ΨBcontours, colors='k', linestyles='solid')

if( np.amax(np.abs(Dz)) <= 1e-10 ):
    cntrf2 = ax2.pcolormesh(yi, xi, Dzi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
else:
    cntrf2 = ax2.pcolormesh(yi, xi, Dzi, vmin=np.amin(Dz), vmax=np.amax(Dz), cmap="Spectral_r")
cntr2 = ax2.contour(yi, xi, ΨDi, ΨDcontours, colors='k', linestyles='solid')

ax1.set_ylabel(r'$x$ (km)',fontsize='16')
ax2.set_ylabel(r'$x$ (km)',fontsize='16')
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
ax1.set_yticks([0,Ly/5,2*Ly/5,3*Ly/5,4*Ly/5],labels=['0','0.2','0.4','0.6','0.8'])
ax2.set_yticks([0,Ly/5,2*Ly/5,3*Ly/5,4*Ly/5],labels=['0','0.2','0.4','0.6','0.8'])
# ax1.set_yticks([0,Ly/5,2*Ly/5,3*Ly/5,4*Ly/5,5*Ly/5],labels=['0','0.2','0.4','0.6','0.8','1'])
# ax2.set_yticks([0,Ly/5,2*Ly/5,3*Ly/5,4*Ly/5,5*Ly/5],labels=['0','0.2','0.4','0.6','0.8','1'])
ax1.tick_params(axis='y',direction='in')
ax2.tick_params(axis='y',direction='in',color='white')
ax1.set_xlabel(r'$y$ (km)',fontsize='16')
ax2.set_xlabel(r'$y$ (km)',fontsize='16')

ax2.set_autoscale_on(False)
x01, y01, w1, h1 = ax1.get_position().bounds
x02, y02, w2, h2 = ax2.get_position().bounds
ax1.set_position([x01, y01, w1, h1])
ax2.set_position([x01+w1, y02, w2, h2])

# Add colorbars
cbar_ycoord = 0.145

cbar_ax1 = fig.add_axes([x01+0.05*w1,cbar_ycoord,0.9*w1,0.04])
cb1 = fig.colorbar(cntrf1, cax=cbar_ax1, orientation="horizontal")

cbar_ax2 = fig.add_axes([x02+0.05*w2,cbar_ycoord,0.9*w2,0.04])
cb2 = fig.colorbar(cntrf2, cax=cbar_ax2, orientation="horizontal")

plt.show()
plt.close()

# DeltaU_B = U_B - U_B[0] #change in magnetic field energy w.r.t. initial value. Should be negative

# #Compute total energy losses to Joule heating and outward Poynting flux. Should be negative
# DeltaU_B[0] = 1e-30 #set first value of DeltaU_B and DeltaJHPF arrays equal to a small number to avoid divide by zero errors at t=0
# DeltaEInt[0] = 1e-30
# # DeltaJHPF = [1e-30]
# # for j in range(1,len(t)):
# #     # if( j == 1 ): #for first point, use 
# #     #     DeltaJHPF.append( 0.5*(JH[0]+PF[0]+JH[1]+PF[1])*(t[1]-t[0]) )
# #     # else:
# #         DeltaJHPF.append( integrate.simpson(JH[0:j+1]+PF[0:j+1],x=t[0:j+1]) )

# # DeltaJHPF = np.asarray(DeltaJHPF) #convert to numpy array
# EConErrorPercent = np.abs(1-DeltaEInt/DeltaU_B)*100

# fig = plt.figure(figsize=(6,4),dpi=150)
# ax = fig.add_subplot(111)
# plt.margins(x=0)
# ax.plot(t*t_0/yr,EConErrorPercent,linewidth='1.5',color='red')
# plt.xlabel(r'$t$ (yr)',fontsize=15)
# plt.ylabel(r'Energy conservation error (%)',fontsize=14,labelpad=1)
# # ax.set_ylim([0.,np.amax(np.abs(1-DeltaU_B/DeltaJHPF)*100)])
# if( np.amax(EConErrorPercent) > 100. ):
#     ax.set_ylim([0.,100.])
# else: ax.set_ylim([0.,np.amax(EConErrorPercent)])
# ax.set_xlim([0.,t[-1]*t_0/yr])
# ax.tick_params(which='both',direction="in",labelsize=14)
# ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1e3))
# # ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(50))
# ax.xaxis.set_ticks_position('both')
# ax.yaxis.set_ticks_position('both')
# # plt.annotate(r'Re$[\omega_g]$',xy=(1,0.85),xytext=(0.15,0.9),fontsize=15,xycoords='axes fraction')
# # plt.hlines(1100, 2.54, 2.98, linestyles="-", linewidth=1.5, color='black')
# # plt.annotate(r'Im$[\omega_g]$',xy=(1,0.85),xytext=(0.30,0.9),fontsize=15,xycoords='axes fraction')
# # plt.hlines(1100, 3.06, 3.49, linestyles="--", linewidth=1.5, color='black')
# # plt.legend(loc=(0.67,0.59),frameon=False,fontsize=14,labelspacing=0.3)
# plt.savefig('EnergyConsPlots/EnergyCons'+OutputFolder+'.pdf',bbox_inches = 'tight',pad_inches = 0.07)

# plt.show()