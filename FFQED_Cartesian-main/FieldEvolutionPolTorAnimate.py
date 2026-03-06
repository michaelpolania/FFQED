#
#       FieldEvolutionPolTorAnimate.py
#       Produces animation of field evolution 
#       showing poloidal field as lines and toroidal field
#       as contours.
#       Also shows local Joule heating rate ~ J^2 as contours
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

from mpl_toolkits.axes_grid1 import make_axes_locatable

import glob
import os.path
import shutil

import h5py

c = 29979245800 #speed of light in cm/s
t_year = 3600*24*365 #1 year in seconds
t_day = 3600*24 #1 day in seconds
L_km = 1e5 #1 km in cm
yr = 3600*24*365 #1 year in seconds
B_0Local = 1e13 #in G

# Auto-detect latest EMHD_Sim_Data_* output folder if present. If you want to
# open a specific run, set OutputFolder to that folder name (e.g. 'EMHD_Sim_Data_LLB').
candidates = sorted(glob.glob('EMHD_Sim_Data_*'))
if len(candidates) == 0:
    raise FileNotFoundError("No EMHD_Sim_Data_* folders found. Run the simulation first or set OutputFolder explicitly in this script.")
OutputFolder = candidates[-1]

# Prefer the ffmpeg binary on PATH or common Homebrew location; fall back to
# ImageMagick or Pillow writer. This avoids hardcoding '/usr/bin/ffmpeg' which
# doesn't exist on many systems (Homebrew installs to /opt/homebrew/bin).
ffmpeg_path = shutil.which('ffmpeg')
if not ffmpeg_path:
    # common Homebrew prefix for Apple Silicon/macOS
    possible = ['/opt/homebrew/bin/ffmpeg', '/usr/local/bin/ffmpeg']
    for p in possible:
        if os.path.exists(p):
            ffmpeg_path = p
            break

if ffmpeg_path:
    mpl.rcParams['animation.ffmpeg_path'] = ffmpeg_path
    print(f"Found ffmpeg at: {ffmpeg_path}")
else:
    print("ffmpeg not found on PATH or common locations. Will try ImageMagick or Pillow writer.")

Nfrac = 5 #use every Nth time step in video

PlotOption = 'vc' # 'J' to plot current density or 'vc' to plot lattice velocity in second panel
vcmax = 3e-5 # maximum value of vc in units of cm/s to plot 

# prefer the canonical file <OutputFolder>/<OutputFolder>.h5 but fall back to any
# .h5 file found inside the folder (this handles cases like EMHD_Sim_Data_LLB/EMHD_Sim_Data_LLB.h5)
h5path = os.path.join(OutputFolder, OutputFolder + '.h5')
if not os.path.exists(h5path):
    h5_files = glob.glob(os.path.join(OutputFolder, '*.h5'))
    if len(h5_files) == 0:
        raise FileNotFoundError(f"No .h5 file found in '{OutputFolder}'. Expected '{h5path}' or any .h5 inside the folder.")
    h5path = h5_files[0]

print(f"Using HDF5 file: {h5path}")

with h5py.File(h5path, mode='r') as file:

    #T_0 = file.attrs['Advection speed']    
    # print(file['tasks'].keys())
    
    Lx = file.attrs['Lx']
    Ly = file.attrs['Ly']
    B_0 = file.attrs['B_0']
    L_0 = file.attrs['L_0']
    t_0 = file.attrs['t_0']
    
    B = file['B'][:]*B_0/B_0Local #in units of B_0 -> B_0Local
    vc = file['vc'][:]*L_0/t_0 #in units of cm/s
    x = file['x'][0,:]*L_0/L_km #in units of L_0->km
    y = file['y'][0,:]*L_0/L_km #in units of L_0->km
    t = file['t'][0,:,0] #in units of t_0
    U_B = file['U_B'][0,:,0] #in units of B_0^2*L_0^2
    JH = file['JH'][0,:,0] #in units of B_0^2*L_0^2/t_0
    PF = file['PF'][0,:,0] #in units of B_0^2*L_0^2/t_0
    DeltaEInt = file['DeltaEInt'][0,:,0] #in units of B_0^2*L_0^2


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

Bxtemp = B[0,0]
Bx = np.empty_like(B[0,0,0:-1])
#For Bx, take averages
for j in range(len(Bxtemp[:,0])-1):
    Bx[j] = ( Bxtemp[j] + Bxtemp[j+1] )/2
Bx = np.flip(np.rot90(Bx,axes=(0,1)),axis=0)
By = np.flip(np.rot90(B[0,1,0:-1],axes=(0,1)),axis=0)
Bz = np.flip(np.rot90(B[0,2,0:-1],axes=(0,1)),axis=0)

Ψ = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
for j in range(0,len(Bx[:,0])): #loop over y-direction, which has been rotated into first index
    for i in range(0,len(Bx[0,:])): #loop over x-direction, which has been rotated into second index
        Ψ[j,i] = np.trapezoid(Bx[0:j+1,i],x=y[0:j+1]) - np.trapezoid(By[j,0:i+1],x=x[0:i+1])
Ψ[0,:] = 0

Ψ = Ψ - np.mean(Ψ)
Ψmax = np.amax(Ψ) #maximum value of Ψ at t=0. Used to enforce correct periodicity on Ψ throughout simulation
Ψmin = np.amin(Ψ) #minimum value of Ψ at t=0. Used to enforce correct periodicity on Ψ throughout simulation

if(PlotOption == 'J'):
    Jx = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
    Jy = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
    Jz_1 = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
    Jz_2 = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)

    for j in range(0,len(Bz[:,0])):
        Jy[j,:] = -np.gradient(Bz[j,:],x[0:-1])
        Jz_1[j,:]= np.gradient(By[j,:],x[0:-1])
    for j in range(1,len(Bz[:,0])-1):
        Jx[j,:] = ( Bz[j+1,:] - Bz[j-1,:] )/(2*deltay_dat)
        Jz_2[j,:]= -( Bx[j+1,:] - Bx[j-1,:] )/(2*deltay_dat)
    Jx[0,:] = ( Bz[1,:] - Bz[-1,:] )/(2*deltay_dat)
    Jz_2[0,:]= -( Bx[1,:] - Bx[-1,:] )/(2*deltay_dat)
    Jx[-1,:] = ( Bz[0,:] - Bz[-2,:] )/(2*deltay_dat)
    Jz_2[-1,:]= -( Bx[0,:] - Bx[-2,:] )/(2*deltay_dat)
        
    J = 1/(4*np.pi)*np.sqrt( Jx**2 + Jy**2 + (Jz_1 + Jz_2)**2 )
elif(PlotOption == 'vc'):
    vcytemp = vc[0,1]
    vcztemp = vc[0,2]
    vc_y = np.empty_like(vc[0,1,0:-1])
    vc_z = np.empty_like(vc[0,2,0:-1])
    #For vc_z, take averages
    for i in range(len(vcztemp[:,0])-1):
        for j in range(len(vcztemp[0,:])-1):  
            vc_z[i,j] = ( vcztemp[i,j] + vcztemp[i,j+1] + vcztemp[i+1,j] + vcztemp[i+1,j+1] )/4
            vc_y[i,j] = (vcytemp[i,j]+vcytemp[i+1,j])/2.
        vc_z[i,-1] = ( vcztemp[i,-1] + vcztemp[i,0] + vcztemp[i+1,-1] + vcztemp[i+1,0] )/4
        vc_y[-1,j] = vcytemp[-1,j]
    vc_y = np.flip(np.rot90(vc_y,axes=(0,1)),axis=0)
    vc_z = np.flip(np.rot90(vc_z,axes=(0,1)),axis=0)
    J = vc_z #np.sqrt( vc_z**2 + vc_y**2 )

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

ΨI = Ψ.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(ΨI,1.))
Ψi_Init = interpolator(Yi,Xi)
Ψi = np.ma.masked_invalid(Ψi_Init)

JI = J.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(JI,1.))
Ji_Init = interpolator(Yi,Xi)
Ji = np.ma.masked_invalid(Ji_Init)

fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize=(10,6.5), dpi=200)
fig.subplots_adjust(left=0.075, bottom=0.3, right=1, top=0.95, wspace=None, hspace=None)

ax1 = ax[0]
ax2 = ax[1]

ax1.set_title(r'$B_{\rm pol}$ (lines) and $B_{z,13}$ (color)',fontsize='16') #-B_{\rm pol})/B_0
if(PlotOption == 'J'): ax2.set_title(r'$J/(cB_0/L_0)$',fontsize='16')
elif(PlotOption == 'vc'): ax2.set_title(r'$v_{c,z}$ (cm/s)',fontsize='16')

# Ψcontours = np.linspace( 1.05*np.amin(Ψi), np.amax(Ψi)-0.05*np.amin(Ψi), 20 )
Ψcontours = np.linspace( np.amin(Ψi), (np.amax(Ψi)-np.amin(Ψi)), 26 )

if( np.amax(np.abs(Bz)) <= 1e-10 ):
    cntrf1 = ax1.pcolormesh(yi, xi, Bzi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
else:
    cntrf1 = ax1.pcolormesh(yi, xi, Bzi, vmin=np.amin(Bz), vmax=np.amax(Bz), cmap="Spectral_r")
cntr1 = ax1.contour(yi, xi, Ψi, Ψcontours, colors='k', linestyles='solid')
if(PlotOption == 'J'):
    if( np.amax(np.abs(Ji)) <= 1e-5 ):
        cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=0, vmax=1, shading='nearest', cmap="inferno")
    else:
        cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=np.amin(Ji), vmax=np.amax(Ji), shading='nearest', cmap="inferno")
elif(PlotOption == 'vc'):
    if( np.amax(np.abs(Ji)) <= 1e-30 ):
        cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=-vcmax, vmax=vcmax, shading='nearest', cmap="inferno")
    else:
        vcmaxtemp = np.amax(np.abs(Ji))
        cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=-vcmaxtemp, vmax=vcmaxtemp, shading='nearest', cmap="inferno")

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
cb1 = fig.colorbar(cntrf1, cax=cbar_ax1, orientation='horizontal')
cb1.ax.xaxis.set_ticks_position('bottom')
cb1.ax.tick_params(rotation=60)

cbar_ax2 = fig.add_axes([x01+1.05*w1,cbar_ycoord,0.9*w1,0.04])
cb2 = fig.colorbar(cntrf2, cax=cbar_ax2, orientation='horizontal')
cb2.ax.xaxis.set_ticks_position('bottom')
cb2.ax.tick_params(rotation=60)   

t_xcoord = 0.025
t_ycoord = 0.215

tannotation = fig.text(t_xcoord, t_ycoord, f'$t=${time/yr:.2f} yr', fontsize=16, rotation=0)

# ax1.axvline(0,0,1,color='black')
# ax2.axvline(0,0,1,color='black')
# ax3.axvline(0,0,1,color='black')
plt.savefig('ContourMagEvo'+OutputFolder+'.pdf',pad_inches = 0.07)

##### Plot B field each component of B on a contour plot #####################################################
def BFieldPolTorPlot(frame,B,t,Nfrac):
    
    global cntrf1, cntr1, cntrf2, fig, vcmax, Ψcontours, Ψmax, Ψmin
    global cb1, cb2, cbar_ax1, cbar_ax2, tannotation
    global triang, yi, xi, Yi, Xi
    
    i_t = frame*Nfrac
    time = t[i_t].item()*t_0
    
    Bxtemp = B[i_t,0]
    Bx = np.empty_like(B[i_t,0,0:-1])
    #For Bx, take averages
    for j in range(len(Bxtemp[:,0])-1):
        Bx[j] = ( Bxtemp[j] + Bxtemp[j+1] )/2
    Bx = np.flip(np.rot90(Bx,axes=(0,1)),axis=0)
    
    # Bx = np.rot90(B[i,0],axes=(0,1))*B_0
    By = np.flip(np.rot90(B[i_t,1,0:-1],axes=(0,1)),axis=0)
    Bz = np.flip(np.rot90(B[i_t,2,0:-1],axes=(0,1)),axis=0)
       
    Ψ = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
    for j in range(0,len(Bx[:,0])): #loop over y-direction, which has been rotated into first index
        for i in range(0,len(Bx[0,:])): #loop over x-direction, which has been rotated into second index
            Ψ[j,i] = np.trapezoid(Bx[0:j+1,i],x=y[0:j+1]) - np.trapezoid(By[j,0:i+1],x=x[0:i+1])
            
    Ψ = Ψ - np.mean(Ψ)
    # Ψ = np.where(Ψ < Ψmin, Ψ + Ψmin, np.where(Ψ > Ψmax, Ψ - 0.8*Ψmax, Ψ))

    if(PlotOption == 'J'):
        Jx = np.flip( np.rot90( np.empty_like(B[i_t,0,0:-1]) ), axis=0)
        Jy = np.flip( np.rot90( np.empty_like(B[i_t,0,0:-1]) ), axis=0)
        Jz_1 = np.flip( np.rot90( np.empty_like(B[i_t,0,0:-1]) ), axis=0)
        Jz_2 = np.flip( np.rot90( np.empty_like(B[i_t,0,0:-1]) ), axis=0)
    
        for j in range(0,len(Bz[:,0])):
            Jy[j,:] = -np.gradient(Bz[j,:],x[0:-1])
            Jz_1[j,:]= np.gradient(By[j,:],x[0:-1])
        for j in range(1,len(Bz[:,0])-1):
            Jx[j,:] = ( Bz[j+1,:] - Bz[j-1,:] )/(2*deltay_dat)
            Jz_2[j,:]= -( Bx[j+1,:] - Bx[j-1,:] )/(2*deltay_dat)
        Jx[0,:] = ( Bz[1,:] - Bz[-1,:] )/(2*deltay_dat)
        Jz_2[0,:]= -( Bx[1,:] - Bx[-1,:] )/(2*deltay_dat)
        Jx[-1,:] = ( Bz[0,:] - Bz[-2,:] )/(2*deltay_dat)
        Jz_2[-1,:]= -( Bx[0,:] - Bx[-2,:] )/(2*deltay_dat)
            
        J = 1/(4*np.pi)*np.sqrt( Jx**2 + Jy**2 + (Jz_1 + Jz_2)**2 )
    elif(PlotOption == 'vc'):
        vcytemp = vc[i_t,1]
        vcztemp = vc[i_t,2]
        vc_y = np.empty_like(vc[i_t,1,0:-1])
        vc_z = np.empty_like(vc[i_t,2,0:-1])
        #For vc_z, take averages
        for i in range(len(vcztemp[:,0])-1):
            for j in range(len(vcztemp[0,:])-1):  
                vc_z[i,j] = ( vcztemp[i,j] + vcztemp[i,j+1] + vcztemp[i+1,j] + vcztemp[i+1,j+1] )/4
                vc_y[i,j] = (vcytemp[i,j]+vcytemp[i+1,j])/2.
            vc_z[i,-1] = ( vcztemp[i,-1] + vcztemp[i,0] + vcztemp[i+1,-1] + vcztemp[i+1,0] )/4
            vc_y[-1,j] = vcytemp[-1,j]
        vc_y = np.flip(np.rot90(vc_y,axes=(0,1)),axis=0)
        vc_z = np.flip(np.rot90(vc_z,axes=(0,1)),axis=0)
        J = vc_z #np.sqrt( vc_z**2 + vc_y**2 )

    # Jx = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
    # Jy = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
    # Jz_1 = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)
    # Jz_2 = np.flip( np.rot90( np.empty_like(B[0,0,0:-1]) ), axis=0)

    # for i in range(0,len(By[:,0])):
    #     Jy[i,:] = -np.gradient(Bz[i,:],x[0:-1])
    #     Jz_1[i,:]= np.gradient(By[i,:],x[0:-1])
    # for j in range(1,len(Bx[0,:])):
    #     Jx[j,:] = ( Bz[j+1,:] - Bz[j-1,:] )/(2*deltay_dat)
    #     Jz_2[j,:]= -( Bx[j+1,:] - Bx[j-1,:] )/(2*deltay_dat)
    # Jx[0,:] = ( Bz[1,:] - Bz[-1,:] )/(2*deltay_dat)
    # Jz_2[0,:]= -( Bx[1,:] - Bx[-1,:] )/(2*deltay_dat)
    # Jx[-1,:] = ( Bz[0,:] - Bz[-2,:] )/(2*deltay_dat)
    # Jz_2[-1,:]= -( Bx[0,:] - Bx[-2,:] )/(2*deltay_dat)
        
    # J = 1/(4*np.pi)*np.sqrt( Jx**2 + Jy**2 + (Jz_1 - Jz_2)**2 )
    #J = vc_z
    
    cntrf1.remove()
    cntr1.remove()
    cntrf2.remove()
    # Don't need to remove color bars- they get updated when the contour plot objects are updated
    
    tannotation.remove()
    
    BzI = Bz.flatten()
    interpolator = tri.LinearTriInterpolator(triang, np.multiply(BzI,1.))
    Bzi_Init = interpolator(Yi,Xi)
    Bzi = np.ma.masked_invalid(Bzi_Init)
    
    ΨI = Ψ.flatten()
    interpolator = tri.LinearTriInterpolator(triang, np.multiply(ΨI,1.))
    Ψi_Init = interpolator(Yi,Xi)
    Ψi = np.ma.masked_invalid(Ψi_Init)
    
    JI = J.flatten()
    interpolator = tri.LinearTriInterpolator(triang, np.multiply(JI,1.))
    Ji_Init = interpolator(Yi,Xi)
    Ji = np.ma.masked_invalid(Ji_Init)
       
    # Ψcontours = np.linspace( 1.05*np.amin(Ψi), np.amax(Ψi)-0.05*np.amin(Ψi), 20 )
    
    if( np.amax(np.abs(Bz)) <= 1e-10 ):
        cntrf1 = ax1.pcolormesh(yi, xi, Bzi, vmin=-0.1, vmax=0.1, cmap="Spectral_r") #"RdBu_r"
    else:
        cntrf1 = ax1.pcolormesh(yi, xi, Bzi, vmin=np.amin(Bz), vmax=np.amax(Bz), cmap="Spectral_r") #"RdBu_r"
    cntr1 = ax1.contour(yi, xi, Ψi, Ψcontours, colors='k', linestyles='solid')
    if(PlotOption == 'J'):
        if( np.amax(np.abs(Ji)) <= 1e-5 ):
            cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=0, vmax=1, shading='nearest', cmap="inferno")
        else:
            cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=np.amin(Ji), vmax=np.amax(Ji), shading='nearest', cmap="inferno")
    elif(PlotOption == 'vc'):
        if( np.amax(np.abs(Ji)) <= 1e-30 ):
            cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=-vcmax, vmax=vcmax, shading='nearest', cmap="inferno")
        else:
            vcmaxtemp = np.amax(np.abs(Ji))
            cntrf2 = ax2.pcolormesh(yi, xi, Ji, vmin=-vcmaxtemp, vmax=vcmaxtemp, shading='nearest', cmap="inferno")
        
    #Delete and recreate color bars. Allows color bar labels to be updated
    del(cb1)
    del(cb2)
    
    cb1 = fig.colorbar(cntrf1, cax=cbar_ax1, orientation='horizontal')
    cb1.ax.xaxis.set_ticks_position('bottom')
    cb1.ax.tick_params(rotation=60)
    
    cb2 = fig.colorbar(cntrf2, cax=cbar_ax2, orientation='horizontal')
    cb2.ax.xaxis.set_ticks_position('bottom')
    cb2.ax.tick_params(rotation=60)   
    
    tannotation = fig.text(t_xcoord, t_ycoord, f'$t=${time/yr:.2f} yr', fontsize=14, rotation=0)
    
    if( frame%10 == 0 ): print(f"Making animation, frame {frame:.0f} complete")

anim = animation.FuncAnimation(fig, partial(BFieldPolTorPlot,B=B,t=t,Nfrac=Nfrac), frames = int(len(t)/Nfrac), interval = 5, blit=False)

#plt.savefig('ContourMagEvo'+outputFolderSuffix+'Dumb.pdf',bbox_inches = 'tight',pad_inches = 0.07)
writerVideo = animation.FFMpegWriter(fps=30)
# Pick a writer dynamically. Prefer ffmpeg, then ImageMagick, then Pillow.
writerVideo = None
save_ext = 'mp4'
if ffmpeg_path:
    try:
        writerVideo = animation.FFMpegWriter(fps=30)
        save_ext = 'mp4'
        print('Using FFMpegWriter')
    except Exception as e:
        print('FFMpegWriter not available:', e)

if writerVideo is None:
    # ImageMagick uses 'convert' or 'magick' on PATH
    imagemagick = shutil.which('convert') or shutil.which('magick')
    if imagemagick:
        try:
            writerVideo = animation.ImageMagickWriter(fps=30)
            save_ext = 'gif'
            print('Using ImageMagickWriter')
        except Exception as e:
            print('ImageMagickWriter not available:', e)

if writerVideo is None:
    try:
        from matplotlib.animation import PillowWriter
        writerVideo = PillowWriter(fps=30)
        save_ext = 'gif'
        print('Using PillowWriter')
    except Exception:
        writerVideo = None

outfile = f'FieldEvoPolTor{OutputFolder}.{save_ext}'
if writerVideo is None:
    print('\nNo suitable animation writer found (ffmpeg, ImageMagick, Pillow).')
    print('To enable MP4 export install ffmpeg (homebrew: brew install ffmpeg) or ImageMagick,')
    print('or install Pillow for GIF export (pip install pillow). Skipping file save; showing animation interactively instead.')
    anim.save = lambda *args, **kwargs: print('Save skipped: no writer available')
else:
    try:
        anim.save(outfile, writer=writerVideo)
        print(f'Saved animation to {outfile}')
    except FileNotFoundError as e:
        print('Save failed because external writer was not found or not executable:', e)
        print('You can install ffmpeg with: brew install ffmpeg')
plt.show()
plt.close()

DeltaU_B = U_B - U_B[0] #change in magnetic field energy w.r.t. initial value. Should be negative

#Compute total energy losses to Joule heating and outward Poynting flux. Should be negative
DeltaU_B[0] = 1e-30 #set first value of DeltaU_B and DeltaJHPF arrays equal to a small number to avoid divide by zero errors at t=0
DeltaEInt[0] = 1e-30
# DeltaJHPF = [1e-30]
# for j in range(1,len(t)):
#     # if( j == 1 ): #for first point, use 
#     #     DeltaJHPF.append( 0.5*(JH[0]+PF[0]+JH[1]+PF[1])*(t[1]-t[0]) )
#     # else:
#         DeltaJHPF.append( integrate.simpson(JH[0:j+1]+PF[0:j+1],x=t[0:j+1]) )

# DeltaJHPF = np.asarray(DeltaJHPF) #convert to numpy array
EConErrorPercent = np.abs(1-DeltaEInt/DeltaU_B)*100

fig = plt.figure(figsize=(6,4),dpi=150)
ax = fig.add_subplot(111)
plt.margins(x=0)
ax.plot(t*t_0/yr,EConErrorPercent,linewidth='1.5',color='red')
plt.xlabel(r'$t$ (yr)',fontsize=15)
plt.ylabel(r'Energy conservation error (%)',fontsize=14,labelpad=1)
# ax.set_ylim([0.,np.amax(np.abs(1-DeltaU_B/DeltaJHPF)*100)])
if( np.amax(EConErrorPercent) > 100. ):
    ax.set_ylim([0.,100.])
else: ax.set_ylim([0.,np.amax(EConErrorPercent)])
ax.set_xlim([0.,t[-1]*t_0/yr])
ax.tick_params(which='both',direction="in",labelsize=14)
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1e3))
# ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(50))
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
# plt.annotate(r'Re$[\omega_g]$',xy=(1,0.85),xytext=(0.15,0.9),fontsize=15,xycoords='axes fraction')
# plt.hlines(1100, 2.54, 2.98, linestyles="-", linewidth=1.5, color='black')
# plt.annotate(r'Im$[\omega_g]$',xy=(1,0.85),xytext=(0.30,0.9),fontsize=15,xycoords='axes fraction')
# plt.hlines(1100, 3.06, 3.49, linestyles="--", linewidth=1.5, color='black')
# plt.legend(loc=(0.67,0.59),frameon=False,fontsize=14,labelspacing=0.3)
plt.savefig('EnergyCons'+OutputFolder+'.pdf',bbox_inches = 'tight',pad_inches = 0.07)

plt.show()