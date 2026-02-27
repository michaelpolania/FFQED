#
#       FieldEvolutionAnimate.py
#       Produces animation of field evolution
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
import shutil
import os

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

from mpl_toolkits.axes_grid1 import make_axes_locatable

import glob
import os.path

import h5py

c = 29979245800 #speed of light in cm/s
t_year = 3600*24*365 #1 year in seconds
t_day = 3600*24 #1 day in seconds
L_km = 1e5 #1 km in cm
yr = 3600*24*365 #1 year in seconds

# Auto-detect latest EMHD_Sim_Data_* output folder if present. If you want to
# open a specific run, set OutputFolder to that folder name (e.g. 'EMHD_Sim_Data_LLB').
candidates = sorted(glob.glob('EMHD_Sim_Data_*'))
if len(candidates) == 0:
    raise FileNotFoundError("No EMHD_Sim_Data_* folders found. Run the simulation first or set OutputFolder explicitly in this script.")
OutputFolder = candidates[-1]

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
    
    B = file['B'][:] #in units of B_0
    x = file['x'][0,:] #in units of L_0
    y = file['y'][0,:] #in units of L_0
    t = file['t'][0,:,0] #in units of t_0
    U_B = file['U_B'][0,:,0] #in units of B_0^2*L_0^2
    JH = file['JH'][0,:,0] #in units of B_0^2*L_0^2/t_0
    PF = file['PF'][0,:,0] #in units of B_0^2*L_0^2/t_0
    
    Lx = file.attrs['Lx']
    Ly = file.attrs['Ly']
    B_0 = file.attrs['B_0']
    L_0 = file.attrs['L_0']
    t_0 = file.attrs['t_0']

# file is closed automatically by the context manager

print(f"Number of timesteps in data set: {len(t):.0f}")
    
##############  Plot B field ##########################################################
def BFieldPlotTorPol(Btor,StreamF,y_array,z_array,t):
    Bx = Btor.evaluate()['g'][0]
    Psi = StreamF.evaluate()['g'][0]
     
    global cntr1
    cntr1.remove()
    
    ngridy = 200
    ngridz = 200
    
    # Create grid values first.
    yi = np.linspace(y_array[0], y_array[-1], ngridy) #y coordinate in L_0
    zi = np.linspace(z_array[0], z_array[-1], ngridz) #z coordinate in L_0
    Yi,Zi = np.meshgrid(yi,zi) #interpolation grid
    
    Bxi_Init = np.empty([len(yi),len(zi)]) #array for interpolated values of Bx
    Byi_Init = np.empty([len(yi),len(zi)]) #array for interpolated values of By
    Bzi_Init = np.empty([len(yi),len(zi)]) #array for interpolated values of Bz
    
    ycoords = np.array([])
    zcoords = np.array([])
    for j in range(len(y_array)):
        for i in range(len(z_array)):
            ycoords = np.append(ycoords,y_array[j])
    for i in range(len(y_array)):
       zcoords = np.concatenate((zcoords,z_array))
    triang = tri.Triangulation(ycoords,zcoords)
    
    BxI = Bx.flatten()
    interpolator = tri.LinearTriInterpolator(triang, np.multiply(BxI,1.))
    Bxi_Init = interpolator(Yi,Zi)
    Bxi = np.ma.masked_invalid(Bxi_Init)
    
    z_Ar, y_Ar = np.meshgrid(z_array,y_array)
    
    fig, ax = plt.subplots(figsize=(9,9), dpi=150)
    
    # ax1.set_title(r'$\delta B_x/B_0$')
    # ax2.set_title(r'$\delta B_y/B_0$')
    # ax3.set_title(r'$\delta B_z/B_0$')
    
    if( np.amax(Bx) <= 1e-10 ):
        levels1 = np.linspace(0,0.1,10)
    else:
        levels1 = np.linspace(np.amin(Bx),np.amax(Bx),10)
    cntr1 = ax.contourf(yi, zi, Bxi, levels1, cmap="Spectral_r") #"RdBu_r"
    ax.contour(y_Ar, z_Ar, Psi, 30, colors='black', linewidths=0.9)
    ax.set_ylabel(r'$z/L_0$')
    ax.tick_params(axis='y',direction='in')
    ax.set_xlabel(r'$y/L_0$')
    
    fig.subplots_adjust(left=0.0)
    fig.subplots_adjust(bottom=0.18)
    x01, y01, w1, h1 = ax.get_position().bounds
    ax.set_position([x01, y01, w1, h1])
    # ax2.set_position([x01+w1, y02, w2, h2])
    # ax3.set_position([x01+w1+w2, y03, w2, h2])
    
    cbar_ax = fig.add_axes([x01,0,0.9*w1,0.04])
    cb = fig.colorbar(cntr1, cax=cbar_ax, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('bottom')
    cb.ax.tick_params(rotation=60)             #Figure out how to rotate the axis label without generating the labels
    
    fig.text(0.0,0.065,f'$t=${t/yr:.2f} yr',fontsize=14,rotation=0)
    plt.show()

# fig, ax = plt.subplots(ncols = 1, nrows = 2, figsize=(3.3,6.5), dpi=150)
# fig.subplots_adjust(left=0.2, bottom=0.1, right=0.96, top=0.96, wspace=None, hspace=None)

# #fig.tight_layout()

mpl.rcParams['mathtext.fontset'] = 'cm' #sets font to LaTeX font

Bxtemp = B[0,0]
Bx = np.empty_like(B[0,0,0:-1])
#For Bx, take averages
for j in range(len(Bxtemp[:,0])-1):
    Bx[j] = ( Bxtemp[j] + Bxtemp[j+1] )/2
Bx = np.flip(np.rot90(Bx,axes=(0,1)),axis=0)

# Bx = np.rot90(B[i,0],axes=(0,1))*B_0
By = np.flip(np.rot90(B[0,1,0:-1],axes=(0,1)),axis=0)
Bz = np.flip(np.rot90(B[0,2,0:-1],axes=(0,1)),axis=0)
time = t[0].item()*t_0

ngridy = 200
ngridx = 200

deltax_dat = np.empty([0])
for i in range(len(x)):
    if i == 0:
        deltax_dat = np.append( deltax_dat, 2*(x[0]-0) )
    else:
        deltax_dat = np.append( deltax_dat, 2*(x[i]-x[i-1])-deltax_dat[-1] )

deltay_dat = ( y[1] - y[0] )

# deltax = ( x[1] - x[0] )*(np.size(x)-1)/ngridx
# deltay = ( y[1] - y[0] )*np.size(y)/ngridy

# Create grid values first.
if( ngridx ):
    xi = np.linspace(x[0]-deltax_dat[0]/2, x[-2]+deltax_dat[-1]/2, ngridx) #x coordinate in L_0. x[-2]+deltax_dat/2 is the upper edge since the entire final cell is outside the domain
    yi = np.linspace(y[0]-deltay_dat/2, y[-1]+deltay_dat/2, ngridy) #y coordinate in L_0    
    
Yi,Xi = np.meshgrid(yi, xi) #interpolation grid

Bxi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of Bx
Byi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of By
Bzi_Init = np.empty([len(yi),len(xi)]) #array for interpolated values of Bz

#Set coordinates for edge grid points to be at edges and not edge cell centers. Makes contour plots look better
y[0] = y[0] - deltay_dat/2
y[-1] = y[-1] + deltay_dat/2
# y = y - deltay_dat/2
# y[-1] = y[-1] + deltay_dat + deltay

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

BxI = Bx.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(BxI,1.))
Bxi_Init = interpolator(Yi,Xi)
Bxi = np.ma.masked_invalid(Bxi_Init)

ByI = By.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(ByI,1.))
Byi_Init = interpolator(Yi,Xi)
Byi = np.ma.masked_invalid(Byi_Init)

BzI = Bz.flatten()
interpolator = tri.LinearTriInterpolator(triang, np.multiply(BzI,1.))
Bzi_Init = interpolator(Yi,Xi)
Bzi = np.ma.masked_invalid(Bzi_Init)

fig, ax = plt.subplots(ncols = 3, nrows = 1, figsize=(11.5,6.5), dpi=150)
fig.subplots_adjust(left=0.05, bottom=0.3, right=1.1, top=0.95, wspace=None, hspace=None)

ax1 = ax[0]
ax2 = ax[1]
ax3 = ax[2]

ax1.set_title(r'$B_x/B_0$',fontsize='14') #-B_{\rm pol})/B_0
ax2.set_title(r'$B_y/B_0$',fontsize='14')
ax3.set_title(r'$B_z/B_0$',fontsize='14')

if( np.amax(np.abs(Bx)) <= 1e-10 ):
    cntrf1 = ax1.pcolormesh(yi, xi, Bxi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
else:
    cntrf1 = ax1.pcolormesh(yi, xi, Bxi, vmin=np.amin(Bz), vmax=np.amax(Bz), cmap="Spectral_r")
if( np.amax(np.abs(By)) <= 1e-10 ):
    cntrf2 = ax2.pcolormesh(yi, xi, Byi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
elif( np.amax(np.abs(By)) - np.amin(np.abs(By)) < 1e-10 ):
    cntrf2 = ax2.pcolormesh(yi, xi, Byi, vmin=0.9*np.amax(np.abs(By)), vmax=1.1*np.amax(np.abs(By)), cmap="Spectral_r")
else:
    cntrf2 = ax2.pcolormesh(yi, xi, Byi, vmin=np.amin(By), vmax=np.amax(By), cmap="Spectral_r")
if( np.amax(np.abs(Bz)) <= 1e-10 ):
    cntrf3 = ax3.pcolormesh(yi, xi, Bzi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
else:
    cntrf3 = ax3.pcolormesh(yi, xi, Bzi, vmin=np.amin(Bz), vmax=np.amax(Bz), cmap="Spectral_r")

ax1.set_ylabel(r'$x/L_0$',fontsize='14')
ax1.set_yticks([0,Ly/5,2*Ly/5,3*Ly/5,4*Ly/5,5*Ly/5],labels=['0','0.2','0.4','0.6','0.8','1'])
ax2.set_yticks([0,Ly/5,2*Ly/5,3*Ly/5,4*Ly/5],labels=['','','','',''])
ax3.set_yticks([0,Ly/5,2*Ly/5,3*Ly/5,4*Ly/5],labels=['','','','',''])
ax1.tick_params(axis='y',direction='in')
ax2.tick_params(axis='y',direction='in')
ax3.tick_params(axis='y',direction='in')
ax1.set_xlabel(r'$y/L_0$',fontsize='14')
ax2.set_xlabel(r'$y/L_0$',fontsize='14')
ax3.set_xlabel(r'$y/L_0$',fontsize='14')

# fig.subplots_adjust(left=0.0, bottom=0.2)
# fig.subplots_adjust(left=0.1, bottom=-0.2)

ax2.set_autoscale_on(False)
ax3.set_autoscale_on(False)
x01, y01, w1, h1 = ax1.get_position().bounds
x02, y02, w2, h2 = ax2.get_position().bounds
x03, y03, w3, h3 = ax3.get_position().bounds
ax1.set_position([x01, y01, w1, h1])
ax2.set_position([x01+w1, y02, w2, h2])
ax3.set_position([x01+w1+w2, y03, w2, h2])

# Add colorbars
cbar_ycoord = 0.145

cbar_ax1 = fig.add_axes([x01,cbar_ycoord,0.9*w1,0.04])
cb1 = fig.colorbar(cntrf1, cax=cbar_ax1, orientation='horizontal')
cb1.ax.xaxis.set_ticks_position('bottom')
cb1.ax.tick_params(rotation=60)

cbar_ax2 = fig.add_axes([x01+w1,cbar_ycoord,0.9*w1,0.04])
cb2 = fig.colorbar(cntrf2, cax=cbar_ax2, orientation='horizontal')
cb2.ax.xaxis.set_ticks_position('bottom')
cb2.ax.tick_params(rotation=60)   

cbar_ax3 = fig.add_axes([x01+w1+w2,cbar_ycoord,0.9*w1,0.04])
cb3 = fig.colorbar(cntrf3, cax=cbar_ax3, orientation='horizontal')
cb3.ax.xaxis.set_ticks_position('bottom')
cb3.ax.tick_params(rotation=60)

t_xcoord = 0.025
t_ycoord = 0.215

tannotation = fig.text(t_xcoord, t_ycoord, f'$t=${time/yr:.2f} yr', fontsize=14, rotation=0)

# ax1.axvline(0,0,1,color='black')
# ax2.axvline(0,0,1,color='black')
# ax3.axvline(0,0,1,color='black')
# plt.savefig('ContourMagEvo'+OutputFolder+'.pdf',pad_inches = 0.07)

##### Plot each component of B on a contour plot ##########################################################
def BFieldPlot(frame,B,t,Nfrac):
    
    global cntrf1, cntrf2, cntrf3, fig
    global cb1, cb2, cb3, cbar_ax1, cbar_ax2, cbar_ax3, tannotation
    global triang, yi, xi, Yi, Xi
    
    i = frame*Nfrac
    time = t[i].item()*t_0
    
    Bxtemp = B[i,0]
    Bx = np.empty_like(B[i,0,0:-1])
    #For Bx, take averages
    for j in range(len(Bxtemp[:,0])-1):
        Bx[j] = ( Bxtemp[j] + Bxtemp[j+1] )/2
    Bx = np.flip(np.rot90(Bx,axes=(0,1)),axis=0)
    
    # Bx = np.rot90(B[i,0],axes=(0,1))*B_0
    By = np.flip(np.rot90(B[i,1,0:-1],axes=(0,1)),axis=0)
    Bz = np.flip(np.rot90(B[i,2,0:-1],axes=(0,1)),axis=0)
    
    cntrf1.remove()
    cntrf2.remove()
    cntrf3.remove()
    # Don't need to remove color bars- they get updated when the contour plot objects are updated
    
    tannotation.remove()
    
    BxI = Bx.flatten()
    interpolator = tri.LinearTriInterpolator(triang, np.multiply(BxI,1.))
    Bxi_Init = interpolator(Yi,Xi)
    Bxi = np.ma.masked_invalid(Bxi_Init)
    
    ByI = By.flatten()
    interpolator = tri.LinearTriInterpolator(triang, np.multiply(ByI,1.))
    Byi_Init = interpolator(Yi,Xi)
    Byi = np.ma.masked_invalid(Byi_Init)
    
    BzI = Bz.flatten()
    interpolator = tri.LinearTriInterpolator(triang, np.multiply(BzI,1.))
    Bzi_Init = interpolator(Yi,Xi)
    Bzi = np.ma.masked_invalid(Bzi_Init)
           
    if( np.amax(np.abs(Bx)) <= 1e-10 ):
        cntrf1 = ax1.pcolormesh(yi, xi, Bxi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
    else:
        cntrf1 = ax1.pcolormesh(yi, xi, Bxi, vmin=np.amin(Bz), vmax=np.amax(Bz), cmap="Spectral_r")
    if( np.amax(np.abs(By)) <= 1e-10 ):
        cntrf2 = ax2.pcolormesh(yi, xi, Byi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
    elif( np.amax(np.abs(By)) - np.amin(np.abs(By)) < 1e-10 ):
        cntrf2 = ax2.pcolormesh(yi, xi, Byi, vmin=0.9*np.amax(np.abs(By)), vmax=1.1*np.amax(np.abs(By)), cmap="Spectral_r")
    else:
        cntrf2 = ax2.pcolormesh(yi, xi, Byi, vmin=np.amin(By), vmax=np.amax(By), cmap="Spectral_r")
    if( np.amax(np.abs(Bz)) <= 1e-10 ):
        cntrf3 = ax3.pcolormesh(yi, xi, Bzi, vmin=-0.1, vmax=0.1, cmap="Spectral_r")
    else:
        cntrf3 = ax3.pcolormesh(yi, xi, Bzi, vmin=np.amin(Bz), vmax=np.amax(Bz), cmap="Spectral_r")
    
    #Delete and recreate color bars. Allows color bar labels to be updated
    del(cb1)
    del(cb2)
    del(cb3)
    
    cb1 = fig.colorbar(cntrf1, cax=cbar_ax1, orientation='horizontal')
    cb1.ax.xaxis.set_ticks_position('bottom')
    cb1.ax.tick_params(rotation=60)
    
    cb2 = fig.colorbar(cntrf2, cax=cbar_ax2, orientation='horizontal')
    cb2.ax.xaxis.set_ticks_position('bottom')
    cb2.ax.tick_params(rotation=60)   
    
    cb3 = fig.colorbar(cntrf3, cax=cbar_ax3, orientation='horizontal')
    cb3.ax.xaxis.set_ticks_position('bottom')
    cb3.ax.tick_params(rotation=60)
    
    tannotation = fig.text(t_xcoord, t_ycoord, f'$t=${time/yr:.2f} yr', fontsize=14, rotation=0)
    
    if( frame%10 == 0 ): print(f"Making animation, frame {frame:.0f} complete")

Nfrac = 5 #use every Nth time step in video

anim = animation.FuncAnimation(fig, partial(BFieldPlot,B=B,t=t,Nfrac=Nfrac), frames = int(len(t)/Nfrac), interval = 5, blit=False)
# anim = animation.FuncAnimation(fig, animate, frames = int(len(t1)/Nfrac), interval = 5, blit=False)

#plt.savefig('ContourMagEvo'+outputFolderSuffix+'Dumb.pdf',bbox_inches = 'tight',pad_inches = 0.07)
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

outfile = f'FieldEvo{OutputFolder}.{save_ext}'
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
DeltaJHPF = [1e-30]
for j in range(1,len(t)):
    # if( j == 1 ): #for first point, use 
    #     DeltaJHPF.append( 0.5*(JH[0]+PF[0]+JH[1]+PF[1])*(t[1]-t[0]) )
    # else:
        DeltaJHPF.append( integrate.simpson(JH[0:j+1]+PF[0:j+1],x=t[0:j+1]) )

DeltaJHPF = np.asarray(DeltaJHPF) #convert to numpy array

fig = plt.figure(figsize=(6,4),dpi=150)
ax = fig.add_subplot(111)
plt.margins(x=0)
ax.plot(t*t_0/yr,np.abs(1-DeltaU_B/DeltaJHPF)*100,linewidth='1.5',label='$g_1$-mode',color='red')
plt.xlabel(r'$t$ (yr)',fontsize=15)
plt.ylabel(r'Energy conservation error (%)',fontsize=14,labelpad=1)
ax.set_ylim([0.,np.amax(np.abs(1-DeltaU_B/DeltaJHPF)*100)])
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