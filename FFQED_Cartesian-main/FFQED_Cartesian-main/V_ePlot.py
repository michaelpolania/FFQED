#
#       V_ePlot.py
#       Produces plot of maximum electron velocity v_e at surface
#
###############################################################################

import re
import numpy as np
import math as mt
import pylab as pylab
from scipy import integrate, optimize
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

c = 29979245800 #speed of light in cm/s
M_sol = 1.98847e33 #solar mass in g
t_year = 3600*24*365 #1 year in seconds
t_day = 3600*24 #1 day in seconds
L_km = 1e5 #1 km in cm
yr = 3600*24*365 #1 year in seconds

VelocityChoice = "Hall" #either "Hall" or "Electron"

OutputFolder = 'EMHD_Sim_Data_8'

with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', mode='r') as file:

    Lx = file.attrs['Lx']
    Ly = file.attrs['Ly']
    B_0 = file.attrs['B_0']
    L_0 = file.attrs['L_0']
    t_0 = file.attrs['t_0']
    
    B = file['B'][:] #in units of B_0
    vc = file['vc'][:]*L_0/(c*t_0) #in units of c
    x = file['x'][0,:] #in units of L_0
    y = file['y'][0,:] #in units of L_0
    t = file['t'][0,:,0] #in units of t_0
    U_B = file['U_B'][0,:,0] #in units of B_0^2*L_0^2
    JH = file['JH'][0,:,0] #in units of B_0^2*L_0^2/t_0
    PF = file['PF'][0,:,0] #in units of B_0^2*L_0^2/t_0

file.close()
    
mpl.rcParams['mathtext.fontset'] = 'cm' #sets font to LaTeX font

########################## Calculating v_e #########################################

deltax_dat = np.empty([0])
for i in range(len(x)):
    if i == 0:
        deltax_dat = np.append( deltax_dat, 2*(x[0]-0) )
    else:
        deltax_dat = np.append( deltax_dat, 2*(x[i]-x[i-1])-deltax_dat[-1] )

deltay = ( y[1] - y[0] ) #since deltay is constant
deltax = deltax_dat[-1] #only need final value for computing v_e at outer surface

########################## Calculating v_e #########################################

unit_e = 4.80320425e-10 #elementary charge in units of statcoulomb
n_e = 3.3e-6 #electron number density at r = 12.4 km in fm^-3

# Arrays of electron velocity at upper surface divided by c
v_ey = np.empty( (len(t),len(y)) )
v_ez = np.empty( (len(t),len(y)) )
v_eperp = np.empty( ( (len(t),len(y)) ) )
v_eperpMean = np.empty( (len(t)) )
v_eperpMax = np.empty( (len(t)) )

#Electron MHD velocities
if(VelocityChoice == "Electron"):
    for n in range(len(t)):
        for j in range(len(y)):
            if(j == len(y)-1):
                v_ez[n,j] = -1./(4.*np.pi*n_e*1e39*unit_e)*B_0/L_0*( ( B[n,1,-1,j] - B[n,1,-2,j] )/deltax - ( B[n,0,-1,1] - B[n,0,-1,j] )/deltay )
            else:
                v_ez[n,j] = -1./(4.*np.pi*n_e*1e39*unit_e)*B_0/L_0*( ( B[n,1,-1,j] - B[n,1,-2,j] )/deltax - ( B[n,0,-1,j+1] - B[n,0,-1,j] )/deltay )
            v_ey[n,j] = -1./(4.*np.pi*n_e*1e39*unit_e)*B_0/L_0*( -( B[n,2,-1,j] - B[n,2,-2,j] )/deltax )
            v_eperp[n,j] = v_ey[n,j]#np.sqrt( v_ey[n,j]**2 + v_ez[n,j]**2 )
        v_eperpMean[n] = np.mean( v_eperp[n,:] )
        v_eperpMax[n] = np.amax( v_eperp[n,:] )
    
#Hall MHD velocities
elif(VelocityChoice == "Hall"):
    for n in range(len(t)):
        for j in range(len(y)):
            if(j == len(y)-1):
                v_ez[n,j] = vc[n,2,-1,j] - 1./(4.*np.pi*n_e*1e39*unit_e)*B_0/L_0*( ( B[n,1,-1,j] - B[n,1,-2,j] )/deltax - ( B[n,0,-1,1] - B[n,0,-1,j] )/deltay )
            else:
                v_ez[n,j] = vc[n,2,-1,j] - 1./(4.*np.pi*n_e*1e39*unit_e)*B_0/L_0*( ( B[n,1,-1,j] - B[n,1,-2,j] )/deltax - ( B[n,0,-1,j+1] - B[n,0,-1,j] )/deltay )
            v_ey[n,j] = vc[n,1,-1,j] - 1./(4.*np.pi*n_e*1e39*unit_e)*B_0/L_0*( -( B[n,2,-1,j] - B[n,2,-2,j] )/deltax )
            v_eperp[n,j] = v_ey[n,j]#np.sqrt( v_ey[n,j]**2 + v_ez[n,j]**2 )
        v_eperpMean[n] = np.mean( v_eperp[n,:] )
        v_eperpMax[n] = np.amax( v_eperp[n,:] )
    
###################################### Plotting  #########################################

mpl.rcParams.update(mpl.rcParamsDefault) #sets font to default...
mpl.rcParams['mathtext.fontset'] = 'cm' #... then sets font to LaTeX font

fig, ax1 = plt.subplots(figsize=(6,5),dpi=100)
plt.margins(x=0.005,y=0.005)   
ax1.plot(t*t_0/yr,v_eperpMean*c*yr,label=r'$\theta=90^{\circ}$',color='black',linewidth=1.5)
ax1.set_xlabel(r'$t$ (years)',fontsize=12)
ax1.set_ylabel(r'$\overline{v}_{L,\perp}(z=L_z)$ (cm/yr)',fontsize=12)
ax1.set_yscale('symlog')
#ax1.set_ylim([0,0.04])
#ax1.set_xlim([0,50])
ax1.legend(loc='best',frameon=False,handlelength=3,fontsize=12)
#ax1.legend(bbox_to_anchor=(1.05, -.1),ncol=2,frameon=False,handlelength=3,fontsize=12)

#ax1.annotate(r'Preliminary',xy=(0.75, -0.13), xycoords='axes fraction',fontsize=14,color='r')
plt.savefig('v_eMeanPlot'+OutputFolder+'.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()
plt.close()

np.savetxt(OutputFolder+'_vL.csv', np.column_stack((t*t_0/yr, v_eperpMean*c*yr)), delimiter=',')

# Timing noise model

nu_0 = 1 #characteristic spin frequency in s^{-1}
nu = np.array( [ 2.579, 1.672, 4.069, 2.469 ] )/nu_0 # spin (angular) frequency in s^-1*nu_0
nudot = np.array( [ -1.1846e-14, -1.7700e-14, -5.8852e-14, -3.65843e-13 ] )*yr/nu_0 # nudot in s^-2*yr*nu_0
R = 1.2e6 #radius in cm
alpha = np.pi/8
sinalpha = np.sin(alpha) #sin of angle between magnetic dipole and rotation axis
Rsintheta = R*1/2**sinalpha #R*sin(theta) in cm

I = 2/5*1.4*M_sol*(1.2e6)**2 #moment of inertia of 1.4 solar mass star with 12 km radius in g*cm^2
Bp_inferred = np.sqrt( -6*I*c**3/R**6*(nudot*nu_0/yr)/(nu*nu_0)**3 )/sinalpha/1e13 #inferred surface dipole magnetic field in 1e12 G

#New model of torque changing

def phi_fit2(t,nu,nudot):
    return nu*t + 0.5*nudot*t**2

def phi_fit3(t,nu,nudot,nudotdot):
    return nu*t + 0.5*nudot*t**2 + 1/6*nudotdot*t**3

t_yr = t*t_0/yr #time in years

#Define timespan to examine
t_st = 10 #in yr
t_end = 50 #in yr
i_st, i_end = np.argmax(t_yr>t_st), np.argmax(t_yr>t_end)
if i_end == 0: i_end = -1 #prevents error if t_end isn't included in t_yr array
if t_yr[-1] < t_end: print("WARNING: t_end is larger than duration of simulation")
t_yr_window = t_yr[i_st:i_end]
v_eperpMean_window = v_eperpMean[i_st:i_end]

nudotdot_t = np.empty( (len(nu),len(t_yr_window)) )
alphadot_t = np.empty( (len(t_yr_window)) )
nudotdotDipole = np.empty( len(nu) )
nu_fit = np.empty( len(nu) )
nudot_fit = np.empty( len(nu) )

phi_fit2_eval = np.empty( (len(nu), len(t_yr_window)) )
phi_fit3_eval = np.empty( (len(nu), len(t_yr_window)) )
residual2 = np.empty( (len(nu), len(t_yr_window)) )
residual3 = np.empty( (len(nu), len(t_yr_window)) )
sigma_TN2 = np.empty( ( len(nu) ) )
sigma_TN3 = np.empty( ( len(nu) ) )

for j in range(len(nu)):
    phi_dat = np.empty( (len(t_yr_window)) )
    for n in range(len(t_yr_window)):
        alphadot = v_eperpMean_window[n]*c/R*yr #alphadot in rad/s*t_0
        alphadot_t[n] = alphadot
        nudotdot = ( 2/np.tan(alpha+alphadot*t_yr_window[n])*alphadot*nudot[j] + 3*nudot[j]**2/nu[j] ) #nu
        nudotdot_t[j,n] = nudotdot #store nudot values
        phi_dat[n] = phi_fit3(t_yr_window[n]-t_yr_window[0],nu[j],nudot[j],nudotdot) #"measurement"
    
    nudotdotDipole[j] = 3*nudot[j]**2/nu[j] #nudotdot if alphadot = 0 i.e. magnetic dipole spindown with no change in alpha
    #Calculate second-order polynomial fit parameters (nu, nudot)
    popt2, pcov2 = optimize.curve_fit(phi_fit2, t_yr_window-t_yr_window[0], phi_dat, p0=(nu[j],nudot[j],), 
                                    bounds=((0.8*nu[j],10*nudot[j]),(1.2*nu[j],0.1*nudot[j])) )
    perr2 = np.sqrt(np.diag(pcov2)) #standard deviation of fit parameters
    
    #Calculate third-order polynomial fit parameters (nu, nudot, nudotdot)
    # popt3, pcov3 = optimize.curve_fit(phi_fit3, t_yr_window-t_yr_window[0], phi_dat, p0=(nu[j],nudot[j],nudotdotDipole[j]), 
    #                                 bounds=((0.9*nu[j],2*nudot[j],-100*nudotdotDipole[j]),(1.1*nu[j],0.5*nudot[j],100*nudotdotDipole[j])) )
    popt3, pcov3 = optimize.curve_fit(phi_fit3, t_yr_window-t_yr_window[0], phi_dat, p0=(nu[j],nudot[j],nudotdotDipole[j]), 
                                    bounds=((0.8*nu[j],10*nudot[j],-np.inf),(1.2*nu[j],0.1*nudot[j],np.inf)) )
    perr3 = np.sqrt(np.diag(pcov3)) #standard deviation of fit parameters
    
    nu_fit[j] = popt2[0] #fitted spin frequency in units of nu_0 using 2nd order polynomial fit
    nudot_fit[j] = popt2[1] #fitted frequency derivative in units of nu_0/yr using 3rd order polynomial fit
    phi_fit2_eval[j] = phi_fit2(t_yr_window-t_yr_window[0],popt2[0],popt2[1])
    phi_fit3_eval[j] = phi_fit3(t_yr_window-t_yr_window[0],popt3[0],popt3[1],popt3[2])

    #Calculate phase residuals
    residual2[j] = phi_dat - phi_fit2_eval[j]
    residual3[j] = phi_dat - phi_fit3_eval[j]

    #Calculate timing noise- standard deviation of residuals
    T = t_yr_window[-1] - t_yr_window[0]
    sigma_TN2[j] = np.sqrt( 1/(T*nu[j]**2)*np.trapezoid( residual2[j]**2, t_yr_window ) )*yr*1e3 #standard deviation of the timing noise (2nd order fit) in ms
    sigma_TN3[j] = np.sqrt( 1/(T*nu[j]**2)*np.trapezoid( residual3[j]**2, t_yr_window ) )*yr*1e3 #standard deviation of the timing noise (3rd order fit) in ms


# Fractional changes of fitted nu and nudot vs. "measured" values
nu_fit_frac = 1-np.abs(nu_fit/nu)
nudot_fit_frac = 1-np.abs(nudot_fit/nudot)

fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
plt.margins(x=0,y=0)  
for j in range(len(nu)):
    ax.plot(t_yr_window,-residual2[j]/nu[j]*yr*1e3,label=rf'$B_p=${Bp_inferred[j]:.1f}$\times 10^{{13}}$ G, $\sigma_{{TN,2}}=${sigma_TN2[j]:.2f} ms')
ax.set_xlabel(r'$t$ (yr)',fontsize=14)
ax.set_ylabel(u'$\mathcal{R}_{t,2}$ (ms)',fontsize=14)
ax.legend(loc='best',ncols=1)
plt.savefig('Residuals2Torque.pdf', bbox_inches='tight')
plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
plt.margins(x=0,y=0)  
for j in range(len(nu)):
    ax.plot(t_yr_window,-residual3[j]/nu[j]*yr*1e3,label=rf'$B_p=${Bp_inferred[j]:.1f}$\times 10^{{13}}$ G, $\sigma_{{TN,2}}=${sigma_TN3[j]:.2f} ms')
ax.set_xlabel(r'$t$ (yr)',fontsize=14)
ax.set_ylabel(u'$\mathcal{R}_{t,3}$ (ms)',fontsize=14)
ax.legend(loc='best',ncols=1)
plt.savefig('Residuals3Torque.pdf', bbox_inches='tight')
plt.show()
plt.close()

names = ['B1642-03','B1818-04']

np.savetxt(OutputFolder+'_Residualst2_'+names[0]+'params.csv', np.column_stack((t_yr_window, -residual2[0]/nu[0]*yr*1e3)), delimiter=',')
np.savetxt(OutputFolder+'_Residualst2_'+names[1]+'params.csv', np.column_stack((t_yr_window, -residual2[1]/nu[1]*yr*1e3)), delimiter=',')

np.savetxt(OutputFolder+'_Residualst3_'+names[0]+'params.csv', np.column_stack((t_yr_window, -residual3[0]/nu[0]*yr*1e3)), delimiter=',')
np.savetxt(OutputFolder+'_Residualst3_'+names[1]+'params.csv', np.column_stack((t_yr_window, -residual3[1]/nu[1]*yr*1e3)), delimiter=',')

# fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
# plt.margins(x=0,y=0)  
# ax.plot(t_yr_window, alphadot_t/yr,label=r'$\dot{\alpha}$')
# ax.set_xlabel(r'$t$ (yr)',fontsize=14)
# ax.set_ylabel(r'$\dot{\alpha}$ (s$^{-1}$)',fontsize=14)
# plt.legend(loc='best')
# plt.show()

# fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
# plt.margins(x=0,y=0)  
# ax.plot(t_yr_window, nudotdot_t[0]/nudotdotDipole[0],label=r'$\dot{\alpha}$')
# ax.set_xlabel(r'$t$ (yr)',fontsize=14)
# ax.set_ylabel(r'$\ddot{\nu}/\ddot{\nu}_0$ (s$^{-1}$)',fontsize=14)
# plt.legend(loc='best')
# plt.show()




# fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
# plt.margins(x=0,y=0)  
# for j in range(len(nu)):
#     ax.plot(t*t_0/yr,nudotdot_t[j]/1e-23,label=rf'$\nu=${nu[j]:.2f} s, $\dot{{\nu}}=${nudot[j]:.2e} s$^{{-1}}$')
#     # ax.plot(t*t_0/yr,nudotdot_t[j]/nudotdot_alphadot0[j],label=rf'$\nu=${nu[j]:.2f} s, $\dot{{\nu}}=${nudot[j]:.2e} s$^{{-1}}$')
# ax.set_xlabel(r'$t$ (yr)')
# # ax.set_ylim([-25,20])
# # ax.set_ylabel(r'$\ddot{\nu}/\ddot{\nu}_0$')
# ax.set_ylabel(r'$\ddot{\nu}$ ($10^{-23}$ s$^{-2}$)')
# ax.legend(loc='best',ncols=2)
# plt.savefig('ResidualsTorque.pdf', bbox_inches='tight')
# plt.show()

# phidot_dat[n] = phidot_fit(t_dim[n],nu[j],nudot[j],nudotdot)
# popt, pcov = optimize.curve_fit(phi_fit_constrained, t*t_0, phi_dat, p0=(nu[j],nudot[j],), bounds=((0.5*nu[j],10*nudot[j]),(1.5*nu[j],0.1*nudot[j])))
# popt, pcov = optimize.curve_fit(phi_fit, t_dim, phi_dat, p0=(nu[j],nudot[j],nudotdotN), 
#                                 bounds=((0.5*nu[j],5*nudot[j],-2*nudotdotN),(1.5*nu[j],0.5*nudot[j],2*nudotdotN)))
# popt, pcov = optimize.curve_fit(lambda t_dim, nudotdot: phi_fit_constrained(t_dim, nudotdot, nu[j], nudot[j]), t_dim, phi_dat, p0=(3*nudot[j]**2/nu[j],) )
# popt2, pcov2 = optimize.curve_fit(phidot2_fit, t_dim, phidot_dat, p0=(nu[j],nudot[j]) )

#Old model of emitting region moving

# def phi_fit(t,nu,nudot):
#     return nu*t + 0.5*nudot*t**2

# for j in range(len(nu)):
#     P = 2.*np.pi/( nu[j] + nudot[j]*t*t_0 ) #model of period in s
#     Pprime = np.empty( (len(t)) )
#     phi_dat = np.empty( (len(t)) )
#     for n in range(len(t)):
#         Pprime[n] = P[n]*( 1 + v_eperpMean[n]*c/(nu[j]*Rsintheta) )
#         phi_dat[n] = 2*np.pi*t_0*np.trapezoid(1./Pprime[0:n+1],t[0:n+1])
        
#     popt, pcov = optimize.curve_fit(phi_fit, t*t_0, phi_dat, p0=(nu[j],nudot[j],))
    
#     phi_fit_eval[j] = phi_fit(t*t_0,popt[0],popt[1])
#     residual[j] = phi_dat - phi_fit_eval[j]
    
#     sigma_TN2[j] = np.sqrt( 1/(t[-1]*t_0*nu[j]**2)*np.trapezoid( residual[j]**2, t )*t_0 )*1e3 #standard deviation of the timing noise in ms

# plt.plot(t*t_0/yr,phi_dat, label=r'$\phi$, data')
# plt.plot(t*t_0/yr,phi_fit_eval, label=r'$\phi$, second-order polynomial fit')
# plt.xlabel(r'$t$ (yr)')
# plt.ylabel(r'$\phi$')
# plt.legend(loc='best')
# plt.show()

# fig, ax = plt.subplots(figsize=(7.5,5),dpi=100)
# plt.margins(x=0,y=0)  
# for j in range(len(nu)):
#     # plt.plot(t*t_0/yr,residual[j]/nu[j]*1e3,label=rf'$\nu=${nu[j]:.2f} s$^{{-1}}$, $B_p=${Bp_inferred[j]:.1f}$\times 10^{{13}}$ G, $\sigma_{{TN,2}}=${sigma_TN2[j]:.2f} ms')
#     ax.plot(t*t_0/yr,-residual[j]/nu[j]*1e3,label=rf'$B_p=${Bp_inferred[j]:.1f}$\times 10^{{13}}$ G, $\sigma_{{TN,2}}=${sigma_TN2[j]:.2f} ms')
# ax.set_xlabel(r'$t$ (yr)')
# ax.set_ylim([-25,20])
# ax.set_ylabel(u'$\mathcal{R}_t$ (ms)')
# ax.legend(loc='best',ncols=2)
# plt.savefig('Residuals.pdf', bbox_inches='tight')
# plt.show()
# plt.close()

# fig, ax1 = plt.subplots(figsize=(6,5),dpi=100)
# plt.margins(x=0.005,y=0.005)   
# ax1.plot(t*t_0/yr,v_eperpMax,label=r'$\theta=90^{\circ}$',color='black',linewidth=1.5)
# #ax1.plot(t_data[1],v_eperpMaxSurface_data[1],label=r'$\theta=80^{\circ}$',color='red',linewidth=1.5)
# #ax1.plot(t_data[0],v_eperpMax_data[0],label=r'$\theta=90^{\circ}$',color='red',linewidth=1.5)
# #ax1.plot(t_data[1],v_eperpMax_data[1],label=r'$\theta=60^{\circ}$',color='blue',linewidth=1.5)
# #ax1.plot(t_data[2],v_eperpMax_data[2],label=r'$\theta=45^{\circ}$',color='red',linewidth=1.5)
# ax1.set_xlabel(r'$t$ (years)',fontsize=12)
# ax1.set_ylabel(r'Max($v_{e,\perp}(z=L_z)/c)$',fontsize=12)
# # ax1.set_yscale('log')
# #ax1.set_ylim([0,0.04])
# #ax1.set_xlim([0,50])
# ax1.legend(loc='best',frameon=False,handlelength=3,fontsize=12)
# #ax1.legend(bbox_to_anchor=(1.05, -.1),ncol=2,frameon=False,handlelength=3,fontsize=12)