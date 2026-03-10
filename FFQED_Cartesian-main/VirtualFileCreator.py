#
#				          VirtualFileCreator.py
#			
#		Combines the h5 data files produced by each processor into
#       a single virtual data file
#
#####################################################################

import re
import numpy as np
import math as mt
import pylab as pyl
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.tri as tri
from matplotlib import figure
from matplotlib.pyplot import *
from matplotlib import rc

import glob
import os.path

import h5py

OutputFolder = 'EMHD_Sim_Data_LLB'

################### Combining parallel-written H5 files into single virtual H5 file ########################

def concatenateVecField(file_names_to_concatenate,key,openkey):
    # entry_key = 'u'  # where the data is inside of the source files.
    shT = h5py.File(file_names_to_concatenate[0],'r')[key].shape[0]
    shComps = h5py.File(file_names_to_concatenate[0],'r')[key].shape[1]
    shX = h5py.File(file_names_to_concatenate[0],'r')[key].shape[2]
    shY = []
    for filename in file_names_to_concatenate:
        shY.append( h5py.File(filename,'r')[key].shape[3] ) # get the first one's shape.
    shYTot = np.sum(shY,axis=0)

    layout = h5py.VirtualLayout(shape=(shT,shComps,shX,shYTot,), dtype=np.float64)
    with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', openkey, libver='latest') as f:
        minIndex = 0
        for i, filename in enumerate(file_names_to_concatenate):
            vsource = h5py.VirtualSource(filename, key, shape=(shT,shComps,shX,shY[i],))
            # layout[i, :, :] = vsource
            layout[:,:,:,minIndex:minIndex+shY[i]] = vsource
            minIndex += shY[i]

        f.create_virtual_dataset(key, layout, fillvalue=0)
        f.close()

# Adds coordinate 'key' to virtual data set. Assumes coordinate is not domain-decomposed.
def concatenateCoords(file_names_to_concatenate,key,openkey):
    # entry_key = 'u'  # where the data is inside of the source files.
    shT = h5py.File(file_names_to_concatenate[0],'r')[key].shape[0]
    shX = h5py.File(file_names_to_concatenate[0],'r')[key].shape[1]

    layout = h5py.VirtualLayout(shape=(shT,shX,), dtype=np.float64)
    with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', openkey, libver='latest') as f:
        vsource = h5py.VirtualSource(file_names_to_concatenate[0], key, shape=(shT,shX,))
        # layout[i, :, :] = vsource
        layout[:,0:shX] = vsource

        f.create_virtual_dataset(key, layout, fillvalue=0)
        f.close()

# Adds coordinate 'key' in domain-decomposed direction to virtual data set
def concatenateCoordsDecomp(file_names_to_concatenate,key,openkey):
    # entry_key = 'u'  # where the data is inside of the source files.
    shT = h5py.File(file_names_to_concatenate[0],'r')[key].shape[0]
    shX = []
    for filename in file_names_to_concatenate:
        shX.append( h5py.File(filename,'r')[key].shape[1] ) # get the first one's shape.
    shXTot = np.sum(shX)

    layout = h5py.VirtualLayout(shape=(shT,shXTot,), dtype=np.float64)
    with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', openkey, libver='latest') as f:
        minIndex = 0
        for i, filename in enumerate(file_names_to_concatenate):
            vsource = h5py.VirtualSource(filename, key, shape=(shT,shX[i],))
            # layout[i, :, :] = vsource
            layout[:,minIndex:minIndex+shX[i]] = vsource
            minIndex += shX[i]

        f.create_virtual_dataset(key, layout, fillvalue=0)
        f.close()

# Adds time series to virtual data set
def add_t(file_names_to_concatenate,key,openkey):
    sh = h5py.File(file_names_to_concatenate[0], 'r')[key].shape  # get the first one's shape.
    layout = h5py.VirtualLayout((1,) + sh, dtype=np.float64)
    vsource = h5py.VirtualSource(file_names_to_concatenate[0], key, shape=sh)
    layout[0] = vsource
    with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', openkey, libver='latest') as f:
        f.create_virtual_dataset(key, layout, fillvalue=0)
        f.close()

#Adds attributes to virtual data set
def add_attrs(file_name,openkey):
    with h5py.File(file_name, 'r', libver='latest') as orig:
        with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', openkey, libver='latest') as f:
            for key, value in orig.attrs.items():
                f.attrs[key] = value
                
# Adds energy conservation date to virtual data set
def add_EC(file_names_to_concatenate,keys,openkey):
    sh = h5py.File(file_names_to_concatenate[0], 'r')[keys[0]].shape  # get the first one's shape.
    with h5py.File(OutputFolder+'/'+OutputFolder+'.h5', openkey, libver='latest') as f:
        for i in range(len(keys)): #adds each entry in keys as an entry to virtual data set
            layout = h5py.VirtualLayout((1,) + sh, dtype=np.float64)
            vsource = h5py.VirtualSource(file_names_to_concatenate[0], keys[i], shape=sh)
            layout[0] = vsource
            f.create_virtual_dataset(keys[i], layout, fillvalue=0)
        f.close()

# If virtual H5 file already exists, do not re-create it. Will need to delete it if new data is generated.
# if os.path.isfile(OutputFolder+'/'+OutputFolder+'.h5') == False:

files = []
Nfiles = len(glob.glob(OutputFolder+'/'+OutputFolder+'/'+'*.h5'))
for n in range(0, Nfiles):
    files.append(OutputFolder+'/'+OutputFolder+'/'+OutputFolder+'_{}.h5'.format(n))
    
concatenateVecField(files,'B','w')
concatenateVecField(files,'H','a')
concatenateVecField(files,'D','a')
concatenateVecField(files,'E','a')
concatenateCoords(files,'x','a')
concatenateCoordsDecomp(files,'y','a')
add_t(files,'t','a')
add_attrs(files[0],'a')
add_EC(files,['U_B','JH','PF','DeltaEInt'],'a')