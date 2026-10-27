import numpy as np
import os

import xarray as xr

from numpy import zeros
import h5py

def initializedOutput(dir,pref):
    from glob import glob
    file_pattern = f"{pref}_G*timestep.nc"  # eg '*_G*_day.nc'
    timesteplist = glob((dir+file_pattern))
    timesteplist.sort()   # not needed, perhaps

    return timesteplist
    
    


def tablegen(name,filelist, mode='LAYERS'):

    grus = []
    if True:
        negGRU=0
        for idx,fn in enumerate(timesteplist):
                
            print(f"{idx}/{len(timesteplist)}: Processing {fn}")

            nc=xr.open_dataset(fn)
            nTime = len(nc['time'])
            nGRU  = len(nc['gru'])
            try:
                dur=nc['time'][-1]-nc['time'][0]
            except:
                continue
    
            numFluxCalls=nc["numFluxCalls"]
            if mode=='LAYERS':
                nSnow = nc["nSnow"]
                veg = nc["nVegState"]

            shp=np.shape(numFluxCalls)
            negFlag = np.sum(numFluxCalls,axis=0)<0
            negGRU+=np.sum(negFlag)
            numFluxCalls[:,negFlag]=0
            if mode=='LAYERS':
                nSnow[:,negFlag] = 0
                veg[:,negFlag] = 0

            numFluxCalls = np.sum(numFluxCalls,axis=0)
            if mode=='LAYERS':
                nSnow = np.sum(nSnow,axis=0)
                veg = np.sum(veg,axis=0)

            gruID = nc['gru']
            for ix in range(len(numFluxCalls)):
                if mode=='LAYERS':
                    grus.append((veg[ix],nSnow[ix],numFluxCalls[ix],gruID[ix]))
                else:
                    grus.append(numFluxCalls[ix],gruID[ix])
        
    grus = sorted(grus)
    grus = [x[-1] for x in grus]
    f = open(name,'w')
    for gru in grus:
        f.write(str(gru)+'\n')

    return grus

import sys
fnames = initializedOutput(sys.argv[1], sys.argv[2])
tablegen(sys.argv[3])
