"""
gribTools
======

Import GRIB2 data (must be converted to netCDF) and return a dictionary of variables.

    gribImport
    gribLevs

"""

import numpy as np
import warnings
warnings.filterwarnings("ignore",category=FutureWarning)
import xarray as xr
from datetime import datetime as dt



def gribImport(gribFile):
    """
    Import select variables from a GRIB2 (netCDF) file,
    convert to masked arrays, and return these data 
    within a dictionary.
    
    Parameters
    ----------
    gribFile : string
        String specifying the path and filename of the GRIB2
        netCDF file.
        
    Returns
    -------
    gribDict : dict
        Dictionary containing all imported GRIB2 data.
    """
    
    gribData = xr.open_dataset(gribFile)
    
    
    u = gribData.UGRD_P0_L100_GLC0.to_masked_array().squeeze()
    v = gribData.VGRD_P0_L100_GLC0.to_masked_array().squeeze()
    lat = gribData.gridlat_0.to_masked_array().squeeze()
    lon = gribData.gridlon_0.to_masked_array().squeeze()
    geoHght = gribData.HGT_P0_L100_GLC0.to_masked_array().squeeze()
    presLev = gribData.lv_ISBL0.to_masked_array()/100
    time = dt.strptime(gribData.grib_source[8:21],'%Y%m%d_%H%M')
    
    gribDict = {'lon': lon, 'lat': lat, 'geoHght': geoHght, 'presLev': presLev,
               'u': u, 'v': v, 'time': time,'allData': gribData}
    
    return gribDict
    
    
def gribLevs(levs,geoHght):
    """
    Search through the geopotential height variable from model data
    to find the indices pertaining to the closest match to the levels
    which were given. Geopotential heights for each vertical level will
    be given as an average of all heights at that level.
    
    Parameters
    ----------
    levs : array-like
        List or 1-D array of levels (km) to be compared against.
    geoHght : array
        3-D array of geopotential heights (m) to search through.
        These are converted to km.
        
    Returns
    -------
    geoHtIx : int array
        1-D array containing vertical level indices which most closely
        match the input array of SAMURAI levels to be plotted.
    """
    
    geoHght = geoHght/1000
   
    geoHtAvg = np.empty([geoHght.shape[0]])
    for ix in range(0,geoHght.shape[0]):
        geoHtAvg[ix] = np.mean(geoHght[ix,:,:])
        
    geoHtIx = np.zeros([len(levs),],dtype=int)
    for iL in range(0,len(levs)):
        levMatch = min(geoHtAvg, key=lambda x: abs(x - levs[iL]))
        geoHtIx[iL] = np.int(np.squeeze(np.where(geoHtAvg == levMatch)))
    
    return geoHtIx