"""
samImport
======

Import SAMURAI data and return a dictionary of variables.

    samImport

"""

import numpy as np
import xarray as xr
import pandas as pd


def samImport(samFile):
    """
    Import select variables from a SAMURAI output file,
    convert to masked arrays, and return these data 
    within a dictionary.
    
    Parameters
    ----------
    samFile : string
        String specifying the path and filename of the SAMURAI
        output netCDF file.
        
    Returns
    -------
    samDict : dict
        Dictionary containing all imported SAMURAI data.
    """
    
    samData = xr.open_dataset(samFile)

    x = samData.x.data[:]
    y = samData.y.data[:]
    lon = samData.longitude.data[:]
    lat = samData.latitude.data[:]
    alt = samData.altitude.data[:]
    dbz = samData.DBZ.to_masked_array().squeeze()
    u = samData.U.to_masked_array().squeeze()
    v = samData.V.to_masked_array().squeeze()
    w = samData.W.to_masked_array().squeeze()
    vort = samData.VORT.to_masked_array().squeeze()
    time = pd.to_datetime(samData.time.data[0])
    
    
    samDict = {'x': x, 'y': y, 'lon': lon, 'lat': lat, 'alt': alt, 'dbz': dbz,
               'u': u, 'v': v, 'w': w, 'vort': vort, 'time': time,
               'allXR': samData}
    
    return samDict