"""
samImport_masked_deprecated
======

Import SAMURAI data from 2-4 separate analyses and determine a new mask where all contain valid (unmasked) data. 
This resultant combo mask is applied to the kinematic variables (u,v,w,vort) of the master (complete) analysis.
A dictionary of all final master variables is returned (just as in samImport).

    samImport_masked

"""

import numpy as np
import numpy.ma as ma
import warnings
warnings.filterwarnings("ignore",category=FutureWarning)
import xarray as xr
import pandas as pd

warnings.filterwarnings('ignore', 'invalid value encountered in less')


def samImport_masked_deprecated(samFile_master,samFile_sub1,samFile_sub2=None,samFile_sub3=None):
    """
    Import select variables from a SAMURAI output file,
    convert to masked arrays, and return these data 
    within a dictionary.
    
    Parameters
    ----------
    samFile_master : string
        String specifying the path and filename of the SAMURAI
        output netCDF file containing the master (complete) analysis using data
        from all sources (i.e., NOAA P-3 TDR data + KUDX 88D data + KABR 88D data).
    samFile_sub1 : string
        String specifying the path and filename of the SAMURAI
        output netCDF file containing the analysis which used a subset of the
        data in the master analysis (i.e., just NOAA P-3 TDR data).
    samFile_sub2 : string, optional
        String specifying the path and filename of a third SAMURAI
        output netCDF file containing the analysis which used a subset of the
        data in the master analysis (i.e., just KUDX 88D data).
    samFile_sub3 : string, optional
        String specifying the path and filename of a fourth SAMURAI
        output netCDF file containing the analysis which used a subset of the
        data in the master analysis (i.e., just KABR 88D data).
        
    Returns
    -------
    samDict : dict
        Dictionary containing appropriately masked versions of the master SAMURAI data.
    """
    
    samData_master = xr.open_dataset(samFile_master)

    x = samData_master.x.data[:]
    y = samData_master.y.data[:]
    lon = samData_master.longitude.data[:]
    lat = samData_master.latitude.data[:]
    alt = samData_master.altitude.data[:]
    time = pd.to_datetime(samData_master.time.data[0])
    dbz = samData_master.DBZ.to_masked_array().squeeze()
    
    u_master = samData_master.U.to_masked_array().squeeze()
    v_master = samData_master.V.to_masked_array().squeeze()
    w_master = samData_master.W.to_masked_array().squeeze()
    vort_master = samData_master.VORT.to_masked_array().squeeze()
    
    master_mask = ma.getmask(u_master)
    
    
    samData_sub1 = xr.open_dataset(samFile_sub1)
    sub1_mask = ma.getmask(samData_sub1.U.to_masked_array().squeeze())

    maskCount = np.add(master_mask.astype(int),sub1_mask.astype(int))

    
    if samFile_sub2:
        samData_sub2 = xr.open_dataset(samFile_sub2)
        
        sub2_mask = ma.getmask(samData_sub2.U.to_masked_array().squeeze())
        maskCount = np.add(maskCount,sub2_mask.astype(int))
        
        if samFile_sub3:
            samData_sub3 = xr.open_dataset(samFile_sub3)
            
            sub3_mask = ma.getmask(samData_sub3.U.to_masked_array().squeeze())
            maskCount = np.add(maskCount,sub3_mask.astype(int))
    
    
    finalMask = (np.zeros_like(master_mask)).ravel()
    for ix, maskC in enumerate(maskCount.ravel()):
        if maskC > 2:
            finalMask[ix] = True
    
    finalMask_reshaped = finalMask.reshape(master_mask.shape)
    
    # Apply new combined mask to kinematic variables before saving into dictionary
    
    u_final = ma.masked_array(u_master,finalMask_reshaped)
    v_final = ma.masked_array(v_master,finalMask_reshaped)
    w_final = ma.masked_array(w_master,finalMask_reshaped)
    vort_final = ma.masked_array(vort_master,finalMask_reshaped)
    
    
    samDict = {'x': x, 'y': y, 'lon': lon, 'lat': lat, 'alt': alt, 
               'dbz': dbz,'time': time,
               'u': u_final, 'v': v_final, 'w': w_final, 'vort': vort_final,
               'u_orig': u_master, 'v_orig': v_master, 'w_orig': w_master, 'vort_orig': vort_master,
               'allXR': samData_master}
    
    return samDict