"""
samImport_masked
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


def samImport_masked(samFile_master,samFile_sub1,samFile_sub2=None,samFile_sub3=None):
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

    samData_sub1 = xr.open_dataset(samFile_sub1)
    sub1_mask = ma.getmask(samData_sub1.U.to_masked_array().squeeze())
    
    spiralMask = sub1_mask

    # Currently assumes that sub1 is a PDD (or other dual-Doppler) analysis, sufficient for
    # valid winds on its own. sub2 and sub3, at least for now, are often 88Ds, and should only
    # contribute to the final wind field (where not already defined by sub1) if they are 
    # overlapping each other
    if samFile_sub2 and samFile_sub3:
        samData_sub2 = xr.open_dataset(samFile_sub2)
        samData_sub3 = xr.open_dataset(samFile_sub3)
        
        sub2_mask = ma.getmask(samData_sub2.U.to_masked_array().squeeze())
        sub3_mask = ma.getmask(samData_sub3.U.to_masked_array().squeeze())
        
        mask88d = np.logical_or(sub2_mask,sub3_mask)
    
        comboMask = np.logical_and(spiralMask,mask88d)
    else:
        comboMask = spiralMask
    
    # Apply new combined mask to kinematic variables before saving into dictionary
    
    u_final = ma.masked_array(u_master,comboMask)
    v_final = ma.masked_array(v_master,comboMask)
    w_final = ma.masked_array(w_master,comboMask)
    vort_final = ma.masked_array(vort_master,comboMask)
    
    
    samDict = {'x': x, 'y': y, 'lon': lon, 'lat': lat, 'alt': alt, 
               'dbz': dbz,'time': time,
               'u': u_final, 'v': v_final, 'w': w_final, 'vort': vort_final,
               'u_orig': u_master, 'v_orig': v_master, 'w_orig': w_master, 'vort_orig': vort_master,
               'allXR': samData_master}
    
    return samDict