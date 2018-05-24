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
    
    u_master_mask = ma.getmask(u_master)
    v_master_mask = ma.getmask(v_master)
    w_master_mask = ma.getmask(w_master)
    vort_master_mask = ma.getmask(vort_master)
    
    
    samData_sub1 = xr.open_dataset(samFile_sub1)
    
    u_sub1_mask = ma.getmask(samData_sub1.U.to_masked_array().squeeze())
    v_sub1_mask = ma.getmask(samData_sub1.V.to_masked_array().squeeze())
    w_sub1_mask = ma.getmask(samData_sub1.W.to_masked_array().squeeze())
    vort_sub1_mask = ma.getmask(samData_sub1.VORT.to_masked_array().squeeze())
    
    u_comboMask = np.logical_or(u_master_mask,u_sub1_mask)
    v_comboMask = np.logical_or(v_master_mask,v_sub1_mask)
    w_comboMask = np.logical_or(w_master_mask,w_sub1_mask)
    vort_comboMask = np.logical_or(vort_master_mask,vort_sub1_mask)
    
    if samFile_sub2:
        samData_sub2 = xr.open_dataset(samFile_sub2)
        
        u_sub2_mask = ma.getmask(samData_sub2.U.to_masked_array().squeeze())
        v_sub2_mask = ma.getmask(samData_sub2.V.to_masked_array().squeeze())
        w_sub2_mask = ma.getmask(samData_sub2.W.to_masked_array().squeeze())
        vort_sub2_mask = ma.getmask(samData_sub2.VORT.to_masked_array().squeeze())
        
        u_comboMask = np.logical_or(u_comboMask,u_sub2_mask)
        v_comboMask = np.logical_or(v_comboMask,v_sub2_mask)
        w_comboMask = np.logical_or(w_comboMask,w_sub2_mask)
        vort_comboMask = np.logical_or(vort_comboMask,vort_sub2_mask)
        
        if samFile_sub3:
            samData_sub3 = xr.open_dataset(samFile_sub3)
            
            u_sub3_mask = ma.getmask(samData_sub3.U.to_masked_array().squeeze())
            v_sub3_mask = ma.getmask(samData_sub3.V.to_masked_array().squeeze())
            w_sub3_mask = ma.getmask(samData_sub3.W.to_masked_array().squeeze())
            vort_sub3_mask = ma.getmask(samData_sub3.VORT.to_masked_array().squeeze())
            
            u_comboMask = np.logical_or(u_comboMask,u_sub3_mask)
            v_comboMask = np.logical_or(v_comboMask,v_sub3_mask)
            w_comboMask = np.logical_or(w_comboMask,w_sub3_mask)
            vort_comboMask = np.logical_or(vort_comboMask,vort_sub3_mask)
    
    
    
    # Apply new combined mask to kinematic variables before saving into dictionary
    
    u_final = ma.masked_array(u_master,u_comboMask)
    v_final = ma.masked_array(v_master,v_comboMask)
    w_final = ma.masked_array(w_master,w_comboMask)
    vort_final = ma.masked_array(vort_master,vort_comboMask)
    
    
    samDict = {'x': x, 'y': y, 'lon': lon, 'lat': lat, 'alt': alt, 
               'dbz': dbz,'time': time,
               'u': u_final, 'v': v_final, 'w': w_final, 'vort': vort_final,
               'u_orig': u_master, 'v_orig': v_master, 'w_orig': w_master, 'vort_orig': vort_master,
               'allXR': samData_master}
    
    return samDict