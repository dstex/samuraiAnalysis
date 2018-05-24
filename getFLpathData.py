import warnings
warnings.filterwarnings("ignore",category=FutureWarning)
import xarray as xr
from datetime import datetime as dt
import numpy as np

warnings.filterwarnings('ignore', 'invalid value encountered in less')

def getFLpathData(flFile,pathStrt,pathEnd,crdsOnly=False):
    """
    Return flight-level data between a path start and end time.
    
    Parameters
    ----------
    flFile : string
        Filename of the flight-level file to be loaded.
    pathStrt,pathEnd : datetimes
        Start/end times of the flight path.
    crdsOnly : bool, optional
        If True, only lat/lon will be extracted and returned (for speed).
    
    Returns
    -------
    allData : dict
        Dictionary containing FL data along the path bounded by the path start/end times.
    """
    
    flData = xr.open_dataset(flFile)

    dtNum = flData.get('datetime_FL').data
    lat = flData.get('lat').to_masked_array()
    lon = flData.get('lon').to_masked_array()
    
    if not crdsOnly:
        alt = flData.get('Alt').to_masked_array()
        ws = flData.get('windSpd').to_masked_array()
        wd = flData.get('windDir').to_masked_array()
        u = flData.get('u').to_masked_array()
        v = flData.get('v').to_masked_array()
        uRel = flData.get('u_relative').to_masked_array()
        vRel = flData.get('v_relative').to_masked_array()
        w = flData.get('w_dpj').to_masked_array()
        staticP = flData.get('staticPres').to_masked_array()
        dynamicP = flData.get('dynamicPres').to_masked_array()
        rh = flData.get('RH_hybrid').to_masked_array()
        temp = flData.get('TA').to_masked_array()
        td = flData.get('TD').to_masked_array()
        
    dtStr = np.char.mod('%d',dtNum)
    flDT = np.asarray([dt.strptime(fDate,'%Y%m%d%H%M%S') for fDate in dtStr])
    
    
    strtIx = np.squeeze(np.where(flDT == pathStrt))
    endIx = np.squeeze(np.where(flDT == pathEnd))
    
    dtPath = flDT[strtIx:endIx+1]
    latPath = lat[strtIx:endIx+1]
    lonPath = lon[strtIx:endIx+1]
    
    if not crdsOnly:
        altPath = alt[strtIx:endIx+1]
        wsPath = ws[strtIx:endIx+1]
        wdPath = wd[strtIx:endIx+1]
        uPath = u[strtIx:endIx+1]
        vPath = v[strtIx:endIx+1]
        uRelPath = uRel[strtIx:endIx+1]
        vRelPath = vRel[strtIx:endIx+1]
        wPath = w[strtIx:endIx+1]
        staticPPath = staticP[strtIx:endIx+1]
        dynamicPPath = dynamicP[strtIx:endIx+1]
        rhPath = rh[strtIx:endIx+1]
        tempPath = temp[strtIx:endIx+1]
        tdPath = td[strtIx:endIx+1]
    
        allData = {'datetime': dtPath, 'lat': latPath, 'lon': lonPath, 'alt': altPath,
                   'windSpd': wsPath, 'windDir': wdPath, 'u': uPath, 'v': vPath,
                   'uRel': uRelPath, 'vRel': vRelPath, 'w': wPath, 'staticPres': staticPPath,
                   'dynamicPres': dynamicPPath, 'rh': rhPath, 'temp': tempPath, 'td': tdPath}
    else:
        allData = {'datetime': dtPath, 'lat': latPath, 'lon': lonPath}
        
    return allData