import xarray as xr
from datetime import datetime as dt
import numpy as np

def getFLpathData(flFile,pathStrt,pathEnd,crdsOnly=False):
    """
    Return flight-level data between a path start and end time. These data
    can then be overlaid on SAMURAI analysis plots if desired.
    
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
    flLat,flLon : 1-D float arrays
        Latitudes and longitudes (degrees) of aircraft between given times.
    flWS : 1-D float array
        Flight-level wind speed (m/s) between given times.
    flWD : 1-D float array
        Flight-level wind direction (degrees from north) between given times.
    flDTPath : 1-D datetime array
        Datetimes between the given times.
    """
    
    flData = xr.open_dataset(flFile,decode_times=False)

    flLat = flData.get('LatGPS.1').to_masked_array()
    flLon = flData.get('LonGPS.1').to_masked_array()
    flHH = (flData.HH.data[:].astype(int)).astype(str)
    flMM = (flData.MM.data[:].astype(int)).astype(str)
    flSS = (flData.SS.data[:].astype(int)).astype(str)
    if not crdsOnly:
        flWS = flData.get('WS.d').to_masked_array()
        flWD = flData.get('WD.d').to_masked_array()
        
    
    
    flDateStr = np.empty(np.shape(flHH),dtype=object)
    for ix in range(0,len(flHH)):
        flDateStr[ix] = dt.strftime(pathStrt, '%Y%m%d') + '-' + flHH[ix] + ':' + flMM[ix] + ':' + flSS[ix]
    flDT = np.asarray([dt.strptime(fDate,'%Y%m%d-%H:%M:%S') for fDate in flDateStr])
    
    flStrtIx = np.squeeze(np.where(flDT == pathStrt))
    flEndIx = np.squeeze(np.where(flDT == pathEnd))
    
    flLatPath = flLat[flStrtIx:flEndIx+1]
    flLonPath = flLon[flStrtIx:flEndIx+1]
    if not crdsOnly:
        flDTPath = flDT[flStrtIx:flEndIx+1]
        flWSPath = flWS[flStrtIx:flEndIx+1]
        flWDPath = flWD[flStrtIx:flEndIx+1]
        return flLatPath,flLonPath,flWSPath,flWDPath,flDTPath
    else:
        return flLatPath,flLonPath
    
    