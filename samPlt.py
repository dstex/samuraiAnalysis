"""
samPlt
======

Common plotting routines for SAMURAI data.

    initMapPlot
    initPlot
    plotContour
    plotVec
"""

import pyart
from matplotlib import pyplot as plt
import numpy as np
from datetime import datetime as dt
import cartopy.crs as ccrs
import cartopy.feature as cf

from getVarLims import getVarLims


def initMapPlot(lon,lat,figsize=(10,10),zoom=False,mapBnds=None):
    """
    Initialize a cartopy map plot for SAMURAI data and return 
    the figure, axes, grid, and map projection handles for further 
    customization.

    Parameters
    ----------
    lon,lat : arrays
        1-D arrays containing the longitude and latitude coordinates of the 
        data to be plotted.
    figsize : tuple, optional
        Tuple indicating the size of the figure in inches (width, height).
    zoom : bool, optional
        False [default] to use the lon,lat bounds for the map extent. When True,
        mapBnds list will be used to specify the map extent.
    mapBnds : list, optional
        List specifying map boundaries, formatted as [minLon,maxLon,minLat,maxLat]. 
        Only used if zoom is True.

    Returns
    -------
    fig,ax,grd,proj : handles
        Figure, axes, and grid handles used for further customization
        of the plot, and proj for transforming data to map coordinates.
    """
    
    if zoom:
        minLon = mapBnds[0]
        maxLon = mapBnds[1]
        minLat = mapBnds[2]
        maxLat = mapBnds[3]
    else:
        minLon = np.min(lon)
        maxLon = np.max(lon)
        minLat = np.min(lat)
        maxLat = np.max(lat)
        
    proj = ccrs.PlateCarree()

    states_provinces = cf.NaturalEarthFeature(
            category='cultural',name='admin_1_states_provinces_lines',
            scale='50m',facecolor='none')
    
    fig = plt.figure(figsize=figsize)
    
    ax = fig.add_subplot(111,projection=proj)
    ax.coastlines()
    ax.add_feature(cf.OCEAN)
    ax.add_feature(cf.LAKES)
    ax.add_feature(cf.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray')
    ax.set_extent([minLon,maxLon,minLat,maxLat])

    grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    grd.xlabels_top = False
    grd.ylabels_right = False
    
    return fig,ax,grd,proj


def initPlot(x,y,figsize=(10,10),zoom=False,xlim=None,ylim=None):
    """
    Initialize an X-Y plot for SAMURAI data and 
    return the figure and axes handles for further customization.

    Parameters
    ----------
    x,y : arrays
        1-D arrays containing the x and y coordinates of the data
        to be plotted.
    figsize : tuple, optional
        Tuple indicating the size of the figure in inches (width, height).
    zoom : bool, optional
        False [default] to use the x,y bounds for the figure extent. When True,
        xlim and/or ylim lists will be used to specify the figure extent.
    xlim,ylim : tuples, optional
        Tuples indicating the (min,max) axis extents. Only used if zoom is True.

    Returns
    -------
    fig,ax : handles
        Figure and axes handles used for further customizing the plot.
    """
    
    fig = plt.figure(figsize=figsize)
    
    ax = fig.add_subplot(111)
    
    if zoom:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    return fig,ax


def plotContour(var,pltLev,x,y,crds,pltFlag,figsize=(10,10),zoom=False,
                   mapBnds=None,xlim=None,ylim=None,vLimMthd='default',
                   vLimLevs='all',vLim=None,savePlt=False,savePath=None,
                   runId='',dT=None,noDisp=False):
    """
    This function contours a given variable at some level.

    Parameters
    ----------
    var : array 
        3-D array (ordered as [level,y/lat,x/lon]) containing the SAMURAI data
        to be plotted. Ideally, this array will be masked.
    pltLev : float
        Height level to plot. This value should be a multiple of the height
        resolution used in the SAMURAI analysis (i.e., if dz=0.5 km, pltLev could
        be 0, 0.5, 1, 1.5, ...).
    x,y : arrays
        1-D arrays containing either the x/y or lon/lat coordinates of the plotted data.
        Chosen variables should coincide appropriately with crds argument.
    crds : {'map', 'xy'}
        Coordinate system to use for plot. 'map' requires x,y arguments be given
        as lon,lat values. 'xy' requires x,y arguments be x,y values.
    pltFlag : {'dbz', 'w', 'u', 'v', 'vort'}
        String specifying what variable parameter set to use. This can be expanded
        to include additional SAMURAI output variables as needed. This is currently used
        only because there is no good way to extract the name of the plotted variable
        as a string (which could then be searched for relevant keywords like 'dbz').
    figsize : tuple, optional
        Tuple indicating the size of the figure in inches (width, height).
    zoom : bool, optional
        False [default] to use the x,y bounds for the figure extent. When True,
        xlim and/or ylim lists will be used to specify the figure extent.
    mapBnds : list, optional
        List specifying map boundaries, formatted as [minLon,maxLon,minLat,maxLat]. 
        Only used if zoom is True, and crds is 'map'.
    xlim,ylim : tuples, optional
        Tuples indicating the (min,max) axis extents. Only used if zoom is True, and crds is 'xy'.
    vLimMthd : {'default', 'tight', 'tightM', 'tightE','tightEM'}, optional
        Variable min/max values for contour plots. If 'default' [default], typical values for given variable
        will be used. If 'tight', rounded max and min values over whole array will be used. If 'tightM', min 
        and max will be equal and of opposite sign. If 'tightE', bordering 3 data points are excluded from
        max/min calculation. If 'tightEM', values are mirrored as in 'tightM', with border exclusion as in 
        'tightE'.
    vLimLevs : int [single value, list, or array], optional
        Integer(s) defining which vertical level to consider in the
        determination of the max and min. Default is 'all' which will use
        all vertical levels.
    vLim : tuple of floats, optional
        Supersedes vLimMthd. Provides max and min contouring bounds. Default is None.
    runId : string, optional
        String indicating a specific run of a SAMURAI analysis. This is used in the title of the figure
        as well as in the figure filename if savePlt is True. In the filename, this string will immediately 
        follow the file datetime string, and the string will be modified to include an underscore at the
        beginning, and all whitespaces will be replaced with hyphens. Examples of suggested use: 'PDD6', 
        'P6 S3 kabd', or 'superHiRes Awesomeness'. 
        Defaults to empty string.
    dT : Timestamp, optional
        Single value Timestamp specifying the date/time of the plotted data. Must be readable by the
        python datetime.datetime.strftime function. If None, 'unknownDT' and
        'Unknown Date/Time' will be used in figure filenames and titles.
    
    

    Returns
    -------
    fig,ax[,grd,proj] : handles
        Figure and axes handles used for further customizing the plot.
        Grid and projection handles are only returned if crds argument is 'map'.
    """
    
    
    pltLevIx = int(pltLev*2)
    
    
    if dT is not None:
        dtTitle = dt.strftime(dT,'%m/%d/%Y %H:%M:%S UTC')
    else:
        dtTitle = 'Unknown Date/Time'
        
        
    titleStr = '{}\n{:g} km AGL {}'.format(dtTitle,pltLev,runId)
        
    
    if pltFlag == 'dbz':
        cmap = pyart.graph.cm.NWSRef
        cmapLbl = 'Reflectivity $(dBZ)$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -8
                vmax = 72
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(var,vLimLevs)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(var,vLimLevs,excldFrame=True)
            # Not including mirrored options as likely not used for reflectivity
        
        
    elif pltFlag == 'w':
        cmap = pyart.graph.cm.GrMg16
        cmapLbl = 'Vertical Velocity $(m s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -15
                vmax = 15
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(var,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(var,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True,excldFrame=True)
        
        
    elif pltFlag == 'u':
        cmap = pyart.graph.cm.GrMg16
        cmapLbl = 'U Wind Component $(m s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -50
                vmax = 50
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(var,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(var,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True,excldFrame=True)
        
        
    elif pltFlag == 'v':
        cmap = pyart.graph.cm.GrMg16
        cmapLbl = 'V Wind Component $(m s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -50
                vmax = 50
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(var,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(var,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True,excldFrame=True)
        
        
    elif pltFlag == 'vort':
        cmap = pyart.graph.cm.GrMg16
        cmapLbl = 'Vertical Vorticity $(s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -700
                vmax = 700
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(var,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(var,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True,excldFrame=True)
        
        
    else:
        raise ValueError('Plotting flag does not match known options')
    
    
    
    if crds == 'map':
        fig,ax,grd,proj = initMapPlot(x,y,figsize=figsize,zoom=zoom,mapBnds=mapBnds)
        plt.pcolormesh(x,y,var[pltLevIx,:,:],cmap=cmap,vmin=vmin,vmax=vmax,transform=proj)
        plt.colorbar(fraction=0.04,pad=0.03,aspect=20,label=cmapLbl)
        plt.xlabel('Longitude ($^{\circ}$)')
        plt.ylabel('Latitude ($^{\circ}$)')
        plt.title(titleStr)

        
    elif crds == 'xy':
        fig,ax = initPlot(x,y,figsize=figsize,zoom=zoom,xlim=xlim,ylim=ylim)
        plt.pcolormesh(x,y,var[pltLevIx,:,:],cmap=cmap,vmin=vmin,vmax=vmax)
        ax.grid(which='both',linewidth=0.5,linestyle='--',color='gray',alpha=0.5)
        plt.colorbar(fraction=0.04,pad=0.03,aspect=20,label=cmapLbl)
        plt.xlabel('Distance from Origin (km)')
        plt.ylabel('Distance from Origin (km)')
        plt.title(titleStr)
        
        
    else:
        raise ValueError('Invalid plot coordinate type chosen')
        
    
    if crds is 'map':
        return fig,ax,grd,proj
    else:
        return fig,ax,None,None
        
        
        
def plotVec(x,y,quivU,quivV,crds,proj=None,quivKeySpd=20,quivKeyUnits='m/s',zoom=False,
            quivKeyX=0.92,quivKeyY=1.05,**kwargs):
            
    """
    Plot wind vectors.
    
    Parameters
    ----------
    x,y : arrays
        Arrays containing either the x/y or lon/lat coordinates of the plotted data.
        Chosen variables should coincide appropriately with crds argument.
    quivU : array
        3-D array containing the U-component of the wind.
    quivV : array
        3-D array containing the V-component of the wind.
    crds : {'map', 'xy'}
        Coordinate system to use for plot. 'map' requires x,y arguments be given
        as lon,lat values. 'xy' requires x,y arguments be x,y values.
    quivKeySpd : float, optional
        Reference value to use when plotting the wind vector key.
    quivKeyUnits : string, optional
        String specifying the units of the plotted wind vectors.
    zoom : bool, optional
        Currently not implemented.
    quivKeyX/Y : float, optional
        Axes relative X/Y coordinates where quiver key will be located.
    **kwargs : optional
        Keyword arguments accepted by the plt.quiver() function are valid.
        
    """
            
    if crds == 'xy':
        quiv = plt.quiver(x,y,quivU,quivV,scale=60,scale_units='inches',**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    elif crds == 'map':
        quiv = plt.quiver(x,y,quivU,quivV,scale=60,scale_units='inches',transform=proj,**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    
    else:
        raise ValueError('Invalid plot coordinate type chosen')