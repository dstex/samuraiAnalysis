"""
samPlt
======

Common plotting routines for SAMURAI data.

    initMapPlot
    initPlot
    plotContour
    plotVec
    plotXS
    xsKDtree
    xsCrdCald
    multXSCrdCalc
"""

import pyart
from matplotlib import pyplot as plt
import numpy as np
import numpy.ma as ma
from datetime import datetime as dt
import cartopy.crs as ccrs
import cartopy.feature as cf
from scipy.spatial import cKDTree

from samuraiAnalysis import getVarLims

getVarLims = getVarLims.getVarLims

def initMapPlot(lon,lat,figsize=(10,10),zoom=False,mapBnds=None,NB=False):
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
    NB : bool, optional
        If True, create a figure without a frame, and with no axis labels. This is particularly
        useful if planning to make a Google Earth KML file from the output figure.

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
    
    if NB:
        fig = plt.figure(figsize=figsize,frameon=False)
    else:
        fig = plt.figure(figsize=figsize)
    
    ax = fig.add_subplot(111,projection=proj)
    ax.coastlines()
    ax.add_feature(cf.OCEAN)
    ax.add_feature(cf.LAKES)
    ax.add_feature(cf.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray')
    ax.set_extent([minLon,maxLon,minLat,maxLat])
    
    if NB:
        grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                          linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    else:
        grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        grd.xlabels_top = False
        grd.ylabels_right = False
    
    return fig,ax,grd,proj


def initPlot(figsize=(10,10),zoom=False,xlim=None,ylim=None):
    """
    Initialize an X-Y plot for SAMURAI data and 
    return the figure and axes handles for further customization.

    Parameters
    ----------
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
                   vLimLevs='all',vLim=None,runId='',dT=None,NB=False,strmRel=False):
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
        String indicating a specific run of a SAMURAI analysis. This is used in the title of the figure.
        Defaults to empty string.
    dT : Timestamp, optional
        Single value Timestamp specifying the date/time of the plotted data. Must be readable by the
        python datetime.datetime.strftime function. If None, 'unknownDT' and
        'Unknown Date/Time' will be used in figure filenames and titles.
    NB : bool, optional
        If True, create a figure without a frame, and with no axis labels. This is particularly
        useful if planning to make a Google Earth KML file from the output figure. 
    strmRel : bool, optional
        If True, label plotted wind variables appropriately.
    

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
#         cmap = pyart.graph.cm.Gray9
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
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'RdBu_r'
        cmapLbl = 'Vertical Velocity $(m s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -7
                vmax = 7
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(var,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(var,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(var,vLimLevs,mirror=True,excldFrame=True)
        
        
    elif pltFlag == 'u':
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'RdBu_r'
        if strmRel:
            cmapLbl = 'Storm-Relative U Wind Component $(m s^{-1})$'
        else:
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
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'RdBu_r'
        if strmRel:
            cmapLbl = 'Storm-Relative V Wind Component $(m s^{-1})$'
        else:
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
                
    elif pltFlag == 'wind':
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'gist_rainbow_r'
        cmapLbl = 'Wind Speed $(s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = 0
                vmax = 60
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        
    else:
        raise ValueError('Plotting flag does not match known options')
    
    
    
    if crds == 'map':
        fig,ax,grd,proj = initMapPlot(x,y,figsize=figsize,zoom=zoom,mapBnds=mapBnds,NB=NB)
        plt.pcolormesh(x,y,var[pltLevIx,:,:],cmap=cmap,vmin=vmin,vmax=vmax,transform=proj)
        if not NB:
            plt.colorbar(fraction=0.04,pad=0.03,aspect=20,label=cmapLbl)
            plt.xlabel('Longitude ($^{\circ}$)')
            plt.ylabel('Latitude ($^{\circ}$)')
            plt.title(titleStr)

        
    elif crds == 'xy':
        fig,ax = initPlot(figsize=figsize,zoom=zoom,xlim=xlim,ylim=ylim)
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
    crds : {'map', 'xy', 'xs'}
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
            
    if crds == 'xs':
        quiv = plt.quiver(x,y,quivU,quivV,scale=45,scale_units='inches',**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    elif crds == 'xy':
        quiv = plt.quiver(x,y,quivU,quivV,scale=60,scale_units='inches',**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    elif crds == 'map':
        quiv = plt.quiver(x,y,quivU,quivV,scale=60,scale_units='inches',transform=proj,**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    
    else:
        raise ValueError('Invalid plot coordinate type chosen')
        
        

def plotXS(pltVar3d,lon1d,lat1d,alt1d,u3d,v3d,w3d,xsStrt,xsEnd,pltFlag,
                   xsAngl=None,xsRes=1.0,leafSz=1,xCrd='lat',figsize=(10,5),vecPlt=False,
                   vLimMthd='default',vLimLevs='all',vLim=None,runId='',dT=None,strmRel=False):
    """
    This function contours a given variable at some level.

    Parameters
    ----------
    pltVar3d : array 
        3-D array (ordered as [level,lat,lon]) containing the SAMURAI data
        to be plotted. Ideally, this array will be masked.
    lon1d,lat1d,alt1d : arrays
        1-D arrays containing the lon/lat/alt coordinates of the plotted data.
    u3d,v3d,w3d : arrays
        3-D arrays (ordered as [level,lat,lon]) containing the SAMURAI u-,v-,w- components 
         of the wind data. Ideally, these arrays will be masked.
    xsStrt : tuple
        Tuple indicating the (latitude, longitude) in degrees of the cross-section start.
    xsEnd : tuple
        Tuple indicating the (latitude, longitude) in degrees of the cross-section end.
    pltFlag : {'dbz', 'w', 'u', 'v', 'vort'}
        String specifying what variable parameter set to use. This can be expanded
        to include additional SAMURAI output variables as needed. This is currently used
        only because there is no good way to extract the name of the plotted variable
        as a string (which could then be searched for relevant keywords like 'dbz').
    xsAngl : float, optional
        Angle of cross-section (from north) by which winds will be rotated. 
        Winds are rotated to place V parallel to the cross-section. Defaults to
        None, in which case the angle is automatically calculated.
    xsRes : float, optional
        Horizontal resolution of the output lat/lon arrays (in km). Default is 1.0 km.
    leafSz : int, optional
        Integer specifying how many neighbors are to be considered in the query of the
        kDtree.
    xCrd : string, optional
        String indicating which variable to use when labeling the x-axis of the cross-section.
        Default is 'lat'.
    figsize : tuple, optional
        Tuple indicating the size of the figure in inches (width, height).
    vecPlt : bool, optional
        If True, vectors parallel to the cross-section will be plotted (rotated V component and W).
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
        String indicating a specific run of a SAMURAI analysis. This is used in the title of the figure.
        Defaults to empty string.
    dT : Timestamp, optional
        Single value Timestamp specifying the date/time of the plotted data. Must be readable by the
        python datetime.datetime.strftime function. If None, 'unknownDT' and
        'Unknown Date/Time' will be used in figure filenames and titles.
    strmRel : bool, optional
        If True, label plotted wind variables appropriately.
    
    

    Returns
    -------
    fig,ax : handles
        Figure and axes handles used for further customizing the plot.
    """
    
    
    
    if dT is not None:
        dtTitle = dt.strftime(dT,'%m/%d/%Y %H:%M:%S UTC')
    else:
        dtTitle = 'Unknown Date/Time'
        
        
    titleStr = '{}\nVertical Cross Section between ({},{}) & ({},{}) -- {}'.format(dtTitle,xsStrt[0],xsStrt[1],
                                                                                xsEnd[0],xsEnd[1],runId)
    
       
        
    
    
    
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
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            # Not including mirrored options as likely not used for reflectivity
        
        
    elif pltFlag == 'w':
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'RdBu_r'
        cmapLbl = 'Vertical Velocity $(m s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -4
                vmax = 4
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        
    elif pltFlag == 'u':
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'RdBu_r'
        if strmRel:
            cmapLbl = 'XS $\parallel$ Storm-Relative Wind Component $(m s^{-1})$'
        else:
            cmapLbl = 'XS $\parallel$ Wind Component $(m s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -30
                vmax = 30
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        
    elif pltFlag == 'v':
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'RdBu_r'
        if strmRel:
            cmapLbl = 'XS $\perp$ Storm-Relative Wind Component $(m s^{-1})$'
        else:
            cmapLbl = 'XS $\perp$ Wind Component $(m s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -30
                vmax = 30
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        
    elif pltFlag == 'vort':
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'RdBu_r'
        cmapLbl = 'Vertical Vorticity $(s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -400
                vmax = 400
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
                
    elif pltFlag == 'wind':
        # cmap = pyart.graph.cm.GrMg16
        cmap = 'gist_rainbow_r'
        cmapLbl = 'Wind Speed $(s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = 0
                vmax = 60
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        
    else:
        raise ValueError('Plotting flag does not match known options')
    
    
    # Flatten 3-D arrays to allow for proper indexing with KdTree
    shape3d = pltVar3d.shape
    pltVar = pltVar3d.reshape(shape3d[0],-1) # Flatten lat/lon into single dimension
    u = u3d.reshape(shape3d[0],-1)
    v = v3d.reshape(shape3d[0],-1)
    w = w3d.reshape(shape3d[0],-1)
    
    # Tile 1-D coordinate arrays
    alt = np.tile(alt1d,(shape3d[1],shape3d[2],1)).transpose((2,0,1))
    lat = np.ravel(np.tile(lat1d,(shape3d[2],1)).transpose())
    lon = np.ravel(np.tile(lon1d,(shape3d[1],1)))
    
    
    matchIx,xsCrds,dRnd,hdng = xsKDtree(lon,lat,xsStrt,xsEnd,xsRes=xsRes,leafSz=leafSz)
    vertLevs = pltVar.shape[0]
    xsPnts = matchIx.shape[0]

    actualCrds = np.array([lon,lat])
    matchdCrds = actualCrds.transpose()[matchIx.ravel(),:]

    if xsAngl is None:
        xsAngl = hdng
    
    # Calculate the rotated wind components then flatten the 3-D arrays
    # rotAngl = np.deg2rad(90-xsAngl) 
    # uRot3d = (u3d * np.sin(rotAngl)) - (v3d * np.cos(rotAngl)) # Rotates V || to XS - and incorrectly sometimes 
    # vRot3d = (u3d * np.cos(rotAngl)) + (v3d * np.sin(rotAngl))
    rotAngl = np.deg2rad(xsAngl) 
    uRot3d = (u3d * np.cos(rotAngl)) + (v3d * np.sin(rotAngl)) # Rotates U || to XS
    vRot3d = (v3d * np.cos(rotAngl)) - (u3d * np.sin(rotAngl))
    
    uRot = uRot3d.reshape(shape3d[0],-1)
    vRot = vRot3d.reshape(shape3d[0],-1)


    xsPltVar = ma.empty((vertLevs,xsPnts))*np.nan
    xsUrot = ma.empty((vertLevs,xsPnts))*np.nan
    xsVrot = ma.empty((vertLevs,xsPnts))*np.nan
    xsW = ma.empty((vertLevs,xsPnts))*np.nan
    
    if leafSz > 1:
        for iz in range(0,vertLevs):
            for ix in range(0,xsPnts):
                xsPltVar[iz,ix] = ma.mean(pltVar[iz,matchIx[ix,:]])
                xsUrot[iz,ix] = ma.mean(uRot[iz,matchIx[ix,:]])
                xsVrot[iz,ix] = ma.mean(vRot[iz,matchIx[ix,:]])
                xsW[iz,ix] = ma.mean(w[iz,matchIx[ix,:]])
    else:
        for iz in range(0,vertLevs):
            for ix in range(0,xsPnts):
                xsPltVar[iz,ix] = pltVar[iz,matchIx[ix]]
                xsUrot[iz,ix] = uRot[iz,matchIx[ix]]
                xsVrot[iz,ix] = vRot[iz,matchIx[ix]]
                xsW[iz,ix] = w[iz,matchIx[ix]]


    if pltFlag == 'u':
        xsPltVar = xsUrot
    if pltFlag == 'v':
        xsPltVar = xsVrot

    if xCrd == 'lat':
        x = xsCrds[1,:]
        xLab = 'Latitude ($^{\circ}$)'
    elif xCrd == 'lon':
        x = xsCrds[0,:]
        xLab = 'Longitude ($^{\circ}$)'
    
    fig,ax = initPlot(figsize=figsize)
    plt.pcolormesh(x,alt1d,xsPltVar,cmap=cmap,vmin=vmin,vmax=vmax)
    ax.grid(which='both',linewidth=0.5,linestyle='--',color='gray',alpha=0.5)
    plt.colorbar(fraction=0.04,pad=0.03,aspect=20,label=cmapLbl)
    plt.xlabel(xLab)
    plt.ylabel('Altitude (km AGL)')
    plt.title(titleStr)

    if vecPlt:
        plotVec(x[1::3],alt1d,xsUrot[:,1::3],xsW[:,1::3],crds='xs',quivKeyX=0.97,quivKeyY=1.07)
        
    return fig,ax
    
def xsKDtree(lon,lat,xsStrt,xsEnd,xsRes=1.0,leafSz=1):
    """
    This function returns an index array for the (leafSz) nearest neighbors for
    every point between given start/end coordinate pairs.

    Parameters
    ----------
    lon,lat : arrays
        3-D arrays containing the lon/lat coordinates of the data for every data point.
    xsStrt : tuple
        Tuple indicating the (latitude, longitude) of the cross-section start.
    xsEnd : tuple
        Tuple indicating the (latitude, longitude) of the cross-section end.
    xsRes : float, optional
        Horizontal resolution of the output lat/lon arrays (in km). Default is 1.0 km.
    leafSz : int, optional
        Integer specifying how many neighbors are to be considered in the query of the
        kDtree.
     
    

    Returns
    -------
    matchIx : array
        2-D array of dimension (nXS_points, leafSz) containing indices of the nearest
        points along the cross-section.
    xsCrds : array
        2-D array of dimension (2, len(lat[lon])) containing coordinates at a given resolution
        along the cross-section. [0,:] gives lon values, [1,:] gives lat values
    dRnd : float
        Distance (rounded to nearest integer value) in km between the start and end XS points.
    hdng : float
        Heading (in degrees) between the start and end cross-section points.
    """
    
    # Create single array of lon/lat then construct a kDTree
    actualCrds = np.array([lon,lat])
    kd = cKDTree(actualCrds.transpose())
    
    # Get lat/lon values between two points and create
    # another single array with lon,lat pairs
    xsLat,xsLon,dRnd,hdng = xsCrdCalc(xsStrt[0],xsStrt[1],xsEnd[0],xsEnd[1],xsRes)
    xsCrds = np.array([xsLon,xsLat])
    
    # Query the k-D tree to determine which points (indices) within our
    # actual lon,lat array correspond most closely with the XS points
    # k will determine how many neighbors to include in the result
    matchIx = kd.query(xsCrds.transpose(),k=leafSz)[1]
    
    return matchIx,xsCrds,dRnd,hdng
    
    
    
def xsCrdCalc(lat1d,lon1d,lat2d,lon2d,xsRes=1.0):
    """
    Calculates arrays of lat/lon at a given resolution between two lat/lon pairs.
    The Haversine formula is used to calculate the great circle distance between the
    two points, and then the coordinates between the start and end points are populated
    into final lat/lon arrays. The process follows that outlined by Chris Veness at
    http://www.movable-type.co.uk/scripts/latlong.html.
    
    Additionally, the heading between the two points is also calculated. This portion
    of the process is outlined at
    http://www.igismap.com/formula-to-find-bearing-or-heading-angle-between-two-points-latitude-longitude/

    Parameters
    ----------
    lat1d,lon1d : float
        Lon/lat coordinates (degrees) of the starting point of the cross-section.
    lat2d,lon2d : float
        Lon/lat coordinates (degrees) of the ending point of the cross-section.
    xsRes : float, optional
        Horizontal resolution of the output lat/lon arrays (in km). Default is 1.0 km.
     
    

    Returns
    -------
    xsLat,xsLon : array
        1-D arrays of length (ceil(GreatCircleDistance/xsRes)+1) containing lat/lon
        values at horizontal resolution xsRes between the start and end points.
    dRnd : float
        Distance (rounded to nearest integer value) in km between the start and end XS points.
    hdng : float
        Heading (in degrees) between the start and end cross-section points.
    """

    lat1tmp = np.deg2rad(lat1d)
    lon1tmp = np.deg2rad(lon1d)
    lat2tmp = np.deg2rad(lat2d)
    lon2tmp = np.deg2rad(lon2d)
    lon1 = lon1tmp
    lat1 = lat1tmp
    lon2 = lon2tmp
    lat2 = lat2tmp
    
    # Ensure that the initial point of the cross-section is the western-most point, and
    #   in the event the line is N-S oriented, set the initial point to the south.
#     if lon1tmp > lon2tmp:
#         print('XS final/initial points swapped to ensure proper wind rotation')
#         lon1 = lon2tmp
#         lat1 = lat2tmp
#         lon2 = lon1tmp
#         lat2 = lat1tmp
#     elif lon1tmp < lon2tmp:
#         lon1 = lon1tmp
#         lat1 = lat1tmp
#         lon2 = lon2tmp
#         lat2 = lat2tmp
#     elif lon1tmp == lon2tmp:
#         if lat1tmp > lat2tmp:
#             print('XS final/initial points swapped to ensure proper wind rotation')
#             lon1 = lon2tmp
#             lat1 = lat2tmp
#             lon2 = lon1tmp
#             lat2 = lat1tmp
#         else:
#             lon1 = lon1tmp
#             lat1 = lat1tmp
#             lon2 = lon2tmp
#             lat2 = lat2tmp
    
    # Calculate the bearing/heading between the start and end points
    # I believe the western-most point must be given as the "initial"
    #    point for the returned heading to work properly with the wind rotation scheme
    # In the event the two longitudes are the same (XS perfectly N-S), the southern-most
    #    point should be given as the initial point.
    tempY = np.sin(lon2-lon1)*np.cos(lat2)
    tempX = (np.cos(lat1)*np.sin(lat2)) - (np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
    hdng = np.rad2deg(np.arctan2(tempX,tempY))

    if hdng < 0:
        hdng += 360

    a = (np.sin((lat2-lat1)/2))**2 + (np.cos(lat1)*np.cos(lat2)*(np.sin((lon2-lon1)/2))**2) #square of half the chord length between the points
    c = 2*np.arctan2(np.sqrt(a),np.sqrt(1-a)) #angular distance in radians
    d = 6371*c #great circle distance in km

    dRnd = np.round(d)

    xsLat = np.zeros((np.int(np.ceil(dRnd/xsRes))+1,))
    xsLon = np.zeros((np.int(np.ceil(dRnd/xsRes))+1,))
    xsLat[0] = np.rad2deg(lat1)
    xsLon[0] = np.rad2deg(lon1)
    xsLat[-1] = np.rad2deg(lat2)
    xsLon[-1] = np.rad2deg(lon2)

    # Calculate lat/lon values between the start and end points
    aIx = 1
    for ix in np.arange(xsRes,dRnd,xsRes):

        frac = ix/dRnd

        ai = np.sin((1-frac)*c)/np.sin(c)
        bi = np.sin(frac*c)/np.sin(c)
        x = ai*np.cos(lat1)*np.cos(lon1) + bi*np.cos(lat2)*np.cos(lon2)
        y = ai*np.cos(lat1)*np.sin(lon1) + bi*np.cos(lat2)*np.sin(lon2)
        z = ai*np.sin(lat1) + bi*np.sin(lat2)
        xsLat[aIx] = np.rad2deg(np.arctan2(z,np.sqrt(x**2 + y**2)))
        xsLon[aIx] = np.rad2deg(np.arctan2(y,x))
        
        aIx += 1
    
    return xsLat, xsLon, dRnd, hdng
    
def multXSCrdCalc(lat1d,lon1d,distKM,hdng):
    """
    Calculates lat/lon pairs at some distance and heading from a given point.
    This is useful for defining consecutive cross-sections spaced at some given distance.

    Parameters
    ----------
    lat1d,lon1d : float
        Lon/lat coordinates (degrees) of the start/end of reference cross-section.
    distKM : float
        Distance from reference point to obtain new point at (in KM)
    hdng : float
        Heading (relative to north) of new point relative to reference point.
     
    

    Returns
    -------
    newLat,newLon : floats
        Lat/lon (in degrees) of new XS start/end point.
    """

    lat1 = np.deg2rad(lat1d)
    lon1 = np.deg2rad(lon1d)
    
    rEarth = 6371
    hdngRad = np.deg2rad(hdng)
    
    lat2 = np.arcsin( np.sin(lat1)*np.cos(distKM/rEarth) + np.cos(lat1)*np.sin(distKM/rEarth)*np.cos(hdngRad) )

    lon2 = lon1 + np.arctan2( np.sin(hdngRad)*np.sin(distKM/rEarth)*np.cos(lat1), np.cos(distKM/rEarth)-np.sin(lat1)*np.sin(lat2) )
    
    newLat = np.rad2deg(lat2)
    newLon = np.rad2deg(lon2)
    
    return newLat, newLon