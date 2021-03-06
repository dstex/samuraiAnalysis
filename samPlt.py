"""
samPlt
======

Common plotting routines for SAMURAI data.

    initMapPlot
    plotContour
    plotVec
    plotXS
    xsKDtree
    xsCrdCald
    multXSCrdCalc
"""

import pyart
import cmocean
import warnings
import matplotlib as mpl
# mpl.use('PDF')
from matplotlib import pyplot as plt
import numpy as np
import numpy.ma as ma
from datetime import datetime as dt
import cartopy.crs as ccrs
import cartopy.feature as cf
from scipy.spatial import cKDTree
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from pyproj import Proj as Proj


from samuraiAnalysis import getVarLims

getVarLims = getVarLims.getVarLims

warnings.filterwarnings('ignore', 'invalid value encountered in less')

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
    pltFlag : {'dbz', 'w', 'u', 'v', 'vort', 'wind'}
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
        cmap = pyart.graph.cm_colorblind.HomeyerRainbow
        cmapLbl = 'Reflectivity $(dBZ)$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -4
                vmax = 60
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(var,vLimLevs)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(var,vLimLevs,excldFrame=True)
            # Not including mirrored options as likely not used for reflectivity
            
        # steps = np.abs(vmin)+np.abs(vmax)+1
        steps = 33
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'w':
        cmap = cmocean.cm.delta
        cmapLbl = 'Vertical Velocity $(m\ s^{-1})$'
        
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
                
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps*5)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
             
    elif pltFlag == 'u':
        cmap = cmocean.cm.delta
        if strmRel:
            cmapLbl = 'Storm-Relative U Wind Component $(m\ s^{-1})$'
        else:
            cmapLbl = 'U Wind Component $(m\ s^{-1})$'
        
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
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'v':
        cmap = cmocean.cm.delta
        if strmRel:
            cmapLbl = 'Storm-Relative V Wind Component $(m\ s^{-1})$'
        else:
            cmapLbl = 'V Wind Component $(m\ s^{-1})$'
        
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
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
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
                
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                
    elif pltFlag == 'wind':
        cmap = pyart.graph.cm_colorblind.HomeyerRainbow
        cmapLbl = 'Wind Speed $(m\ s^{-1})$'
        
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
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    else:
        raise ValueError('Plotting flag does not match known options')
    
    
    if crds == 'map':
        fig,ax,grd,proj = initMapPlot(x,y,figsize=figsize,zoom=zoom,mapBnds=mapBnds,NB=NB)
        plt.pcolormesh(x,y,var[pltLevIx,:,:],cmap=cmap,norm=norm,vmin=vmin,vmax=vmax,transform=proj)
        if not NB:
            cb = plt.colorbar(shrink=0.7, pad = 0.01, aspect=25)
            cb.set_label(cmapLbl,size=15)
            cb.ax.tick_params(labelsize=14)
            plt.xlabel('Longitude ($^{\circ}$)',size=14)
            plt.ylabel('Latitude ($^{\circ}$)',size=14)
            plt.title(titleStr)

        
    elif crds == 'xy':
        if xlim is None:
            xlim = (np.min(x),np.max(x))
        if ylim is None:
            ylim = (np.min(y),np.max(y))
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        
        plt.pcolormesh(x,y,var[pltLevIx,:,:],cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(which='major',linewidth=0.5,linestyle='--',color='gray',alpha=0.7)
        ax.set_aspect('equal')
        ax.tick_params(labelsize=14)
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        cb = plt.colorbar(shrink=0.7, pad = 0.01, aspect=25)
        cb.set_label(cmapLbl,size=15)
        cb.ax.tick_params(labelsize=14)
        plt.xlabel('Distance from Origin (km)',size=15)
        plt.ylabel('Distance from Origin (km)',size=15)
        plt.title(titleStr,size=18)
        
        
    else:
        raise ValueError('Invalid plot coordinate type chosen')
        
    
    if crds is 'map':
        return fig,ax,grd,proj
    else:
        return fig,ax
        
        
        
def plotVec(x,y,quivU,quivV,crds,proj=None,scale=None,quivKeySpd=40,quivKeyUnits='m/s',zoom=False,
            quivKeyX=0.92,quivKeyY=1.05,scaleVmag=None,**kwargs):
            
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
    proj : map projection handle, required if crds = 'map'
        Map projection handle for tranforming coordinates. Defaults to None.
    scale : integer, optional
        Defaults to None in which case an auto scaling algorithm in quiver scales vectors.
        See quiver documentation for more details.
    quivKeySpd : float, optional
        Reference value to use when plotting the wind vector key.
    quivKeyUnits : string, optional
        String specifying the units of the plotted wind vectors.
    zoom : bool, optional
        Currently not implemented.
    quivKeyX/Y : float, optional
        Axes relative X/Y coordinates where quiver key will be located.
    scaleVmag : float, optional
        Scales the magnitude of the V wind component by this value. Only implemented 
        if crds = 'xs' currently.
    **kwargs : optional
        Keyword arguments accepted by the plt.quiver() function are valid.
        
    Returns
    -------
    quiv,quivKey: handles
        quiver and quiver key handles used for further customizing the plot.
        
    """
            
    if crds == 'xs':
        if scaleVmag is not None:
            quivV *= scaleVmag
        quiv = plt.quiver(x,y,quivU,quivV,scale=scale,scale_units='inches',color='k',width=0.0025,**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    elif crds == 'xy':
        quiv = plt.quiver(x,y,quivU,quivV,scale=scale,scale_units='inches',color='k',width=0.004,**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    elif crds == 'map':
        quiv = plt.quiver(x,y,quivU,quivV,scale=scale,scale_units='inches',transform=proj,color='#3e3e3e',**kwargs)
        quivKey = plt.quiverkey(quiv, quivKeyX, quivKeyY, quivKeySpd, 
                                repr(quivKeySpd) + ' ' + quivKeyUnits, labelpos='S')
    
    else:
        raise ValueError('Invalid plot coordinate type chosen')
        
    return quiv,quivKey
        
        
def plotBarbs(x,y,barbU,barbV,crds,length=6.25,proj=None,scaleVmag=None,**kwargs):
            
    """
    Plot wind barbs.
    
    Parameters
    ----------
    x,y : arrays
        Arrays containing either the x/y, lon/lat, or x/z coordinates of the plotted data.
        Chosen variables should coincide appropriately with crds argument.
    barbU : array
        3-D array containing the U-component of the wind.
    barbV : array
        3-D array containing the V-component of the wind.
    crds : {'map', 'xy', 'xs'}
        Coordinate system to use for plot. 'map' requires x,y arguments be given
        as lon,lat values. 'xy' requires x,y arguments be x,y values. 'xs' requires
        x,y arguments be x,z values.
    length : float, optional
        Length of wind barbs. Stock default is 7; using 6.25 as default here as it is less crowded.
    proj : map projection handle, required if crds = 'map'
        Map projection handle for tranforming coordinates. Defaults to None.
    scaleVmag : float, optional
        Only implemented if crds = 'xs'. The 'barbs' function in matplotlib requires a special modification
        for this to work. "Scales" the wind barb such that the staff angle changes appropriately, but the 
        magnitude indicated by the barbs/flags is the combo of the true U and W values.
    **kwargs : optional
        Keyword arguments accepted by PolyCollection are valid.
        
    """
            
    if crds == 'xs':
        if scaleVmag is not None:
            brbs = plt.barbs(x,y,barbU,barbV*scaleVmag,length=length,barbcolor='#3e3e3e',scaleVmag=scaleVmag,**kwargs)
        else:
            brbs = plt.barbs(x,y,barbU,barbV,length=length,barbcolor='#3e3e3e',**kwargs)
        
    elif crds == 'xy':
        brbs = plt.barbs(x,y,barbU,barbV,length=length,barbcolor='k',**kwargs)
        
    elif crds == 'map':
        brbs = plt.barbs(x,y,barbU,barbV,transform=proj,length=length,color='#3e3e3e',**kwargs)
    
    else:
        raise ValueError('Invalid plot coordinate type chosen')
        
        
def plotStream(x,y,U,V,crds,proj=None,strmRel=False,**kwargs):
            
    """
    Plot streamlines.
    
    Parameters
    ----------
    x,y : arrays
        Arrays containing either the x/y, lon/lat, or x/z coordinates of the plotted data.
        Chosen variables should coincide appropriately with crds argument.
    U : array
        3-D array containing the U-component of the wind.
    V : array
        3-D array containing the V-component of the wind.
    crds : {'map', 'xy', 'xs'}
        Coordinate system to use for plot. 'map' requires x,y arguments be given
        as lon,lat values. 'xy' requires x,y arguments be x,y values. 'xs' requires
        x,y arguments be x,z values.
    proj : map projection handle, required if crds = 'map'
        Map projection handle for tranforming coordinates. Defaults to None.
    strmRel : bool, optional
        If True, modifies the streamline linewidth scaling factor to emphasize the often
        weaker SR winds.
    **kwargs : optional
        Keyword arguments accepted by PolyCollection are valid.
        
    """
    speed = np.sqrt(U*U + V*V)
    if strmRel:
        lw = 4*speed / speed.max()
    else:
        lw = 3*speed / speed.max()
    
    if crds == 'xs':
        strms = plt.streamplot(x,y,U,V,color='k',density=2,arrowsize=1.5,linewidth=lw, **kwargs)
        
    elif crds == 'xy':
        strms = plt.streamplot(x,y,U,V,color='k',density=2,arrowsize=1.5,linewidth=lw, **kwargs)
        
    elif crds == 'map':
        strms = plt.streamplot(x,y,U,V,transform=proj,color='k',density=1,arrowsize=2, **kwargs)
    
    else:
        raise ValueError('Invalid plot coordinate type chosen')
        
        

def plotXS_mapCrds(pltVar3d,lon1d,lat1d,alt1d,u3d,v3d,w3d,xsStrt,xsEnd,pltFlag,
                   xsAngl=None,xsRes=1.0,leafSz=1,xCrd='lat',figsize=(10,5),vecType=None,
                   vLimMthd='default',vLimLevs='all',vLim=None,runId='',dT=None,strmRel=False):
    """
    This function contours a given variable within a cross section.

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
    vecType : {'vec','vecW','barb','stream'}, optional
        Determines what type of wind representation is used. If None (default), no winds are plotted.
        If 'vec', vectors parallel to the cross-section will be plotted (rotated U component and W).
        If 'vecW', only the W component is plotted (0*U and W).
        If 'barb', barbs parallel to the cross-section will be plotted (rotated U component and W).
        If 'stream', streamlines parallel to the cross-section will be plotted (rotated U component and W).
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
        cmap = pyart.graph.cm_colorblind.HomeyerRainbow
        cmapLbl = 'Reflectivity $(dBZ)$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -4
                vmax = 60
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            # Not including mirrored options as likely not used for reflectivity
        
        # steps = np.abs(vmin)+np.abs(vmax)+1
        steps = 33
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
        
    elif pltFlag == 'w':
        # cmap = pyart.graph.cm.GrMg16
        cmap = plt.get_cmap('RdBu_r')
        cmapLbl = 'Vertical Velocity $(m\ s^{-1})$'
        
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
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps*5)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'u':
        # cmap = pyart.graph.cm.GrMg16
        cmap = plt.get_cmap('RdBu_r')
        if strmRel:
            cmapLbl = 'XS $\parallel$ Storm-Relative Wind Component $(m\ s^{-1})$'
        else:
            cmapLbl = 'XS $\parallel$ Wind Component $(m\ s^{-1})$'
        
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
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'v':
        # cmap = pyart.graph.cm.GrMg16
        cmap = plt.get_cmap('RdBu_r')
        if strmRel:
            cmapLbl = 'XS $\perp$ Storm-Relative Wind Component $(m\ s^{-1})$'
        else:
            cmapLbl = 'XS $\perp$ Wind Component $(m\ s^{-1})$'
        
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
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'vort':
        # cmap = pyart.graph.cm.GrMg16
        cmap = plt.get_cmap('RdBu_r')
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
                
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                
    elif pltFlag == 'wind':
        # cmap = pyart.graph.cm.GrMg16
        # cmap = 'gist_rainbow_r'
        cmap = pyart.graph.cm_colorblind.HomeyerRainbow
        cmapLbl = 'Wind Speed $(m\ s^{-1})$'
        
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
                
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
        
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
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    plt.pcolormesh(x,alt1d,xsPltVar,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    ax.grid(which='both',linewidth=0.5,linestyle='--',color='gray',alpha=0.5)
    cb = plt.colorbar(shrink=0.7, pad = 0.01, aspect=25)
    cb.set_label(cmapLbl,size=15)
    cb.ax.tick_params(labelsize=14)
    plt.xlabel(xLab)
    plt.ylabel('Altitude (km AGL)')
    plt.title(titleStr)

    if vecType == 'vec':
        plotVec(x[1::3],alt1d,xsUrot[:,1::3],xsW[:,1::3],crds='xs',quivKeyX=0.97,quivKeyY=1.07)
    if vecType == 'barb':
        plotBarbs(x[1::4],alt1d[1::2],xsUrot[1::2,1::4],xsW[1::2,1::4],crds='xs')
    if vecType == 'stream':
        plotStream(xsX,alt1d,xsUrot,xsW,crds='xs')

        
    return fig,ax


def plotXS(pltVar3d,x1d,y1d,alt1d,u3d,v3d,w3d,xsStrt,xsEnd,pltFlag,
           xsRes=1,figsize=(10,5),vecType=None,xStrd=4,zStrd=2,doPCM=False,
           vLimMthd='default',vLimLevs='all',vLim=None,runId='',xsIX='',dT=None,
           strmRel=False,scaleVmag=None,doUcont=False,doUVoutline=False):
    """
    This function contours a given variable within a cross section.

    Parameters
    ----------
    pltVar3d : array 
        3-D array (ordered as [level,y,x]) containing the SAMURAI data
        to be plotted. Ideally, this array will be masked.
    x1d,y1d,alt1d : arrays
        1-D arrays containing the x/y/alt coordinates of the plotted data.
    u3d,v3d,w3d : arrays
        3-D arrays (ordered as [level,y,x]) containing the SAMURAI u-,v-,w- components 
         of the wind data. Ideally, these arrays will be masked.
    xsStrt : tuple
        Tuple indicating the (x, y) in cartesian coordinates of the cross-section start.
    xsEnd : tuple
        Tuple indicating the (x, y) in cartesian coordinates of the cross-section end.
    pltFlag : {'dbz', 'w', 'u', 'v', 'vort'}
        String specifying what variable parameter set to use. This can be expanded
        to include additional SAMURAI output variables as needed. This is currently used
        only because there is no good way to extract the name of the plotted variable
        as a string (which could then be searched for relevant keywords like 'dbz').
    xsRes : int, optional
        Horizontal resolution of the output x/y arrays. Default is 1, which for the default
        grid spacing in SAMURAI of 1.0 km, this yields one interpolated x/y pair per 1 km. If
        this were to be changed to 10, for example, there would be a point every 100 meters along
        the cross-section.
    figsize : tuple, optional
        Tuple indicating the size of the figure in inches (width, height).
    vecType : {'vec','vecW','barb','stream'}, optional
        Determines what type of wind representation is used. If None (default), no winds are plotted.
        If 'vec', vectors parallel to the cross-section will be plotted (rotated U component and W).
        If 'vecW', only the W component is plotted (0*U and W).
        If 'barb', barbs parallel to the cross-section will be plotted (rotated U component and W).
        If 'stream', streamlines parallel to the cross-section will be plotted (rotated U component and W).
    doPCM : bool, optional
        If True, pcolormesh is used. Otherwise, contourf is used (default).
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
    
    if pltFlag == 'dbz':
        cmap = pyart.graph.cm_colorblind.HomeyerRainbow
        cmapLbl = 'Reflectivity $(dBZ)$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -4
                vmax = 60
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            # Not including mirrored options as likely not used for reflectivity
            
        # steps = np.abs(vmin)+np.abs(vmax)+1
        steps = 33
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'w':
        cmap = cmocean.cm.delta
        cmapLbl = 'Vertical Velocity $(m\ s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -40
                vmax = 40
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps*5)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'u':
        cmap = cmocean.cm.delta
        if strmRel:
            cmapLbl = 'XS $\parallel$ Storm-Relative Wind Component $(m\ s^{-1})$'
        else:
            cmapLbl = 'XS $\parallel$ Wind Component $(m\ s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -50
                vmax = 50
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'v':
        cmap = cmocean.cm.delta
        if strmRel:
            cmapLbl = 'XS $\perp$ Storm-Relative Wind Component $(m\ s^{-1})$'
        else:
            cmapLbl = 'XS $\perp$ Wind Component $(m\ s^{-1})$'
        
        if vLim:
            vmin = vLim[0]
            vmax = vLim[1]
        else:
            if vLimMthd is 'default':
                vmin = -50
                vmax = 50
            elif vLimMthd is 'tight':
                vmin,vmax = getVarLims(pltVar,vLimLevs)
            elif vLimMthd is 'tightM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True)
            elif vLimMthd is 'tightE':
                vmin,vmax = getVarLims(pltVar,vLimLevs,excldFrame=True)
            elif vLimMthd is 'tightEM':
                vmin,vmax = getVarLims(pltVar,vLimLevs,mirror=True,excldFrame=True)
        
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    elif pltFlag == 'vort':
        cmap = plt.get_cmap('RdBu_r')
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
                
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                
    elif pltFlag == 'wind':
        cmap = pyart.graph.cm_colorblind.HomeyerRainbow
        cmapLbl = 'Wind Speed $(m\ s^{-1})$'
        
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
                
        steps = np.abs(vmin)+np.abs(vmax)+1
        bounds = np.linspace(vmin,vmax,steps)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
          
    else:
        raise ValueError('Plotting flag does not match known options')
    
    
    x1,y1 = xsStrt
    x2,y2 = xsEnd
    
    dx = x2-x1
    dy = y2-y1
    xsAngle = np.arctan2(dy,dx)
    
    xsLen_km = np.int(np.round(np.hypot(np.abs(dx),np.abs(dy))))
    xsNumPts = xsLen_km*xsRes
    xsX = np.arange(0,xsNumPts)

    if dT is not None:
        dtTitle = dt.strftime(dT,'%m/%d/%Y %H:%M:%S UTC')
    else:
        dtTitle = 'Unknown Date/Time'
        
    titleStr = '{}\nCross Section between ({},{}) & ({},{}) -- XS-{} -- {}'.format(dtTitle,x1,y1,x2,y2,xsIX,runId)
    
    
    uRot3d = (u3d * np.cos(xsAngle)) + (v3d * np.sin(xsAngle)) # Rotates U || to XS
    vRot3d = (v3d * np.cos(xsAngle)) - (u3d * np.sin(xsAngle))
    
    
    vertLevs = alt1d.shape[0]
    xsPltVar = ma.empty((vertLevs,xsNumPts))*np.nan
    xsUrot = ma.empty((vertLevs,xsNumPts))*np.nan
    xsVrot = ma.empty((vertLevs,xsNumPts))*np.nan
    xsW = ma.empty((vertLevs,xsNumPts))*np.nan
    
    x_XY = np.array([x1,x2])
    y_XY = np.array([y1,y2])
    # Convert from x/y coords to pixel coordinates
    x_pxl = x1d.shape[0] * (x_XY - x1d.min()) / x1d.ptp()
    y_pxl = y1d.shape[0] * (y_XY - y1d.min()) / y1d.ptp()
    x_pxlXS, y_pxlXS = np.linspace(x_pxl[0], x_pxl[1], xsNumPts), np.linspace(y_pxl[0], y_pxl[1], xsNumPts)

    xs_x1d = x1d[x_pxlXS.astype(int)]
    xs_y1d = y1d[y_pxlXS.astype(int)]

    for iz in range(vertLevs):
        xsPltVar[iz,:] = pltVar3d[iz,y_pxlXS.astype(int), x_pxlXS.astype(int)]
        xsUrot[iz,:] = uRot3d[iz,y_pxlXS.astype(int), x_pxlXS.astype(int)]
        xsVrot[iz,:] = vRot3d[iz,y_pxlXS.astype(int), x_pxlXS.astype(int)]
        xsW[iz,:] = w3d[iz,y_pxlXS.astype(int), x_pxlXS.astype(int)]
        
    xStrd = xStrd*xsRes
    
    if pltFlag == 'u':
        xsPltVar = xsUrot
    if pltFlag == 'v':
        xsPltVar = xsVrot

    xPntRes = 1000/xsRes
    
    if xPntRes == 1000:
        xLab = 'Distance from ({},{}) (km)'.format(x1,y1)
    elif xPntRes < 1000:
        xLab = 'Distance from ({},{}) (each pt = {} m)'.format(x1,y1,xPntRes)
    elif xPntRes > 1000:
        xLab = 'Distance from ({},{}) (each pt = {:.2f} km)'.format(x1,y1,xPntRes/1000)
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    if doPCM:
        contF = plt.pcolormesh(xsX,alt1d,xsPltVar,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    else:
        contF = plt.contourf(xsX,alt1d,xsPltVar,bounds,cmap=cmap,norm=norm,extend='both')
    
    if doUcont:
        contU = plt.contour(xsX,alt1d,xsUrot,[20],colors='k',linewidths=3)
        
    if doUVoutline:
        if pltFlag == 'u':
            uOutln = plt.contour(xsX,alt1d,xsUrot,bounds[::5],colors='k')
            if (uOutln.levels[0]-5)%10 != 0.0:
                lblBnds = np.arange(uOutln.levels[0],uOutln.levels[-1],10.)
            else:
                lblBnds = np.arange(uOutln.levels[1],uOutln.levels[-1],10.)
            plt.clabel(uOutln,lblBnds,fmt='%.0f')
            
        if pltFlag == 'v':
            vOutln = plt.contour(xsX,alt1d,xsVrot,bounds[::5],colors='k')
            if (vOutln.levels[0]-5)%10 != 0.0:
                lblBnds = np.arange(vOutln.levels[0],vOutln.levels[-1],10.)
            else:
                lblBnds = np.arange(vOutln.levels[1],vOutln.levels[-1],10.)
            plt.clabel(vOutln,lblBnds,fmt='%.0f')
        if pltFlag == 'wind':
            wndOutln = plt.contour(xsX,alt1d,xsPltVar,bounds[::5],colors='k')
            lblBnds = np.arange(wndOutln.levels[1],wndOutln.levels[-1],5.)
            plt.clabel(wndOutln,lblBnds,fmt='%.0f')
        
        
    ax.grid(which='major',linewidth=0.5,linestyle='--',color='gray',alpha=0.5)
    ax.tick_params(labelsize=14)
    ax.xaxis.set_major_locator(MultipleLocator(20*xsRes))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(5*xsRes))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    cb = plt.colorbar(contF,shrink=0.8, pad = 0.01, aspect=25)
    cb.set_label(cmapLbl,size=15)
    cb.ax.tick_params(labelsize=14)
    plt.xlabel(xLab,size=15)
    plt.ylabel('Altitude (km AGL)',size=15)
    plt.title(titleStr,size=18)

    if vecType == 'vec':
        plotVec(xsX[2::xStrd],alt1d[2::zStrd],xsUrot[2::zStrd,2::xStrd],xsW[2::zStrd,2::xStrd],crds='xs',
                scaleVmag=scaleVmag,quivKeyX=0.97,quivKeyY=1.07,scale=110)
    if vecType == 'vecW':
        xStrd = 3
        zStrd = 3
        plotVec(xsX[2::xStrd],alt1d[2::zStrd],xsUrot[2::zStrd,2::xStrd]*0,xsW[2::zStrd,2::xStrd],crds='xs',
                quivKeyX=0.97,quivKeyY=1.07,quivKeySpd=5,scale=15)
    if vecType == 'barb':
        plotBarbs(xsX[2::xStrd],alt1d[2::zStrd],xsUrot[2::zStrd,2::xStrd],xsW[2::zStrd,2::xStrd],crds='xs',
                  scaleVmag=scaleVmag)
    if vecType == 'stream':
        plotStream(xsX,alt1d,xsUrot,xsW,crds='xs',strmRel=strmRel)
        
    return fig,ax,np.rad2deg(xsAngle),xsLen_km,xsNumPts,xs_x1d,xs_y1d

    
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
    
def map2xy(lon,lat,lon_0,lat_0):
    """
    This function converts given map coordinates in degrees lon/lat
    into x/y cartesian coordinates about an origin specified by (lon_0,lat_0).

    Parameters
    ----------
    lon/lat : arrays 
        Arrays containing lon/lat data in degrees.
    lon_0/lat_0 : floats
        Values specifying the origin point (0,0) of the converted coordinates.
        
    Returns
    -------
    x,y : arrays 
        Arrays of same size as input lon/lat containing the converted map coordinates.
    """
    p = Proj(proj='laea',ellps='WGS84',lat_0=lat_0,lon_0=lon_0)
    
    x,y = p(lon,lat)
    x /= 1000
    y /= 1000
    
    return x,y
    
    
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