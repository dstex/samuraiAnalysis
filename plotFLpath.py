from matplotlib import pyplot as plt

def plotFLpath(proj,lon,lat,ax=None,fig=None,
                dubLine=False,**kwargs):
    """
    Plot flight path on top of SAMURAI data (only when plotted atop a map).
    
    Parameters
    ----------
    proj : cartopy projection handle
        This parameter is a returned value from the contourSamurai function, and specifies
        map projection parameters to be used when transforming the FL lat/lon.
    lon,lat : 1-D float arrays
        Latitudes and longitudes (degrees) of aircraft between given times.
    ax : axis handle, optional
        Defaults to the current axis.
    fig : figure handle, optional
        Defaults to the current figure.
    dubLine : bool, optional
        If True, two lines will be plotted, one think and black,
        and the other thin and white. This improves visibility against
        most background fields. If False [default], a single line is plotted
        using arguments supplied by the user.
    **kwargs : optional
        Keyword arguments accepted by the plt.plot() function are valid.
        
    """    
    
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
        
    if dubLine:
        ax.plot(lon,lat,transform=proj,color='k',linewidth=3,linestyle='-')
        ax.plot(lon,lat,transform=proj,color='w',linewidth=1.5,linestyle='-')
    else:
        ax.plot(lon,lat,transform=proj,**kwargs)