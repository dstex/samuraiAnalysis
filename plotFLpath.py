from matplotlib import pyplot as plt
import numpy as np

def plotFLpath(x,y,crds,proj=None,ax=None,fig=None,
                dubLine=False,**kwargs):
    """
    Plot flight path on top of SAMURAI data (only when plotted atop a map).
    
    Parameters
    ----------
    x,y : 1-D float arrays
        x/y coordinates (if crds = 'xy') or longitudes/latitudes (degrees) of aircraft between given times.
    crds : string
        Options are 'map', which uses the native lat/lon path values and will transform them
        according to the provided projection; 'xy', which will plot plane coordinates provided
        in cartesian coordinates.
    proj : cartopy projection handle, required when 'crds' = 'map'
        This parameter is a returned value from the samPlt.plotContour function, and specifies
        map projection parameters to be used when transforming the FL lat/lon.
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
    
    if crds == 'map':
        if proj is not None:
            if dubLine:
                ax.plot(x,y,transform=proj,color='k',linewidth=3,linestyle='-')
                ax.plot(x,y,transform=proj,color='w',linewidth=1.5,linestyle='-')
            else:
                ax.plot(x,y,transform=proj,**kwargs)
        else:
            print('No map projection was provided for transforming FL path... Skipping FL plot.')
    
    if crds == 'xy':
        if dubLine:
            ax.plot(x,y,color='k',linewidth=3,linestyle='-')
            ax.plot(x,y,color='w',linewidth=1.5,linestyle='-')
        else:
            ax.plot(x,y,**kwargs)