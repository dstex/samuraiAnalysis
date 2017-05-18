import numpy as np

def getVarLims(var,levs='all',excldFrame=False,rnd=True,mirror=False):
    """
    Return the max and min values of a given variable. Useful
    for determing the contour limits needed for comparing multiple
    levels of data.

    Parameters
    ----------
    var : array
        3-D array containing the variable in question
    levs : int [single value, list, or array], optional
        Integer(s) defining which vertical level to consider in the
        determination of the max and min. Default is 'all' which will use
        all vertical levels.
    excldFrame : bool, optional
        If True, input variable is sliced to exlude a frame of 3 data points
        from the max/min determination. Useful in many instances (e.g., vorticity)
        where boundaries are poorly constrained leading to unrealistic values.
    rnd : bool, optional
        If True [default], results will be np.floor(varMin) and np.ceil(varMax).
    mirror : bool, optional
        If True, the max absolute value between varMin and varMax will be used for
        both varMin and varMax, resulting in contouring limits centered at 0. This
        is very useful when using a colormap with white in the middle (i.e., for
        vertical velocity).
    

    Returns
    -------
    min,max : floats
        Masked min/max values converted to floats (default return type of np.ma.max[min]imum()
        is a 1-D array).
    """

    if levs is not 'all':
        var = var[levs,:,:]
    
    if excldFrame:
        var = var[:,2:-3,2:-3]
        
    if rnd:
        varMax = np.ceil(np.float(np.ma.maximum(var).data))
        varMin = np.floor(np.float(np.ma.minimum(var).data))
    else:
        varMax = np.float(np.ma.maximum(var).data)
        varMin = np.float(np.ma.minimum(var).data)
        
    if mirror:
        varMin = varMax = max(abs(varMin),abs(varMax))
        varMin = -1*varMin

    
    return varMin,varMax