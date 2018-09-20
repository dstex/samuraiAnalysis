import os
import shutil
import warnings
from datetime import datetime as dt
import matplotlib as mpl
# mpl.use('PDF')
from matplotlib import pyplot as plt
import numpy as np
import yaml
import argparse

from samuraiAnalysis import plotFLpath,samImport_masked,samImport,samPlt,makeKML,getFLpathData

getFLpathData = getFLpathData.getFLpathData
plotFLpath = plotFLpath.plotFLpath
samImport_masked = samImport_masked.samImport_masked
samImport = samImport.samImport
plotContour = samPlt.plotContour
plotVec = samPlt.plotVec
plotBarbs = samPlt.plotBarbs
plotStream = samPlt.plotStream
plotXS = samPlt.plotXS
map2xy = samPlt.map2xy

warnings.filterwarnings('ignore', 'invalid value encountered in less')
warnings.filterwarnings('ignore', 'invalid value encountered in sqrt')
warnings.filterwarnings('ignore', 'invalid value encountered in greater_equal')
warnings.filterwarnings('ignore', 'tight_layout : falling back to Agg renderer')
warnings.filterwarnings('ignore', 'Warning: converting a masked element to nan.')

### Get YAML parameter file name/path
parser = argparse.ArgumentParser()
parser.add_argument("pFile", help="Path to PlotSamuraiXS YAML parameters file")
args = parser.parse_args()
pFile = args.pFile


#### `planContour` ####
# Takes a given variable and variable label and using other parameters defined 
# in this script, creates plan view contour plots.
def planContour(pltVar,pltVarLbl):

    fig,ax = plotContour(pltVar,lev,x1d,y1d,'xy',pltVarLbl,
                                  figsize=(20,16),dT=dT,runId=runId,strmRel=strmRel)
    
    
    if vecType == 'vec':
        plotVec(x1d[1::7],y1d[1::7],u[int(lev*2),1::7,1::7],v[int(lev*2),1::7,1::7],'xy',scale=150)
    if vecType == 'barb':
        plotBarbs(x1d[1::7],y1d[1::7],u[int(lev*2),1::7,1::7],v[int(lev*2),1::7,1::7],'xy',length=7)
    if vecType == 'stream':
        plotStream(x1d,y1d,u[int(lev*2),:,:],v[int(lev*2),:,:],'xy',strmRel=strmRel)
    
    
    if pltFT:
        fl_x,fl_y = map2xy(flLon,flLat,lon_0,lat_0)
        plotFLpath(fl_x,fl_y,crds,dubLine=True)
        
    if plotXSloc:
        if xsLocs is not None:
            locIter = xsLocs
        else:
            locIter = range(0,len(xsStrt))
        for iz in locIter:
            x1_tmp = xsStrt[iz][0]
            y1_tmp = xsStrt[iz][1]
            x2_tmp = xsEnd[iz][0]
            y2_tmp = xsEnd[iz][1]
            
            ax.plot([x1_tmp,x2_tmp],[y1_tmp,y2_tmp],linestyle='-',linewidth=2,color='k')
            
            # Plot markers on XS line - need to adjust the increment depending on XS length (need to automate this)
#             xsLen = 170 #S2-diag
#             xsLen = 175 #S2-WE
#             xsLen = 200 #S5-diag
#             xsLen = 175 #S5-WE
#             xsLen = 185 #S7-both
#             spc = (xsLen/5)+1
#             ax.scatter(np.linspace(x1_tmp,x2_tmp,spc),np.linspace(y1_tmp,y2_tmp,spc),s=85,color='b',edgecolor='w')
            
            # Alternate labeling start/end of XS lines
#             if iz%2 == 0:
#                 xyP = (x1_tmp,y1_tmp) 
#                 xyT = (-25,-25)
#             else:
#                 xyP = (x2_tmp,y2_tmp) 
#                 xyT = (-25,25)
#             ax.annotate(str(iz+1),xy=xyP,xytext=xyT,zorder=10,
#                         textcoords='offset points',fontsize=13,color='k',
#                         bbox=dict(boxstyle="round", fc="w",alpha=0.7,pad=0.04),
#                         arrowprops=dict(fc='w',ec='k',alpha=0.7,shrink=0.05,width=2,headwidth=8,headlength=8))
    

    if strmRel:
        srStr = '-SR'
    else:
        srStr = ''
        
    if vecType is not None:
        vecStr = '-{}'.format(vecType)
    else:
        vecStr = ''
        
    saveStr = '{}{}{}_{}_{}{}{}_{:.1f}km.{}'.format(savePath,dtSave,runIdSv,'xy',pltVarLbl,vecStr,srStr,lev,fType)  
    fig.savefig(saveStr,bbox_inches='tight')

    plt.close()


#### `xsContour` ####
# Takes a given variable, cross-section start/end points, 
# and using other parameters defined in this script, creates cross-section contour plots.
def xsContour(pltVar,pltVarLbl,xsStrtPt,xsEndPt,xsIX):
    fig,ax,xsAngle,xsLen_km,xsNumPts,xs_x1d,xs_y1d = plotXS(pltVar,x1d,y1d,alt,u,v,w,xsStrtPt,xsEndPt,
                                                            pltFlag=pltVarLbl,dT=dT,runId=runId,xsIX=xsIX+1,
                                                            vecType=vecTypeXS,figsize=(27,9),
                                                            strmRel=strmRel,scaleVmag=scaleVmag,
                                                            doUcont=doUcont,doUVoutline=doUVoutline)
    
    print('\t\txsAngle = {}'.format(xsAngle))
    print('\t\txsLen_km = {}'.format(xsLen_km))
    print('\t\txsNumPts = {}'.format(xsNumPts))
    #print('xs_x1d:\n\t{}\n'.format(np.array_str(xs_x1d, precision=2)))
    #print('xs_y1d:\n\t{}\n'.format(np.array_str(xs_y1d, precision=2)))

    xsStrng = '[{},{}]-[{},{}]'.format(xsStrtPt[0],xsStrtPt[1],xsEndPt[0],xsEndPt[1])
    
    if strmRel:
        srStr = '-SR'
    else:
        srStr = ''
        
    if vecTypeXS is not None:
        vecStr = '-{}'.format(vecTypeXS)
    else:
        vecStr = ''
    
    saveStr = '{}{}{}_{}{}{}_{}-{}_{}.{}'.format(savePath,dtSave,runIdSv,pltVarLbl,vecStr,srStr,'XS',xsIX+1,xsStrng,fType)
    fig.savefig(saveStr,bbox_inches='tight')
    
    plt.close()



with open(pFile, 'r') as pltSamPrms:
    prmsIn = yaml.load(pltSamPrms)

for key,val in prmsIn.items():
        exec(key + '=val')


print('\nWorking on case: {}'.format(outPrefix_master))
print('\tUsing SAMURAI parameter file: {}\n'.format(pFile))

### Import data and initialize save path ###
samFileMaster = samPrefix + flight + '/' + outPrefix_master + '/output/samurai_XYZ_analysis.nc'
if maskCombo:
    if outPrefix_sub1 and outPrefix_sub2 and outPrefix_sub3:
        print('\tMasking velocities on 3 additional analyses...\n')
        samFileSub1 = samPrefix + flight + '/' + outPrefix_sub1 + '/output/samurai_XYZ_analysis.nc'
        samFileSub2 = samPrefix + flight + '/' + outPrefix_sub2 + '/output/samurai_XYZ_analysis.nc'
        samFileSub3 = samPrefix + flight + '/' + outPrefix_sub3 + '/output/samurai_XYZ_analysis.nc'
        samDataMaster = samImport_masked(samFileMaster,samFileSub1,samFile_sub2=samFileSub2,samFile_sub3=samFileSub3)
    elif outPrefix_sub1 and outPrefix_sub2:
        print('\tMasking velocities on 1 additional analyses...\n')
        samFileSub1 = samPrefix + flight + '/' + outPrefix_sub1 + '/output/samurai_XYZ_analysis.nc'
        samFileSub2 = samPrefix + flight + '/' + outPrefix_sub2 + '/output/samurai_XYZ_analysis.nc'
        samDataMaster = samImport_masked(samFileMaster,samFileSub1,samFile_sub2=samFileSub2)
    elif outPrefix_sub1:
        print('\tMasking velocities on 1 additional analysis...\n')
        samFileSub1 = samPrefix + flight + '/' + outPrefix_sub1 + '/output/samurai_XYZ_analysis.nc'
        samDataMaster = samImport_masked(samFileMaster,samFileSub1)
    else:
        print('\tmaskCombo is True, but no outPrefixes for masking analyses were given. Running normal data import...\n')
        samDataMaster = samImport(samFileMaster)
else:
    samDataMaster = samImport(samFileMaster)
    


if dtDir:
    dtNow = dt.strftime(dt.now(),'%Y%m%d_%H%M')
    savePath = samPrefix + flight + '/' + outPrefix_master + '/figs/' + dirPost + dtNow + '/'
else:
    savePath = samPrefix + flight + '/' + outPrefix_master + '/figs/' + dirPost + '/'



### Import Data / Initialize File Saving Parameters
dT = samDataMaster['time']
dbz = samDataMaster['dbz']
u = samDataMaster['u']
v = samDataMaster['v']
w = samDataMaster['w']
vort = samDataMaster['vort']
lat = samDataMaster['lat']
lon = samDataMaster['lon']
x1d = samDataMaster['x']
y1d = samDataMaster['y']
alt = samDataMaster['alt']

wndSpd = np.sqrt(np.square(u)+np.square(v))

if strmRel:
    u = u-strmU
    v = v-strmV


# Determine if we're plotting every or every other level
# (if neither -1 nor -2, then pltLevs is given a list of levels in the params file)
if pltLevs[0] == -1:
    pltLevs = samDataMaster['alt']
elif pltLevs[0] == -2:
    pltLevs = samDataMaster['alt'][2::2]


# Import FL data if plotting the flight track
if pltFT and plotPlanViews:
    flFile = flPrefix + flight + '_FltLvl_Processed.nc'
    flS = dt.strptime(flSs,'%Y%m%d-%H%M%S')
    flE = dt.strptime(flEs,'%Y%m%d-%H%M%S')
    flData = getFLpathData(flFile,pathStrt=flS,pathEnd=flE,crdsOnly=True)
    flLat = flData['lat']
    flLon = flData['lon']


# Make the runID and datetime filename-friendly
runId = outPrefix_master
runIdSv = '_' + runId.replace('_','-')


if dT is not None:
    dtSave = dt.strftime(dT,'%Y%m%d_%H%M')
else:
    dtSave = 'unknownDT'


# Make the output figure directory if it doesn't exist
if not os.path.exists(savePath):
    os.makedirs(savePath)

# Save a copy of the YAML paramaters file to the output directory
pFile_cpApnd = savePath[(savePath.rindex('figs/')+5):-1]
cpPfile_name = '{}{}_{}.yml'.format(savePath,pFile[(pFile.rindex('/')+1):-4],pFile_cpApnd)
shutil.copy(pFile,cpPfile_name)




### Generate all desired plan view plots
if plotPlanViews:

    for ix in range(0,len(pltLevs)):
        lev = pltLevs[ix]

        print('\tNow plotting plan views at {:.2f} km ({} of {})...'.format(lev,ix+1,len(pltLevs)))

        if pltDBZ:
            planContour(dbz,'dbz')
        if pltU:
            planContour(u,'u')
        if pltV:
            planContour(v,'v')
        if pltW:
            planContour(w,'w')
        if pltVort:
            planContour(vort,'vort')
        if pltWndSpd:
            planContour(wndSpd,'wind')


### Generate all desired cross-section plots
if plotXScts:
    if xsLocs is not None:
        locIter = xsLocs
    else:
        locIter = range(0,len(xsStrt))                                                                                    
    
    for ix in locIter:
        xsStrtTmp = xsStrt[ix]
        xsEndTmp = xsEnd[ix]
        
        print('\tNow plotting cross-section between ({},{}) and ({},{}) ({} of {})...'.format(xsStrtTmp[0],xsStrtTmp[1],xsEndTmp[0],
                                                                                          xsEndTmp[1],ix+1,len(xsStrt)))

        if pltDBZx:
            xsContour(dbz,'dbz',xsStrtTmp,xsEndTmp,ix)
        if pltUx:
            xsContour(u,'u',xsStrtTmp,xsEndTmp,ix)
        if pltVx:
            xsContour(v,'v',xsStrtTmp,xsEndTmp,ix)
        if pltWx:
            xsContour(w,'w',xsStrtTmp,xsEndTmp,ix)
        if pltVortX:
            xsContour(vort,'vort',xsStrtTmp,xsEndTmp,ix)
        if pltWndSpdX:
            xsContour(wndSpd,'wind',xsStrtTmp,xsEndTmp,ix)