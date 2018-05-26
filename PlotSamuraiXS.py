import os
import shutil
import warnings
from datetime import datetime as dt
import matplotlib as mpl
# mpl.use('PDF')
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import yaml
from glob import glob
import argparse

from samuraiAnalysis import getVarLims,gribTools,plotFLpath,samImport_masked,samImport,samPlt,makeKML,getFLpathData

getFLpathData = getFLpathData.getFLpathData
getVarLims = getVarLims.getVarLims
gribImport = gribTools.gribImport
gribLevs = gribTools.gribLevs
plotFLpath = plotFLpath.plotFLpath
samImport_masked = samImport_masked.samImport_masked
samImport = samImport.samImport
plotContour = samPlt.plotContour
plotVec = samPlt.plotVec
plotBarbs = samPlt.plotBarbs
plotXS_mapCrds = samPlt.plotXS_mapCrds
plotXS = samPlt.plotXS
xsCalc = samPlt.xsCrdCalc
multXSCrdCalc = samPlt.multXSCrdCalc
makeKML = makeKML.makeKML

warnings.filterwarnings('ignore', 'invalid value encountered in less')
warnings.filterwarnings('ignore', 'tight_layout : falling back to Agg renderer')

### Get YAML parameter file name/path
parser = argparse.ArgumentParser()
parser.add_argument("pFile", help="Path to PlotSamuraiXS YAML parameters file")
args = parser.parse_args()
pFile = args.pFile


#### `planContour` ####
# Takes a given variable and variable label and using other parameters defined 
# in this script, creates plan view contour plots.
def planContour(pltVar,pltVarLbl,crds):
    
    if crds == 'map':
        fig,ax,grd,proj = plotContour(pltVar,lev,lon,lat,'map',pltVarLbl,
                                      figsize=(10,8),dT=dT,runId=runId,
                                      zoom=zoom,mapBnds=mapBnds,NB=noBorder,strmRel=strmRel)
    
        if pltFT:
            plotFLpath(proj,flLon,flLat,dubLine=True)

        if pltVec:
            plotVec(lon[1::5],lat[1::5],u[int(lev*2),1::5,1::5],
                    v[int(lev*2),1::5,1::5],'map',proj=proj)
        if pltBarb:
            plotBarbs(x1d[1::5],y1d[1::5],u[int(lev*2),1::5,1::5],
                      v[int(lev*2),1::5,1::5],'map',proj=proj)
        
        if pltRAP:
                plotVec(gribData['lon'],gribData['lat'],gribData['u'][rapLevIx[ix],:,:],
                        gribData['v'][rapLevIx[ix],:,:],crds='map',proj=proj,color='white',
                        linewidth=1.5,edgecolor='white',quivKeyX=0.08)
                plotVec(gribData['lon'],gribData['lat'],gribData['u'][rapLevIx[ix],:,:],
                        gribData['v'][rapLevIx[ix],:,:],crds='map',proj=proj,color='blue',
                        edgecolor=None,quivKeyX=0.08)

        if pltGRlocs:
            ax.scatter(lonGR,latGR,marker='d',color='w',s=200,zorder=7,transform=proj,clip_on=True)
            ax.scatter(lonGR,latGR,marker='d',color='k',s=100,zorder=8,transform=proj,clip_on=True)

            for labelGR, c, d in zip(grTxt, lonGR, latGR):
                ax.annotate(labelGR, xy = (c, d), xytext = (-8, 5),zorder=10,
                             textcoords = 'offset points',fontsize=13,color='w',
                             bbox=dict(boxstyle="round", fc="b",alpha=0.6,pad=0.01),
                             transform=proj,clip_on=True)
            
        if plotXSloc:
            if len(xsStrt_map) > 1:
                lbls = []
                for iz in range(0,len(xsStrt_map)):
                    xsLat,xsLon,dRnd,hdng = xsCalc(xsStrt_map[iz][0],xsStrt_map[iz][1],xsEnd_map[iz][0],xsEnd_map[iz][1])
                    ax.plot(xsLon,xsLat,transform=proj,linestyle='-',linewidth=2,color='k')
                    if iz%2 == 0:
                        xyP = (xsStrt_map[iz][1],xsStrt_map[iz][0]) # Label the start of the cross section
                        xyT = (-25,-25)
                    else:
                        xyP = (xsEnd_map[iz][1],xsEnd_map[iz][0]) # Label the end of the cross section
                        xyT = (-25,25)
                    ax.annotate(str(iz+1),xy=xyP,xytext=xyT,zorder=10,transform=proj,
                                textcoords='offset points',fontsize=13,color='k',
                                bbox=dict(boxstyle="round", fc="w",alpha=0.7,pad=0.04),
                                arrowprops=dict(fc='w',ec='k',alpha=0.7,shrink=0.05,width=2,headwidth=8,headlength=8))
                    if iz == len(xsStrt_map)-1:
                        lbls.append('{} - ({:.2f},{:.2f})-({:.2f},{:.2f})'.format(iz+1,xsStrt_map[iz][0],xsStrt_map[iz][1],
                                                                                  xsEnd_map[iz][0],xsEnd_map[iz][1]))
                    else:
                        lbls.append('{} - ({:.2f},{:.2f})-({:.2f},{:.2f})\n'.format(iz+1,xsStrt_map[iz][0],xsStrt_map[iz][1],
                                                                                    xsEnd_map[iz][0],xsEnd_map[iz][1]))

                lblStr = ''.join(lbls)
                # ax.text(0.99,0.99,lblStr,ha='right',va='top',transform=ax.transAxes,
                #          bbox=dict(boxstyle='square',fc='w',ec='k',pad=0.15,alpha=0.75))
                xsStrng = 'multXS'
            else:
                xsLat,xsLon,dRnd,hdng = xsCalc(xsStrt_map[0][0],xsStrt_map[0][1],xsEnd_map[0][0],xsEnd_map[0][1])
                ax.plot(xsLon,xsLat,transform=proj,linestyle='-',linewidth=2,color='k')
                xsStrng = '{:.0f}{:.0f}-{:.0f}{:.0f}'.format(xsStrt_map[0][0]*10,xsStrt_map[0][1]*-10,
                                                           xsEnd_map[0][0]*10,xsEnd_map[0][1]*-10)
        
            if noBorder:
                if strmRel:
                    saveStr = '{}NB_{}{}_{}_{}-SR_{:.1f}km_{}.{}'.format(savePath,dtSave,runIdSv,'map',pltVarLbl,lev,xsStrng,fType)
                    fig.savefig(saveStr,bbox_inches='tight',pad_inches=0)
                if strmRel:
                    saveStr = '{}NB_{}{}_{}_{}_{:.1f}km_{}.{}'.format(savePath,dtSave,runIdSv,'map',pltVarLbl,lev,xsStrng,fType)
                    fig.savefig(saveStr,bbox_inches='tight',pad_inches=0)
            else:
                if strmRel:
                    saveStr = '{}{}{}_{}_{}-SR_{:.1f}km_{}.{}'.format(savePath,dtSave,runIdSv,'map',pltVarLbl,lev,
                                                                                           xsStrng,fType)
                else:
                    saveStr = '{}{}{}_{}_{}_{:.1f}km_{}.{}'.format(savePath,dtSave,runIdSv,'map',pltVarLbl,lev,
                                                                                           xsStrng,fType)
                fig.savefig(saveStr,bbox_inches='tight')
        
        else:
            fig.tight_layout()
            if noBorder:
                saveStr = '{}NB_{}{}_{}_{}_{:.1f}km.{}'.format(savePath,dtSave,runIdSv,'map',pltVarLbl,lev,fType)
                fig.savefig(saveStr,bbox_inches='tight',pad_inches=0)
    
            else:
                saveStr = '{}{}{}_{}_{}_{:.1f}km.{}'.format(savePath,dtSave,runIdSv,'map',pltVarLbl,lev,fType)
                fig.savefig(saveStr,bbox_inches='tight')
    
    
        if saveKML:
            kmlSaveStr = '{}NB_{}{}_{}_{}_{:.1f}km.kmz'.format(savePath,dtSave,runIdSv,'map',pltVarLbl,lev)
            makeKML(lon.min(),lat.min(),lon.max(),lat.max(),figs=[saveStr],kmzfile=kmlSaveStr)
            os.remove(saveStr)


    if crds == 'xy':
        fig,ax = plotContour(pltVar,lev,x1d,y1d,'xy',pltVarLbl,
                                      figsize=(10,8),dT=dT,runId=runId,strmRel=strmRel)
    

        if pltVec:
            plotVec(x1d[1::7],y1d[1::7],u[int(lev*2),1::7,1::7],v[int(lev*2),1::7,1::7],'xy')
        if pltBarb:
            plotBarbs(x1d[1::7],y1d[1::7],u[int(lev*2),1::7,1::7],v[int(lev*2),1::7,1::7],'xy')
        
            
        if plotXSloc:
            if len(xsStrt) > 1:
                lbls = []
                for iz in range(0,len(xsStrt)):
                    x1_tmp = xsStrt[iz][0]
                    y1_tmp = xsStrt[iz][1]
                    x2_tmp = xsEnd[iz][0]
                    y2_tmp = xsEnd[iz][1]
                    ax.plot([x1_tmp,x2_tmp],[y1_tmp,y2_tmp],linestyle='-',linewidth=2,color='k')
                    if iz%2 == 0:
                        xyP = (x1_tmp,y1_tmp) # Label the start of the cross section
                        xyT = (-25,-25)
                    else:
                        xyP = (x2_tmp,y2_tmp) # Label the end of the cross section
                        xyT = (-25,25)
                    ax.annotate(str(iz+1),xy=xyP,xytext=xyT,zorder=10,
                                textcoords='offset points',fontsize=13,color='k',
                                bbox=dict(boxstyle="round", fc="w",alpha=0.7,pad=0.04),
                                arrowprops=dict(fc='w',ec='k',alpha=0.7,shrink=0.05,width=2,headwidth=8,headlength=8))
                    if iz == len(xsStrt)-1:
                        lbls.append('{} - ({},{})-({},{})'.format(iz+1,x1_tmp,y1_tmp,x2_tmp,y2_tmp))
                    else:
                        lbls.append('{} - ({},{})-({},{})\n'.format(iz+1,x1_tmp,y1_tmp,x2_tmp,y2_tmp))

                lblStr = ''.join(lbls)

                xsStrng = 'multXS'
            else:
                x1_tmp = xsStrt[0][0]
                y1_tmp = xsStrt[0][1]
                x2_tmp = xsEnd[0][0]
                y2_tmp = xsEnd[0][1]
                ax.plot([x1_tmp,x2_tmp],[y1_tmp,y2_tmp],linestyle='-',linewidth=2,color='k')
                xsStrng = '[{},{}]-[{},{}]'.format(x1_tmp,y1_tmp,x2_tmp,y2_tmp)

            if strmRel:
                saveStr = '{}{}{}_{}_{}-SR_{:.1f}km_{}.{}'.format(savePath,dtSave,runIdSv,'xy',pltVarLbl,lev,xsStrng,fType)
            else:
                saveStr = '{}{}{}_{}_{}_{:.1f}km_{}.{}'.format(savePath,dtSave,runIdSv,'xy',pltVarLbl,lev,xsStrng,fType)
            
            fig.savefig(saveStr,bbox_inches='tight')
        
        else:
            # fig.tight_layout()
            if strmRel:
                saveStr = '{}{}{}_{}_{}-SR_{:.1f}km.{}'.format(savePath,dtSave,runIdSv,'xy',pltVarLbl,lev,fType)
            else:
                saveStr = '{}{}{}_{}_{}_{:.1f}km.{}'.format(savePath,dtSave,runIdSv,'xy',pltVarLbl,lev,fType)
            
            fig.savefig(saveStr,bbox_inches='tight')


    plt.close()


#### `xsContour` ####
# Takes a given variable, cross-section start/end points, 
# and using other parameters defined in this script, creates cross-section contour plots.
def xsContour(pltVar,pltVarLbl,xsStrtPt,xsEndPt,xsIX,crds):
    
    if crds == 'map':
        lonDiff = np.abs(xsEndPt[1]-xsStrtPt[1])
        latDiff = np.abs(xsEndPt[0]-xsStrtPt[0])
    
        if latDiff >= lonDiff:
            xCrd = 'lat'
        else:
            xCrd = 'lon'
    
        fig,ax = plotXS_mapCrds(pltVar,lon,lat,alt,u,v,w,xsStrtPt,xsEndPt,
                      pltFlag=pltVarLbl,dT=dT,leafSz=leafSz,runId=runId,
                      vecPlt=pltVecXS,barbPlt=pltBarbXS,figsize=(15,5),xCrd=xCrd,strmRel=strmRel)
    
        if strmRel:
            saveStr = '{}{}{}_{}-SR_{}-{}_{:.0f}{:.0f}-{:.0f}{:.0f}.{}'.format(savePath,dtSave,runIdSv,pltVarLbl,'XS',xsIX+1,
                                                                          xsStrtPt[0]*10,xsStrtPt[1]*-10,
                                                                          xsEndPt[0]*10,xsEndPt[1]*-10,
                                                                           fType)
        else:
            saveStr = '{}{}{}_{}_{}-{}_{:.0f}{:.0f}-{:.0f}{:.0f}.{}'.format(savePath,dtSave,runIdSv,pltVarLbl,'XS',xsIX+1,
                                                                          xsStrtPt[0]*10,xsStrtPt[1]*-10,
                                                                          xsEndPt[0]*10,xsEndPt[1]*-10,
                                                                           fType)
        fig.savefig(saveStr,bbox_inches='tight')
        
    if crds == 'xy':
        fig,ax = plotXS(pltVar,x1d,y1d,alt,u,v,w,xsStrtPt,xsEndPt,
                      pltFlag=pltVarLbl,dT=dT,runId=runId,
                      vecPlt=pltVecXS,barbPlt=pltBarbXS,figsize=(15,5),strmRel=strmRel)
    
        x1_tmp = xsStrtPt[0]
        y1_tmp = xsStrtPt[1]
        x2_tmp = xsEndPt[0]
        y2_tmp = xsEndPt[1]
        xsStrng = '[{},{}]-[{},{}]'.format(x1_tmp,y1_tmp,x2_tmp,y2_tmp)
        
        if strmRel:
            saveStr = '{}{}{}_{}-SR_{}-{}_{}.{}'.format(savePath,dtSave,runIdSv,pltVarLbl,'XS',xsIX+1,xsStrng,fType)
        else:
            saveStr = '{}{}{}_{}_{}-{}_{}.{}'.format(savePath,dtSave,runIdSv,pltVarLbl,'XS',xsIX+1,xsStrng,fType)
        
        fig.savefig(saveStr,bbox_inches='tight')
    
    plt.close()




# prmFiles = sorted(glob('/Users/danstechman/GoogleDrive/PECAN-Data/samurai/20150620/p3s2_udx/*.yml'))
# for pFile in prmFiles:
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
    
     
flFile = flPrefix + flight + '_FltLvl_Processed.nc'

if dtDir:
    dtNow = dt.strftime(dt.now(),'%Y%m%d_%H%M')
    savePath = samPrefix + flight + '/' + outPrefix_master + '/figs/' + dtNow + dirPost + '/'
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


# Various settings for map plots
if crds == 'map':
    # Define set of parallel cross-sections if `multXS` is True
    if plotXScts:
        if multXS:
            # Generate a number of parallel cross-sections spaced at a given interval
            # starting at some initial set of XS start/end points
            dists = np.arange(1,numXS)

            xsStrt_map = []
            xsEnd_map = []

            for d in dists:
                newLatStrt,newLonStrt = multXSCrdCalc(multXSinitLatStrt,multXSinitLonStrt,d,multXShdng)
                xsStrt_map.append((newLatStrt,newLonStrt))

                newLatEnd,newLonEnd = multXSCrdCalc(multXSinitLatEnd,multXSinitLonEnd,d,multXShdng)
                xsEnd_map.append((newLatEnd,newLonEnd))
        if not avgXS:
            leafSz = 1
            
    if unifyBnds:
        zoom = True
        mapBnds = unifiedMapBnds
    else:
        zoom = False
        mapBnds = None

    # Check to make sure we're creating borderless figures if saving KML files
    if saveKML and not noBorder:
        print('\t\'No border\' must be True if saving KML. Setting noBorder equal to True...')
        noBorder = True
        
    # Import FL data if plotting the flight track
    if pltFT and plotPlanViews:
        flS = dt.strptime(flSs,'%Y%m%d-%H%M%S')
        flE = dt.strptime(flEs,'%Y%m%d-%H%M%S')
        flData = getFLpathData(flFile,pathStrt=flS,pathEnd=flE,crdsOnly=True)
        flLat = flData['lat']
        flLon = flData['lon']


# Make the runID and datetime filename-friendly
if not pltRAP:
    runId = outPrefix_master
    # runId = outPrefix_master[5:]
elif crds == 'map':
    runId = outPrefix_master + '_rapOvr'

    # Import model (RAP) data if desired and match levels
    gribData = gribImport(gribFile)
    rapLevIx = gribLevs(pltLevs,gribData['geoHght'])


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
            planContour(dbz,'dbz',crds)
        if pltU:
            planContour(u,'u',crds)
        if pltV:
            planContour(v,'v',crds)
        if pltW:
            planContour(w,'w',crds)
        if pltVort:
            planContour(vort,'vort',crds)
        if pltWndSpd:
            planContour(wndSpd,'wind',crds)


### Generate all desired cross-section plots
if plotXScts:

    if crds == 'map':
        
        for ix in range(0,len(xsStrt_map)):

            xsStrtTmp = xsStrt_map[ix]
            xsEndTmp = xsEnd_map[ix]
            
            print('\tNow plotting cross-section between ({:.2f},{:.2f}) and ({:.2f},{:.2f}) ({} of {})...'.format(xsStrtTmp[0],xsStrtTmp[1],
                                                                                                              xsEndTmp[0],xsEndTmp[1],
                                                                                                              ix+1,len(xsStrt_map)))
        
            

            if pltDBZx:
                xsContour(dbz,'dbz',xsStrtTmp,xsEndTmp,ix,crds)
            if pltUx:
                xsContour(u,'u',xsStrtTmp,xsEndTmp,ix,crds)
            if pltVx:
                xsContour(v,'v',xsStrtTmp,xsEndTmp,ix,crds)
            if pltWx:
                xsContour(w,'w',xsStrtTmp,xsEndTmp,ix,crds)
            if pltVortX:
                xsContour(vort,'vort',xsStrtTmp,xsEndTmp,ix,crds)
            if pltWndSpdX:
                xsContour(wndSpd,'wind',xsStrtTmp,xsEndTmp,ix,crds)
                
    elif crds == 'xy':
        
                                                                                              
        for ix in range(0,len(xsStrt)):

            xsStrtTmp = xsStrt[ix]
            xsEndTmp = xsEnd[ix]
            
            print('\tNow plotting cross-section between ({},{}) and ({},{}) ({} of {})...'.format(xsStrtTmp[0],xsStrtTmp[1],xsEndTmp[0],
                                                                                              xsEndTmp[1],ix+1,len(xsStrt)))
        
            

            if pltDBZx:
                xsContour(dbz,'dbz',xsStrtTmp,xsEndTmp,ix,crds)
            if pltUx:
                xsContour(u,'u',xsStrtTmp,xsEndTmp,ix,crds)
            if pltVx:
                xsContour(v,'v',xsStrtTmp,xsEndTmp,ix,crds)
            if pltWx:
                xsContour(w,'w',xsStrtTmp,xsEndTmp,ix,crds)
            if pltVortX:
                xsContour(vort,'vort',xsStrtTmp,xsEndTmp,ix,crds)
            if pltWndSpdX:
                xsContour(wndSpd,'wind',xsStrtTmp,xsEndTmp,ix,crds)