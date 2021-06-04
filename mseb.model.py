#   ---------------------------------------------------------------------------------
#   The moist static energy equation for first baroclinic mode anomalies (MSEB) model
#   ---------------------------------------------------------------------------------
# 
#      Authors: Chen-Shuo Fan and Dietmar Dommenget
# 
#
#   References: "A Diagnostic Model for The Large-Scale Tropical Circulation Based on Moist Static Energy Balance"
#               by Chen-Shuo Fan and Dietmar Dommenget, 
#               submitted to Climate Dynamics (2021)
#
#
#   input data: The MSEB model needs the following fields to be defined before running it,
# 
#               temp(level,latitude,longitude): air temperature                 [K]
#               shum(level,latitude,longitude): specific humidity               [kg/kg]
#               geop(level,latitude,longitude): geopotential                    [m^2/s^2]
#               uwnd(level,latitude,longitude): eastward wind                   [m/s]
#               vwnd(level,latitude,longitude): northward wind                  [m/s]
#                     slhf(latitude,longitude): surface latent heat flux        [J/m^2]
#                     sshf(latitude,longitude): surface sensible heat flux      [J/m^2]
#                     snsr(latitude,longitude): surface net solar radiation     [J/m^2]
#                     sntr(latitude,longitude): surface net thermal radiation   [J/m^2]
#                     tnsr(latitude,longitude): top net solar radiation         [J/m^2]
#                     tntr(latitude,longitude): top net thermal radiation       [J/m^2]
#
#
#   experiments: log_exp=0 vertical velocity [Pa/s] of each month (default)
#                log_exp=1 vertical velocity [Pa/s] of climatology (annual mean)
#                log_exp=2 vertical velocity [Pa/s] of seasonal cycle (JJA minus DJF)   
#                log_exp=3 vertical velocity [Pa/s] of El Nino year (El Nino case minus climatology)

from geopy.distance import geodesic
from scipy.interpolate import CubicSpline
from scipy.integrate import simps
from metpy.units import units
from metpy import calc
import xarray as xr
import numpy as np



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def deep_convec_mode(temp):
    '''This function returns 2 outputs, omega(level): deep convection mode / ptop: tropopause pressure level'''
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # tropical-mean temperature profile
    temp_tropic = np.mean(temp,axis=(1,2))
    # interpolate to higher vertical resolution temperature profile
    temp_spline = CubicSpline(temp.level, temp_tropic); splinePts = 1000
    pHiRes=units.Quantity(np.linspace(temp.level[0], temp.level[-1], num=splinePts),"hPa")
    zHiRes=calc.pressure_to_height_std(pHiRes).magnitude

    # define tropopause as the lowest level at which the lapse rate >= -2 K/km in high-resolution vertical temperature profile
    i=0
    while ( float ( ( (temp_spline(pHiRes[splinePts-2-i-1]) - temp_spline(pHiRes[splinePts-2-i]))/ (zHiRes[splinePts-2-i-1]-zHiRes[splinePts-2-i]) ) ) ) <= (-2.0):
        ptop=pHiRes[splinePts-2-i]
        i+=1
            
    # replace low-resolution tropopause pressure level with high-resolution tropopause 
    pIdx=len(temp.level)
    while temp.level[pIdx-1] > ptop.magnitude:
        pIdx-=1 
    temp_tropic[pIdx-1] = temp_spline(ptop)
    plev_new = temp.level[pIdx-1::].values
    plev_new[0] = ptop.magnitude

    # saturation specific humidity (q_sat) of eq.(A2) - Wills et al. (2017)
    e_sat = 0.611*np.exp( ( 17.3*(temp_tropic[pIdx-1::]-273.15) )/( (temp_tropic[pIdx-1::]-273.15)+237.3 ) )
    q_sat = 0.622*( e_sat/( plev_new*10**(-1) ) )

    # eq.(A2) and eq.(A3) - Wills et al. (2017)
    q_t = ( (Lv*Lv*q_sat) / (Cp*Rv*temp_tropic[pIdx-1::]*temp_tropic[pIdx-1::]) )*( 1+(Rv/Rd)*q_sat-q_sat )
    z_t = (-1.0)*(Rd/Cp)*np.log10(plev_new/plev_new[-1])

    # eq.(A2) and eq.(A3) at LCL - Wills et al. (2017)
    q_b = ( ( Lv*Lv*q_sat[-3] ) / ( Cp*Rv*temp_tropic[-3]*temp_tropic[-3] ) )*( 1+(Rv/Rd)*q_sat-q_sat )[-3] 
    z_b = (-1.0)*(Rd/Cp)*np.log10(plev_new[-3]/plev_new[-1])

    # eq.(A5) - Wills et al. (2017)
    A1 = (1+q_b+z_b)/(1+q_t+z_t)

    # eq.(A8) - Wills et al. (2017)
    A1_plus=np.zeros( len(plev_new) )
    for i in range(0,len(plev_new)):
        A1_plus[i] = simps(np.flip(A1)[0:i+1],np.flip(plev_new)[0:i+1])			
    A1_plus=A1_plus/np.flip(plev_new)	
            
    # eq.(A7) - Wills et al. (2017)
    V1=A1_plus-( simps(A1_plus,np.flip(plev_new))/( np.flip(plev_new)[-1]-np.flip(plev_new)[0] ) )	
            
    # eq.(A9) - Wills et al. (2017)
    omega=np.zeros( len(plev_new) )
    for i in range(0,len(plev_new)):
        omega[i] = simps(V1[0:i+1],np.flip(plev_new)[0:i+1])		
    omega=abs(omega/np.flip(plev_new)[0])	

    return omega, ptop.magnitude



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def mseb_model(slhf,sshf,snsr,sntr,tnsr,tntr,temp,shum,geop,uwnd,vwnd,dmode,ptop):
    '''This function returns 1 output, wMSEB(latitude,longitude): vertical velocity at 500hPa [Pa/s]'''
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # eq.(4) - Fan and Dommenget (2021)
    ftop = (tnsr+tntr)/(24*60*60) # convert unit to [W/m^2]
    # eq.(5) - Fan and Dommenget (2021)
    fsur = (snsr+sntr+slhf+sshf)/(24*60*60) # convert unit to [W/m^2]
    # eq.(3) - Fan and Dommenget (2021)
    fnet = ftop - fsur

    # replace low-resolution tropopause pressure level with high-resolution tropopause 
    pIdx=0
    while temp.level[pIdx] > ptop:
        pIdx+=1 
    plev = temp.level[0:pIdx+1].values
    plev[pIdx] = ptop
            
    # interpolate data points to high-resolution tropopause pressure level 
    temp[pIdx,:,:]=((temp[pIdx,:,:]-temp[pIdx-1,:,:])/(temp.level[pIdx]-temp.level[pIdx-1]))*(ptop-temp.level[pIdx-1])+temp[pIdx-1,:,:]
    shum[pIdx,:,:]=((shum[pIdx,:,:]-shum[pIdx-1,:,:])/(shum.level[pIdx]-shum.level[pIdx-1]))*(ptop-shum.level[pIdx-1])+shum[pIdx-1,:,:]
    geop[pIdx,:,:]=((geop[pIdx,:,:]-geop[pIdx-1,:,:])/(geop.level[pIdx]-geop.level[pIdx-1]))*(ptop-geop.level[pIdx-1])+geop[pIdx-1,:,:]
    uwnd[pIdx,:,:]=((uwnd[pIdx,:,:]-uwnd[pIdx-1,:,:])/(uwnd.level[pIdx]-uwnd.level[pIdx-1]))*(ptop-uwnd.level[pIdx-1])+uwnd[pIdx-1,:,:]
    vwnd[pIdx,:,:]=((vwnd[pIdx,:,:]-vwnd[pIdx-1,:,:])/(vwnd.level[pIdx]-vwnd.level[pIdx-1]))*(ptop-vwnd.level[pIdx-1])+vwnd[pIdx-1,:,:]

    # horizontal gradient of temperature/specific humidity
    dTdx,dQdx=(np.zeros((len(temp.level),len(temp.latitude),len(temp.longitude))) for i in range(2)); gridSize=abs(temp.latitude[1]-temp.latitude[0])
    for i in range( len(temp.latitude) ):
        dTdx[:,i,:] = np.gradient(temp[:,i,:],geodesic( (temp.latitude[i],0),(temp.latitude[i],gridSize) ).meters,axis=1)
        dQdx[:,i,:] = np.gradient(shum[:,i,:],geodesic( (temp.latitude[i],0),(temp.latitude[i],gridSize) ).meters,axis=1)
    dTdy=np.gradient(temp,geodesic( (0,0),(gridSize,0) ).meters,axis=1)
    dQdy=np.gradient(shum,geodesic( (0,0),(gridSize,0) ).meters,axis=1)

    # horizontal advection of eq.(10) - Fan and Dommenget (2021)
    advp = (uwnd*dTdx+vwnd*dTdy)*Cp+(uwnd*dQdx+vwnd*dQdy)*Lv
    adv = simps(advp[0:pIdx+1,:,:],plev[0:pIdx+1]*100.0,axis=0)*(-1.0/g)
            
    # moist static enrgy (MSE)
    mse = temp*Cp + shum*Lv + geop

    # gross moist stability (M) of eq.(10) - Fan and Dommenget (2021)
    gmsp = np.zeros((len(temp.level),len(temp.latitude),len(temp.longitude)))
    for z in range(1,pIdx):
        gmsp[z,:,:] = dmode[z]*((mse[z+1,:,:]-mse[z-1,:,:])/(plev[z+1]*100.0-plev[z-1]*100.0))
    gms = simps(gmsp[0:pIdx+1,:,:],plev[0:pIdx+1]*100.0,axis=0)*(-1.0/g)

    # eq.(10) - Fan and Dommenget (2021)
    wMSEB = dmode[5]*((fnet-adv)/gms) # dmode[5]: deep convction mode at 500hPa

    return wMSEB



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# test MSEB model with different input data sets as used in Fan and Dommenget (2021)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# experiment switch: 0. monthly / 1. climatology / 2. seasonal cycle / 3. El Nino 
log_exp = 0

# read monthly pressure levels/surface data
fn = 'era5prs.nc'; prs = xr.open_dataset(fn) 
fn = 'era5sfc.nc'; sfc = xr.open_dataset(fn)

#(!!! read data with ascendant pressure level/latitude)
slhf = sfc['slhf'] [:,::-1,:]; sshf = sfc['sshf'] [:,::-1,:]; snsr = sfc['ssr']  [:,::-1,:]
sntr = sfc['str']  [:,::-1,:]; tnsr = sfc['tsr']  [:,::-1,:]; tntr = sfc['ttr']  [:,::-1,:] 
temp = prs['t']  [:,:,::-1,:]; shum = prs['q']  [:,:,::-1,:]; geop = prs['z']  [:,:,::-1,:] 
uwnd = prs['u']  [:,:,::-1,:]; vwnd = prs['v']  [:,:,::-1,:]

# begin/end data time in year
yrBegin = 1979; yrEnd = 2018 

# physical parameter (natural constants)
Lv = 2.5e6  # specific latent heat of vaporization          [J/kg]
Cp = 1004.0 # specific heat capacity at constant pressure   [J/kg/K]
Rv = 461.52	# specific gas constants of water vapor         [J/kg/K]
Rd = 287.0  # specific gas constants of dry air             [J/kg/K]
g  = 9.8    # gravitational acceleration                    [m/s^2]

if log_exp == 0:
    print("exp.0: running default (monthly) experiment...")

    wMon = []
    for yr in range(yrEnd-yrBegin+1):

        for mon in range(12):
            # compute deep convection mode (!!! read data with descendant pressure level)
            dmode,ptop=deep_convec_mode(temp[yr*12+mon,::-1,:,:])
            # compute 500hPa vertical velocity 
            wMSEB=mseb_model(slhf.isel(time=yr*12+mon),sshf.isel(time=yr*12+mon),snsr.isel(time=yr*12+mon), \
                             sntr.isel(time=yr*12+mon),tnsr.isel(time=yr*12+mon),tntr.isel(time=yr*12+mon), \
                             temp.isel(time=yr*12+mon),shum.isel(time=yr*12+mon),geop.isel(time=yr*12+mon), \
                             uwnd.isel(time=yr*12+mon),vwnd.isel(time=yr*12+mon),dmode,ptop)
            wMon.append(wMSEB)

    # write out MSEB model output
    ds = xr.Dataset({"w": (("time","latitude","longitude"), wMon)},coords={"time": temp.time,"latitude": temp.latitude,"longitude": temp.longitude,},)
    ds.to_netcdf("mseb_output.exp0.nc"); print("exp.0: write results to mseb_output.exp0.nc"); print("exp.0: done")

if log_exp == 1:
    print("exp.1: running climatology experiment...")

    wMon = []
    for yr in range(yrEnd-yrBegin+1):

        for mon in range(12):
            # compute deep convection mode (!!! read data with descendant pressure level)
            dmode,ptop=deep_convec_mode(temp[yr*12+mon,::-1,:,:])
            # compute 500hPa vertical velocity 
            wMSEB=mseb_model(slhf.isel(time=yr*12+mon),sshf.isel(time=yr*12+mon),snsr.isel(time=yr*12+mon), \
                             sntr.isel(time=yr*12+mon),tnsr.isel(time=yr*12+mon),tntr.isel(time=yr*12+mon), \
                             temp.isel(time=yr*12+mon),shum.isel(time=yr*12+mon),geop.isel(time=yr*12+mon), \
                             uwnd.isel(time=yr*12+mon),vwnd.isel(time=yr*12+mon),dmode,ptop)
            wMon.append(wMSEB)

    wClim = np.mean(wMon,axis=0)
    # write out MSEB model output
    ds = xr.Dataset({"w": (("latitude","longitude"), wClim)},coords={"latitude": temp.latitude,"longitude": temp.longitude,},)
    ds.to_netcdf("mseb_output.exp1.nc"); print("exp.1: write results to mseb_output.exp1.nc"); print("exp.1: done")

if log_exp == 2:
    print("exp.2: running seasonal cycle experiment...")

    wDJF = []; wJJA = []
    for yr in range(yrEnd-yrBegin+1):

        for mon in ([0,1,5,6,7,11]):
            # compute deep convection mode (!!! read data with descendant pressure level)
            dmode,ptop=deep_convec_mode(temp[yr*12+mon,::-1,:,:])
            # compute 500hPa vertical velocity 
            wMSEB=mseb_model(slhf.isel(time=yr*12+mon),sshf.isel(time=yr*12+mon),snsr.isel(time=yr*12+mon), \
                             sntr.isel(time=yr*12+mon),tnsr.isel(time=yr*12+mon),tntr.isel(time=yr*12+mon), \
                             temp.isel(time=yr*12+mon),shum.isel(time=yr*12+mon),geop.isel(time=yr*12+mon), \
                             uwnd.isel(time=yr*12+mon),vwnd.isel(time=yr*12+mon),dmode,ptop)
            if mon in ([0,1,11]): wDJF.append(wMSEB)
            if mon in ([5,6,7]) : wJJA.append(wMSEB)

    wSeason = (np.mean(wJJA,axis=0) - np.mean(wDJF,axis=0))/2.0
    # write out MSEB model output
    ds = xr.Dataset({"w": (("latitude","longitude"), wSeason)},coords={"latitude": temp.latitude,"longitude": temp.longitude,},)
    ds.to_netcdf("mseb_output.exp2.nc"); print("exp.2: write results to mseb_output.exp2.nc"); print("exp.2: done")

if log_exp == 3:
    print("exp.3: running El Nino experiment...")

    ninoYr = [1982,1991,1997,2002,2009,2015]; wClim = []; wNino = []
    for yr in range(yrEnd-yrBegin+1):

        for mon in ([9,10,11]): 
            # compute deep convection mode (!!! read data with descendant pressure level)
            dmode,ptop=deep_convec_mode(temp[yr*12+mon,::-1,:,:])
            # compute 500hPa vertical velocity 
            wMSEB=mseb_model(slhf.isel(time=yr*12+mon),sshf.isel(time=yr*12+mon),snsr.isel(time=yr*12+mon), \
                             sntr.isel(time=yr*12+mon),tnsr.isel(time=yr*12+mon),tntr.isel(time=yr*12+mon), \
                             temp.isel(time=yr*12+mon),shum.isel(time=yr*12+mon),geop.isel(time=yr*12+mon), \
                             uwnd.isel(time=yr*12+mon),vwnd.isel(time=yr*12+mon),dmode,ptop)
            wClim.append(wMSEB)
            if (yr+yrBegin) in ninoYr: wNino.append(wMSEB)

    wEnso = np.mean(wNino,axis=0) - np.mean(wClim,axis=0)
    # write out MSEB model output
    ds = xr.Dataset({"w": (("latitude","longitude"), wEnso)},coords={"latitude": temp.latitude,"longitude": temp.longitude,},)
    ds.to_netcdf("mseb_output.exp3.nc"); print("exp.3: write results to mseb_output.exp3.nc"); print("exp.3: done")