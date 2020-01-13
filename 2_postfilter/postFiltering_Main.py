import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

from IPython.display import display, clear_output, SVG, Image, display_svg
import time

import os
import xarray as xr
import floater
from floater.generators import FloatSet
from floater import rclv
from skimage.feature import peak_local_max
from skimage.measure import find_contours
from skimage.filters import rank
from skimage.morphology import disk
from skimage.filters import sobel
from skimage.measure import label, regionprops

import csv
import pandas as pd
from scipy import interpolate
import datetime as dt

import shapely.geometry as sgeom
import warnings; warnings.simplefilter('ignore')

import shapely.geometry as sgeom
import multiprocessing
from functools import partial 
import gc
import scipy.ndimage as ndimage      

import math

# ---------------------------
def lat_lon_largest_local_min(data,loc):
    x = np.zeros(loc.shape[0])
    for a in np.arange(0,loc.shape[0],1):
        x[a] = data.load().data[loc[a][0],loc[a][1]]
    
    posi = np.argmin(x)
    localmin = x[posi]/100.0
    lat_min = data[loc[posi][0],loc[posi][1]].coords['lat'].values
    lon_min = data[loc[posi][0],loc[posi][1]].coords['lon'].values
    
    return localmin,lat_min,lon_min

# ---------------------------
def relative_location(lat,lon, xlat,xlon, ylat,ylon, Txlat,Txlon):
    # xlat, xlon: latitude and longitude of actually vortex max
    # ylat, ylon: latitude and longitude of actually slp min center
    # Txlat, Txlon: coordinates of xlat, xlon in floater
    i_xlat = np.abs(lat-xlat).argmin()
    i_xlon = np.abs(lon-xlon).argmin()
    
    i_ylat = np.abs(lat-ylat).argmin()
    i_ylon = np.abs(lon-ylon).argmin()
    
    delta_lat = i_xlat - i_ylat
    delta_lon = i_xlon - i_ylon
    
    #print(delta_lat)
    #print(delta_lon)
    #print(Txlat-delta_lat)
    #print(Txlon-delta_lon)
    
    Exlat = Txlat-delta_lat
    Exlon = Txlon-delta_lon
    
    if (Exlat<=0)|(Exlon<=0):
        #raise ValueError("the vortex center is not in the contour!")
        return None,None
    else:
        return Exlat, Exlon
    
# ---------------------------    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# ---------------------------    
def distance_on_unit_sphere(lat1, long1, lat2, long2):

    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
        
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
        
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
        
    # Compute spherical distance from spherical coordinates.
        
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
    
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos(cos)

    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc*6371.00

# -----------------------------------
def filtering(track,year):
    print("----------- year "+str(year)+" START -----------")
    
    track_iyr = track[track.year==year]
    reserve   = []
    plot      = False
    
    # -----------------------------------
    # loop through each LPS track
    for iLPS in track_iyr.index.unique():  
        # -----------------------------------
        # print the progress of filtering
        if iLPS%200==0:
            print("year: "+str(year)+"       |    progress:"\
             +str( np.around(100*(iLPS-track_iyr.index.unique()[0])/len(track_iyr.index.unique())) ) )
        
        # -----------------------------------
        # genesis and terminate location
        track_lon_b, track_lat_b, track_yr, track_mm, track_dd, track_hh, j,k,l = track.loc[iLPS].values[0]
        track_lon_e, track_lat_e, track_yr, track_mm, track_dd, track_hh, j,k,l = track.loc[iLPS].values[len(track.loc[iLPS].datetime)-1]
        
        # to the west of 100E
        if (track_lon_b<100)|(track_lat_b<18):continue
        # form over land
        landsea = lsm.sel(lat=track_lat_b,lon=track_lon_b,method='nearest')
        if (landsea<0.7):continue
        # form below surface
        #orography_of_vortex = orography.sel(lat=track_lat_b,lon=track_lon_b,method='nearest')
        #if orography_of_vortex==1:continue

        # if the system formed poleward of the monsoon domain and moves southward (to filter out mid-to-high latitidue systems)
        if (track_lat_b>bdry_lat[find_nearest(bdry_lon,track_lon_b)])&(track_lat_b>track_lat_e):continue 
        
        # cyclones must have four time steps equatorward of the northern boundary of the monsoon domain
        monsoon_step = 0
        for it in np.arange(0,len(track.loc[iLPS].datetime),1):
            track_lon, track_lat, track_yr, track_mm, track_dd, track_hh, j,k,l = track.loc[iLPS].values[it]
            if (track_lat<=bdry_lat[find_nearest(bdry_lon,track_lon)]):
                monsoon_step = monsoon_step+1
        if monsoon_step<2:continue

        # -----------------------------------
        # four time steps over land
        t_land = 0.0
        for it in np.arange(0,len(track.loc[iLPS].datetime),1):  
            track_lon, track_lat, track_yr, track_mm, track_dd, track_hh, j,k,l = track.loc[iLPS].values[it]
            landsea = lsm.sel(lat=track_lat,lon=track_lon,method='nearest')     
            if landsea>=0.7: 
                t_land = t_land+1.0
            else:break
        if (t_land<4):continue
        
        # path distance: sum of all great-circule distance between track positions >= 1000km
        travel_distance = 0
        for it in np.arange(0,len(track.loc[iLPS].datetime)-1,1):            
            track_lon, track_lat, track_yr, track_mm, track_dd, track_hh, j,k,l = track.loc[iLPS].values[it]
            track_lon0, track_lat0, track_yr0, track_mm0, track_dd0, track_hh0, j,k,l = track.loc[iLPS].values[it+1]
            travel_distance = travel_distance + \
                            distance_on_unit_sphere(track_lat,track_lon,track_lat0,track_lon0)
        if travel_distance<1000.00:continue

        # -----------------------------------
        # 3) Loop through each time step
        oro = 0
        contour = None 
        
        for it in np.arange(0,len(track.loc[iLPS].datetime),1):
            if (contour!=None):break
                
            distance = 99999.9
            slp_min = -1000.0        
            exp = None
              
            track_lon, track_lat, track_yr, track_mm, track_dd, track_hh, j,k,l = track.loc[iLPS].values[it]

            # ---------------------------------
            # 1) check monsoon domain and orography
            monsoon_domain = monsoon.sel(lat=track_lat,lon=track_lon,method='nearest')
            if monsoon_domain==0:
                continue

            orography_of_vortex = orography.sel(lat=track_lat,lon=track_lon,method='nearest')
            if orography_of_vortex==1:
                oro = oro+1
                if oro<=(len(track.loc[iLPS].datetime)/10.0):
                    continue
                else:
                    contour = None
                    break

            # ---------------------------------
            # 2) read in the SLPA data (relative to 21-day running average)
            time_step = str(track_yr)+"-"+str(track_mm).zfill(2)+"-"+str(track_dd).zfill(2)+"T"+str(track_hh)+":00:00"      

            file = "/home/yujia/Data/ERA-Interim/6hourly/slp/0.75/r_21d/T63/"+str(track_yr)+".Asia.nc"
            ds = xr.open_dataset(file)
            
            slp = ds.slp.sel(time=time_step)
            slp = slp[::-1,:]/100.0
            slp_iday = slp.load().data

            lat = slp.lat
            lon = slp.lon
            lon2d, lat2d = np.meshgrid(lon,lat)

            ilat = np.abs(lat-track_lat).argmin()
            ilon = np.abs(lon-track_lon).argmin()
            loc = [ilat.values,ilon.values]
            threshold = slp_iday[ilat,ilon]

            # --------------------------------------
            # 3) find a contour surrounding the vortex center and SLPA minimum within the contour <=-2hPa
            try:
                exp = list(rclv.find_contour_around_maximum(-slp_iday,loc,-threshold-0.05,\
                            max_footprint=784,periodic=(True,True)))
                area0, hull0, cd0 = rclv.contour_area(exp[0])
                if area0>400:continue         
            except ValueError as err:
                continue
            
            # --------------------------------------
            # is_inside2: double check if the vortex center falls within the contour
            is_inside2 = rclv.point_in_contour(exp[0],[exp[5],exp[4]])
            # --------------------------------------
            # xy_min: SLPA local min <= -?1? hPa (can be multiple minimums)
            xy_min = peak_local_max(exp[1],min_distance=1,threshold_abs=1,exclude_border=0,indices=True)

            # --------------------------------------
            # find the largest slpa minimum within the contour (if multiple SLPA extrame present)
            for i in xy_min:
                is_inside1 = rclv.point_in_contour(exp[0],i)
                if (is_inside1)&(is_inside2)&(slp_min<=exp[1][i[0],i[1]]):
                    slp_min     = exp[1][i[0],i[1]]
                    slp_min_loc = i
                    slp_min_lat = lat[ilat - (exp[4] - i[0])].values
                    slp_min_lon = lon[ilon - (exp[5] - i[1])].values
                    distance    = distance_on_unit_sphere(track_lat,track_lon,slp_min_lat,slp_min_lon)

            # --------------------------------------
            # 4) the distance between vortex center and SLPA min must be less than 500km
            if (distance<=500.0): 
                ilat = np.abs(lat-slp_min_lat).argmin()
                ilon = np.abs(lon-slp_min_lon).argmin()
                loc  = [ilat.values,ilon.values]
                localmin = slp_iday[ilat,ilon]
                orography_of_slp = orography.sel(lat=slp_min_lat,lon=slp_min_lon,method='nearest')
                if orography_of_slp==1:continue

                track.loc[iLPS].intensity.values[it] = -localmin
                
                # find convex contour -> DEFINE the boundary of cyclones
                upper_bound = -localmin
                low_bound = np.amax([-localmin+threshold,1.0]) 
                
                max_cd   = 0.05
                min_area = 0.0
                min_grad = 0.0
                exp      = None
                # find convex contour with the largest mean slpa gradient within the contour
                for deltaP in np.arange(low_bound,upper_bound+0.1,0.1):
                    try:
                        exp_tmp = list(rclv.find_contour_around_maximum(-slp_iday,loc,-localmin-deltaP,\
                                       max_footprint=400,periodic=(True,True)))

                        area, hull, cd = rclv.contour_area(exp_tmp[0])
                        xy_min         = np.amax(exp_tmp[1])
                        
                        r_mask = np.zeros_like(exp_tmp[1], dtype='bool')
                        r_mask[np.round(exp_tmp[0][:, 0]).astype('int'), np.round(exp_tmp[0][:, 1]).astype('int')] = 1
                        r_mask = ndimage.binary_fill_holes(r_mask)
                        r_mask = ~r_mask
                        exp_tmp[1][r_mask] = None

                        gradient = np.gradient(exp_tmp[1]) 
                        mag_grad = np.sqrt(gradient[0]**2+gradient[1]**2)
                        ave_grad = np.nanmean(mag_grad)
                       
                        if (cd>max_cd)&(exp!=None):
                            track.loc[iLPS].intensity.values[it] = deltaP
                            break
                        if (-localmin==xy_min)&(cd<=max_cd)&(area>=min_area)&(ave_grad>=min_grad):
                            exp = exp_tmp
                            min_grad = ave_grad
                            track.loc[iLPS].intensity.values[it] = ave_grad
                        #if (-localmin==xy_min)&(cd<=max_cd)&(area>=min_area):    
                        #    track.loc[iLPS].intensity.values[it] = deltaP
                            
                    except ValueError as err:
                        continue
                        
                # --------------------------------------
                # 5) if find a contour satisfying the criteria -> compute the averaged gradient within the convex contour
                if exp!=None:  
                    area, hull, cd = rclv.contour_area(exp[0])

                    r_mask = np.zeros_like(exp[1], dtype='bool')
                    r_mask[np.round(exp[0][:, 0]).astype('int'), np.round(exp[0][:, 1]).astype('int')] = 1
                    r_mask = ndimage.binary_fill_holes(r_mask)

                    # --------------------------------------
                    r_mask = ~r_mask
                    exp[1][r_mask] = None
                    gradient = np.gradient(exp[1]) 
                    mag_grad = np.sqrt(gradient[0]**2+gradient[1]**2)
                    ave_grad = np.nanmean(mag_grad)
                    
                    # --------------------------------------
                    # 6) averaged gradient >= 2hPa/6degree (0.25)
                    if (ave_grad>=0.30):
                        contour = exp

                        if plot:
                            area, hull, cd = rclv.contour_area(contour[0])

                            fig,ax = plt.subplots(figsize=(8,3),ncols=2)
                            plot = ax[0].pcolormesh(contour[1],cmap='Spectral')
                            ax[0].plot(contour[5],contour[4],marker="X",color="red",markersize=20.0)
                            ax[0].plot(contour[0][:,1], contour[0][:,0],'r',linewidth=1.5)

                            ax[1].pcolormesh(mag_grad,cmap='viridis')
                            ax[1].plot(contour[0][:,1], contour[0][:,0],'r',linewidth=1.5)
                            ax[1].plot(contour[5],contour[4],marker="X",color="red",markersize=20.0)

        # ------------------------------------------
        if (contour!= None): 
            reserve.append([iLPS])

            track.loc[iLPS].closed = 1
            
            if plot:
                fig, ax = plt.subplots(figsize=(4,3),subplot_kw={'projection': ccrs.PlateCarree(central_longitude=150)})
                extent = [55, 185, -5, 55]    
                lat = track.loc[iLPS].lat.copy()
                lon = track.loc[iLPS].lon.copy()
                if lon.max()-lon.min()>180:
                    transform = np.where(lon>180)[0]
                    lon.iloc[transform] = lon.iloc[transform]-360
                    track_istorm = sgeom.LineString(zip(lon.values, lat.values))
                    ax.add_geometries([track_istorm], ccrs.PlateCarree(),facecolor='none',edgecolor='r',linewidth=1.5)
                else:
                    track_istorm = sgeom.LineString(zip(lon.values, lat.values))
                    ax.add_geometries([track_istorm], ccrs.PlateCarree(),facecolor='none',edgecolor='r',linewidth=1.5)

                    ax.coastlines()
                    ax.stock_img()
                    ax.set_extent(extent)
                    ax.set_title(str(iLPS))


    print("----------- year "+str(year)+" END -----------")
    return reserve

# ===========
f_track   = '../post_processing/trx_Asian.dat'
f_genesis = '../post_processing/geninits_Asian.dat'
f_out     = './trx_Asian.IC.dat'

# -----------------------------------
# read in tracks
track = pd.read_csv(f_track,sep=' ',\
                    header=None,names=['lon','lat','year','month','day','hour'],\
                    usecols=np.arange(0,6))
genesis = pd.read_csv(f_genesis,sep=' ',header=None,names=['nrow'])
track['datetime'] = pd.to_datetime(track['year']*1000000+track['month']*10000+track['day']*100+track['hour'],format='%Y%m%d%H')
track['closed'] = track['datetime']
track['closed'] = -1
track['intensity'] = track['lat']
track['intensity'] = -999.9

# Generate Key for each LPS using genesis information
LPS_id = np.zeros(track['datetime'].shape,dtype=int)
for i in range(genesis.size):
    if i<=genesis.size-2:
        idStr = genesis.values[i][0]-1
        idEnd = genesis.values[i+1][0]-2
        LPS_id[idStr:idEnd+1] = i+1
    elif i==genesis.size-1:
        LPS_id[idStr:] = i+1
        
# Set LPS_id to be the index, hence the single key represents a LPS
track['id'] = LPS_id
track = track.set_index('id')

# ---------------------------
# Orography
orography = xr.open_dataset("/home/yujia/Floater/Orography/orography.mask.850hPa.nc").orography_mask
orography = orography.rename({'longitude': 'lon', 'latitude': 'lat'})

# ---------------------------
# monsoon domain
monsoon = xr.open_dataset("/home/yujia/Floater/MonsoonDomain/MonsoonDomain.nc").monsoon
bdry = np.loadtxt("/home/yujia/Floater/MonsoonDomain/MonsoonBoundary.dat")
bdry_lon = bdry[:,0]
bdry_lat = bdry[:,1]

# ---------------------------
lsm = xr.open_dataset("/home/yujia/Floater/MonsoonDomain/land_sea_mask.1.0.nc").lsm.sel(time="1989-01-01")
lsm = lsm.rename({'longitude': 'lon', 'latitude': 'lat'})

# ---------------------------
# call filtering 
for yr in range(1979,2019):
    output = filtering(track,yr)
    print(yr,len(output))
    print("----")
    print(output)
    #for i in output:
    #    track.loc[i[0]].closed = 1

# ---------------------------
# output
track = track.drop(columns=['datetime'])
track.to_csv(f_out, sep='\t',index=False,header=False)

# ---------------------------
track_PostFiltering = track[track.closed==1]
track_PostFiltering = track_PostFiltering.drop(columns=['closed'])
track_PostFiltering.to_excel("track_PostFiltering.IC.xlsx")
