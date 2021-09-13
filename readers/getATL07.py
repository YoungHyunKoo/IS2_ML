import warnings
import h5py
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import h5py
import scipy
# from astropy.time import Time

def convert_time(delta_time):
    times = []
    years = []
    months = []
    days = []
    hours = []
    minutes = []
    seconds =[]
    
    for i in range(0, len(delta_time)):
        times.append(dt.datetime(1980, 1, 6) + dt.timedelta(seconds = delta_time[i]))
        years.append(times[i].year)
        months.append(times[i].month)
        days.append(times[i].day)
        hours.append(times[i].hour)
        minutes.append(times[i].minute)
        seconds.append(times[i].second)
    
    temp = pd.DataFrame({'time':times, 'year': years, 'month': months, 'day': days,
                         'hour': hours, 'minute': minutes, 'second': seconds
                        })
    return temp

# Function to read ATL03 data (.h5 format)
def getATL07(fname, beam_number, bbox, maxh = 1000):
    # 0, 2, 4 = Strong beam; 1, 3, 5 = weak beam
    # sfype: surface type (0=land, 1=ocean , 2=sea ice, 3=land ice, 4=inland water)
    
    f = h5py.File(fname, 'r')
    
    orient = f['orbit_info']['sc_orient'][:]  # orientation - 0: backward, 1: forward, 2: transition
    
    if len(orient) > 1:
        print('Transitioning, do not use for science!')
        return [[] for i in beamlist]
    elif (orient == 0):
        beams=['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']                
    elif (orient == 1):
        beams=['gt3r', 'gt3l', 'gt2r', 'gt2l', 'gt1r', 'gt1l']
    # (strong, weak, strong, weak, strong, weak)
    # (beam1, beam2, beam3, beam4, beam5, beam6)

    beam = beams[beam_number]    
    
    # height of each received photon, relative to the WGS-84 ellipsoid
    # (with some, not all corrections applied, see background info above)
    heights=f[beam]['sea_ice_segments']['heights']['height_segment_height'][:]
    # latitude (decimal degrees) of each received photon
    lats=f[beam]['sea_ice_segments']['latitude'][:]
    # longitude (decimal degrees) of each received photon
    lons=f[beam]['sea_ice_segments']['longitude'][:]
    # longitude (decimal degrees) of each received photon
    x_atc=f[beam]['sea_ice_segments']['seg_dist_x'][:]
    # seconds from ATLAS Standard Data Product Epoch. use the epoch parameter to convert to gps time
    deltatime=f[beam]['sea_ice_segments']['delta_time'][:]
    
    if bbox != None:
        if bbox[0] < bbox[2]:
            valid = (lats >= bbox[1]) & (lats <= bbox[3]) & (lons >= bbox[0]) & (lons <= bbox[2])
        else:
            valid = (lats >= bbox[1]) & (lats <= bbox[3]) & ((lons >= bbox[0]) | (lons <= bbox[2]))
    
    if len(heights[valid]) == 0:
        return pd.DataFrame({})
    
    else:    
        # width of best fit gaussian
        width=f[beam]['sea_ice_segments']['heights']['height_segment_w_gaussian'][:]
        # RMS difference between sea ice modeled and observed photon height distribution
        rms=f[beam]['sea_ice_segments']['heights']['height_segment_w_gaussian'][:]
        # number of laser pulses
        n_pulse=f[beam]['sea_ice_segments']['heights']['height_segment_n_pulse_seg'][:]
        
        # Length of each height segment
        seg_len=f[beam]['sea_ice_segments']['heights']['height_segment_length_seg'][:]        
        # Beam azimuth
        b_azi = f[beam]['sea_ice_segments/geolocation/beam_azimuth'][:]
        # Beam elevation
        b_ele = f[beam]['sea_ice_segments/geolocation/beam_coelev'][:]
        # Solar azimuth
        s_azi = f[beam]['sea_ice_segments/geolocation/solar_azimuth'][:]
        # Solar elevation
        s_ele = f[beam]['sea_ice_segments/geolocation/solar_elevation'][:]
        
        # Default surface type of ATL07 product
        stype=f[beam]['sea_ice_segments']['heights']['height_segment_type'][:]

        # Calculated background count rate based on sun angle, surface slope, unit reflectance
        bck_cal=f[beam]['sea_ice_segments']['stats']['backgr_calc'][:]
        # Background count rate, averaged over the segment based on 200 hz atmosphere
        bck_r200=f[beam]['sea_ice_segments']['stats']['backgr_r_200'][:]
        # Background count rate, averaged over the segment based on 25 hz atmosphere
        bck_r25=f[beam]['sea_ice_segments']['stats']['backgr_r_25'][:]

        # Background rate normalized to a fixed solar elevation angle
        bck_norm=f[beam]['sea_ice_segments']['stats']['background_r_norm'][:]

        # Segment histogram width estimate
        hist_w=f[beam]['sea_ice_segments']['stats']['hist_w'][:]
        # photon count rate, averaged over segment
        photon_rate=f[beam]['sea_ice_segments']['stats']['photon_rate'][:]
        
        # mean height of histogram
        h_mean=f[beam]['sea_ice_segments']['stats']['hist_mean_h'][:]
        # median height of histogram
        h_median=f[beam]['sea_ice_segments']['stats']['hist_median_h'][:]
        # Difference between mean and median
        h_diff=h_mean - h_median
        
        # Sea ice concentration
        sic=f[beam]['sea_ice_segments']['stats']['ice_conc'][:]

        # Delta time to gps seconds
        atlas_epoch=f[beam]['/ancillary_data/atlas_sdp_gps_epoch'][0]
        temp = convert_time(deltatime + atlas_epoch)    

        df07=pd.DataFrame({'beam': beam, 'lat':lats, 'lon':lons, 'x': x_atc, 'deltatime':deltatime,
                           'height':heights, 'h_mean': h_mean, 'h_median': h_median, 'h_diff': h_diff,
                           'width': width, 'rms': rms, 'n_pulse': n_pulse,
                           'bck_cal': bck_cal, 'bck_r200': bck_r200, 'bck_r25': bck_r25,
                           'bck_norm': bck_norm, 'hist_w': hist_w, 'ph_rate': photon_rate,
                           'sic': sic, 'stype': stype, 'seg_len': seg_len,
                           'b_azi': b_azi, 'b_ele': b_ele, 's_azi': s_azi, 's_ele': s_ele
                          })

        # Concatenate ATL03 dataframe and time dataframe
        df07 = pd.concat([df07, temp], axis=1).reset_index(drop = True)
        
        df07 = df07[df07['height'] <= maxh]

        if bbox != None:
            df07 = df07[df07['lat'] >= bbox[1]].reset_index(drop = True)
            df07 = df07[df07['lat'] <= bbox[3]].reset_index(drop = True)
            if bbox[0] < bbox[2]:
                df07 = df07[df07['lon'] >= bbox[0]].reset_index(drop = True)
                df07 = df07[df07['lon'] <= bbox[2]].reset_index(drop = True)
            else:
                df07 = df07[(df07['lon'] >= bbox[0]) | (df07['lon'] <= bbox[2])].reset_index(drop = True)

        return df07
    
    