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

# Function to calculate the photon rate (photons / shots(pulses))
def count_pid(pid, conf, conf_th = 3):
    idcount = np.zeros(len(pid))
    idcount2 = np.zeros(len(pid))
    firstid = pid[0]
    
    count = 0
    count2 = 0
    
    for i in range(0, len(pid)):
        if pid[i] == firstid:
            count += 1
            if conf[i] >= conf_th:
                count2 += 1
        else:
            idcount[i-count:i] = count
            idcount2[i-count:i] = count2
            firstid = pid[i]
            count = 1
            if conf[i] >= conf_th:
                count2 = 1
            else:
                count2 = 0
        
    return idcount, idcount2

# Function to calculate the background count and background rate
def cal_bckgrd(pid, frame_cnt, bck_pce0, bck_cnt0, bck_rate0):
    bck_cnt1 = np.zeros(len(pid))
    bck_rate1 = np.zeros(len(pid))

    frames, idx = np.unique(bck_pce0, return_index=True)
    
    for i in range(0, len(pid)):      
        
        idx0 = idx[np.argmin(abs(frames - frame_cnt[i]))]
    
        v = (pid[i]-1) // 50
        bck_cnt1[i] = bck_cnt0[idx0 + v]
        bck_rate1[i] = bck_rate0[idx0 + v]

#     for frame in frames:
#         idx0 = np.where(frame_cnt == frame)[0][0]
#         pid_part = pid[frame_cnt == frame]
#         bcnt_part = bck_cnt[bck_pce == frame]
#         brate_part = bck_rate[bck_pce == frame]

#         n = 0
#         for j in range(0, 4):
#             idx1 = np.where((pid_part >= j*50+1) & (pid_part <= (j+1)*50))[0]
#             if len(idx1) > 0:
#                 bck_cnt1[idx0+idx1[0]:idx0+idx1[-1]+1] = bcnt_part[min(n, len(bcnt_part)-1)]
#                 bck_rate1[idx0+idx1[0]:idx0+idx1[-1]+1] = brate_part[min(n, len(bcnt_part)-1)]
#                 n+=1
                
    return bck_cnt1, bck_rate1

# Function to read ATL03 data (.h5 format)
def getATL03(fname, beam_number, bbox, stype = 2):
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
    
    try:
        # height of each received photon, relative to the WGS-84 ellipsoid
        # (with some, not all corrections applied, see background info above)
        heights=f[beam]['heights']['h_ph'][:]
        # latitude (decimal degrees) of each received photon
        lats=f[beam]['heights']['lat_ph'][:]
        # longitude (decimal degrees) of each received photon
        lons=f[beam]['heights']['lon_ph'][:]
        # seconds from ATLAS Standard Data Product Epoch. use the epoch parameter to convert to gps time
        deltatime=f[beam]['heights']['delta_time'][:]
    except:
        return pd.DataFrame({})   
    
    if bbox != None:
        if bbox[0] < bbox[2]:
            valid = (lats >= bbox[1]) & (lats <= bbox[3]) & (lons >= bbox[0]) & (lons <= bbox[2])
        else:
            valid = (lats >= bbox[1]) & (lats <= bbox[3]) & ((lons >= bbox[0]) | (lons <= bbox[2]))
    
    if len(heights[valid]) == 0:
        return pd.DataFrame({})    
    else:
        # Surface types for signal classification confidence
        # 0=Land; 1=Ocean; 2=SeaIce; 3=LandIce; 4=InlandWater    
        conf=f[beam]['heights']['signal_conf_ph'][:, stype] #choose column 2 for confidence of sea ice photons
        # confidence level associated with each photon event
        # -2: TEP
        # -1: Events not associated with a specific surface type
        #  0: noise
        #  1: buffer but algorithm classifies as background
        #  2: low
        #  3: medium
        #  4: high
        # number of ATL03 20m segments
        n_seg, = f[beam]['geolocation']['segment_id'].shape
        # first photon in the segment (convert to 0-based indexing)
        Segment_Index_begin = f[beam]['geolocation']['ph_index_beg'][:] - 1
        # number of photon events in the segment
        Segment_PE_count = f[beam]['geolocation']['segment_ph_cnt'][:]
        # along-track distance for each ATL03 segment
        Segment_Distance = f[beam]['geolocation']['segment_dist_x'][:]
        # along-track distance (x) for photon events
        x_atc = np.array(f[beam]['heights']['dist_ph_along'][:])
        # cross-track distance (y) for photon events
        y_atc = np.array(f[beam]['heights']['dist_ph_across'][:])
        # The photon event counter is part of photon ID and counts from 1 for each channel until reset by laser pulse counter
        pulse_id = f[beam]['heights']['ph_id_pulse'][:]
        frame_cnt = f[beam]['heights']['pce_mframe_cnt'][:]
        # pcount: number of the returned photon that has the same pulse id
        pulse_cnt, pulse_cnth = count_pid(pulse_id, conf)
        # pulse_cnt: photon rate (photons/pulse); pulse_cnth: photon rate for high confidence photons

        ph_id_count = f[beam]['heights']['ph_id_count'][:]

        # Calculate background factors
        bck_pce0 = f[beam]['bckgrd_atlas']['pce_mframe_cnt'][:]
        bck_cnt0 = f[beam]['bckgrd_atlas']['bckgrd_counts'][:]
        bck_rate0 = f[beam]['bckgrd_atlas']['bckgrd_rate'][:]    
        bck_cnt, bck_rate = cal_bckgrd(pulse_id, frame_cnt, bck_pce0, bck_cnt0, bck_rate0)
        del bck_pce0, bck_cnt0, bck_rate0

        # Remove the uneffective reference photons (no geo-correction parameters)
        mask_ind = (Segment_Index_begin >= 0)
        Segment_Index_begin = Segment_Index_begin[mask_ind]
        Segment_PE_count = Segment_PE_count[mask_ind]
        n_seg = len(Segment_PE_count)

        # Geographical correction parameters (refer to the ATL03 documents)
        seg_lat = f[beam]['geolocation/reference_photon_lat'][mask_ind]

        dac0 = f[beam]['geophys_corr/dac'][mask_ind]
        dac0[dac0 > 3e+38] = np.nan
        geoid0 = f[beam]['geophys_corr/geoid'][mask_ind]
        geoid0[geoid0 > 3e+38] = np.nan

        tide_earth = f[beam]['geophys_corr/tide_earth'][mask_ind]
        tide_earth[tide_earth > 3e+38] = np.nan
        tide_load = f[beam]['geophys_corr/tide_load'][mask_ind]
        tide_load[tide_load > 3e+38] = np.nan
        tide_oc = f[beam]['geophys_corr/tide_ocean'][mask_ind]
        tide_oc[tide_oc > 3e+38] = np.nan
        tide_pole = f[beam]['geophys_corr/tide_pole'][mask_ind]
        tide_pole[tide_pole > 3e+38] = np.nan
        tide_oc_pole = f[beam]['geophys_corr/tide_oc_pole'][mask_ind]
        tide_oc_pole[tide_oc_pole > 3e+38] = np.nan
        tide_eq = f[beam]['geophys_corr/tide_equilibrium'][mask_ind]
        tide_eq[tide_eq > 3e+38] = np.nan
        tide0 = tide_earth + tide_load + tide_oc + tide_pole + tide_oc_pole + tide_eq

        sea0 = f[beam]['geolocation/solar_elevation'][mask_ind] # Solar elevation angle
        sea0[sea0 > 3e+38] = np.nan
        saa0 = f[beam]['geolocation/solar_azimuth'][mask_ind] # Solar azimuth angle
        saa0[saa0 > 3e+38] = np.nan

        # Remove unaffective geo-correction values
        dac0 = dac0[tide0 != np.nan]
        geoid0 = geoid0[tide0 != np.nan]
        sea0 = sea0[tide0 != np.nan]
        saa0 = saa0[tide0 != np.nan]
        tide0 = tide0[tide0 != np.nan]
        seg_lat = seg_lat[tide0 != np.nan]

        # Since the number of reference points are less than the original photons,
        # reference heights of all photons should be interpolated from the existing reference points
        dac = np.interp(lats, seg_lat, dac0)
        geoid = np.interp(lats, seg_lat, geoid0)
        tide = np.interp(lats, seg_lat, tide0)
        sea = np.interp(lats, seg_lat, sea0)
        saa = np.interp(lats, seg_lat, saa0)

        del dac0, geoid0, tide0, sea0, saa0
        
        # Along track distance = x_atc
        for j in range(n_seg):
            # index for 20m segment j
            idx = Segment_Index_begin[j]
            # number of photons in 20m segment
            cnt = Segment_PE_count[j]
            # add segment distance to along-track coordinates
            x_atc[idx:idx+cnt] += Segment_Distance[j]   

        # Delta time to gps seconds
        atlas_epoch=f[beam]['/ancillary_data/atlas_sdp_gps_epoch'][0]
        temp = convert_time(deltatime + atlas_epoch)

        df03=pd.DataFrame({'beam': beam, 'lat':lats, 'lon':lons, 'x':x_atc, 'y':y_atc,
                           'height':heights, 'dac': dac, 'geoid': geoid, 'tide': tide,
                           's_azi': saa, 's_ele': sea, 
                           'deltatime':deltatime, 'conf':conf,
                           'pid': pulse_id, 'pcnt': pulse_cnt, 'pcnth': pulse_cnth,
                           'bcnt': bck_cnt, 'brate': bck_rate
                          })

        # Concatenate ATL03 dataframe and time dataframe
        df03 = pd.concat([df03, temp], axis=1).reset_index(drop = True)

        if bbox != None:
            df03 = df03[df03['lat'] >= bbox[1]].reset_index(drop = True)
            df03 = df03[df03['lat'] <= bbox[3]].reset_index(drop = True)
            if bbox[0] < bbox[2]:
                df03 = df03[df03['lon'] >= bbox[0]].reset_index(drop = True)
                df03 = df03[df03['lon'] <= bbox[2]].reset_index(drop = True)
            else:
                df03 = df03[(df03['lon'] >= bbox[0]) | (df03['lon'] <= bbox[2])].reset_index(drop = True)

        return df03