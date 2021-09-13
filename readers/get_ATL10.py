import numpy as np
import pandas as pd
import h5py
import datetime as dt

# Adapted from a notebook by Tyler Sutterly 6/14/2910

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

## Read h5 ATL10 files ========================================================
def get_ATL10data(fileT, maxFreeboard, bounding_box, beamlist=None):
    # Pandas/numpy ATL10 reader
        
    f1 = h5py.File(fileT, 'r')

    orient = f1['orbit_info']['sc_orient'][:]  # orientation - 0: backward, 1: forward, 2: transition
    
    if len(orient) > 1:
        print('Transitioning, do not use for science!')
        return [[] for i in beamlist]
    elif (orient == 0):
        beams=['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']                
    elif (orient == 1):
        beams=['gt3r', 'gt3l', 'gt2r', 'gt2l', 'gt1r', 'gt1l']
    # (strong, weak, strong, weak, strong, weak)
    # (beam1, beam2, beam3, beam4, beam5, beam6)

    if beamlist == None:
        beams = [ beams[i] for i in [0, 2, 4]]
    else:
        beams = [ beams[i] for i in beamlist ]
    # use only strong beams

    dL = []

    for beam in beams:
        if beam in list(f1.keys()):
            freeboard=f1[beam]['freeboard_beam_segment']['beam_freeboard']['beam_fb_height'][:]

            freeboard_confidence=f1[beam]['freeboard_beam_segment']['beam_freeboard']['beam_fb_confidence'][:]
            freeboard_quality=f1[beam]['freeboard_beam_segment']['beam_freeboard']['beam_fb_quality_flag'][:]
            atlas_epoch=f1['/ancillary_data/atlas_sdp_gps_epoch'][0]

            # Delta time in gps seconds
            delta_time = f1[beam]['freeboard_beam_segment']['beam_freeboard']['delta_time'][:]
            temp = convert_time(delta_time + atlas_epoch)

            # Height segment ID (10 km segments)
            height_segment_id=f1[beam]['freeboard_beam_segment']['beam_freeboard']['height_segment_id'][:]

            lons=f1[beam]['freeboard_beam_segment']['beam_freeboard']['longitude'][:]
            lats=f1[beam]['freeboard_beam_segment']['beam_freeboard']['latitude'][:]

            dF = pd.DataFrame({'freeboard':freeboard, 'beam':beam, 
                               'lon':lons, 'lat':lats, 'delta_time':delta_time, 
                               'height_segment_id':height_segment_id
                              })
            dF = pd.concat([dF, temp], axis=1)
            
        else:
            dF = pd.DataFrame(columns=['freeboard','beam','lon','lat','delta_time',
                            'height_segment_id', 'time', 'year', 'month', 'day',
                            'hour', 'minute', 'second'])

        if len(dF) > 0:

            #dF['months'] = pd.Series(months, index=dF.index)
            dF = dF[(dF['freeboard']>0)]
            dF = dF[(dF['freeboard']<maxFreeboard)]
            dF = dF[(dF['lat']>=bounding_box[1])]
            dF = dF[(dF['lat']<=bounding_box[3])]
            # dF = dF[(dF['lon']>=bounding_box[0]) | (dF['lon']<=bounding_box[2])]
            dF = dF[(dF['lon']>=bounding_box[0])]
            dF = dF[(dF['lon']<=bounding_box[2])]

            # Reset row indexing
            dF=dF.reset_index(drop=True)

            dL.append(dF)
        else:
            dL.append([])              
        
    return dL

def get_ATL10lead(fileT, maxFreeboard, bounding_box, beamlist=None):
    # Pandas/numpy ATL10 reader
        
    f1 = h5py.File(fileT, 'r')

    orient = f1['orbit_info']['sc_orient']  # orientation - 0: backward, 1: forward, 2: transition

    if orient == 1: # forward
        beams = ['gt3r', 'gt3l', 'gt2r', 'gt2l', 'gt1r', 'gt1l']
    else: # backward
        beams = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
    # (strong, weak, strong, weak, strong, weak)
    # (beam1, beam2, beam3, beam4, beam5, beam6)

    if beamlist == None:
        beams = [ beams[i] for i in [0, 2, 4]]
    else:
        beams = [ beams[i] for i in beamlist ]
    # use only strong beams

    dL = []

    for beam in beams:
        if beam in list(f1.keys()):
            lead_height=f1[beam]['leads']['lead_height'][:]
            lead_length=f1[beam]['leads']['lead_length'][:]
            lead_sigma=f1[beam]['leads']['lead_sigma'][:]

            atlas_epoch=f1['/ancillary_data/atlas_sdp_gps_epoch'][0]

            # Delta time in gps seconds
            delta_time = f1[beam]['leads']['delta_time'][:]
            temp = convert_time(delta_time + atlas_epoch)

            # Height segment ID (10 km segments)
            ssh_n=f1[beam]['leads']['ssh_n'][:]
            ssh_ndx=f1[beam]['leads']['ssh_ndx'][:]

            lons=f1[beam]['leads']['longitude'][:]
            lats=f1[beam]['leads']['latitude'][:]

            dF = pd.DataFrame({'beam':beam, 'height':lead_height, 'length': lead_length, 'sigma': lead_sigma,
                               'lon':lons, 'lat':lats, 'delta_time':delta_time, 
                               'ssh_n':ssh_n, 'ssh_ndx': ssh_ndx
                              })
            dF = pd.concat([dF, temp], axis=1)
            
        else:
            dF = pd.DataFrame(columns=['beam', 'height','length', 'sigma', 'delta_time',
                            'ssh_n', 'ssh_ndx', 'year', 'month', 'day',
                            'hour', 'minute', 'second'])

        if len(dF) > 0:

            dF = dF[(dF['lat']>=bounding_box[1])]
            dF = dF[(dF['lat']<=bounding_box[3])]
            # dF = dF[(dF['lon']>=bounding_box[0]) | (dF['lon']<=bounding_box[2])]
            dF = dF[(dF['lon']>=bounding_box[0])]
            dF = dF[(dF['lon']<=bounding_box[2])]

            # Reset row indexing
            dF=dF.reset_index(drop=True)

            dL.append(dF)
        else:
            dL.append([])              
        
    return dL
