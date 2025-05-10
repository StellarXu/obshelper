#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import ephem as ep
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import Angle

from .utils.util import line_set

def bj_time(s):
    """
    Convert BJ time string to Time object with 8 hours added.
    
    Parameters:
    s (str): BJ time string in format 'YYYY-MM-DD HH:MM:SS' or 'YYYY/MM/DD HH:MM:SS'
    
    Returns:
    Time: Converted Time object
    """
    t = Time(str(s).replace('/','-'))
    return t + 8*u.hour

def ephem_coord(coo):
    """
    Convert coordinate string to ephem format.
    
    Parameters:
    coo (str): Coordinate string in format 'hh:mm:ss.ss dd:mm:ss.ss'
    
    Returns:
    list: Coordinate in ephem format
    """
    coord = coo.split(' ')
    for s in ['h','m','s','d']:
        coord[0] = coord[0].replace(s,' ') 
        coord[1] = coord[1].replace(s,' ') 
    return coord

def plot_za(BJ_time, coord):
    """
    Plot zenith angle of an object over a day.
    
    Parameters:
    BJ_time (str): BJ time string in format 'YYYY-MM-DD HH:MM:SS'
    coord (str): Coordinate string in format 'hh:mm:ss.ss dd:mm:ss.ss'
    
    Returns:
    tuple: Minimum zenith angle and corresponding time
    """
    t0 = Time(BJ_time)
    ts = t0 + np.arange(-12, 12, 0.1)*u.hour

    obj = ep.FixedBody()
    obj._ra, obj._dec = ephem_coord(coord)
    obj._epoch = '2000'

    zas = []
    for i in range(len(ts)):
        t = ts[i]
        ut = t - 8*u.hour

        fast = ep.Observer()
        fast.lon = '106.8566667'    # East positive, deg
        fast.lat = '25.65294444'    # North positive, deg
        fast.elevation = 1110.0288      # Altitude, m 
        fast.date = ut.value

        obj.compute(fast)

        alt = Angle(str(obj.alt),unit = u.deg)
        za = 90*u.deg - alt
        zas += [za.value,]
        
        if np.abs((t-t0).to_value(u.hour)) < 0.01:
            print("#### Object ####")
            print("Previous Transit time:",bj_time(fast.previous_transit(obj)))
            print("Next Transit time:",bj_time(fast.next_transit(obj)))
            print("#### Sun ####")
            sun_rise = bj_time(fast.next_rising(ep.Sun()))
            sun_set = bj_time(fast.previous_setting(ep.Sun()))
            # at night
            if np.abs((sun_rise - t0).to_value(u.hour)) < 12:
                print("Next Rise time:",sun_rise)
                print("Previous Set time:",sun_set)
            else:
                sun_rise = bj_time(fast.previous_rising(ep.Sun()))
                sun_set = bj_time(fast.next_setting(ep.Sun()))
                print("Previous Rise time:",sun_rise)
                print("Next Set time:",sun_set)
    zas = np.array(zas)

    unrise = (zas > 90)
    zas[unrise] = np.nan

    fig, ax = plt.subplots(figsize = (10,5))
    ax.plot_date(ts.plot_date, zas, fmt = '-', lw = 3, label = f"obj:{str(obj._ra)}\n    {str(obj._dec)}")  
    plt.gcf().autofmt_xdate()  # orient date labels at a slant  
    [ax.axhline(y = ii,c='C1',ls='--') for ii in [20,30,40]]
    ax.axhline(y = 26.4, c='C3',ls='--')
    ax.invert_yaxis()
    ax.grid()
    
    if sun_rise > sun_set:
        s1 = sun_set
        s2 = sun_rise
        at_night = True
    else:
        s1 = sun_rise
        s2 = sun_set
        at_night = False
    
    print("input time at night?:", at_night)
    
    if at_night:  
        ax.fill_betweenx(y = [np.nanmin(zas),np.nanmax(zas)], x1 = s1.plot_date,x2 = s2.plot_date,alpha=.1,label='night',color = 'k')
    else:
        ax.fill_betweenx(y = [np.nanmin(zas),np.nanmax(zas)], x1 = ts[0].plot_date,x2 = s1.plot_date,alpha=.1,label='night',color = 'k')
        ax.fill_betweenx(y = [np.nanmin(zas),np.nanmax(zas)], x1 = s2.plot_date,x2 = ts[-1].plot_date,alpha=.1,color = 'k')
    line_set(ax, xlabel = 'time (date-hour)', ylabel = 'ZA (deg)',)
    plt.show()
    
    za_min = np.nanmin(zas)
    ts_min = ts[zas == za_min]
    
    return za_min, ts_min


def plot_za_in_1day(self, BJ_time, sample_time, gap = 20):
    """
    Plot zenith angle of an object in one day.
    
    Parameters:
    self: Instance of the class
    BJ_time (str): BJ time string in format 'YYYY/MM/DD HH:MM:SS'
    sample_time (float): Sample time in seconds
    gap (int, optional): Gap between points. Default is 20
    
    Returns:
    tuple: Minimum zenith angle and corresponding time
    """
    ra = self.ra1
    dec = self.dec1
    
    coo = SkyCoord(ra, dec, unit = (u.deg, u.deg)).to_string('hmsdms')
    coord = [ephem_coord(coo[j]) for j in range(len(coo))]
    
    t0 = Time(BJ_time)
    ts = t0 + np.arange(0, len(ra)) * sample_time * u.second

    zas = []
    for i in np.arange(0, len(ra), gap):
        obj = ep.FixedBody()
        obj._ra, obj._dec = coord[i]
        obj._epoch = '2000'
    
        t = ts[i]
        ut = t - 8*u.hour

        fast = ep.Observer()
        fast.lon = '106.8566667'    # East positive, deg
        fast.lat = '25.65294444'    # North positive, deg
        fast.elevation = 1110.0288      # Altitude, m 
        fast.date = ut.value

        obj.compute(fast)
        alt = Angle(str(obj.alt),unit = u.deg)

        za = 90*u.deg - alt
        zas += [za.value,]
        
    zas = np.array(zas)
    
    unrise = (zas > 90)
    zas[unrise] = np.nan

    ts = ts[::gap]
    fig, ax = plt.subplots(figsize = (10,5))
    ax.plot_date(ts.plot_date, zas, fmt = '-', lw = 3, label = f"obj:{str(obj._ra)}\n    {str(obj._dec)}")  
    plt.gcf().autofmt_xdate()  # orient date labels at a slant  
    [ax.axhline(y = ii,c='C1',ls='--') for ii in [20,30,40]]
    ax.axhline(y = 26.4, c='C3',ls='--')
    ax.invert_yaxis()
    ax.grid()
    
    line_set(ax, xlabel = 'time (date-hour)', ylabel = 'ZA (deg)',)
    plt.show()
    
    za_min = np.nanmin(zas)
    ts_min = ts[zas == za_min]
    
    return za_min, ts_min


def plot_sun_AngleDist(BJ_time, coord):
    """
    Plot angular distance between an object and the sun over a year.
    
    Parameters:
    BJ_time (str): BJ time string in format 'YYYY/MM/DD HH:MM:SS'
    coord (str): Coordinate string in format 'hh:mm:ss.ss dd:mm:ss.ss'
    """
    t0 = Time(BJ_time)
    ts = t0 + np.arange(0, 365, 1) * u.day

    obj = ep.FixedBody()
    obj._ra, obj._dec = ephem_coord(coord)
    obj._epoch = '2000'
    
    angs = []
    for i in range(len(ts)):
        t = ts[i]
        ut = t - 8*u.hour

        fast = ep.Observer()
        fast.lon = '106.8566667'    # East positive, deg
        fast.lat = '25.65294444'    # North positive, deg
        fast.elevation = 1110.0288      # Altitude, m 
        fast.date = ut.value

        obj.compute(fast)

        sun = ep.Sun(ut.value)

        angle = ep.separation(obj, sun)
        angs += [np.rad2deg(angle),]

    angs = np.array(angs)
    
    fig, ax = plt.subplots(figsize = (10,5))
    ax.plot_date(ts.plot_date, angs, fmt = '-', lw = 3, label = f"obj:{str(obj._ra)}\n    {str(obj._dec)}")  
    plt.gcf().autofmt_xdate()  # orient date labels at a slant  
    ax.grid()

    line_set(ax, xlabel = 'time (year-month)', ylabel = 'Angular Distance to Sun (deg)',)
    plt.show()
    
#     return ts, angs