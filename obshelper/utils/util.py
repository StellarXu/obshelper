#!/usr/bin/env python
# coding: utf-8

import numpy as np
from astropy import  log

restfreq = 1420.405751768 # Mhz
c = 299792.458 # km/s


def area(ra1, ra2, dec1, dec2, cross_zero = False):
    """
    Calculate the area on the celestial sphere given a range of right ascension and declination.

    Args:
        ra1 (float): Starting right ascension in degrees.
        ra2 (float): Ending right ascension in degrees.
        dec1 (float): Starting declination in degrees.
        dec2 (float): Ending declination in degrees.
        cross_zero (bool, optional): Whether the right ascension crosses zero. Defaults to False.

    Returns:
        astropy.units.Quantity: Area in square degrees.
    """
    S2 = 2 * np.pi * (1 - np.sin(np.deg2rad(dec2)))
    S1 = 2 * np.pi * (1 - np.sin(np.deg2rad(dec1)))
    S = S1 - S2
    if not cross_zero:
        S *= (ra2 - ra1) / 360
    else:
        S *= (ra1 + 360 - ra2) / 360
    from astropy import units as u
    return (S * u.rad**2).to(u.deg**2)

def vopt2vrad(vopt):
    """
    Convert optical velocity to radio velocity.

    Args:
        vopt (float): Optical velocity in km/s.

    Returns:
        float: Radio velocity in km/s.
    """
    vrad = c - c**2 / (c + vopt)
    return vrad

def vrad2vopt(vrad):
    """
    Convert radio velocity to optical velocity.

    Args:
        vrad (float): Radio velocity in km/s.

    Returns:
        float: Optical velocity in km/s.
    """
    vopt = c**2 / (c - vrad) - c
    return vopt

def vopt2freq(vopt):
    """
    Convert optical velocity to frequency.

    Args:
        vopt (float): Optical velocity in km/s.

    Returns:
        float: Frequency in MHz.
    """
    freq = restfreq / (vopt / c + 1)
    
    return freq

def vrad2freq(vrad):
    """
    Convert radio velocity to frequency.

    Args:
        vrad (float): Radio velocity in km/s.

    Returns:
        float: Frequency in MHz.
    """
    freq = restfreq * (1 - vrad / c)
    
    return freq


def freq2vopt(freq):
    """
    Convert frequency to optical velocity.

    Args:
        freq (float): Frequency in MHz.

    Returns:
        float: Optical velocity in km/s.
    """
    vopt = (restfreq / freq - 1) * c
    
    return vopt

def redshift(v , relative = False):
    """
    Calculate redshift given velocity.

    Args:
        v (float): Velocity in km/s.
        relative (bool, optional): Whether to use relativistic formula. Defaults to False.

    Returns:
        float: Redshift.
    """
    beta = v / c
    
    if relative:
        g = 1 / np.sqrt(1 - beta ** 2)
        z = (1 + v / c) * g - 1
        return z
    else:
        return beta
    
def redshift2vel(z, relative = False):
    """
    Convert redshift to velocity.

    Args:
        z (float): Redshift.
        relative (bool, optional): Whether to use relativistic formula. Defaults to False.

    Returns:
        float: Velocity in km/s.
    """
    if not relative:
        return c * z
    else:
        z1 = z + 1
        return c * (z1 ** 2 - 1) / (z1 ** 2 + 1)
    
def line_set(ax,xlabel,ylabel,direction='in',xlim=None,ylim=None,legend=True,
             title=None,loc='best', size = None,frameon=True, bottom=True, top=True,left=True,right=True):
    """
    Set plot style and formatting.
    """
    import matplotlib as mpl            
    mpl.rcParams['figure.facecolor']=(1, 1, 1, 1)

    if size is None:
        size = {'mz': 1,   # Set thickness of the tick marks
                'lz': 3,   # Set length of the tick marks
                'lbz': 16,  # Set label size
                'tkz': 14,  # Set tick size
                }
    mz = size['mz']; lz = size['lz']
    lbz = size['lbz']; tkz = size['tkz']
    # Make tick lines thicker
    for l in ax.get_xticklines():
        l.set_markersize(lz)
        l.set_markeredgewidth(mz)
    for l in ax.get_yticklines():
        l.set_markersize(lz)
        l.set_markeredgewidth(mz)

    # Make figure box thicker
    for s in ax.spines.values():
        s.set_linewidth(mz)
    ax.minorticks_on()
    ax.tick_params("both",which = 'both',direction=direction,labelsize=tkz,)
#                   bottom=bottom, top=top,left=left,right=right)
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    ax.set_xlabel(xlabel,fontsize=lbz)
    ax.set_ylabel(ylabel,fontsize=lbz)
    if legend: ax.legend(loc = loc,fontsize=tkz, frameon=frameon)
    if title is not None: ax.set_title(title,fontsize=lbz)