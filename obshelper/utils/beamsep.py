#!/usr/bin/env python
# coding: utf-8

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

# no rotation beam position
Xc = [0, 5.74, 2.86, -2.87, -5.74, -2.87, 2.88, 11.5, 8.6, 5.72, -0.00866, -5.74, -8.61, -11.5, -8.61, -5.74, 0.00869, 5.76, 8.62]
Yc = [0, -0.00783, -4.97, -4.97, 0.00218, 4.97, 4.97, -0.0213, -4.98, -9.95, -9.94, -9.94, -4.97, -0.00129, 4.97, 9.94, 9.94, 9.93, 4.96]

def beam_pos(dec = 0, rotate_angle = 23.4):
    """
    Calculate beam positions after rotation.

    Parameters:
    dec (float): declination in degree (affects beam position calculation)
    rotate_angle (float or astropy Angle): rotation angle in degree or as an astropy Angle object

    Returns:
    tuple: two numpy arrays representing x and y beam positions in arcmin
    """
    if hasattr(rotate_angle, 'unit'):
        theta = rotate_angle.to(u.rad).value
    else:
        theta = np.deg2rad(rotate_angle)

    from math import sin, cos
    rot_mat = np.array([[cos(theta), -sin(theta)], 
                        [sin(theta), cos(theta)]])

    x, y = np.dot(rot_mat, np.vstack([Xc, Yc]))
    x /= cos(np.deg2rad(dec))
    
    return x, y

def beam_sep2center(center, beam):
    """
    Calculate the position of a specific beam relative to the center.

    Parameters:
    center (SkyCoord): center position
    beam (str): beam name (e.g., 'M01', 'M02', etc.)

    Returns:
    SkyCoord: new position after moving to the beam's position
    """
    print(f"M01 centered on {center.to_string('hmsdms')}")
    
    keys = np.arange(1,20).astype(str)

    beam_sep = {}
    for i in range(19):
        beam_sep["M"+ keys[i].zfill(2)] = [Xc[i] - Xc[0], Yc[i] - Yc[0]]

    ra = center.ra + beam_sep[beam][0] * u.arcmin / np.cos(center.dec)
    dec = center.dec + beam_sep[beam][1] * u.arcmin
    
   
    coff = SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg))
    print(f"Moved M01 to {beam}'s position as source off on {coff.to_string('hmsdms')}")
    print("Remember the beam as source on is on the other side!")
    return coff