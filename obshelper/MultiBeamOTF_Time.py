import numpy as np
import matplotlib.pyplot as plt
import hiviewer as hv
import os
from astropy import log

from astropy import coordinates as coords
import astropy.units as u

# Function to format RA and Dec into string format
def format_radec(pos):
    ra1 = pos.ra.to_string(unit = u.hourangle, sep = ':', precision=2, pad = True)
    dec1 = pos.dec.to_string(unit = u.degree, sep = ':', precision=1, alwayssign = True)
    return ra1, dec1

class MultiOTFcalculator(object):
    """
    A class used to calculate OTF observation parameters and footprints
    """
    def __init__(self, ra, dec, dist = None, diff_ra = None, diff_dec = None, start = None, end = None):
        """
        Initialize the MultiOTFcalculator
        
        Parameters:
        ra, dec: center position in degree
        diff_ra: width in arcmin
        diff_dec: width in arcmin
        start, end: optional Coord
        """
        self.ra = ra
        self.dec = dec
        
        if diff_ra is not None and diff_dec is not None:
            self.diff_ra = diff_ra * u.arcmin
            self.diff_dec = diff_dec * u.arcmin
        self.start = start
        self.end = end
        self.center = coords.SkyCoord(ra = self.ra * u.degree, dec = self.dec * u.degree, 
                               frame='icrs')
    
    def get_startend(self):
        """
        Calculate start and end positions if not provided
        """
        if self.start is None or self.end is None:
            self.start = coords.SkyCoord(ra = self.center.ra - self.diff_ra/2, dec = self.center.dec - self.diff_dec/2, frame='icrs')
            self.end = coords.SkyCoord(ra = self.center.ra + self.diff_ra/2, dec = self.center.dec + self.diff_dec/2, frame='icrs')
        else:
            self.start = coords.SkyCoord(self.start, unit=(u.hourangle, u.deg), frame='icrs')
            self.end = coords.SkyCoord(self.end, unit=(u.hourangle, u.deg), frame='icrs')
            
            self.diff_ra = np.abs(self.start.ra - self.end.ra).to(u.arcmin)
            self.diff_dec = np.abs(self.start.dec - self.end.dec).to(u.arcmin)
        
    def input_OTF_para(self, direction = '-', scan_gap = 21.66, scan_speed = 15, multibeamOTF = True):
        """
        Input OTF observation parameters
        
        Parameters:
        direction: scan direction, '-' or '|'
        scan_gap: scan gap in arcmin
        scan_speed: scan speed in arcsec/second
        multibeamOTF: whether it's multibeam OTF
        """
        self.direction = direction
        self.scan_gap = scan_gap * u.arcmin
        self.scan_speed = scan_speed * u.arcsec / u.second
        
        if direction == '-':
            self.switch_time = 90 * u.second if multibeamOTF else 18 * u.second
        elif direction == '|':
            self.switch_time = 54 * u.second if multibeamOTF else np.round(12 * scan_gap) * u.second
            
        print("MultiBeam OTF, direction: ", self.direction)
        print("Scan gap:", self.scan_gap)
        print("Scan speed:", self.scan_speed)
            
    def calculate_time(self, Print = True):
        """
        Calculate total observation time
        """
        
        if self.direction == '-':
            self.onePathTime = np.round((self.diff_ra / self.scan_speed).to(u.second))
            self.scanTimes = np.round(self.diff_dec / self.scan_gap + 1.0)
            self.switchTimes = np.round(self.diff_dec / self.scan_gap)
        elif self.direction == '|':
            self.onePathTime = np.round((self.diff_dec / self.scan_speed).to(u.second))
            self.scanTimes = np.round(self.diff_ra / self.scan_gap + 1.0)
            self.switchTimes = np.round(self.diff_ra / self.scan_gap)
            
        self.tot_time = self.onePathTime * self.scanTimes + self.switch_time * self.switchTimes
        if Print:
            print(f"Scan {self.scanTimes} times along {self.direction}, switch {self.switchTimes} times.")
            print(f"Need total {self.tot_time} = {self.tot_time.to(u.minute)}.")
        
        
    def sort_footprint(self, arr, first_ascending = True):
        """
        Sort the footprint array
        
        Parameters:
        arr: ra or dec array
        first_ascending: sorting order
        
        Returns:
        sorted array
        """
        if first_ascending:
            arr[1::2] = arr[1::2][:, ::-1]
        else:
            arr[0::2] = arr[0::2][:, ::-1]

        return arr
        
    def footprints_M01(self, sample_time, rotate_angle = None, head = 1, Print = True):
        """
        Calculate footprints for M01 mode
        
        Parameters:
        sample_time: sample time in second
        rotate_angle: rotation angle in degree
        head: scan start from +1 or -1 direction
        Print: whether to print information
        """
        self.sample_time = sample_time * u.second
        
        if self.direction == '-':
            rotate_angle = 23.4 if rotate_angle is None else rotate_angle
        elif self.direction == '|':
            rotate_angle = 53.4 if rotate_angle is None else rotate_angle
        
        self.rotate_angle = rotate_angle * u.deg
        if Print: 
            print("Sample time:", self.sample_time)
            print("Rotation angle:", self.rotate_angle)
        
        self.calculate_time(Print = False)
        
        onePathPoints = self.onePathTime / self.sample_time
        onePathStep = self.onePathTime / onePathPoints * self.scan_speed
        if self.direction == '-':
            ra1 = self.start.ra + head * np.arange(0, onePathStep.value * (onePathPoints + 1), onePathStep.value) * onePathStep.unit 
            dec1 = self.start.dec + head * np.arange(0, self.scan_gap.value * (self.switchTimes + 1), self.scan_gap.value) * self.scan_gap.unit
        elif self.direction == '|':
            dec1 = self.start.dec + head * np.arange(0, onePathStep.value * (onePathPoints + 1), onePathStep.value) * onePathStep.unit 
            ra1 = self.start.ra + head * np.arange(0, self.scan_gap.value * (self.switchTimes + 1), self.scan_gap.value) * self.scan_gap.unit  / np.cos(self.dec * u.degree)
        
        self.mesh1 = np.meshgrid(ra1.value, dec1.value)
        
        # sort direction
        if self.direction == '-':
            self.mesh1 = np.meshgrid(ra1.value, dec1.value)
            r = 0
            d = 1
        elif self.direction == '|':
            self.mesh1 = np.meshgrid(dec1.value, ra1.value)
            r = 1
            d = 0

        self.mesh1[0] = self.sort_footprint(self.mesh1[0], first_ascending = True)

        self.ra1 = self.mesh1[r].flatten()
        self.dec1 = self.mesh1[d].flatten()
        
        if self.direction == '|':
            self.renew_vertical_scan()
        
    def renew_vertical_scan(self, ):
        print("renew obs time and center")
        self.start = coords.SkyCoord(self.ra1.min(), self.dec1.min(), unit = (u.deg, u.deg))
        self.end = coords.SkyCoord(self.ra1.max(), self.dec1.max(), unit = (u.deg, u.deg))
        self.center = coords.SkyCoord((self.ra1.min() + self.ra1.max()) / 2, 
                          (self.dec1.min() + self.dec1.max()) / 2, 
                          unit = (u.deg, u.deg))
        
        self.diff_ra = np.abs(self.start.ra - self.end.ra).to(u.arcmin)
        self.diff_dec = np.abs(self.start.dec - self.end.dec).to(u.arcmin)
        self.calculate_time()
    
        
    def footprints_all(self,):
        """
        Calculate all footprints
        """
        from .utils.beamsep import beam_pos
            
        ras, decs = beam_pos(self.dec, rotate_angle = self.rotate_angle) # arcmin
        ras, decs = ras / 60, decs / 60
        
        self.ra_all = np.full(((19,len(self.ra1))), self.ra1) + np.full((len(self.ra1),19),ras).T
        self.dec_all = np.full(((19,len(self.dec1))), self.dec1) + np.full((len(self.dec1),19),decs).T
        
    def show_footprints(self, ax, opt, gap = 0, color = 'C0', ms = 1, alpha = 1, **kwargs):
        """
        Show footprints on the plot
        
        Parameters:
        ax: axes object
        opt: options
        gap: gap between points
        color: color of points
        ms: marker size
        alpha: alpha value
        """
        self.footprints_M01(**kwargs)
        self.footprints_all()
        
        for i in range(19):
            c = 'r' if i == 0 else color
            if gap > 0:
                px, py = opt.deg2pix(self.ra_all[i][::gap], self.dec_all[i][::gap],around = False)
            else:
                px, py = opt.deg2pix(self.ra_all[i], self.dec_all[i],around = False)

            ax.plot(px, py, '.', c = c, ms = ms, alpha = alpha)
        ax.set_xlabel('ra')    
        ax.set_ylabel('dec')  
        
        
    def plot_ZA_in_1day(self, BJ_time, sample_time, gap = 20):
        """
        Plot zenith angle in one day
        
        Parameters:
        BJ_time: Beijing time
        sample_time: sample time
        gap: gap between points
        
        Returns:
        za_min: minimum zenith angle
        ts_min: corresponding time
        """
        from .calculate_ZA import plot_za_in_1day
        
        za_min, ts_min = plot_za_in_1day(self, BJ_time, sample_time, gap)
        return za_min, ts_min