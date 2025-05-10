import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import hiviewer as hv
import os
from astropy import log

from astroquery.skyview import SkyView
from astropy import coordinates as coords
import astropy.units as u
from astropy.table import Table

class OptSource(object):
    """
    A class used to represent an optical source and provide related functionalities.
    """
    def __init__(self, name, ra, dec, vel):
        """
        Initialize the OptSource instance.

        Parameters:
        name (str): Name of the source.
        ra (float): Right Ascension in degree.
        dec (float): Declination in degree.
        vel (float): Velocity in km/s.
        """
        self.name = name
        self.ra = float(ra)
        self.dec = float(dec)
        self.vel = vel
        self.get_pos()
    
    @property
    def center(self):
        """
        Calculate the center position of the source.

        Returns:
        coords.SkyCoord: Center position in ICRS frame.
        """
        return coords.SkyCoord(ra = self.ra * u.degree, dec = self.dec * u.degree, frame='icrs')
    
    def get_pos(self):  
        """
        Calculate the position string of the source.
        """
        pos = self.center
        self.pos = pos.to_string()
        self.pos_hmsdms = pos.to_string('hmsdms')
        
    def get_opt_image(self, surveys = ["DSS2 Blue", "SDSSdr7g"], height = None, width = None, write = True):
        """
        Retrieve and save optical images of the source from SkyView.

        Parameters:
        surveys (list, optional): List of surveys to retrieve images from. Default is ["DSS2 Blue", "SDSSdr7g"].
        height (float, optional): Height of the image in degree. Default is None.
        width (float, optional): Width of the image in degree. Default is None.
        write (bool, optional): Whether to write the images to files. Default is True.

        Returns:
        list: List of file paths of the saved images.
        """
        self.surveys = surveys
        if write:
            sv = SkyView()
            print("Retrieving...")
            paths = sv.get_images(
                    position=f"{self.ra}, {self.dec}",
                    survey=surveys,
                    coordinates='J2000',
                    width = width,
                    height = height,
                    cache = False)
            
            optnames = []
            for j in range(len(surveys)):
                for i in range(len(paths[j])):
                    hdu_opt = paths[j][i]
                    if not os.path.exists('./tmp'):
                        os.mkdir('./tmp')
                    optname = './tmp/' + self.name + '_' +surveys[j] +f'_{i}-opt.fits'
                    hdu_opt.writeto(optname, overwrite = True)
                    print(f"Saved {optname}")
                    optnames += [optname,]
        else:
            from glob import glob
            paths = glob('./tmp/' + self.name + f'_*-opt.fits')
            if len(paths) > 0:
                print("Loading...")
                optnames = paths
            else:
                raise ValueError("Cannot find opt files in tmp, try 'write = True'")
        print("Finish!")
        self.optnames = optnames    
        
    def show_img(self,survey_index = 0,**kwargs):
        """
        Show the optical image of the source.

        Parameters:
        survey_index (int, optional): Index of the survey. Default is 0.

        Returns:
        ax: Matplotlib axes object.
        """
        j = survey_index
        bk = hv.FitsPic(self.optnames[j])
        ax = bk.plot_slice(**kwargs)
        return ax
    

    def plot_ZA(self, BJ_time = '2022-08-15 03:00:00'):
        """
        Plot the zenith angle of the source over a day.

        Parameters:
        BJ_time (str, optional): Beijing time in format 'YYYY-MM-DD HH:MM:SS'. Default is '2022-08-15 03:00:00'.

        Returns:
        tuple: Minimum zenith angle and corresponding time.
        """
        from .calculate_ZA import plot_za
        
        plot_za(BJ_time, self.pos_hmsdms)
        
    
    def plot_Sun_AngleDist(self, BJ_time = '2022-08-15 03:00:00'):
        """
        Plot the angular distance between the source and the sun over a year.

        Parameters:
        BJ_time (str, optional): Beijing time in format 'YYYY-MM-DD HH:MM:SS'. Default is '2022-08-15 03:00:00'.
        """
        from .calculate_ZA import plot_sun_AngleDist
        
        plot_sun_AngleDist(BJ_time, self.pos_hmsdms)
    
        
    def iESASky_url(self, fov = 5/60, survey= 'DSS2+color'):
        """
        Generate the ESASky URL for the source.

        Parameters:
        fov (float, optional): Field of view in degree. Default is 5/60.
        survey (str, optional): Survey name. Default is 'DSS2+color'.

        Returns:
        str: ESASky URL.
        """
        url = f"https://sky.esa.int/?target={self.ra}%20{self.dec}&hips={survey}&fov={fov}&cooframe=J2000&sci=true&lang=zh"
        print("Click this url to jump:", url)
        
    def iNED_cone_url(self, r = 3):
        """
        Generate the NED cone search URL for the source.

        Parameters:
        r (float, optional): Radius in arcmin. Default is 3.

        Returns:
        str: NED cone search URL.
        """
        ra, dec = self.pos_hmsdms.split(' ')
        url = f"https://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&in_csys=Equatorial&in_equinox=J2000&ra={ra}&dec={dec}&radius={r}"
        print("Click this url to jump:", url)

    def iNED_box_url(self, sizex = 3, sizey = 3):
        """
        Generate the NED box search URL for the source.

        Parameters:
        sizex (float, optional): Box half size x in arcmin. Default is 3.
        sizey (float, optional): Box half size y in arcmin. Default is 3.

        Returns:
        str: NED box search URL.
        """
        x = sizex * u.arcmin
        y = sizey * u.arcmin

        pos = self.center
        ra1 = str(np.around((pos.ra - x).to(u.hourangle).value, 4)) + 'h'
        ra2 = str(np.around((pos.ra + x).to(u.hourangle).value, 4)) + 'h'

        dec1 = str(np.around((pos.dec - y).to(u.degree).value, 4)) + 'd'
        dec2 = str(np.around((pos.dec + y).to(u.degree).value, 4)) + 'd'

        print('ra range', ra1, ra2)
        print('dec range', dec1, dec2)
        url = "https://ned.ipac.caltech.edu/byparams"
        print("Please copy the ranges above and then click this url to jump:", url)
        
        
    def ialadin(self, fov = 5/60, survey= 'DSS2+color',width='100%'):
        """
        Display the source in Aladin Lite.

        Parameters:
        fov (float, optional): Field of view in degree. Default is 5/60.
        survey (str, optional): Survey name. Default is 'DSS2+color'.
        width (str, optional): Width of the Aladin Lite widget. Default is '100%'.

        Returns:
        ipyaladin.Aladin: Aladin Lite widget.
        """
        import ipyaladin as ipyal
        from ipywidgets import Layout, Box, widgets

        self.aladin = ipyal.Aladin(target=self.pos, fov = fov, survey = survey,layout = Layout(width=width))
        
        return self.aladin      
    
    def cut_table(self, tab, r = None, v = None, center = None, keys = ['RA','DEC', 'vel']):
        """
        Cut a table based on radius and velocity.

        Parameters:
        tab: Table to cut.
        r (float, optional): Radius in arcmin. Default is None.
        v (float, optional): Velocity in km/s. Default is None.
        center: Center position. Default is None.
        keys (list, optional): Keys of RA, DEC, and velocity. Default is ['RA','DEC', 'vel'].

        Returns:
        Table: Cutted table.
        """
        if r is not None:
            ra = tab[keys[0]].value * u.degree
            dec = tab[keys[1]].value * u.degree
            c = coords.SkyCoord(ra, dec , frame='icrs')
            if center is None:
                cc = self.center
            else:
                cc = coords.SkyCoord(ra = center[0] * u.degree, dec = center[1] * u.degree,  
                               frame='icrs')
            use = (cc.separation(c) < r * u.arcmin)

        if (v is not None) & (len(keys)==3):
            velo = tab[keys[2]]
            vuse = (np.abs(velo - self.vel) < v) & (~velo.mask)
            use &= vuse
            
        tab = tab[use]

        return tab
    
    def select_types(self, tab):
        """
        Select specific types of objects from the table.

        Parameters:
        tab: Table to select from.

        Returns:
        Table: Selected table.
        """
        obj_type = tab['Type']
        plot_types = ['GPair', 'GTrpl', 'GGroup', 'GClstr', 'QSO', 'QGroup', 'G_Lens', 'Q_Lens', 'AbLS', 'EmLS',
                    'RadioS', 'SmmS', 'MCld']
        G_use = (obj_type == 'G')
        for ty in plot_types: 
            G_use |= (obj_type == ty)
        return tab[G_use]
    
    def retrive_table(self, r=3, v = None, survey = None, tablename = None, keys = None, center = None, **kwargs):
        """
        Retrieve and process a table from a survey.

        Parameters:
        r (float, optional): Radius in arcmin. Default is 3.
        v (float, optional): Velocity in km/s. Default is None.
        survey (str, optional): Survey name. Default is None.
        tablename (str, optional): Table name. Default is None.
        keys (list, optional): Keys of RA, DEC, and velocity. Default is None.
        center: Center position. Default is None.
        **kwargs: Additional keyword arguments.

        Returns:
        Table: Processed table.
        """
        print("retrive table", survey)
        if survey == 'sdss':
            from astroquery.sdss import SDSS
            table = SDSS.query_region(self.pos_hmsdms, radius = r * u.arcmin,verbose=True)
            if keys is None: keys = ['ra', 'dec']
            table = self.cut_table(table, r, v, center, keys)
        elif survey == 'ned':
            from astroquery.ipac.ned import Ned
            table = Ned.query_region(self.pos_hmsdms, radius = r * u.arcmin,verbose=True)
            if keys is None: keys = ['RA','DEC', 'Velocity']
            table = self.cut_table(table, r, v, center, keys)
            # limit galaxies
            table = self.select_types(table)
        elif survey == 'a100':
            from ..conf.common import conf_common
            a100catalog = conf_common['Common']['a100catalog']
            table = Table.read(a100catalog)
            if keys is None: keys = ['RAdeg_HI','DECdeg_HI', 'Vhelio']
            table = self.cut_table(table, r, v, center, keys)
        elif survey == 'bridges':
            if tablename is None:
                tablename = self.io.csvname
            table = Table.read(tablename)
            table = self.cut_table(table, r, v, center, keys)
        else:
            table = Table.read(tablename)
            table = self.cut_table(table, r, v, center, keys)
        
        table = self.get_sep(table, keys)
            
        if ('ned' in survey) or ('NED' in survey):
            # limit galaxies
            table = self.select_types(table)
        
        self.aladin.add_table(table, **kwargs)
        print("added table.")
        
        if not hasattr(self,'tables'):
            self.tables = {}
        self.tables[survey] = table
    
    def get_sep(self, table, keys = ['RA','DEC', 'Velocity']):
        """
        Calculate the separation between the source and objects in the table.

        Parameters:
        table: Table to calculate separation for.
        keys (list, optional): Keys of RA, DEC, and velocity. Default is ['RA','DEC', 'Velocity'].

        Returns:
        Table: Table with added separation columns.
        """
        ra, dec = self.ra, self.dec
        c1 = coords.SkyCoord(ra = ra * u.degree, dec = dec * u.degree, frame='icrs')
        c2 = coords.SkyCoord(ra = table[keys[0]].value * u.degree, dec = table[keys[1]].value * u.degree, frame='icrs')
        table.add_column(c1.separation(c2).to(u.arcmin), name = 'cood_sep', index = np.where(np.array(table.keys()) == keys[1])[0][0]+1)
        if len(keys) == 3:
            v1 = self.vel * u.km / u.s
            v2 = table[keys[2]].value * u.km / u.s

            table.add_column(v1 - v2,name = 'vel_sep', index = np.where(np.array(table.keys()) == keys[2])[0][0])
        
        table.sort('cood_sep')
        table.add_column(np.arange(len(table))+1,name = 'ID', index = 0)

        return table