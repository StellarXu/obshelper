import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u

def cal_rms1(Sigma_HI = 1 * u.M_sun * u.pc**(-2), theta = 15 * u.arcsec, W = 20 * u.km/u.s, dv = 5.6 * u.km/u.s):
    """
    Calculate the root mean square (rms) sensitivity based on the HI surface density and observational parameters.

    Parameters:
    Sigma_HI (Quantity): HI surface density (default: 1 M_sun/pc²)
    theta (Quantity): Angular resolution (default: 15 arcsec)
    W (Quantity): HI line width (default: 20 km/s)
    dv (Quantity): Velocity resolution (default: 5.6 km/s)

    Returns:
    Quantity: Calculated rms sensitivity in mJy
    """
    print(Sigma_HI, theta, W, dv)
    Sigma_HI = Sigma_HI.to(u.M_sun * u.Mpc**(-2)).value
    theta = theta.to(u.rad).value
    W = W.to(u.km/u.s).value
    dv = dv.to(u.km/u.s).value
    
    sigma_rms = (Sigma_HI * 1.133 * theta**2 ) / (2.35e5  * np.sqrt(W * dv)) *u.Jy
    return sigma_rms.to(u.mJy)

def cal_rms2(Sigma_HI = 1 * u.M_sun * u.pc**(-2), theta = 15 * u.arcsec, W = 20 * u.km/u.s, dv = 5.6 * u.km/u.s):
    """
    Alternative calculation of the root mean square (rms) sensitivity based on the HI column density and observational parameters.

    Parameters:
    Sigma_HI (Quantity): HI surface density (default: 1 M_sun/pc²)
    theta (Quantity): Angular resolution (default: 15 arcsec)
    W (Quantity): HI line width (default: 20 km/s)
    dv (Quantity): Velocity resolution (default: 5.6 km/s)

    Returns:
    Quantity: Calculated rms sensitivity in mJy
    """
    print(Sigma_HI, theta, W, dv)
    from astropy.constants import m_p
    m_H = m_p.to(u.g)
    Sigma_HI = Sigma_HI.to(u.g * u.cm**(-2))
    print("Column density", Sigma_HI / m_p.to(u.g))
    
    theta = theta.to(u.arcsec).value
    W = W.to(u.km/u.s).value
    dv = dv.to(u.km/u.s).value
   
    sigma_rms = (Sigma_HI.value * theta**2 ) / (1.1e24  * np.sqrt(W * dv) * m_H.value) * u.Jy
    return sigma_rms.to(u.mJy)